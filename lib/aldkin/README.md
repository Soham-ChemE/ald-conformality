# aldkin: multiscale ALD wall kinetics as a standalone library

`libaldkin` is the physics of the two OpenFOAM boundary conditions in `reactor/lib/aldWallKinetics` as a C++17
library with a C ABI, so that any CFD solver (OpenFOAM, Fluent, Star-CCM+, an in-house code) or a Python script can
put the same self-limiting wall on its wafer faces. Every wall face carries

- a flat Langmuir surface with coverage `thetaFlat`, and
- optionally a quasi-steady 1D Knudsen-diffusion model of a cylindrical hole (diameter `d`, aspect ratio `AR`,
  `nCells` cells, Langmuir sink on the wall; Ylilammi et al. 2018) whose mouth sees the face gas density `c0 T_face`.

Net uptake per unit wafer area is `j = K_eff c0 T_face` with

    K_eff = (1 - phi) (vbar/4) [ s (1 - thetaFlat) + r ]  +  phi * j_hole / (c0 T_face),     phi = holeDensity pi d^2/4,

`j_hole` the Knudsen flux into the hole mouth. The host solver gets `K_eff` (m/s) per face and uses it as a Robin
(mixed) wall condition; the library integrates the coverages. A face with `holeDensity 0` or `AR 0` is a flat
Langmuir surface. With `flatOnly = 1` the hole model is absent and a first-order desorption rate `kdes` may be set
(the `langmuirSticking` model: `dtheta/dt = (vbar/4) s (1 - theta) c0 T / N0 - kdes theta`).

The arithmetic is transcribed expression by expression from the boundary conditions (now in
`reactor/lib/aldWallKinetics/legacy/`), in the same evaluation order and with the same floating-point contract
(`-ffp-contract=off`, as wmake). A host that calls `aldkin_step` once per time step with the previous step's face value
reproduces them **bit for bit**: the regression case `reactor/showerhead_pulse_layout_regress` rerun through the
adapters gives `T` files byte-identical to the pre-library run at all 201 write times (see `paper/wp/wp8_library.md`).

## Build

    cmake -S lib/aldkin -B lib/aldkin/build -DCMAKE_BUILD_TYPE=Release
    cmake --build lib/aldkin/build
    ctest --test-dir lib/aldkin/build --output-on-failure

Products in `lib/aldkin/build`: `libaldkin.dylib` / `libaldkin.so` / `aldkin.dll` (shared, versioned 1.1.0) and
`libaldkin.a` (static, position independent; the OpenFOAM adapter links this one so the cases need no extra library
path). `include/aldkin.h` is the whole API. Options `-DALDKIN_BUILD_TESTS=OFF`, `-DALDKIN_BUILD_EXAMPLES=OFF`;
`cmake --install` puts headers and libraries under the usual prefix. Windows: the shared build defines
`ALDKIN_BUILD_DLL` (dllexport); a consumer of the DLL defines `ALDKIN_DLL` (dllimport); the static library needs
neither. Requirements: CMake >= 3.15, a C++17 compiler (tested: Apple clang 21, GCC on Ubuntu 24.04). `reproduce.sh
--build` runs this build and ctest before compiling the OpenFOAM adapters.

## API (include/aldkin.h)

| Function | Purpose |
|---|---|
| `aldkin_version()` | `"1.1.0"` |
| `aldkin_params_default(&p)` | the boundary conditions' defaults (s 1, r 0, kdes 0, vbar 410, N0 4e18, c0 1, d 1e-7, AR 50, holeDensity 0, nCells 100, flatOnly 0, scheme Euler; two-species fields nu, vbarB, sB, kB, kEtch, cB0 all 0) |
| `aldkin_create(&p)` / `aldkin_clone(h)` / `aldkin_destroy(h)` | handle lifetime; `create` returns NULL on invalid parameters (`aldkin_last_error(NULL)` says why) |
| `aldkin_set_faces(h, n, radius_m)` | define `n` faces by their centre radius (m); `radius_m` may be NULL when no profiles are used; coverages reset to 0 |
| `aldkin_set_profile(h, which, n, r, v)` | radial table for `ALDKIN_PROFILE_AR` or `ALDKIN_PROFILE_DENSITY` (independent radii, as the OpenFOAM `ARProfile` / `holeDensityProfile`); linear between points, clamped outside, repeated radius = step, negatives clamp to 0; `n = 0` removes it |
| `aldkin_set_profiles(h, n, r, AR, density)` | both tables on one radius vector (NULL keeps the scalar) |
| `aldkin_step(h, T_face, dt, flux, keff)` | advance every face by `dt` with `T_face[i] = n_i / c0`, then evaluate; outputs the net uptake flux (molecules/(m2 s)) and `K_eff` (m/s), either may be NULL |
| `aldkin_evaluate(h, T_face, flux, keff)` | re-evaluate `K_eff` for new face values without advancing (inner iterations) |
| `aldkin_mixed_bc(h, D_over_delta, refValue, refGrad, valueFraction)` | the mixed-BC coefficients, see below |
| `aldkin_theta(h, which, out)` / `aldkin_set_theta(h, which, in)` | coverages: `ALDKIN_THETA_FLAT`, `_MID` (cell nCells/2), `_BOTTOM` (cell nCells-1) per face, `_HOLE` the whole profile (n_faces x nCells, cell index fastest); `_BLOCKED_FLAT`, `_BLOCKED_BOTTOM`, `_BLOCKED_HOLE` the sites blocked by the byproduct (two-species option); set FLAT, HOLE, BLOCKED_FLAT or BLOCKED_HOLE for restarts |
| `aldkin_set_byproduct_mouth(h, nB_face)`, `aldkin_byproduct_flux(h, out)`, `aldkin_byproduct_balance(h, face, out5)`, `aldkin_hole_byproduct_gas(h, face, out)` | two-species option: per-face byproduct density at the face and hole mouth (1/m3, NULL = the scalar `cB0`); net byproduct flux into the chamber per unit wafer area; the hole's bookkeeping rates produced, desorbed, escaped, re-adsorbed, etched per unit mouth area; the byproduct density along the hole |
| `aldkin_keff`, `aldkin_desorption_flux`, `aldkin_hole_gas(h, face, nrel)`, `aldkin_face_pattern`, `aldkin_n_faces`, `aldkin_n_cells`, `aldkin_get_params` | read-back of the last evaluation and the setup |
| `aldkin_last_error(h)` | message of the last failed call on `h` |

Return codes: `ALDKIN_OK 0`, `ALDKIN_ERR_ARG 1` (NULL pointer, bad size/index, non-finite input), `ALDKIN_ERR_PARAM 2`,
`ALDKIN_ERR_STATE 3` (step before `set_faces`), `ALDKIN_ERR_ALLOC 4`, `ALDKIN_ERR_INTERNAL 5`. No exception crosses the
ABI. Handles are independent; one handle must not be used from two threads at once (the OpenFOAM adapter owns one per
patch, a Fluent compute node one per partition).

### Parameters (`aldkin_params`) and units

| Field | Meaning | Unit |
|---|---|---|
| `s` | sticking coefficient on an open site | - |
| `r` | loss probability on every collision regardless of coverage (radical recombination; 0 for thermal ALD) | - |
| `kdes` | first-order desorption rate of the adsorbed layer; `flatOnly = 1` only (the patterned model has no desorption term and refuses `kdes != 0`) | 1/s |
| `vbar` | mean thermal speed sqrt(8 k T / pi m) of the precursor | m/s |
| `N0` | adsorption sites filled by one saturated half-cycle | 1/m2 |
| `c0` | number density at `T_face = 1`; pass `c0 = 1` to work with absolute densities in 1/m3 | 1/m3 |
| `d`, `AR`, `holeDensity` | hole diameter, depth/diameter, holes per m2 of wafer (wafer-wide defaults; the profiles override them face by face) | m, -, 1/m2 |
| `nCells` | cells along the hole (the chamber runs use 200) | - |
| `flatOnly` | 1: flat Langmuir surface only | - |
| `scheme` | `ALDKIN_SCHEME_EULER` (the boundary conditions: forward Euler once per step, clamped to [0, 1]) or `ALDKIN_SCHEME_EXPONENTIAL` (the Python twin `aldconf.solver`: `theta <- 1 - (1 - theta) exp(-a dt)`, exact at frozen gas, stable for any `dt`) | - |
| `nu` | two-species option: byproduct molecules released per adsorbed precursor (3 CH4 per TMA over an Al2O3 cycle; 4 HCl per TiCl4 over a TiN cycle, about 0.5 of them in the TiCl4 half-step) | - |
| `vbarB` | mean thermal speed of the byproduct (sets its Knudsen diffusivity `d vbarB/3`); required > 0 when `sB > 0` or `kEtch > 0` | m/s |
| `sB` | byproduct sticking on a free site (reversible site blocking); 0 (default) with `kEtch 0` keeps the single-species path bit for bit | - |
| `kB` | desorption rate of a blocked site (0 = irreversible blocking) | 1/s |
| `kEtch` | probability that a byproduct collision removes a precursor-covered site (0 = no etching) | - |
| `cB0` | byproduct number density at the face and hole mouth, absolute (0 = the chamber sweeps it away; `aldkin_set_byproduct_mouth` overrides it per face, for a host that transports the byproduct) | 1/m3 |

Two-species option (version 1.1, `aldconf.multispecies` is the Python twin): coverages `theta` (precursor, irreversible)
and `thetaB` (blocked) per cell with the free fraction `1 - theta - thetaB` in the precursor sink and in both
adsorption terms; the byproduct gas is quasi-steady in the hole with the source `(4/d)[nu x precursor uptake +
kB N0 thetaB]` and the sink `(4/d)(vbarB/4)[sB free + kEtch theta]`, mouth density `cB0`. Exponential scheme: exact
closed-form solution of the linear 2 x 2 coverage system at frozen gas; Euler scheme: one forward-Euler step of both,
clamped. The boundary conditions do not read these fields yet (next step: dictionary keywords `nu`, `vbarB`, `sB`,
`kB`, `kEtch`, `cB0` and a transported byproduct field for `aldkin_set_byproduct_mouth`).

Gas densities go in and come out relative to `c0` (`T_face = n/c0`); negative face values are treated as 0. `dt` in s.

### Call sequence per solver time step

    aldkin_params p; aldkin_params_default(&p); p.s = 5e-3; ... ;
    aldkin_t* h = aldkin_create(&p);
    aldkin_set_faces(h, n_faces, radius);                 // once; radii only matter with profiles
    aldkin_set_profile(h, ALDKIN_PROFILE_AR, nAR, rAR, AR);   // optional

    for each time step {
        aldkin_step(h, T_face_previous_step, dt, flux, keff);        // first call of the step: advance + evaluate
        for each further inner iteration of the host (optional):
            aldkin_evaluate(h, T_face_current, flux, keff);          // same coverages, new face values
        aldkin_mixed_bc(h, D_over_delta, refValue, refGrad, valueFraction);   // or use keff directly
        ... host applies the wall condition, solves the transport, writes aldkin_theta(...) if wanted
    }
    aldkin_destroy(h);

Order inside `aldkin_step`, as in the boundary conditions: the hole gas is solved with the coverages before the
update; the coverages are updated with that gas; `K_eff` is formed from the updated flat coverage and the pre-update
hole flux. The flat coverage uses the face value directly, the hole cells the face value times the relative hole gas
profile. The host's own `dt` is the sub-step; a host with a coarse step can call `aldkin_step` several times with
`dt/k` between its own steps (with the exponential scheme any `dt` is stable).

### Mixed-BC mapping (what the OpenFOAM boundary conditions use)

The wall balance for the transported scalar `T` (density over `c0`) is

    K_eff T_face - J_des/c0 = D (T_cell - T_face) / delta

with `J_des = kdes theta N0` (flat-only; 0 otherwise). In OpenFOAM's `mixed` form
`T_face = f refValue + (1 - f) (T_cell + refGrad/deltaCoeffs)`:

    refValue      = (J_des/c0) / max(K_eff, 1e-15)        (0 without desorption)
    refGrad       = 0
    valueFraction = K_eff / (K_eff + D/delta)             with D/delta = D * deltaCoeffs (the face's D/delta per face)

`aldkin_mixed_bc` returns exactly these three arrays for a `D_over_delta` array; the adapter copies them into
`refValue()`, `refGrad()`, `valueFraction()`. A solver without a mixed/Robin condition can use the equivalent flux
form: wall flux `= K_eff c_face - J_des` with `c_face` from the previous iteration (explicit, first order), or
linearise it as `K_eff` times the unknown face value in its own matrix.

## Hooking it into a solver

- **OpenFOAM**: `reactor/lib/aldWallKinetics` is the reference adapter. `patternedSurfaceFvPatchScalarField` and
  `langmuirStickingFvPatchScalarField` own an `aldkinHandle` (RAII wrapper in `aldkinHandle.H`), pass
  `patch().Cf()` radii at construction, call `aldkin_step` on the first `updateCoeffs()` of a time index and
  `aldkin_evaluate` on the later ones, and `aldkin_mixed_bc` for the coefficients. `Make/options` adds
  `-I$(ALDKIN_ROOT)/include` and links `$(ALDKIN_ROOT)/build/libaldkin.a` (`ALDKIN_ROOT` defaults to
  `../../../lib/aldkin`; `reproduce.sh --build` runs cmake first, then wmake). The cases are unchanged
  (`libs ("libaldWallKinetics.dylib")`, every dictionary keyword and written field as before).
- **Fluent**: `adapters/fluent/ald_wall_udf.c`: `DEFINE_INIT` (create, faces from `F_CENTROID`), `DEFINE_ADJUST`
  (`aldkin_step` on the first iteration of a time step from `N_TIME`, `aldkin_evaluate` otherwise, face density
  from `F_R` and `F_YI`), `DEFINE_SR_RATE` (surface reaction rate `K_eff rho Y / MW` in kmol/(m2 s)) or
  `DEFINE_PROFILE` (species flux), `DEFINE_EXECUTE_AT_END` (coverages into UDM 0..3). Compiled and linked here against
  `udf_stub.h`; not yet run in Fluent. Link `libaldkin.a` into `libudf`.
- **Any other solver (C/C++/Fortran/Python)**: `examples/plain_c_example.c` steps one patterned face through the
  long chamber run's recorded centre-face density (`examples/centre_face_drive.csv`, exported by
  `examples/export_centre_face_drive.py`) with the solver's 1 ms step and reproduces the recorded hole-bottom t90
  (14.01 s from the 50 ms-sampled drive vs 14.05 s recorded). Fortran: `bind(C)` interfaces to the same symbols.
  Python: `src/aldconf/aldkin.py` (ctypes; `AldKin`, `Params`, `drive`), tested in `tests/test_aldkin.py`.
- **Star-CCM+**: a user-code library (C/Fortran `uclib`) registering a boundary profile function; see
  `adapters/README.md` for what it needs.

## Tests

`ctest` (16 tests, `tests/test_aldkin.cpp` plus the C example):

| test | checks |
|---|---|
| `flat_langmuir_analytic` | exponential scheme vs `1 - exp(-s vbar/4 n t/N0)`: max difference 2.4e-14 over 10 s at 1 ms steps (bound 1e-9) |
| `flat_langmuir_euler_convergence` | the Euler scheme's error halves with dt (ratios 2.000, 2.000) |
| `flat_langmuir_desorption` | one Euler step equals the `langmuirSticking` line bit for bit; steady state `a/(a + kdes)`; `refValue = J_des/K`; `kdes` refused on a patterned handle |
| `hole_constant_density_t90` | AR-100 hole at `T_face = 0.1` (0.3 Pa at 573 K), exponential scheme, 1 ms steps: bottom t90 15.616 s and `thetaBottom(15.616 s) = 0.900093212121`, the numbers of `aldconf.solver.simulate_driven` (to 1e-6); Euler 15.612 s (within 1 %) |
| `hole_euler_literal_step` | one step of a 4-cell hole equals the boundary condition's lines evaluated by hand (gas, coverages, `K_eff`, flux, mixed-BC coefficients), bit for bit |
| `profiles_interpolate` | AR 80 at r = 0.0875 m from the Night-3 table, clamping, step at a repeated radius, negatives to 0, independent tables, density 0 gives a flat face |
| `zero_density_and_zero_holes_are_safe` | T = 0 and T < 0, `holeDensity 0`, `AR 0`, `nCells 1`, zero faces, `dt = 0` |
| `clone_and_restart`, `error_codes`, `version_and_defaults`, `plain_c_example` | state copy and restart continue bit-identically; every error path returns a code and leaves the handle usable; the C example runs through the recorded drive |
| `two_species_defaults_are_bit_identical` | `sB = kEtch = 0` with any `nu`, `vbarB`, `kB`, `cB0`: coverages, `K_eff` and flux `==` the handle without the fields, both models, both schemes |
| `two_species_fast_desorption_limit` | `kB = 1e14`: precursor coverages within 1e-9 of the single-species handle, blocked coverage below 1e-9 |
| `two_species_mass_balance` | produced + desorbed = escaped + re-adsorbed + etched after every step to 1e-9 (with and without etching, with and without a mouth byproduct density); the per-face byproduct flux is the open-area blend; the per-face mouth override |
| `two_species_blocking_slows_the_bottom` | TiCl4/HCl-like hole at 0.3 Pa: re-adsorption lowers the bottom coverage at 12 s, etching lowers it further, Euler and exponential schemes agree within 0.01 |
| `two_species_errors` | `vbarB = 0` with `sB > 0`, out-of-range `sB`, `kEtch`, `nu`, `cB0`, mouth override before `set_faces`, bad face indices |

Python (`tests/test_aldkin.py`, 9 tests, skipped with the build command when the library is absent): the analytic
flat face; profiles; the hole model against `aldconf.solver.simulate_driven` on the long case's centre-face drive
(whole coverage history within 1e-6, bottom t90 identical at 14.05 s); the Euler scheme against the coverages the
OpenFOAM run recorded; the mixed-BC mapping; zero density; clone/restart; errors as exceptions.
