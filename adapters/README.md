# Solver adapters for the aldkin wall kinetics

The multiscale ALD wall (flat Langmuir surface plus a quasi-steady Knudsen/Langmuir hole per face, coverages
integrated per time step, `K_eff` returned for a Robin wall condition) lives in `lib/aldkin` as a C-ABI library.
Three hooks exist; a fourth is described.

| Hook | Where | Status |
|---|---|---|
| OpenFOAM boundary conditions `patternedSurface`, `langmuirSticking` | `reactor/lib/aldWallKinetics` (`Make/options` links `lib/aldkin/build/libaldkin.a`; `reproduce.sh --build`) | reference adapter; bit-identical to the pre-library boundary conditions on the regression case (201 write times, all fields, difference 0) and on the first 3 s of the long AR-100 case |
| Fluent UDF | `adapters/fluent/ald_wall_udf.c`, `adapters/fluent/udf_stub.h` | compiles and links here against the stub and `libaldkin` (`cc -std=c99 -Wall -Wextra -pedantic -I. -I../../lib/aldkin/include -c ald_wall_udf.c`; shared-object link check in the `library` CI job). **Not run inside Fluent.** |
| Plain C (any solver) | `lib/aldkin/examples/plain_c_example.c` | runs in `ctest`; reproduces the long case's recorded hole-bottom t90 from the 50 ms-sampled centre-face drive (14.01 vs 14.05 s) |

## OpenFOAM

The boundary conditions keep every dictionary keyword (`s r kdes vbar N0 c0 D d AR holeDensity nCells ARProfile
holeDensityProfile thetaFlat thetaHole theta value`) and every written field (`ARField holeDensityField thetaFlat
thetaBottom thetaMid Keff thetaHole value`; `theta value` for `langmuirSticking`). Each patch owns one handle
(`aldkinHandle.H`, RAII; copying a patch field clones the kinetics state). Per `updateCoeffs()`: on a new time index
`aldkin_step(h, T_face, deltaT)` with the face values of the previous step, otherwise `aldkin_evaluate`; then
`aldkin_mixed_bc(h, D*deltaCoeffs, refValue, refGrad, valueFraction)`. Face radii for the profiles are
`sqrt(Cf.y^2 + Cf.z^2)` (wedge axis x), as before. Build: `cmake` for `lib/aldkin` (static library, `-fPIC`,
`-ffp-contract=off` as wmake), then `wmake libso` with `ALDKIN_ROOT` pointing at `lib/aldkin` (default
`../../../lib/aldkin` relative to the adapter directory).

## Fluent (`adapters/fluent/ald_wall_udf.c`)

Hooks used and why:

- `DEFINE_INIT(ald_wall_init, d)`: creates the handle with `c0 = 1` (absolute number densities), loops the wall
  zone's faces (`Lookup_Thread(d, ALD_WALL_ZONE_ID)`, `begin_f_loop`) to collect the face-centre radii from
  `F_CENTROID` and calls `aldkin_set_faces`; optional radial tables through `aldkin_set_profiles`.
- `DEFINE_ADJUST(ald_wall_adjust, d)`: runs every iteration. It reads the precursor number density on each face
  (`F_R(f,t) * F_YI(f,t,species) / MW * 1000 N_A`), and on the first iteration of a new time step (`N_TIME` changed)
  calls `aldkin_step(h, n_face, CURRENT_TIMESTEP, NULL, keff)`; on the later iterations `aldkin_evaluate`. This is
  the same "advance once per step with the previous step's face value, re-evaluate within the step" sequence as the
  OpenFOAM adapter.
- `DEFINE_SR_RATE(ald_wall_sr_rate, f, t, r, mw, yi, rr)`: the wall sink as a surface reaction, `rr = K_eff * C_face`
  in kmol/(m2 s) with `C_face = rho Y / MW`. Fluent's species solver treats the wall surface reaction implicitly, which
  is the Robin condition the OpenFOAM mixed BC implements. Define one surface reaction `precursor -> site species`
  (or a dummy product) on the wafer wall in the Reactions panel and point its rate at this UDF.
- `DEFINE_PROFILE(ald_wall_flux, t, i)`: alternative for cases without the surface-reaction model: a specified
  species mass flux `-K_eff rho Y` kg/(m2 s) on the wall (explicit in the face value, first order in the step).
- `DEFINE_EXECUTE_AT_END(ald_wall_record)`: writes `thetaFlat`, `thetaMid`, `thetaBottom`, `K_eff` into user-defined
  face memory 0..3 (allocate 4 UDMs) for contours and monitors.
- `DEFINE_ON_DEMAND(ald_wall_free)`: releases the handle.

Edit the `ALD_*` defines at the top for the case (zone id, species index, molar mass, kinetics). Linking: add
`lib/aldkin/build/libaldkin.a` to the `libudf` makefile's link line (or `libaldkin.so/.dll` on the loader path with
`-laldkin`) and `-I<repo>/lib/aldkin/include` to the include path. Parallel: each compute node steps its own faces
(`begin_f_loop` over the local partition); the host does nothing (`#if !RP_HOST`). The face index `f` from
`begin_f_loop` is assumed to be `0..THREAD_N_ELEMENTS_INT-1` and stable between iterations; if a case violates that,
build an index map in `ald_wall_init`.

What has been verified: the file compiles with `-std=c99 -Wall -Wextra -pedantic` and links into a shared object
against `libaldkin` using the stub macros. What has not: anything inside Fluent (macro arguments, species indexing,
the surface-reaction panel wiring, parallel face order).

## Star-CCM+ (not written)

A user-code library (`uclib`, C or Fortran, `uclib.c` registering functions with `ucfunc`/`ucarg`) would need:

1. a function registered as a boundary profile (`"Boundary"` dimension, per-face result) returning either the species
   wall flux `K_eff c_face` or the coefficients of a Robin condition if the version exposes one; its `ucarg` inputs are
   the face centroids (for the radii), the face species density and the time step;
2. a hook that advances the coverages once per time step, for instance a second registered function evaluated
   first in each step (Star-CCM+ calls user functions on demand, so the advance must be guarded by the time-step
   index as `ald_wall_adjust` does with `N_TIME`);
3. the same `aldkin_create` / `aldkin_set_faces` setup at the first call, keyed by the boundary's face count, and
   `aldkin_theta` to fill user field functions for post-processing.

Everything else (units, parameter set, call order, restart via `aldkin_set_theta`) is identical to the Fluent adapter.
