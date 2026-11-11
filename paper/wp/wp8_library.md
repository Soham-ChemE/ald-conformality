# Work package 8: the wall kinetics as a standalone library (lib/aldkin)

Repository state: no commits made; everything below is in the working tree for the integrator. No file outside the
list below was edited; README.md, CLAIMS.md, paper/assumptions.md, paper/draft.md, main.tex and docs/ are untouched
(their text is proposed at the end of this document).

## What was built

The physics of the two OpenFOAM boundary conditions (`patternedSurface`: per-face quasi-steady 1D Knudsen-diffusion
hole with a Langmuir sink, coverage update once per time step, open-area blending of the flat and hole fluxes into
`K_eff`, thetaFlat/thetaMid/thetaBottom, radial ARProfile/holeDensityProfile tables, holeDensity 0 and AR 0 as flat
faces; `langmuirSticking`: flat Langmuir with recombination `r` and desorption `kdes`) now lives in `lib/aldkin`, a
C++17 library with a C ABI. The boundary conditions are thin adapters around it and give the same numbers bit for
bit. A Fluent UDF, a plain-C example, a ctypes binding and tests use the same library.

| File | Purpose |
|---|---|
| `lib/aldkin/include/aldkin.h` | the API: opaque `aldkin_t`, `aldkin_params` (s, r, kdes, vbar, N0, c0, d, AR, holeDensity, nCells, flatOnly, scheme), create/clone/destroy, `aldkin_set_faces`, `aldkin_set_profile`/`aldkin_set_profiles`, `aldkin_step`, `aldkin_evaluate`, `aldkin_mixed_bc`, `aldkin_theta`/`aldkin_set_theta`, `aldkin_keff`, `aldkin_desorption_flux`, `aldkin_hole_gas`, `aldkin_face_pattern`, `aldkin_last_error`, `aldkin_version`; return codes; Windows export macros (`ALDKIN_BUILD_DLL`/`ALDKIN_DLL`), `visibility("default")` elsewhere |
| `lib/aldkin/src/aldkin.cpp` | the physics, transcribed expression by expression from the boundary-condition sources (comments marked `[BC]` quote the original lines): same tridiagonal discretisation and Thomas solve of the hole, same forward-Euler update once per step with the same clamps, same blending `(1 - phi)(vbar/4)[s(1 - thetaFlat) + r] + phi j_hole`, same `kdes` term and `refValue = J_des/max(K, 1e-15)`, same handling of holeDensity 0 / AR 0 (no hole solved, phi 0, hole coverages stay 0), same profile interpolation and clamping. Adds a second update scheme (`ALDKIN_SCHEME_EXPONENTIAL`, the exact relaxation of `aldconf.solver`) selectable per handle; the default is the boundary conditions' Euler scheme. No exception crosses the ABI (`try/catch` at every entry), one thread-local message for a failed create |
| `lib/aldkin/CMakeLists.txt` | shared (`libaldkin.dylib`/`.so`/`aldkin.dll`, versioned 1.0.0) and static (`libaldkin.a`, `-fPIC`) libraries, `-ffp-contract=off` as wmake, hidden visibility, install rules, 11 ctest tests, the C example also registered as a test |
| `lib/aldkin/tests/test_aldkin.cpp` | the ctest suite (no framework; one test per command-line name): flat Langmuir vs the analytic solution (2.4e-14 over 10 s), Euler first-order convergence (ratios 2.000/2.000), desorption (literal BC step bit-identical, steady state `a/(a + kdes)`, `refValue`), AR-100 hole at constant density vs the Python solver's t90 (15.616 s, theta 0.900093212121, see below), literal 4-cell Euler step by hand (gas, coverages, `K_eff`, flux, mixed-BC coefficients all `==`), profile interpolation (AR 80 at 0.0875 m from the Night-3 table, clamps, step, negatives, independent tables, density 0 = flat face), zero density / zero holeDensity / AR 0 / nCells 1 / zero faces / dt 0, clone and restart, error codes |
| `lib/aldkin/examples/export_centre_face_drive.py` | exports the long case's innermost-face `value`, `thetaBottom`, `thetaFlat` at the 441 write times with `coupling/patterned_lib.load_patterned` into `centre_face_drive.csv` (parameter header line) |
| `lib/aldkin/examples/centre_face_drive.csv` | that export: `reactor/showerhead_pulse_patterned_AR100_long`, face 0 (r = 0.67 mm), 0..22 s every 50 ms, 20 s pulse; the BC's own thetaBottom reaches 0.9 at 14.05 s |
| `lib/aldkin/examples/plain_c_example.c` | plain C: one patterned face stepped at the solver's 1 ms step through the linearly interpolated 50 ms drive, prints coverages, `K_eff` and flux every second, then its hole-bottom t90 against the recorded one (14.008 vs 14.050 s; exits 1 if they differ by more than 0.1 s) |
| `lib/aldkin/README.md` | API table, parameters and units, the call sequence per solver step, the mixed-BC mapping, the three hooks, the test list |
| `reactor/lib/aldWallKinetics/aldkinHandle.H` | RAII owner of an `aldkin_t*` for the adapters: copy clones the state, `check()` turns a return code into a FatalError with `aldkin_last_error` |
| `reactor/lib/aldWallKinetics/patternedSurfaceFvPatchScalarField.{H,C}` | rewritten as the adapter: same four constructors, same dictionary keywords, same written entries in the same order; the state is in the library (face radii from `patch().Cf()`, profiles via `aldkin_set_profile`, restart via `aldkin_set_theta`); `updateCoeffs` calls `aldkin_step` on a new time index and `aldkin_evaluate` otherwise, then `aldkin_mixed_bc` with `D*deltaCoeffs` |
| `reactor/lib/aldWallKinetics/langmuirStickingFvPatchScalarField.{H,C}` | the same for the flat model (`flatOnly = 1`, `kdes`); `theta()` accessor kept through a mirrored field |
| `reactor/lib/aldWallKinetics/Make/options` | `ALDKIN_ROOT ?= ../../../lib/aldkin`; `-I$(ALDKIN_ROOT)/include`; links `$(ALDKIN_ROOT)/build/libaldkin.a` statically (the cases keep `libs ("libaldWallKinetics.dylib")` and need no extra library path). Verified that wmake resolves the default relative path without the variable |
| `reactor/lib/aldWallKinetics/legacy/*.orig`, `legacy/README.txt` | the original sources for the record (suffix `.orig` so wmake's lnInclude step does not link duplicates of the live names); not compiled |
| `adapters/fluent/ald_wall_udf.c` | Fluent UDF: `DEFINE_INIT` (create, faces from `F_CENTROID`), `DEFINE_ADJUST` (step on a new `N_TIME`, evaluate otherwise, face density from `F_R`/`F_YI`), `DEFINE_SR_RATE` (`K_eff rho Y/MW` kmol/(m2 s)), `DEFINE_PROFILE` (species flux alternative), `DEFINE_EXECUTE_AT_END` (coverages to UDM 0..3), `DEFINE_ON_DEMAND` (free). Picks `udf.h` when present, else the stub |
| `adapters/fluent/udf_stub.h` | the Fluent types and macros the UDF uses, inert, so `cc -std=c99 -Wall -Wextra -pedantic -I. -I../../lib/aldkin/include -c ald_wall_udf.c` compiles without warnings and `cc -shared ... -laldkin` links |
| `adapters/README.md` | the three hooks, the Fluent hook rationale and its unverified points, what a Star-CCM+ user-code hook needs |
| `src/aldconf/aldkin.py` | ctypes binding: `load_library` (ALDKIN_LIBRARY, ALDKIN_BUILD_DIR, `lib/aldkin/build[/Release]`, system loader), `is_available`, `Params` dataclass, `AldKin` class (set_faces, set_profile(s), step, evaluate, theta, set_theta, keff, mixed_bc, desorption_flux, hole_gas, face_pattern, clone, context manager), `drive()` mirroring `simulate_driven`'s sub-stepping |
| `tests/test_aldkin.py` | 9 tests; module-level skip with the build command when the library is absent (verified with `ALDKIN_BUILD_DIR=/nonexistent`) |
| `reproduce.sh` | `build_kinlib()`: cmake configure + build (`lib/aldkin/build`, log `lib/aldkin/log.cmake`) + ctest, before wmake; cmake from `ALDCONF_CMAKE`, PATH or `~/.local/bin`; `ALDKIN_ROOT` exported for wmake. `--build`, `--regress`, `--quick`, `--full` all go through it |
| `.github/workflows/ci.yml` | new `library` job (ubuntu-latest: apt cmake g++, cmake build, ctest, the Fluent stub compile and link check, `pip install -e '.[test]'`, `pytest tests/test_aldkin.py -rs`); the `openfoam` job installs cmake and g++ and builds the library through `reproduce.sh --build` before wmake, uploads `lib/aldkin/log.cmake`; the `python` job's test step prints skip reasons |
| `Dockerfile`, `.dockerignore`, `.gitignore` | `cmake g++` in the apt line (the build step `reproduce.sh --build` then builds the library); `lib/aldkin/build` and `log.cmake` excluded/ignored |

Toolchain: `cmake` was absent on the Mac; CMake 4.4.5 (Kitware's `cmake-4.4.5-macos-universal.tar.gz`, SHA-256
verified against Kitware's `cmake-4.4.5-SHA-256.txt`) was unpacked into `~/.local/opt` and linked as
`~/.local/bin/cmake` and `~/.local/bin/ctest`; `reproduce.sh` falls back to that location when `cmake` is not on
PATH (the launcher's PATH includes it here). Compiler: Apple clang 21 (`/usr/bin/c++`), the same compiler and the
same `-O3 -ffp-contract=off` as wmake.

## API in ten lines

```c
aldkin_params p; aldkin_params_default(&p);                 // s r kdes vbar N0 c0 d AR holeDensity nCells flatOnly scheme
aldkin_t* h = aldkin_create(&p);                            // NULL + aldkin_last_error(NULL) on bad parameters
aldkin_set_faces(h, n, radius_m);                           // faces by centre radius (NULL if no profiles)
aldkin_set_profile(h, ALDKIN_PROFILE_AR, nt, r, AR);        // optional radial tables (also _DENSITY; or aldkin_set_profiles)
aldkin_step(h, T_face, dt, flux, keff);                     // once per solver step: advance coverages, K_eff [m/s], flux [1/(m2 s)]
aldkin_evaluate(h, T_face, flux, keff);                     // further inner iterations: re-evaluate, no advance
aldkin_mixed_bc(h, D_over_delta, refValue, refGrad, vf);    // refValue = J_des/max(K,1e-15) (0 without kdes), refGrad 0, vf = K/(K + D/delta)
aldkin_theta(h, ALDKIN_THETA_BOTTOM, out);                  // FLAT | MID | BOTTOM per face, HOLE = n_faces x nCells; aldkin_set_theta for restarts
aldkin_hole_gas(h, face, nrel); aldkin_keff(h, out);        // relative hole gas profile; last K_eff
aldkin_destroy(h);                                          // all calls return ALDKIN_OK or an error code; aldkin_last_error(h) explains
```

`T_face` is the gas density over `c0` (`c0 = 1` for absolute densities); `scheme` is `ALDKIN_SCHEME_EULER` (the
boundary conditions) or `ALDKIN_SCHEME_EXPONENTIAL` (the `aldconf` solver's update).

## Bit-identity

Regression case `reactor/showerhead_pulse_layout_regress` (AR 50, holeDensity 1, 1 s pulse, 2 s, 150 wafer faces,
201 write times), run through `./reproduce.sh --regress` in the launcher (steady flow 2000 iterations, 20 s; pulse
~30 s; "3 passed in 2.48s"), compared with the run of the same chain made with the pre-library boundary conditions
earlier the same day (its 201 time directories were copied aside before the rerun; the steady `U` of both chains is
byte-identical):

| quantity | compared over | max abs difference |
|---|---|---|
| thetaFlat, thetaBottom, thetaMid, Keff, value (parsed with `coupling/patterned_lib.py`) | 201 times x 150 faces | 0 for every field |
| the `T` files themselves | all 201 time directories | byte-identical (0 of 201 differ) |
| thetaFlat, thetaBottom at t = 2 vs `tests/reference/regress_t2.csv` | 150 faces | 0 |
| value at t = 2 vs the reference CSV | 150 faces | 1.0e-9 (the CSV's 8 significant digits; identical to the pre-library run's deviation) |
| thetaBottom at t = 1, first 20 faces vs `regress_t1_thetaBottom_first20.csv` | 20 faces | 0 |

Long case `reactor/showerhead_pulse_patterned_AR100_long` (AR 100, 1.78e13 /m2, 20 s pulse): its `0/`, `constant/`
and `system/` were copied to a scratch directory, `endTime` set to 3 s, run with the adapter library (3000 steps,
61 write times) and compared with the existing time directories 0..3 s: thetaFlat, thetaBottom, thetaMid, Keff and
value differ by **0** at all 61 times and 150 faces (thetaBottom at 3 s: 8.76e-5 to 1.51e-4). The `T` files are not
byte-identical only because the existing run predates Night 3 and lacks the `ARField`/`holeDensityField` entries
the boundary condition has written since (the only lines `diff` reports).

Why it is exact: the adapter passes the same `T_face`, `deltaT`, `c0` and `N0` to expressions with the same operand
order, compiled by the same compiler with the same contract (`-ffp-contract=off`), and the per-step sequence (hole gas
from the old coverages, update, `K_eff` from the new flat coverage and the old hole flux; evaluate-only on the further
`updateCoeffs` calls of a time index) is preserved. The C++ test `hole_euler_literal_step` checks this sequence
against the lines evaluated by hand with `==`.

## Library against the Python twin

`aldconf.solver.simulate_driven` shares the hole discretisation but uses the exact exponential coverage update with
sub-steps; the library's `ALDKIN_SCHEME_EXPONENTIAL` is that update. Driven identically (5 sub-steps per 50 ms
sample, mid-interval interpolation, `T` chosen so `mean_speed(T, 72.09 amu) = 410 m/s`), the whole bottom and mouth
coverage histories agree within 1e-6 and the bottom t90 is 14.05 s in both (`tests/test_aldkin.py`). At constant
`T_face = 0.1` (3.79e19 /m3, 0.3 Pa at 573 K), 1 ms steps: Python t90 15.616 s with thetaBottom(15.616) =
0.900093212121; the library gives the same time and the coverage to 1e-6 (ctest); the Euler scheme gives 15.612 s
(0.03 % earlier). The boundary conditions' Euler scheme fed the 50 ms-sampled drive lands on the recorded bottom t90
within 0.05 s (14.008 vs 14.050 s; the recorded run saw the true 1 ms face history) and on the recorded coverages
within the sampling error (bottom within 0.05, flat within 0.02 at the write times).

## Fluent adapter status

`adapters/fluent/ald_wall_udf.c` compiles with `cc -std=c99 -Wall -Wextra -pedantic` without warnings against
`udf_stub.h` and links into a shared object with `libaldkin` (both also run in the `library` CI job). It has not been
loaded into Fluent: the zone id, species index, molar mass, the surface-reaction panel wiring, the face-index
assumption of `begin_f_loop` and the parallel partitioning are documented as the points to check
(`adapters/README.md`). The hook choice: `DEFINE_ADJUST` for the once-per-step advance (guarded by `N_TIME`) and
`DEFINE_SR_RATE` for the sink, because Fluent's wall surface reaction is the implicit Robin condition that the
OpenFOAM mixed BC implements; `DEFINE_PROFILE` is the explicit flux alternative.

## Tests

ctest (lib/aldkin/build): `100% tests passed, 0 tests failed out of 11` (Total Test time 0.19 s).

pytest, `PYTHONPATH=src .venv/bin/python -m pytest -q tests` with the regression case present and the library built:
`95 passed in 12.74s` (86 before this work package + 9 in `tests/test_aldkin.py`). With `ALDKIN_BUILD_DIR=/nonexistent`
the 9 skip with "libaldkin not built: cmake -S lib/aldkin -B lib/aldkin/build ... (or set ALDKIN_LIBRARY /
ALDKIN_BUILD_DIR)". `./reproduce.sh --regress` through the launcher: cmake + ctest (11 passed), wmake, steady flow,
regression pulse, "3 passed".

## What the integrator must check after the first CI run

1. `library` job: `apt-get install cmake g++` on ubuntu-latest and the cmake build; `ctest` output; the Fluent stub
   compile with GCC (clang was used here; `-pedantic` with `__has_include` in a `#if defined(__has_include)` guard is
   fine on GCC >= 5); the shared-object link `-L../../lib/aldkin/build -laldkin`.
2. `openfoam` job (container `opencfd/openfoam-default:2512`, root): `apt-get install ... cmake g++` succeeds in the
   image; `reproduce.sh --build` finds `cmake` on PATH; the static `libaldkin.a` (built with `-fPIC`) links into
   `libaldWallKinetics.so` under wmake's GCC; `reproduce.sh --regress` reaches "3 passed". GCC with `-std=c++17` does
   not contract floating-point operations and x86-64 has no FMA by default, so the Linux numbers should equal the
   reference as before; the 1e-6 tolerance of the regression test covers any residual.
3. The `python` job now prints skip reasons (`-rs`): expect 3 regression + 9 library skips (83 passed).
4. Docker: `docker build` runs `reproduce.sh --build` and therefore cmake, ctest and wmake in the image.
5. If the `openfoam` job's wmake cannot find the library, set `ALDKIN_ROOT` explicitly in the step (reproduce.sh
   already exports `$ROOT/lib/aldkin`).

## README text to insert (section "Use the kinetics in your own solver")

```
## Use the kinetics in your own solver

The wall model (flat Langmuir surface plus a Knudsen/Langmuir hole per face, coverages integrated per time step,
`K_eff` returned for a Robin wall condition) is a standalone C-ABI library, `lib/aldkin` (`cmake -S lib/aldkin -B
lib/aldkin/build && cmake --build lib/aldkin/build`; `ctest` runs 11 checks). Three hooks exist:
- OpenFOAM: `reactor/lib/aldWallKinetics` are thin adapters around it (`patternedSurface`, `langmuirSticking`;
  bit-identical to the earlier boundary conditions on the regression case).
- Fluent: `adapters/fluent/ald_wall_udf.c` (DEFINE_INIT / DEFINE_ADJUST / DEFINE_SR_RATE); compiled against a stub,
  not yet run in Fluent.
- Any other code: `lib/aldkin/examples/plain_c_example.c` steps one face through a recorded chamber density history;
  Python through `aldconf.aldkin` (ctypes). API, units and the per-step call sequence: `lib/aldkin/README.md`.
```

## CLAIMS rows (CLAIMS.md format)

```
| Boundary conditions vs standalone library lib/aldkin (regression case, 201 write times x 150 faces) | thetaFlat, thetaBottom, thetaMid, Keff, value: max difference 0; T files byte-identical | paper/wp/wp8_library.md | ./reproduce.sh --regress (launcher), tests/test_regression_reference.py |
| Library vs the earlier run of the long AR-100 case, first 3 s | all five wafer fields identical at 61 write times (difference 0) | paper/wp/wp8_library.md | copy of reactor/showerhead_pulse_patterned_AR100_long with endTime 3 |
| Library hole model vs aldconf.solver.simulate_driven on the centre-face drive | coverage history within 1e-6; bottom t90 14.05 s in both | tests/test_aldkin.py | pytest tests/test_aldkin.py (library built) |
| Library flat Langmuir vs 1 - exp(-s vbar/4 n t/N0) | max difference 2.4e-14 over 10 s (exponential scheme); Euler scheme first order (error ratios 2.000) | lib/aldkin/tests/test_aldkin.cpp | ctest --test-dir lib/aldkin/build |
| AR-100 hole at constant 0.3 Pa (T_face 0.1), 1 ms steps | bottom t90 15.616 s (library and Python twin), 15.612 s with the Euler scheme | lib/aldkin/tests/test_aldkin.cpp | same |
| Plain-C example on the 50 ms-sampled centre-face drive | bottom t90 14.008 s vs 14.05 s recorded by the boundary condition | lib/aldkin/examples/plain_c_example.c | ctest (aldkin.plain_c_example) |
```

## Run log for paper/assumptions.md

```
## Night 4, part 1: standalone kinetics library (lib/aldkin; reactor/lib/aldWallKinetics as adapters; adapters/fluent; src/aldconf/aldkin.py)
The wall kinetics of patternedSurface and langmuirSticking were moved into lib/aldkin, a C++17 library with a C ABI
(include/aldkin.h: aldkin_create with the parameter struct s, r, kdes, vbar, N0, c0, d, AR, holeDensity, nCells,
flatOnly, scheme; aldkin_set_faces by face-centre radius; aldkin_set_profile for the ARProfile/holeDensityProfile
tables; aldkin_step once per solver step with the face values and dt, returning K_eff and the uptake flux;
aldkin_evaluate for the further iterations of a step; aldkin_mixed_bc giving refValue = J_des/max(K, 1e-15), refGrad 0,
valueFraction = K_eff/(K_eff + D/delta); aldkin_theta for flat/mid/bottom/whole hole; errors as return codes). The
arithmetic was transcribed line by line in the same operand order; the boundary conditions became adapters that keep
every keyword, written field and constructor, link lib/aldkin/build/libaldkin.a statically (Make/options; reproduce.sh
--build runs cmake and ctest before wmake), and the old sources sit in reactor/lib/aldWallKinetics/legacy/ (.orig).
Verification: showerhead_pulse_layout_regress rerun through ./reproduce.sh --regress gives thetaFlat, thetaBottom,
thetaMid, Keff and value identical to the pre-library run at all 201 write times and 150 faces (max difference 0; the
201 T files are byte-identical) and passes tests/test_regression_reference.py (3 passed); the first 3 s of
showerhead_pulse_patterned_AR100_long rerun in a scratch copy give all five fields identical to the existing time
directories at the 61 write times (difference 0). The library's second update scheme (exponential, the aldconf
solver's) reproduces aldconf.simulate_driven on the long case's centre-face drive within 1e-6 over the whole coverage
history (bottom t90 14.05 s in both) and the Python t90 of 15.616 s for an AR-100 hole at constant 0.3 Pa; the
boundary conditions' Euler scheme gives 15.612 s there. A flat face reproduces 1 - exp(-s vbar/4 n t/N0) to 2.4e-14.
A Fluent UDF (DEFINE_INIT/DEFINE_ADJUST/DEFINE_SR_RATE/DEFINE_PROFILE/DEFINE_EXECUTE_AT_END) compiles and links
against a stub header and libaldkin but has not been run in Fluent. A plain-C example steps one face through the
recorded 50 ms centre-face drive and reaches bottom t90 at 14.008 s (14.05 s recorded). ctest: 11 passed; pytest:
95 passed (9 new in tests/test_aldkin.py, which skip with the build command when the library is absent). CI: new
library job (cmake, ctest, Fluent stub check, binding tests); the openfoam job and the Dockerfile build the library
before wmake. CMake 4.4.5 installed from Kitware's macOS binary into ~/.local/bin (not on the Mac before).
Caveats: kdes is available in the flat-only model only (the patterned model had no desorption term and none was
invented); the library's step uses the host's dt as the coverage sub-step, so a host with a coarse step should call
aldkin_step several times per step or use the exponential scheme; Fluent and Star-CCM+ hooks are untested in those
codes.
```
