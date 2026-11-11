# Work package 9: a 3D voxel Monte Carlo feature solver (lib/aldmc)

Repository state: no commits made; everything below is in the working tree for the integrator. README.md, CLAIMS.md,
paper/assumptions.md, paper/draft.md, main.tex and docs/ are untouched (their text is proposed at the end of this
document).

## What was built

A second, independent feature-scale solver: a 3D voxel grid with test-particle Monte Carlo in the molecular-flow
regime, written as a C++17 library with a C ABI (`lib/aldmc`), parallel with `std::thread`, with a ctypes binding and
geometry builders in `src/aldconf/voxel.py` and a driver that mirrors `aldconf.simulate_driven`. It is the ballistic
counterpart of the repository's 1D Knudsen-diffusion/Langmuir solver: no diffusion equation, no hydraulic-diameter
assumption, arbitrary 3D shapes from an occupancy array (straight and tapered holes, slits, the word-line recess
structure), cross-checked against the 1D solver, against Argonne's Machball (ballistic view factors) and, driven by
the chamber's centre-face density history, against the 14.05 s of the long AR-100 case.

### Physics (as stated in lib/aldmc/README.md and include/aldmc.h)

Molecular flow: particles move in straight lines between wall collisions and never collide with each other. They
enter through the feature mouth (the open plane z = 0) with a cosine (Lambertian) angular distribution and a uniform
position over the open mouth voxels. On hitting a solid voxel face they stick with probability s (1 − θ_face)
(irreversible Langmuir adsorption on the face's own site coverage), are lost with probability r (radical
recombination, default 0), and otherwise re-emit diffusely (cosine law about the face normal). Particles leaving
through the mouth are lost to the gas phase ("escaped"). Coverage lives on the exposed voxel faces (site density N0
times the face area = sites per face). Film growth does not narrow the feature (not implemented; the 1D solver's
`narrowing` has no counterpart here).

### Discretisation, weighting and the coverage update

- Grid: nx × ny × nz voxels of size h, occupancy 1 = solid, 0 = gas, z increasing into the feature, mouth plane z = 0.
  Lateral directions can be periodic (slits). Exposed faces are the gas-voxel faces adjacent to a solid voxel; a
  floor face has direction +z. Flights use an exact voxel traversal (Amanatides–Woo) accelerated by a Chebyshev
  distance field: a particle in a voxel at distance D from the nearest solid or boundary flies D − 1 voxel lengths
  without checks.
- Time stepping and weight: the mouth flux is n(t) vbar/4 per unit mouth area; one `aldmc_step(h, n_mouth, dt, M)`
  launches M particles that represent n vbar/4 · A_mouth · dt molecules, each carrying w = n vbar/4 A_mouth dt / M
  molecules.
- Update (the choice): the coverage field is frozen while the M particles of a step fly (their fates are sampled with
  the coverage at the start of the step) and the wall collisions are counted per face. The default update applies
  the exact Langmuir relaxation of the impinged dose, θ_f ← 1 − (1 − θ_f) exp(−s · hits_f · w /(N0 h²)), which is the
  1D solver's `theta <- 1 - (1 - theta) exp(-s vbar/4 n dt / N0)` with the Monte Carlo impingement in place of
  n vbar/4 dt. It uses every wall collision, so its variance is about 1/s (two hundred times at s = 5e-3) smaller than
  counting sticking events. The literal alternative is also implemented (`update = "deposit"`: every sampled sticking
  event deposits w molecules on its face, θ capped at 1). Both are first order in dt: the frozen start-of-step
  coverage makes too many particles stick near the mouth and too few reach the depth, so t90 grows linearly with dt
  (table b); the headline runs use dt = t90/200, where the residual is about 0.6 %. Mass is bookkept per particle
  (every particle ends escaped, stuck, lost or left; the counts add up to M in every run).
- Staircase area correction: a voxel surface is a staircase; the voxel faces of a cylinder have 4/π times the area of
  the smooth wall, so a plain voxel solver carries 27 % too many sites and saturates 26 % too late (the `plainfaces`
  run below). A smooth surface normal n is estimated for each exposed face from the first moment of the occupancy
  around the face centre (cone kernel 1 − r/R, R = 5 voxels; a flat-topped kernel over-counts the stair corners by
  3–7 %), and the face gets the weight w_face = |n · e_face| (its projected area fraction): it carries N0 h² w_face
  sites, and sticking and loss per collision are s w_face (1 − θ) and r w_face. In an equilibrium gas every face is hit
  at n vbar/4 h² whatever its orientation, so the uptake per unit smooth area and the rate per site are the smooth
  surface's. Re-emission stays diffuse about the actual facet normal, which is what keeps those hit rates (detailed
  balance); re-emitting about the smooth normal instead was tried and rejected: it under-counts the uptake by 10 % at
  θ = 0 and makes the AR-20 hole 11 % slower. The staircase does not slow the axial transport: the Clausing
  transmission of open voxel tubes is within 1 % of Berman's values (L/d 10: 0.108 vs 0.1093; L/d 20: 0.059 vs 0.0595,
  at 16 and 32 voxels across alike), and the bare-surface sticking fraction of the AR-50 voxel hole, 0.1450 at both
  resolutions, equals Machball's Markov solve at Nz = 400 (0.1457). The weighted site inventory of the voxel
  cylinders is within 0.3 % of the smooth one (`inventory_ratio` in the CSV: 1.0000 to 1.0026 for the holes), and
  `voxel.hole` picks exactly round(πr²) gas voxels per slice so the mouth area that sets the entering flux is the
  smooth one (a plain radius threshold gives digital discs 1–4 % off in area, which shifted t90 by 3–4 %).
- Threads and reproducibility: the M particles are processed in fixed blocks of 64; block b of step k has its own
  xoshiro256** stream seeded from (seed, k, b), and each thread accumulates integer collision counts per face that are
  merged by integer addition. The result is bitwise identical for any thread count and for a given seed
  (`ctest aldmc.determinism_and_thread_invariance`, `pytest test_same_seed_is_deterministic_and_thread_invariant`).
  Blocks are handed to threads dynamically (an atomic counter) so the heavy tail of long-walking particles does not
  idle the other threads (with blocks of 1024 the AR-100 run used three of eight threads). `-ffp-contract=off` as in
  lib/aldkin.

| File | Purpose |
|---|---|
| `lib/aldmc/include/aldmc.h` | the API: opaque `aldmc_t`, `aldmc_params` (s, r, N0, vbar, seed, threads, periodicX/Y, areaCorrection, update, normalRadius, emission, maxSegments), `aldmc_create(nx, ny, nz, h, occupancy, params)`, `aldmc_step(h, n_mouth, dt, M)`, `aldmc_reset`, getters for per-face coverage (`aldmc_face_coverage`/`aldmc_set_face_coverage`), per-face hits and sticks of the last step, face geometry (`aldmc_face_info`: voxel, direction, weight, centre), per-slice mean coverage along z (`aldmc_slice_coverage`), floor coverage (`aldmc_bottom_coverage`), particle counts escaped/stuck/lost/left (`aldmc_counts`, `aldmc_last_counts`), segments traversed, mouth area, site area, time; return codes; Windows export macros (`ALDMC_BUILD_DLL`/`ALDMC_DLL`), `visibility("default")` elsewhere; the physics and weighting stated in the header comment |
| `lib/aldmc/src/aldmc.cpp` | the solver: BFS Chebyshev distance field, face enumeration, cone-kernel normal estimate and face weights, exact voxel traversal with free-flight skips, per-block xoshiro256** streams, `std::thread` workers with per-thread integer hit/stick counters, the two coverage updates, the two emission modes; `try/catch` at every entry, thread-local message for a failed create |
| `lib/aldmc/CMakeLists.txt` | shared (`libaldmc.dylib`/`.so`/`aldmc.dll`, 1.0.0) and static (`libaldmc.a`, `-fPIC`) targets, `Threads::Threads`, `-ffp-contract=off`, hidden visibility, install rules, 10 ctest tests, the libm link rule for glibc |
| `lib/aldmc/tests/test_aldmc.cpp` | the ctest suite (no framework): version and defaults; open box loses every particle; flat surface vs 1 − exp(−s n vbar/4 t/N0) to 0.005 (exponential update) and the forward-Euler value 1 − 0.9¹⁰ (deposit update); periodic slit conserves the particle count; r = 1 loses every particle; same seed bitwise identical at 1 and 4 threads; cylinder inventory within 3 % of smooth and plain faces at 1.24×; AR-5 hole saturates monotonically with coverage falling with depth; error codes, reset, set coverage, zero density/zero dt |
| `lib/aldmc/README.md` | physics, weighting, area correction, emission modes, threads, API table, test list |
| `src/aldconf/voxel.py` | ctypes binding: `load_library` (ALDMC_LIBRARY, ALDMC_BUILD_DIR, `lib/aldmc/build[/Release]`, system loader), `is_available`, `Params`, `VoxelGeometry`, builders `hole(d, L, voxel)`, `tapered_hole(d, L, bottom_ratio, voxel)`, `slit(w, L, voxel)` (periodic in x), `recessed_slit(w, L, g_r, L_r, p_r, voxel)` (the word-line structure of `aldconf.recess`, periodic along the slit, g_r and p_r rounded to whole voxels and recorded in `meta`), `from_array(occupancy, voxel, periodic)`, `open_box`, `flat_surface`; `VoxelMC` (step, reset, face/slice/bottom coverage, faces, counts, rate); `simulate_driven_voxel(geometry, params, t_series, n0_series, M, substeps_per_sample, record_depths, resample_dt, snapshot_times, ...)` returning a `VoxelDrivenResult` that unpacks like `simulate_driven`'s (x, theta, hist_t, hist) and also carries the floor history, snapshots, counts and timings; `simulate_constant_voxel`; `time_to_level`. Default threads min(8, cpu count) or `ALDMC_THREADS` |
| `tests/test_voxel_mc.py` | 12 tests, module-level skip with the build command when the library is absent: AR-20 hole vs the 1D t90 and vs Machball's t95 (skipped without machball), flat Langmuir, open box, periodic slit conservation, recombination, determinism and thread invariance, area correction, geometry builders, error handling, the driven result's shape |
| `coupling/voxel_mc_crosscheck.py` | the cross-checks (parts `holes`, `convergence`, `tapered`, `recess`, `driven`, `slit`, `machball`; results cached in `figures/voxel_mc_crosscheck_runs.json` with the histories; `--figure` assembles the outputs; `--threads`) |
| `figures/voxel_mc_crosscheck.csv`, `.png`, `figures/voxel_mc_profiles.csv`, `figures/voxel_mc_crosscheck_runs.json` | one row per voxel run (24) with its 1D / Machball / nested-model references, particle counts, wall times and rates; the six-panel figure; the chamber-driven depth profiles at 5, 10 and 13 s; the cache |
| `reproduce.sh` | `build_kinlib()` now builds both libraries through `cmake_lib()` (lib/aldkin then lib/aldmc; configure, build, ctest; logs `lib/<name>/log.cmake`); `--build`, `--regress`, `--quick`, `--full` go through it (verified here by running the function; the OpenFOAM part was not rerun) |
| `.github/workflows/ci.yml` | the `library` job builds and ctests lib/aldmc after lib/aldkin and runs `pytest tests/test_aldkin.py tests/test_voxel_mc.py` with `ALDMC_THREADS=2` (installs `.[test,crosscheck]` so the Machball test runs there too); the openfoam job's `reproduce.sh --build` builds both and uploads `lib/aldmc/log.cmake` |
| `.gitignore`, `.dockerignore` | `lib/aldmc/build`, `lib/aldmc/log.cmake` |

Toolchain: CMake 4.4.5 (`~/.local/bin`), Apple clang 21 (no OpenMP; `std::thread`), Python 3.12.14, Machball from the
venv. Machine: Apple M4 (10 cores); 7–8 worker threads. ctest: 10 passed (lib/aldmc; lib/aldkin's 11 still pass).
pytest: 118 passed in 21 s (`PYTHONPATH=src .venv/bin/python -m pytest -q tests`; 12 new, none skipped here).

## Cross-checks

Conditions of `tests/test_machball_crosscheck.py` throughout unless stated: s = 5e-3, N0 = 4e18 m⁻², 573 K, 0.3 Pa TMA
(n = 3.79e19 m⁻³, vbar = 410 m/s), 100 nm holes. Times are to 90 % (t90) and 95 % (t95) coverage of the hole floor
(voxel: area-weighted floor faces; 1D: the deepest cell; Machball: its base element, the last to saturate), linearly
interpolated between steps. Voxel runs: 16 voxels across the diameter (h = 6.25 nm) unless stated, seed 1.

### Machball's own discretisation (found on the way)

The repository's Machball numbers (CLAIMS row "Hole-model vs Machball", draft 2.3, assumptions.md) come from
`Via(AR, Nz = 2 AR)`. Machball's times are not converged at that Nz: t90 / t95 in s

| via | Nz = 2 AR | 4 AR | 8 AR | 16 AR (32 AR for AR 20) |
|---|---|---|---|---|
| AR 20 | 1.123 / 1.289 | | 1.162 / 1.327 | 1.158 / 1.334 |
| AR 50 | 4.247 / 4.465 | | 4.461 / 4.689 | 4.476 / 4.704 |
| AR 100 | 15.279 / 15.590 | 15.968 / 16.277 | 16.151 / 16.460 | 16.191 / 16.511 |
| trench AR 50 | 1.114 / 1.290 | | 1.117 / 1.293 | |

The converged values are 3–6 % longer than the quoted ones (the bare-surface sticking fraction of the AR-50 via is
0.1500 at Nz = 100 and 0.1457 at Nz = 400). The tables below use the converged Machball values (Nz = 8 AR, 16 AR at
AR 100; `refs` in the cache, part `machball`). Against them the 1D Knudsen solver is 11 % fast at AR 20, 5.5 % fast at
AR 50 and 1.6 % fast at AR 100, so the two-method statement in draft 2.3 should be re-worded (proposal below).

### (a) Straight holes AR 20 / 50 / 100 (figure panel 1)

| AR | voxel dt, M per step | voxel t90 / t95 (s) | 1D Knudsen t90 / t95 | Machball converged t90 / t95 | voxel vs Machball | voxel vs 1D | particles, wall, rate |
|---|---|---|---|---|---|---|---|
| 20 | 0.005 s, 25 000 | 1.148 / 1.316 | 1.035 / 1.204 | 1.158 / 1.334 | −0.9 % / −1.3 % | +10.9 % / +9.3 % | 6.9e6, 28 s, 30 800 p/s/thread |
| 50 | 0.0125 s, 50 000 | 4.513 / 4.735 | 4.230 / 4.451 | 4.476 / 4.704 | +0.8 % / +0.7 % | +6.7 % / +6.4 % | 1.9e7, 138 s, 17 400 p/s/thread |
| 100 | 0.025 s, 50 000 | 16.454 / 16.784 | 15.931 / 16.249 | 16.191 / 16.511 | +1.6 % / +1.7 % | +3.3 % / +3.3 % | 3.4e7, 344 s, 14 000 p/s/thread |
| 20 | 0.02 s, 100 000 | 1.153 / 1.320 | | | −0.4 % / −1.0 % | | 6.9e6, 18 s |
| 50 | 0.05 s, 200 000 | 4.558 / 4.788 | | | +1.8 % / +1.8 % | | 1.9e7, 97 s |
| 100 | 0.1 s, 200 000 | 16.789 / 17.127 | | | +3.7 % / +3.7 % | | 3.5e7, 286 s |

The two ballistic codes agree within 2 % at every aspect ratio (within 1 % after the dt → 0 extrapolation of table b,
16.35 s at AR 100); the target of 5 % is met with the converged Machball values and missed at AR 100 only against
the unconverged Nz = 2 AR value (16.45 vs 15.28 s, +7.7 %). The continuum 1D model is the outlier at low AR, where its
Dirichlet mouth condition and the diffusion approximation are least justified; the three methods converge at AR 100
(15.93 / 16.19 / 16.45 s).

### (b) Convergence (AR 100, constant 0.3 Pa; figure panel 2)

| variation | runs (t90 in s) | reading |
|---|---|---|
| particles per step M (16 across, dt 0.1 s) | 25 000: 16.851; 50 000: 16.740; 100 000: 16.789; 200 000: 16.789; 400 000: 16.780 | flat from 25 000 on; scatter ±0.3 % |
| seed (16 across, dt 0.1 s, M 200 000) | seed 1: 16.789; seed 2: 16.751 | 0.2 % |
| time step (16 across, equal total particles) | 0.025 s: 16.454; 0.05 s: 16.556; 0.1 s: 16.789; 0.2 s: 17.129 | linear in dt (3.9 s per s of dt): the frozen-coverage bias; dt → 0 gives 16.35 s; dt = t90/200 is within 0.6 % |
| voxel resolution (dt 0.1 s, M 200 000) | 16 across: 16.789; 24: 16.719; 32: 16.612 | 1.1 % over a doubling of the resolution; site inventory 1.0016 / 1.0004 / 1.0000 of smooth |
| update scheme (16 across, dt 0.1 s, M 200 000) | impingement (default): 16.789; deposit: 17.034 | the deposit update's own first-order bias, +1.5 % at this dt |
| area correction off (plain voxel faces) | 21.138 | +25.9 %: the 4/π staircase excess (inventory 1.273), the reason for the correction |

### (c) Tapered hole, bottom 60 % of the 100 nm mouth (constant 0.3 Pa; figure panel 3)

Voxel (16 across at the mouth, dt 0.1 s, M 200 000, 4.1e7 particles, 342 s): floor t90 20.02 s, t95 20.32 s. 1D frustum
model (`TaperedHole(100 nm, 10 µm, 0.6)`, `dose_for_coverage`, figures/tapered_hole.csv): 19.12 / 19.42 s. Voxel
+4.7 %, the same offset as the straight hole at the same dt (+5.4 %); the taper's own effect agrees to 0.7 %:
tapered / straight = 1.192 (voxel) against 1.200 (1D). Weighted inventory 1.0026 of the frustum's lateral area plus
its floor.

### (d) Word-line slit w100_Lr365 (constant 0.3 Pa; figure panel 4)

Geometry `recessed_slit(100 nm, 5 µm, 29 nm, 365 nm, 51 nm, h = 5 nm)`: g_r rounds to 30 nm and p_r to 50 nm (L_r 365 nm
and w 100 nm are exact), 98 recesses per wall, 4 × 168 × 1001 voxels, 124 880 faces, periodic along the slit; dt 0.1 s,
M 100 000, 1.3e7 particles, 299 s (6 000 p/s/thread: 5 800 voxel segments per particle, the recesses are expensive).

| quantity | voxel | nested model (`simulate_recessed_constant`), nominal 29/51 nm | nested, rounded 30/50 nm |
|---|---|---|---|
| slit floor / deepest slit-wall slice t90 | 10.70 / 10.70 s | 14.53 s (slit bottom) | 14.78 s |
| deepest recess end wall t90 | 10.88 s | 14.68 s | 14.92 s |
| recess at the mouth / mid-depth t90 | 0.71 / 7.39 s | 0.56 / 6.32 s | 0.56 / 6.39 s |
| recess lag behind the slit bottom | 0.18 s | 0.15 s | 0.14 s |

The voxel solver saturates the structure 27 % sooner than the nested model, with the same ordering (slit bottom, then
the deepest recess, lag 0.15–0.18 s) and the same early recesses. The cause is the slit transport, not the recesses:
part (f) below shows the 1D slit model 24 % slow against both ballistic codes on a plain AR-50 slit, and the nested
model is that slit model with recess sub-problems of the same kind. (The dotted "slowest recess end" curve in the
figure is the pointwise minimum over 98 noisy recess histories, 12.06 s; it is biased late by the noise and the deepest
recess is the honest number.) Inventory 0.984 of the smooth structure (the recess mouths' corner blending).

### (e) Chamber-driven AR-100 hole (figure panels 5 and 6)

Drive: the long case's innermost-face density (`lib/aldkin/examples/centre_face_drive.csv`, 0.05 s samples, n0 =
T_face × 3.79e20 m⁻³, 20 s pulse), the same series `docs/build_assets.py` feeds to `simulate_driven`. 1D
(`simulate_driven`, 5 substeps per sample): bottom t90 14.025 s (14.05 s is the boundary condition's own figure;
14.007 s in the recorded thetaBottom column of the CSV). Voxel runs stop at 16 s.

| voxel run | t90 / t95 (s) | vs 1D 14.025 s | particles, wall |
|---|---|---|---|
| 16 across, history subsampled to 0.25 s, M 400 000 | 14.553 / 14.657 | +3.8 % | 2.6e7, 198 s |
| 16 across, native 0.05 s samples, M 80 000 | 14.319 / 14.385 | +2.1 % | 2.6e7, 211 s |
| 24 across, 0.25 s, M 400 000 | 14.540 / 14.650 | +3.7 % | 2.6e7, 252 s |

The 0.25 s subsampling (as planned) costs 1.7 % through the frozen-coverage bias; the 0.05 s run is within 2 % of
the 1D solver and within the 1D-vs-ballistic offset seen at constant density (+3.3 % at AR 100 in table a). Coverage
along the depth at 5, 10 and 13 s (`figures/voxel_mc_profiles.csv`, 41 relative depths, recomputed from aldconf, not
from the docs data): the voxel front sits slightly ahead of the 1D front early and behind it late (bottom coverage at
13 s: 0.28 voxel vs 0.46 1D; mid-depth at 10 s: 0.97 vs 0.99); the largest pointwise difference is 0.10 / 0.18 / 0.23
at 5 / 10 / 13 s with the 0.25 s drive and 0.06 / 0.11 / 0.14 with the 0.05 s drive, rms 0.04–0.10. The two solvers
agree on the saturated and empty parts of the hole and differ in the width of the front, as expected between a
ballistic and a diffusive description.

### (f) Slit AR 50 (added: constant 0.3 Pa, gap 100 nm, periodic, no staircase)

Voxel (16 across, dt 0.0125 s, M 25 000; and dt 0.05 s, M 100 000; 11 s each): t90 1.083 / 1.109 s, t95 1.251 / 1.283 s.
Machball `Trench(50, 400)`: 1.117 / 1.293 s. 1D slit model (`Slit`, D_K = 2 w vbar/3): 1.381 / 1.556 s. The two
ballistic codes agree to 1–3 %; the 1D slit model is 24 % slow. This is the known logarithmic correction of the
parallel-plate Knudsen diffusivity that `transport.py` mentions (`knudsen_factor`, 1 by default): the cylinder's
d vbar/3 is exact, the slit's 2 w vbar/3 is not. It affects the slit-based results of the repository (the word-line
recess times in figures/wordline_recess.csv are long by roughly this factor; the PillarHall fit absorbed it into the
fitted s). The hole results, which carry the chamber conclusions, are unaffected.

### Timings

24 voxel runs, 6.1e8 particles, 93 min of wall time in all on 7–8 threads; no single run above 10 min (the longest:
32 across, 514 s; M 400 000, 585 s). Rates per thread: 31 000–49 000 particles/s for AR 20 (720 voxel segments per
particle), 17 000–25 000 for AR 50 (1 150), 14 000–17 000 for AR 100 at 16 across (1 950), 12 100 at 24 across (2 290),
9 500 at 32 across (2 510), 35 000 for the slit (1 020), 6 000 for the recessed slit (5 800). Early in a run (empty
hole) the rate is 80 000–90 000 particles/s per thread; the long walkers of the saturated hole bring the average down.
The 1D solver's equivalent answers take under a second each; Machball's converged runs 1–131 s.

## Limits

- No film narrowing: the opening does not shrink with the deposited film (the 1D solver's `narrowing` has no voxel
  counterpart). For the single-half-cycle saturation times compared here it does not matter (the 1D runs were made
  with `narrowing=False` too); multi-cycle growth would need the occupancy to be rebuilt from the thickness.
- No gas-phase collisions by construction (molecular flow; Kn ≫ 1 inside the features at 300 Pa, mean free path
  about 20 µm against 100 nm openings). The transition-regime (Bosanquet) option of the 1D solver has no counterpart.
- Statistical noise: the floor t90 scatters by 0.2–0.4 % between seeds and particle budgets at the budgets used;
  per-face coverages are noisier (a 6 nm face holds about 150 sites and sees about 30 collisions per 0.1 s step at
  0.3 Pa), and minima over many faces or recesses are biased. The impingement-based update makes the noise
  independent of s.
- Discretisation: the frozen coverage within a step leaves a first-order time-step bias (3.9 s per s of dt at
  AR 100; 0.6 % at dt = t90/200); the staircase correction leaves 1 % over a doubling of the resolution; the
  geometry builders round slit and recess dimensions to whole voxels (recorded in `meta`).
- Cost: 10⁴ voxel segments per particle and 10⁷–10⁸ particles per saturation curve: minutes to ten minutes per
  feature on 8 threads against sub-second 1D solves. The voxel solver is the check and the tool for shapes the 1D
  model cannot describe, not a replacement for the 1D model inside the chamber coupling.

## Proposed paper text (paper/draft.md, new subsection after 2.5)

```
### 2.6 A 3D voxel Monte Carlo feature solver

To check the one-dimensional Knudsen-diffusion model against a description that makes none of its assumptions, the
same features were solved with a test-particle Monte Carlo on a voxel grid (lib/aldmc). In the molecular-flow regime
that holds inside the features (Kn ~ 200), particles move in straight lines between wall collisions and never
collide with each other. They enter through the feature mouth with a cosine angular distribution and a uniform
position over the opening; on hitting a solid voxel face they stick with probability s (1 - theta) on that face's
own coverage, are lost with probability r, and otherwise re-emit diffusely; particles leaving through the mouth are
lost. A time step launches M particles that represent the mouth flux n vbar/4 times the mouth area and the step, and
the coverage of every face is advanced by the Langmuir relaxation of the collisions it received, the same update as
the 1D solver with the Monte Carlo impingement in place of n vbar/4 dt. Because a voxel surface is a staircase with
4/pi times the area of a smooth cylinder wall, each face carries sites and sticking in proportion to its projected
area on the estimated smooth surface; without this correction the hole saturates 26 % late. The solver is
deterministic for a given seed at any thread count, needs no diffusion coefficient or hydraulic diameter, and takes
arbitrary shapes from an occupancy array (straight and tapered holes, slits, the word-line recess structure).

For the straight 100 nm holes at 0.3 Pa TMA and 573 K the voxel solver and Machball agree within 2 % at AR 20, 50
and 100 once Machball's wall discretisation is converged (t90 1.15 / 4.51 / 16.45 s against 1.16 / 4.48 / 16.19 s;
Machball's values move by 3-6 % between Nz = 2 AR and 8 AR), while the 1D Knudsen solver is 11 % fast at AR 20, 7 %
at AR 50 and 3 % at AR 100 (1.04 / 4.23 / 15.93 s): the continuum description converges to the ballistic one as the
aspect ratio grows. The particle budget (25 000 to 400 000 per step), the voxel size (16 to 32 across the diameter)
and the seed change t90 by at most 1 %; the time step enters linearly through the coverage frozen within a step
(0.6 % at dt = t90/200). The frustum with a 60 % bottom lengthens t90 by the factor 1.19 in the voxel solver and 1.20
in the 1D model. Driven by the chamber's centre-face density history of the long AR-100 case, the voxel hole reaches
90 % bottom coverage at 14.32 s (14.55 s with the history subsampled to 0.25 s) against 14.03 s in the 1D solver and
14.05 s in the boundary condition, with the same saturated and empty regions and a somewhat sharper front at 5, 10
and 13 s. For slits the two ballistic codes again agree (AR 50: 1.08-1.11 s voxel, 1.12 s Machball), but the 1D slit
model with D_K = 2 w vbar/3 is 24 % slow, the known logarithmic correction of the parallel-plate Knudsen
diffusivity; the word-line recess structure therefore saturates 27 % sooner in the voxel solver (10.7 s against
14.5 s at 0.3 Pa) with the same ordering of slit bottom and recess ends. The chamber results rest on the hole model,
whose Knudsen diffusivity d vbar/3 is exact; the slit-based numbers of Section 3.11 are conservative by this factor.
```

Draft 2.3's Machball bullet should become: "Two-code check: Machball (ballistic, cosine re-emission) vs this
Knudsen-continuum solver, 95 % bottom coverage time at 0.3 Pa TMA, 573 K: 1.20 vs 1.33 s (AR 20), 4.45 vs 4.70 s
(AR 50), 16.2 vs 16.5 s (AR 100) with Machball's wall discretisation converged (Nz = 8-16 AR; its Nz = 2 AR values are
1.29 / 4.46 / 15.6 s); Section 2.6 adds a third, voxel Monte Carlo method that agrees with Machball within 2 %."

## README text

Row for the "What is in the box" table (after the "Feature solver" row):

```
| 3D feature solver | `lib/aldmc`, `src/aldconf/voxel.py` | C++17 voxel test-particle Monte Carlo (molecular flow, cosine entry and re-emission, Langmuir sticking on every voxel face, std::thread, C ABI) for arbitrary shapes from an occupancy array; agrees with Machball within 2 % and checks the 1D solver (`coupling/voxel_mc_crosscheck.py`, `figures/voxel_mc_crosscheck.png`) |
```

Paragraph for the Reproduce section (after the kinetics-library lines):

```
# 3D voxel Monte Carlo feature solver (C++17, std::thread; the Python binding finds lib/aldmc/build)
cmake -S lib/aldmc -B lib/aldmc/build -DCMAKE_BUILD_TYPE=Release && cmake --build lib/aldmc/build
ctest --test-dir lib/aldmc/build                                       # 10 checks
PYTHONPATH=src .venv/bin/python -m pytest -q tests/test_voxel_mc.py    # 12 tests, skipped until the library is built
PYTHONPATH=src .venv/bin/python coupling/voxel_mc_crosscheck.py --part holes,convergence,tapered,recess,driven,slit,machball
                                                                       # the cross-checks, about 95 min on 8 threads
```

And in "Continuous integration and Docker": "the `library` job builds and tests both C-ABI libraries, lib/aldkin
(kinetics) and lib/aldmc (voxel Monte Carlo), and runs their Python binding tests".

## CLAIMS rows (CLAIMS.md format)

```
| Voxel Monte Carlo vs Machball (converged Nz) vs 1D Knudsen, straight 100 nm holes at 0.3 Pa TMA, 573 K, s 5e-3: t90 | AR 20: 1.15 / 1.16 / 1.04 s; AR 50: 4.51 / 4.48 / 4.23 s; AR 100: 16.45 / 16.19 / 15.93 s (voxel 16 across, dt t90/200; voxel vs Machball within 2 %, 1D 11 / 7 / 3 % fast) | figures/voxel_mc_crosscheck.csv, figures/voxel_mc_crosscheck.png | PYTHONPATH=src .venv/bin/python coupling/voxel_mc_crosscheck.py --part holes,machball |
| Machball's wall discretisation | t95 1.29 / 4.46 / 15.6 s at Nz = 2 AR become 1.33 / 4.70 / 16.5 s at Nz = 8-16 AR (3-6 % longer); the bare-surface sticking fraction of the AR-50 via 0.1500 (Nz 100) vs 0.1457 (Nz 400), the voxel hole 0.1450 | figures/voxel_mc_crosscheck_runs.json (refs.mb_nz) | same, --part machball |
| Voxel solver convergence, AR 100 at 0.3 Pa | M 25 000-400 000 per step and seed: t90 within 0.4 % (16.74-16.85 s at dt 0.1 s); resolution 16 / 24 / 32 across: 16.79 / 16.72 / 16.61 s; time step 0.025 / 0.05 / 0.1 / 0.2 s: 16.45 / 16.56 / 16.79 / 17.13 s (linear, dt -> 0: 16.35 s); deposit update 17.03 s; plain voxel faces (no area correction) 21.14 s (+26 %, inventory 1.273) | figures/voxel_mc_crosscheck.csv | same, --part convergence |
| Voxel tapered AR-100 hole, bottom 60 %, constant 0.3 Pa | t90 20.02 s vs 1D frustum 19.12 s; tapered/straight 1.192 (voxel) vs 1.200 (1D) | figures/voxel_mc_crosscheck.csv, figures/tapered_hole.csv | same, --part tapered |
| Voxel word-line slit w100_Lr365 at constant 0.3 Pa (h 5 nm: g_r 30, p_r 50 nm) | slit floor t90 10.70 s, deepest recess end 10.88 s (lag 0.18 s) vs nested model 14.53 / 14.68 s (lag 0.15 s); the 27 % is the 1D slit transport | figures/voxel_mc_crosscheck.csv, figures/wordline_recess.csv | same, --part recess |
| Slit AR 50 at 0.3 Pa: voxel vs Machball Trench vs 1D slit model | t90 1.08-1.11 / 1.12 / 1.38 s: the ballistic codes agree, the 1D slit (D_K = 2 w vbar/3) is 24 % slow | figures/voxel_mc_crosscheck.csv | same, --part slit |
| Voxel AR-100 hole driven by the long case's centre-face density | bottom t90 14.32 s (native 0.05 s history) / 14.55 s (subsampled to 0.25 s) / 14.54 s (24 across) vs 14.03 s simulate_driven and 14.05 s in the boundary condition; depth profiles at 5 / 10 / 13 s within 0.06 / 0.11 / 0.14 pointwise (0.05 s drive) | figures/voxel_mc_crosscheck.csv, figures/voxel_mc_profiles.csv | same, --part driven |
| Voxel solver verification | flat surface follows 1 - exp(-s n vbar/4 t/N0) to 5e-3; open box loses every particle; periodic slit conserves the count; same seed bitwise identical at 1 and 4 threads; Clausing transmission of open voxel tubes 0.108 / 0.059 vs 0.1093 / 0.0595 (L/d 10 / 20) | lib/aldmc/tests/test_aldmc.cpp, tests/test_voxel_mc.py, paper/wp/wp9_voxel_mc.md | ctest --test-dir lib/aldmc/build; pytest tests/test_voxel_mc.py |
```

## Run log for paper/assumptions.md

```
## Night 4, part 2: 3D voxel Monte Carlo feature solver (lib/aldmc; src/aldconf/voxel.py; coupling/voxel_mc_crosscheck.py; figures/voxel_mc_crosscheck.png, .csv)
A second feature solver on a voxel grid: test-particle Monte Carlo in the molecular-flow regime (straight flights,
no gas-phase collisions; cosine entry over the open mouth; sticking s (1 - theta) on the face's own coverage, loss
r, diffuse re-emission about the facet normal; escape through the mouth), C++17 with a C ABI (include/aldmc.h:
aldmc_create from nx, ny, nz, voxel size, uint8 occupancy, params s, r, N0, vbar, seed, threads, periodicX/Y,
areaCorrection, update, emission; aldmc_step(h, n_mouth, dt, M); per-face, per-slice and floor coverage; counts
escaped / stuck / lost / left), std::thread (Apple clang has no OpenMP), per-block xoshiro256** streams and integer
hit counts merged per thread, so a seed gives bitwise the same result at any thread count. A step of M particles
carries n vbar/4 A_mouth dt molecules; the coverage update is the Langmuir relaxation of the collisions counted per
face (the 1D solver's update with the Monte Carlo impingement), which uses every collision and so has 1/s times the
variance of counting sticking events; the literal deposit-per-stick update is kept as an option (+1.5 % at dt 0.1 s).
Faces carry N0 h^2 |n . e| sites with n the smooth normal estimated from the occupancy (cone kernel, 5 voxels), and
stick and lose in the same proportion, which removes the 4/pi staircase excess (plain faces: +26 % on t90). The
digital disc of a hole takes exactly round(pi r^2) voxels per slice (a radius threshold gave mouth areas 1-4 % off
and t90 3-4 % off). Open voxel tubes transmit 0.108 / 0.059 at L/d 10 / 20 (Berman 0.1093 / 0.0595); the bare
AR-50 hole sticks 0.1450 of the entering particles at 16 and 32 across (Machball 0.1457 at Nz 400). Re-emission about
the smooth normal was tried and rejected (10 % under-counted uptake). No film narrowing.
Cross-checks (s 5e-3, N0 4e18, 573 K, 0.3 Pa, 100 nm; 16 voxels across; t90 of the floor): straight holes AR 20 / 50 /
100: voxel 1.148 / 4.513 / 16.454 s (dt t90/200), 1D 1.035 / 4.230 / 15.931 s, Machball 1.158 / 4.476 / 16.191 s. The
last are Machball's converged values (Nz = 8-16 AR); at the Nz = 2 AR of tests/test_machball_crosscheck.py Machball
gives 1.123 / 4.247 / 15.279 s (t95 1.289 / 4.465 / 15.590 s, the numbers quoted so far), 3-6 % short. Voxel vs
Machball within 2 % everywhere (within 1 % at dt -> 0); the 1D continuum model 11 / 7 / 3 % fast, converging to the
ballistic answer with AR. Convergence at AR 100: M 25 000-400 000 and seed within 0.4 %; 16 / 24 / 32 across 16.79 /
16.72 / 16.61 s; dt 0.025-0.2 s linear (16.45-17.13 s, 3.9 s per s of dt). Frustum 60 %: 20.02 vs 19.12 s 1D
(tapered/straight 1.192 vs 1.200). Word-line slit w100_Lr365 at 0.3 Pa (h 5 nm, g_r 30 / p_r 50 nm): slit floor
10.70 s, deepest recess end 10.88 s vs the nested model's 14.53 / 14.68 s; the 27 % is the 1D slit transport: on a plain
AR-50 slit voxel 1.08-1.11 s and Machball Trench 1.12 s agree while the 1D slit model (D_K = 2 w vbar/3) gives 1.38 s
(the logarithmic parallel-plate correction transport.py mentions; the hole's d vbar/3 is exact). Chamber-driven AR-100
hole on the long case's centre-face history: 14.32 s with the native 0.05 s samples, 14.55 s subsampled to 0.25 s,
14.54 s at 24 across, vs 14.03 s simulate_driven / 14.05 s boundary condition; depth profiles at 5 / 10 / 13 s agree on
the saturated and empty parts and differ in the front width (pointwise up to 0.06 / 0.11 / 0.14 with the 0.05 s drive).
Cost: 24 runs, 6.1e8 particles, 93 min on 7-8 threads of an M4 (6 000-49 000 particles/s per thread; 700-5 800
voxel segments per particle), no run above 10 min. ctest 10 passed; pytest 118 passed in 21 s (12 new, skipped with
the build command when the library is absent). CI: the library job builds and ctests lib/aldmc and runs its binding
tests; reproduce.sh --build cmakes both libraries. Caveats: no narrowing, no gas-phase collisions by construction,
first-order time-step bias (use dt <= t90/100), slit and recess dimensions rounded to whole voxels, statistical
noise on single faces and on minima over many faces.
```

## What the integrator must check after the first CI run

1. `library` job: `cmake -S lib/aldmc` with GCC on ubuntu-latest (`Threads::Threads`, `-ffp-contract=off`), the
   10 ctests (the determinism test runs 4 threads on a 2-core runner, fine), and `pytest tests/test_aldkin.py
   tests/test_voxel_mc.py -rs` with `ALDMC_THREADS=2`: the AR-20 test takes about 10 s on 4 threads here, expect
   20–40 s there; the Machball test needs `.[test,crosscheck]`, which the step now installs.
2. `openfoam` job and the Dockerfile: `reproduce.sh --build` now cmakes lib/aldmc after lib/aldkin (same
   `cmake g++` apt line; the image's GCC needs pthreads, linked through `Threads::Threads`); `lib/aldmc/log.cmake` is
   uploaded with the logs.
3. `python` job: 12 more skips (`libaldmc not built: ...`) in the `-rs` summary.
4. The AR-20 pytest against the 1D t90 is at 15 %, not the 10 % first intended: a ballistic solver sits 10–11 %
   above the continuum one at AR 20 (Machball's converged t90 is 12 % above it), so 10 % is not a property of the
   physics; the tight check is the 5 % against Machball's t95, which passes at 1.3 %.
5. The existing claims that quote Machball at Nz = 2 AR (CLAIMS row "Hole-model vs Machball", draft 2.3,
   assumptions.md "Cross-check of the hole model against Machball") should be updated with the converged values
   above; `tests/test_machball_crosscheck.py` still passes (its tolerances 12 % / 8 % cover both Nz choices) but a
   `(100, 200, 0.08)` case would fail at Nz = 2 AR and pass at Nz = 8 AR (16.25 vs 16.46 s).
6. The 1D slit model's 24 % transport deficit (part f) touches the word-line recess numbers (Section 3.11,
   figures/wordline_recess.csv) and the PillarHall fit's absolute s; the hole-based chamber results are unaffected.
   `knudsen_factor` in `Process`/`simulate_recessed` is the existing hook if a corrected slit diffusivity is wanted.
7. README.md as it stands has an empty code block under "Use the kinetics in your own solver" and a stray line
   "Also update the test count in the existing Reproduce block ..." at the end of the CI section, left over from
   the WP8 merge; the test count there still says 77 (now 118 with both libraries built).
