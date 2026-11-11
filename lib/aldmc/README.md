# aldmc: 3D voxel Monte Carlo feature solver (C ABI)

A feature-scale ALD solver on a voxel grid: test-particle Monte Carlo in the molecular-flow (Knudsen) regime,
C++17 with a C ABI, parallel with `std::thread`. The Python binding with geometry builders is `src/aldconf/voxel.py`;
the cross-checks against the repository's 1D solver and Machball are in `coupling/voxel_mc_crosscheck.py` and
`paper/wp/wp9_voxel_mc.md`.

```
cmake -S lib/aldmc -B lib/aldmc/build -DCMAKE_BUILD_TYPE=Release
cmake --build lib/aldmc/build
ctest --test-dir lib/aldmc/build --output-on-failure        # 10 checks
```

Targets: `aldmc` (shared, `libaldmc.dylib` / `.so` / `aldmc.dll`, versioned 1.0.0) and `aldmc_static` (`libaldmc.a`,
position independent). No OpenMP: Apple clang builds it as is; on Linux `Threads::Threads` is linked.

## Physics

Molecular flow: particles move in straight lines between wall collisions and never collide with each other. They
enter through the feature mouth (the open plane z = 0) with a cosine (Lambertian) angular distribution and a uniform
position over the open mouth voxels. On hitting a solid voxel face they stick with probability s (1 − θ_face)
(irreversible Langmuir adsorption on the face's own site coverage), are lost with probability r (radical
recombination, default 0), and otherwise re-emit diffusely (cosine law about the face normal). Particles leaving
through the mouth are lost to the gas phase ("escaped"). Coverage lives on the exposed voxel faces: a face carries
N0 h² w_face sites (h the voxel size, w_face the projected-area weight below). Film growth does not narrow the
feature (not implemented; the 1D solver's `narrowing` has no counterpart here).

## Time stepping and particle weight

The mouth flux is n(t) vbar/4 per unit mouth area. One `aldmc_step(h, n_mouth, dt, M)` launches M particles that
represent n vbar/4 · A_mouth · dt molecules, each carrying w = n vbar/4 A_mouth dt / M molecules. The coverage field
is frozen while the M particles fly (sticking is sampled with the coverage at the start of the step) and the wall
collisions of the step are counted per face. The default update (`ALDMC_UPDATE_EXPONENTIAL`) applies the exact
Langmuir relaxation of the impinged dose,

    θ_f ← 1 − (1 − θ_f) exp(− s · hits_f · w / (N0 h²)),

the 1D solver's update `theta <- 1 - (1 - theta) exp(-s vbar/4 n dt / N0)` with the Monte Carlo impingement in place
of n vbar/4 dt. It uses every wall collision, so its variance is about 1/s times smaller than counting sticking
events. `ALDMC_UPDATE_DEPOSIT` is the literal alternative (each sampled sticking event deposits w molecules on its
face, θ capped at 1). Both are first order in dt and agree as M and 1/dt grow (`coupling/voxel_mc_crosscheck.py`,
part `convergence`). Mass is bookkept per particle: every particle ends as escaped, stuck, lost or left.

## Area correction (staircase)

A voxel surface is a staircase: the voxel faces of a cylinder have 4/π times the area of the smooth wall, and a
plain voxel solver would carry 27 % too many sites and saturate 27 % too late. With `areaCorrection = 1` (default) a
smooth surface normal n is estimated for every exposed face from the first moment of the occupancy around the face
centre (cone kernel 1 − r/R, R = `normalRadius`, default 5 voxels), and the face is given the weight
w_face = |n · e_face|, its projected area fraction: it carries N0 h² w_face sites, and sticking and loss per
collision are s w_face (1 − θ) and r w_face. In an equilibrium gas every face is hit at n vbar/4 h² whatever its
orientation, so the uptake per unit smooth area and the rate per site are the smooth surface's, and the staircase's
extra collisions (diffuse re-emission about the actual face normal keeps detailed balance) leave the transport that
of the smooth wall. For cylinders of 12 to 32 voxels across the weighted site inventory is within 2 % of the smooth
one (`aldmc_site_area`); `areaCorrection = 0` gives plain voxel faces (w = 1) for comparison. Re-emission stays
diffuse about the actual facet normal (`emission = ALDMC_EMISSION_FACET`, default): that is what keeps the facet hit
rates at their equilibrium value n vbar/4 h², which the weights rely on. Re-emitting about the estimated smooth
normal instead (`ALDMC_EMISSION_SMOOTH`, kept as an option) removes the staircase's grazing re-collisions, under-counts
the uptake by about 10 % at θ = 0 and makes an AR-20 hole saturate 11 % later; the Clausing transmission of an open
voxel tube is within 1 % of Berman's values with either mode (L/d 10: 0.108 vs 0.1093; L/d 20: 0.059 vs 0.0595), so the
axial transport of the staircase is that of the smooth tube already. The digital cross-section of `voxel.hole` is
chosen with exactly round(πr²) gas voxels per slice, so the mouth area that sets the entering flux equals the smooth
one (a radius threshold alone gives discs 1–4 % off in area).

## Threads and reproducibility

The M particles of a step are processed in fixed blocks of 64; block b of step k uses its own xoshiro256** stream
seeded from (seed, k, b), and each thread accumulates integer collision counts per face which are merged by integer
addition. The result is therefore bitwise identical for any thread count and reproducible for a given seed
(`ctest aldmc.determinism_and_thread_invariance`). `threads = 0` uses the hardware concurrency; the Python binding
defaults to min(8, cpu count) or `ALDMC_THREADS`. Flights use a Chebyshev distance field (a particle in a voxel at
distance D from the nearest solid flies D − 1 voxel lengths without checks) and an exact voxel traversal for the
last voxels.

## API

| Call | Purpose |
|---|---|
| `aldmc_version`, `aldmc_params_default` | version string; defaults (s 1, r 0, N0 4e18, vbar 410, seed 1, threads 0, no periodicity, area correction on, normal radius 5, facet emission, exponential update) |
| `aldmc_create(nx, ny, nz, h, occupancy, params)` | solver on an nx·ny·nz grid of voxel size h; `occupancy[(k ny + j) nx + i]` 1 solid, 0 gas; mouth plane z = 0; `periodicX/Y` wrap the lateral directions. Computes the exposed faces, the mouth voxels, the face weights and the distance field |
| `aldmc_step(h, n_mouth, dt, M)` | M particles for mouth density n_mouth (1/m³) over dt; advances the coverage. n_mouth = 0 or M = 0 only advances time |
| `aldmc_reset` | coverage, counters and time to 0 (the random streams restart) |
| `aldmc_face_coverage`, `aldmc_set_face_coverage`, `aldmc_face_hits`, `aldmc_face_info` | per-face coverage (get/set), collisions and sticking events of the last step, gas voxel, direction (0 +x, 1 −x, 2 +y, 3 −y, 4 +z, 5 −z: from the gas voxel to the solid), weight and centre coordinates |
| `aldmc_slice_coverage` | weighted mean coverage of the faces in every z slice (and the weighted face area per slice) |
| `aldmc_bottom_coverage` | weighted mean coverage of the floor faces (direction +z) in the deepest slice that has any |
| `aldmc_counts`, `aldmc_last_counts`, `aldmc_last_segments` | particles escaped (mouth), stuck, lost (recombination), left (other boundary or the segment cap), cumulative and per step; voxel segments traversed in the last step |
| `aldmc_n_faces`, `aldmc_n_mouth_voxels`, `aldmc_mouth_area`, `aldmc_site_area`, `aldmc_time`, `aldmc_get_params` | sizes, open mouth area, weighted face area (sites / N0), accumulated time, parameters in use |
| `aldmc_last_error` | message of the last failed call on the handle (NULL: of the last failed create in this thread) |

Errors are return codes (`ALDMC_OK` 0, `ALDMC_ERR_ARG`, `ALDMC_ERR_PARAM`, `ALDMC_ERR_STATE`, `ALDMC_ERR_ALLOC`,
`ALDMC_ERR_INTERNAL`); no exception crosses the ABI. Windows export macros `ALDMC_BUILD_DLL` / `ALDMC_DLL`.

## Tests (ctest)

`version_and_defaults`; `open_box_loses_all` (no walls: every particle leaves, no faces); `flat_langmuir_analytic`
(floor under the mouth: mean coverage 1 − e⁻¹ after ten steps of 0.1/rate to 0.005); `flat_langmuir_deposit_scheme`
(the deposit update gives the forward-Euler value 1 − 0.9¹⁰ of the same rate); `periodic_slit_conserves` (escaped +
stuck + lost + left = launched, nothing leaves laterally); `recombination_loses_all` (r = 1, s = 0); `determinism_and_
thread_invariance` (same seed: bitwise identical with 1 and 4 threads, different seed differs); `area_weights_cylinder`
(weighted inventory within 3 % of the smooth cylinder, plain faces 1.24 ×); `hole_saturates` (AR-5 hole: floor
coverage monotone to above 0.9, coverage falls with depth mid-way); `error_codes`.
