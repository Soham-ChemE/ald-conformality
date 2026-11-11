# ald-conformality: complete project record

What was built, in which software, with which settings, and what came out. Written 2026-11-11 at commit 19f47ec
(paper draft 0.9.5). Every number below is taken from the repository's own run log (`paper/assumptions.md`), the
claims audit (`CLAIMS.md`) and the figure CSVs; nothing is from memory. Where a number is an assumption or a
literature range rather than a result, the text says so.

Repository: https://github.com/Soham-ChemE/ald-conformality · Live page: https://soham-cheme.github.io/ald-conformality/
· Paper: `paper/main.pdf` (42 pages, 25 figures, 17 tables) · Author: Soham Kavathekar, University of Pennsylvania.

---

## 1. The question and the answer

ALD recipes are tuned on blanket (flat) wafers. A 3D-NAND wafer carries tens of times the surface area of the wafer
under it, all of which has to saturate in one precursor pulse. The project built a 300 mm showerhead ALD chamber in
OpenFOAM in which every face of the wafer runs its own feature-scale model of a high-aspect-ratio hole, two-way
coupled to the chamber gas, and asked how long the pulse must be.

**Result.** The pulse that saturates every hole bottom on the wafer is

    t_pulse ≈ 1.2 × (1 + A_wall/A_wafer) × N0 / (U·c0)

1.2 times the time the showerhead needs to deliver one monolayer of precursor to the whole patterned area
(N0 sites per m², U·c0 precursor molecules per m² per second at the showerhead). It held within 25 % across aspect
ratios 20 to 100 and a 20× range of precursor feed (12 patterned chamber runs: 1.05 to 1.27 × the supply limit).
Blanket wafer: 0.345 s. Reference 3D-NAND-like wafer (100 nm holes, AR 100, 1.78e13 /m², wall area 55.9 × wafer
area): 14.05 s. A pulse sized on a blanket wafer is 40.7× too short for that wafer. Pressure and gap do not enter the
rule; carrier flow and precursor fraction enter only through U·c0.

**Checks against measurements nobody here controlled.** Argonne's fluidized-bed Al₂O₃ run (940 m² of silica gel)
saturates at 1.12 × its supply limit by the authors' own numbers, inside the model's 1.05 to 1.27 band; Argonne's
rotary-drum utilization data order with loading the same way; the 1D hole model fits Aalto's public PillarHall
thickness profiles (Zenodo 3974438); an axisymmetric copy of the Fiji F200 reactor of Gakis et al. reproduces their
CFD's delivered exposure to 12 to 22 % and the measured saturation state. The 3D voxel Monte Carlo solver written here
agrees with Argonne's Machball ballistic code within 2 %; the 1D continuum hole model is 3 to 11 % fast against both.

---

## 2. Timeline

| When | What |
|---|---|
| 2026-11-09 evening (night 1) | Chamber model, Danckwerts inlet, finite sticking and Langmuir wafer, operating-point sweep, one-way then two-way chamber-to-feature coupling, H₂O half-cycle, 16-run pulse-time map, the rule, plasma-ALD recombination ceiling, temporal-vs-spatial comparison, PillarHall fit, Machball cross-check. 57 OpenFOAM cases. |
| 2026-11-10 (night 2) | Validation against the fluidized bed and the rotary drum, edge-ring geometry, reversible water adsorption, temperature and hole-size sweeps, 25-run central composite DOE with JMP response surface and MATLAB Gaussian-process surrogate, Monte Carlo uncertainty. Repository made public, live page built, profile repo updated, LaTeX paper started. 86 cases. |
| 2026-11-10 night (night 3) | Reviewer list B: spatially varying layouts (ARProfile/holeDensityProfile), sticking sensitivity (9 runs), fab units (step coverage, purge, cycle time, wafers per hour), chamber-scale validation (literature exposures, Gakis Fiji F200), full 3D chamber with a one-sided pump (106,000 cells), tapered holes, word-line recesses, TDMAH/HfO₂ chamber run, TiN and W by the rule, layout-to-recipe tool on the page, CI, Docker, tests. |
| 2026-11-11 early (night 4) | Reviewer list C: standalone C-ABI kinetics library (aldkin) with OpenFOAM and Fluent adapters, two-species "byproduct" solver (aldkin 1.1), 3D voxel Monte Carlo feature solver (aldmc), slit transport correction and corrected PillarHall refit, GAA nanosheet and 3D-NAND word-line cells as voxel geometries, WebGPU solver on the page, GAA/NAND cells as live 3D scenes. Paper 0.9.4. |
| 2026-11-11 afternoon | External review checked line by line: carrier-flow sentence corrected, sticking range extended to s = 2e-2 (2 runs), headline reworded, production context, test count and CI badge; CI made independent of the Ubuntu mirrors. Paper 0.9.5. |

Wall clock from first commit to last: 2026-11-09 22:40 to 2026-11-11 15:30. 52 commits.

---

## 3. Software, versions and machine

| Component | Version / source | Used for |
|---|---|---|
| OpenFOAM | v2606, gerlero/openfoam-app on Apple silicon (`~/Applications/OpenFOAM-v2606.app`), launcher `.../Contents/Resources/etc/openfoam -c '...'` | chamber flow (simpleFoam), precursor transport (scalarTransportFoam), blockMesh, topoSet/createPatch, custom boundary conditions compiled with wmake |
| OpenFOAM in CI | `opencfd/openfoam-default:2512` container (newest tag on Docker Hub; the cases use nothing newer) | the regression chain on every push |
| Python | 3.12 (venv `.venv`), numpy, scipy, pandas, matplotlib, scikit-learn, pytest, Playwright, pymupdf | feature solvers, post-processing, surrogates, page QA |
| MATLAB | R2025-series, `fitrgp` | Gaussian-process surrogate, grid optimum, Monte Carlo on the rule (`analysis/doe_gp.m`) |
| JMP | response-surface model, script `analysis/doe_rsm.jsl` | DOE fit (`figures/doe_rsm_report.txt`) |
| C++17 / CMake | Apple clang; CMake 4.4.5 (Kitware binary in `~/.local/bin`); ctest | `lib/aldkin` (kinetics, C ABI), `lib/aldmc` (voxel Monte Carlo, std::thread) |
| Machball | Argonne, Yanguas-Gil and Elam (ballistic view-factor model, Python) | cross-check of the hole and slit models |
| three.js 0.186.1, GSAP 3.15, Lenis 1.3.27, WebGPU/WGSL | live page | 26 interactive 3D scenes, browser GP calculator, layout-to-recipe tool, in-browser voxel solver |
| Tectonic | `~/.local/bin/tectonic` | LaTeX build of `paper/main.tex` |
| Node 24 + Lighthouse 13.6 | `~/.local/node` | page audit (100/100/100/100) |
| GitHub Actions, Docker | `.github/workflows/ci.yml`, `Dockerfile`, `reproduce.sh` | four CI jobs; one-command reproduction |
| Machine | one Apple-silicon Mac (M4-class, 7 to 8 threads used) | everything |

---

## 4. The chamber model (reactor/)

**Geometry.** Axisymmetric single-wafer showerhead chamber as a 5° wedge: r = 0 to 200 mm, gap 15 mm, wafer on the
floor at r < 150 mm, non-absorbing floor from 150 to 200 mm, annular outlet at r = 200 mm over the full gap,
showerhead face over r < 150 mm as a uniform downward inflow. blockMesh, 6,000 orthogonal cells, 0.25 mm first cell
at the wafer. Variants: 25 mm edge ring with a pumping slot 25 mm below the wafer plane (`reactor/edgering_*`); gap
10 and 25 mm; a full 3D cylinder (section 7.10); an axisymmetric Fiji F200 (section 7.9).

**Operating point (baseline).** N₂ carrier, 300 Pa, 573 K (300 °C), 300 sccm (Q = 3.54e-3 m³/s at chamber
conditions), showerhead face velocity U = 0.050 m/s, 1 % TMA (c0 = 3.79e20 /m³, 3 Pa). ρ = 1.76e-3 kg/m³, μ =
2.9e-5 Pa s, ν = 1.65e-2 m²/s. Re = U g/ν = 0.045 (creeping flow), residence time V/Q = 0.53 s. TMA diffusivity in
N₂ D = 8.4e-3 m²/s (Chapman-Enskog scaling of 8e-6 m²/s at 1 atm, 300 K; a parameter). Gap Péclet number U g/D =
0.09: transport across the gap is diffusive. H₂O: v̄ = 821 m/s, D = 2.5e-2 m²/s. TDMAH: M = 354.7, v̄ = 185 m/s,
D = 3.62e-3 m²/s.

**Equations.**
1. Steady laminar incompressible N₂ flow, `simpleFoam`, frozen for the species step (the precursor fraction is small).
   Mass balance closes to machine precision (4.90e-5 m³/s per wedge in and out).
2. Precursor as a passive scalar, `scalarTransportFoam`, dimensionless c/c0. Pulse as a time table on the inlet
   (0.1 s to 24 s depending on the case), dt = 1 ms, written every 5 to 100 ms.

**Boundary conditions.**
- Showerhead: Danckwerts total-flux inlet, U (c_face − c0(t)) = D (c_face − c_cell)/δ, as `uniformMixed` with value
  fraction f = U/(U + D/δ) = 1.49e-3. A fixed-concentration inlet is wrong here: with D/L = 0.56 m/s against U =
  0.05 m/s it injects 11× more precursor by diffusion than the gas carries (the first run "absorbed" 1179 % of the
  injected TMA).
- Wafer: the custom conditions of section 5.
- Floor, top wall: zero flux. Outlet: `inletOutlet`. Isothermal.
- Wafer flux = D c_cell/δ from the wafer-adjacent cell row (δ = 0.25 mm).

**Time to saturation** is read as t90, the first write at which the face's coverage reaches 0.9 (interpolated
between writes in the later studies); "worst" means the slowest wafer face.

---

## 5. Surface kinetics and the two-way coupling (reactor/lib/aldWallKinetics, lib/aldkin)

- `langmuirSticking` (blanket wafer): per-face coverage θ integrated in time, uptake (v̄/4) s (1 − θ) c_face, optional
  radical loss r and a first-order desorption rate k_des (reversible water), implemented as a mixed condition.
- `patternedSurface` (patterned wafer): every wafer face carries a flat Langmuir surface plus a 200-cell
  Knudsen/Langmuir model of a hole (Ylilammi 2018 form: quasi-steady Knudsen diffusion D_K = d v̄/3 with a Langmuir
  sink, coverage advanced per cell per time step). The face's net uptake K_eff c_face = (1 − φ) flat uptake + φ ×
  hole-mouth flux, with φ the open-area fraction, enters OpenFOAM as a Robin (mixed) condition: two-way coupling.
  Optional `ARProfile` and `holeDensityProfile` tables of (radius, value) pairs make depth and density functions of
  radius.
- Parameters (baseline): s_TMA = 5e-3 (literature range 5e-4 to 2e-2), s_H2O = 1e-4, N0 = 4e18 sites/m², hole
  d = 100 nm at 180 nm hexagonal pitch, half the wafer patterned (1.78e13 holes/m², 14 % open area), AR 20 / 50 / 100
  (wall-area ratio 11.2 / 28 / 55.9).
- Verification: with negligible hole density the condition's hole-bottom coverage matches the standalone Python
  solver driven by the same face density within 1 % at every time; the night-4 library adapters reproduce the earlier
  boundary conditions bit for bit (201 written fields identical) on the regression case.
- From night 4 the arithmetic lives in `lib/aldkin` (C++17, C ABI: `aldkin_create`, `aldkin_set_faces`,
  `aldkin_set_profile`, `aldkin_step`, `aldkin_mixed_bc`, `aldkin_theta`; 16 ctest cases). The OpenFOAM conditions are
  thin adapters; a Fluent UDF (`adapters/fluent/ald_wall_udf.c`, DEFINE_INIT/ADJUST/SR_RATE) compiles against a stub but
  has not been run in Fluent; Python binding `aldconf.aldkin` (ctypes).

---

## 6. Feature-scale solvers (src/aldconf, lib/aldmc)

| Solver | What it is | Checks |
|---|---|---|
| `aldconf.solver` (1D hole) | quasi-steady Knudsen diffusion + irreversible Langmuir, 200 cells, film narrowing optional, `simulate` (constant mouth density) and `simulate_driven` (recorded chamber history) | Machball: t90 1.035 / 4.230 / 15.93 s at AR 20 / 50 / 100 against converged Machball 1.158 / 4.476 / 16.19 s (3 to 11 % fast; the earlier 0.3 to 7 % used unconverged Machball values at Nz = 2 AR) |
| `TaperedHole` | conservative (1/A) d/dx (A D dn/dx) form on a frustum | bottom_ratio 1 is bitwise the cylinder; zero sink gives uniform density to 1e-9 |
| `Slit(knudsen_correction="plates")` | slit diffusivity × (3/8) ln(L/w) (Clausing/Berman); default "none" keeps earlier numbers | Machball Trench: −3.6 / −3.8 / −2.5 / −0.3 % at AR 10 / 20 / 50 / 100 (uncorrected was 24 % slow at AR 50) |
| `aldconf.recess` `RecessedSlit` | slit (200 cells) feeding word-line recesses (50 cells each), vectorised Thomas solves | AR-1 recesses equal an area-equivalent slit within 3 % |
| `aldconf.multispecies` | precursor + byproduct, two tridiagonal solves per substep, re-adsorption s_B, desorption k_B, etch k_etch, exact 2×2 coverage advance | s_B = 0 reproduces the single-species solver exactly; byproduct balance closes to 1e-14 |
| `lib/aldmc` + `aldconf.voxel` (3D voxel Monte Carlo) | test particles in molecular flow, cosine entry, Langmuir sticking on the face's own coverage, diffuse re-emission, staircase-area correction with smoothed normals, xoshiro256** per block, bitwise reproducible at any thread count | Berman tube transmission 0.108 / 0.059 at L/d 10 / 20 (0.1093 / 0.0595); Machball within 2 % at AR 20 to 100; chamber-driven AR-100 hole 14.32 s vs 14.05 s boundary condition |
| `aldmc-gpu.wgsl` (browser) | the same Monte Carlo in WebGPU with a CPU fallback | GAA preset 0.558 s on Metal, 0.555 s CPU, against 0.558 s in the paper; AR-20 hole 1.144 vs 1.148 s |

---

## 7. Studies and results

### 7.1 Blanket wafer and operating-point sweep (night 1)
Absorbing wafer: the centre flux plateaus at U c0, so the wafer takes everything the showerhead delivers (supply-
limited). Finite sticking: capture 96.5 % (s = 1) to 65 % (s = 1e-3); lower s gives a flatter dose and a longer gas
tail. Langmuir wafer, 1 % TMA: one monolayer of supply every 0.211 s; t90 0.345 s centre / 0.385 s edge. One-factor
sweep (100 to 1,000 Pa, 150 to 600 sccm, 10 to 25 mm, 0.5 to 3 % TMA): pulse time follows supply (doubling fraction or
flow halves it), pressure barely moves it at fixed sccm but sets the edge lag (5 % at 100 Pa, 14 % at 1,000 Pa), the
10 mm gap beats 25 mm on both.

### 7.2 Coupling and loading (night 1)
One-way: the flat wafer "shields" its holes until it saturates (AR 50 bottom 0.92 s against 4.45 s at constant
density). Two-way with the 3D-NAND-like array: the surface density collapses to 0.019 c0, the wafer captures 94 % of
the feed, and a 1 s pulse delivers 8 % of the monolayer-equivalents the wafer needs. Long pulses: hole bottoms reach
90 % at 14.05 s (1 % TMA) and 4.75 s (3 %); the flat regions between holes take 1.95 s instead of 0.35 s.

### 7.3 Pulse-time map and the rule (night 1)
16 two-way runs: blanket and AR 20 / 50 / 100 wafers at 0.5, 1, 3, 10 % TMA. Worst-point t90 (s):

| TMA | blanket | AR 20 | AR 50 | AR 100 |
|---|---|---|---|---|
| 0.5 % | 0.65 | 5.37 | 13.4 | 28.1 |
| 1 % | 0.38 | 2.71 | 6.73 | 14.1 |
| 3 % | 0.18 | 0.95 | 2.32 | 4.79 |
| 10 % | 0.09 | 0.33 | 0.75 | 1.51 |

Patterned points sit at 1.05 to 1.27 × the supply limit; blanket wafers at 1.5 to 4.3 × theirs (kinetics-limited).

### 7.4 H₂O half-cycle, plasma ceiling, temporal vs spatial (night 1)
Water (s = 1e-4): blanket 1.74 s (surface-limited, 5× slower than TMA), patterned 13.4 s (supply-limited, same as
TMA). Plasma ALD: with loss probability r per wall collision, AR 100 costs 10 to 25 % extra at r ~ 6e-5 (O atoms on
oxides) and becomes unreachable between r = 3e-4 and 3e-3. Spatial ALD with a zone flux k = 30 × the showerhead's beats
temporal ALD by a factor that grows linearly with the wall-area ratio (crossover at ratio 5 to 10).

### 7.5 Validation against measurements (nights 1 and 2)
- PillarHall (Zenodo 3974438, 8 profiles, 500 cycles, 103 pm/cycle): s = 8.5e-3 (×/÷ 1.5) with the hydraulic slit
  diffusivity, exposure 7.5 Pa s per pulse, rms 0.038 of the plateau; with the parallel-plate correction s = 1.96e-2 at
  32.3 Pa, same rms, so the profile alone does not pin s.
- Argonne fluidized bed (Lu et al. 2022, 10 g silica gel, 940 m²): τ_sat = 24 min / 21.5 min = 1.12 (0.92 to 1.31);
  the model's patterned cases give 1.05 to 1.27. Not tuned.
- Argonne rotary drum (Coile et al. 2020): utilization 24.7 % at 2 g to 70.7 % at 75 g, falling τ_sat with area, the
  same trend as the wafers.
- Machball: see section 6.

### 7.6 Geometry, water, temperature, hole size (night 2)
Edge ring: saturation times move by at most one write interval; the pulse results do not depend on the chamber
outboard of the wafer at this Péclet number. Reversible water: k_des = 0.1 /s holds the exhaust above 1 % for 8 s
after the pulse; above ~0.3 /s at 1 % H₂O the surface never reaches 90 %. Temperature 200 / 300 / 400 °C at fixed sccm
and fraction: identical pulse times (U c0 is temperature-independent). Holes 50 / 100 / 200 nm at fixed open area and
AR: identical within 1 % (the wall-area ratio is the only geometric input).

### 7.7 DOE, response surface, surrogate, uncertainty (night 2)
Face-centred central composite design, 25 two-way runs on the AR-100 wafer: TMA 0.5 to 3 % (log), N₂ 150 to 600 sccm
(log), pressure 100 to 1,000 Pa (log), gap 10 to 25 mm. JMP: log₁₀ t90 scales as x^−1.0 Q^−0.8 p^−0.12, R² 0.993, gap
not significant; utilization 51 to 97 %, up with pressure and down with carrier flow. MATLAB GP (ARD squared-
exponential, linear basis), leave-one-out R² 0.999 / 0.988 / 0.78 (log t90 / utilization / edge lag); scikit-learn GP
in Python 0.999 / 0.986 / 0.746. Same-supply optimum: 2.5 %, 150 sccm, 631 Pa, 10 mm gives 9.83 s at 94 % utilization
against 11.6 s at 85 % at the centre point. Monte Carlo with N0 lognormal (factor 1.25): 11.9 s, 90 % band 8.2 to
17.2 s; the site density is the dominant absolute uncertainty, the ratios do not depend on it. Pareto front cuts TMA per
wafer pass by 12 % at the centre pulse time.

### 7.8 Spatially varying layouts, sticking sensitivity, fab units (night 3)
- Layouts: a dense AR-100 array at r < 75 mm with bare wafer beyond saturates in 11.85 s instead of 14.05 s and the
  bare edge in 0.40 to 0.50 s instead of 0.35; an AR ramp 100 → 60 gives 12.5 s centre, 7.6 s edge. The rule must use
  the local wall-area ratio of the most demanding region (conservative by 15 to 21 %); the wafer average under-doses.
- Sticking: eleven runs, s = 5e-4 to 2e-2. Patterned bottoms 14.00 to 14.20 s throughout (1.17 to 1.18 × supply);
  blanket 1.01 → 0.28 s (about s^−0.4); flat surface between holes 6.00 → 1.05 s.
- Fab units: the rule's pulse gives 99 % step coverage (front 0.3 s wide). A 10 µm hole purges in 14 µs; the chamber
  reaches 1 % of the pulse-end level in 1.9 s; desorbing water holds it for 8.3 s at k_des = 0.1. Al₂O₃ cycle 12.2 s
  blanket, 37.6 s AR-100, 29.4 s at the same-supply optimum; 5 nm (49 cycles at 0.103 nm) is 5.5 against 1.9 wafers per
  hour with 60 s handling.

### 7.9 Chamber-scale validation (night 3)
Blanket surface exposure to 90 % coverage 0.146 Pa s (1,090 L) sits within a factor 1.8 of the three published 300 °C
kinetic values. Fiji F200 of Gakis et al. (geometry from their figure, 13,020 cells, 180 sccm Ar at 9.6 Pa, 0.63 mg
TMA per 25 ms pulse with their pulse shape, s = 3.44e-3, N0 = 4.1e18): delivered exposure 0.67 Pa s against their
CFD's 0.55 to 0.60 (12 to 22 %), saturation state reproduced (25 ms at the knee, 60 ms past it), their 2.4 % non-
uniformity not reproduced (0.0007 %): it comes from asymmetries and thermal gradients an axisymmetric isothermal
wedge cannot carry.

### 7.10 Full 3D chamber with a one-sided pump (night 3)
blockMesh O-grid, 105,984 hexahedra, 12 cells across the gap, 60° pump sector on the side wall, simpleFoam converged
in 93 iterations (18 s), scalarTransportFoam 16 s pulse at dt 1 ms, 1,939 s wall serial, 3.4 GB of output. The gas
moves (outward speed 0.85 m/s port side, 0.01 far side), but hole-bottom t90 spreads only 0.96 % azimuthally at
r = 140 mm; worst face 14.17 s against 14.05 s in the wedge (+0.8 %). A one-sided pump costs about 1 % of pulse time.

### 7.11 Device geometry and other chemistries (night 3)
Tapered AR-100 holes (bottom 100 → 60 nm): +1 to 4 % on the supply-limited wafer. Word-line recesses (29 nm gap,
51 nm pitch, 5 µm slit, 250 or 100 nm slit width, 365 nm or 1 µm recess length): the last recess lags the slit bottom
by 0.5 to 9 %; an area-equivalent slit within 4 %. TDMAH/HfO₂ chamber run (v̄ 185 m/s, s 7e-3, N0 2.77e18): 11.60 s
against a rule of 9.96 s (ratio 1.165). By the rule alone: TiN/TiCl₄ 4.44 s, W/WF₆ 56.7 s; per nanometre of film HfO₂
100, Al₂O₃ 144, TiN 148, W 227 s/nm.

### 7.12 Two-species kinetics (night 4)
At Kn = 123 to 1,233 in a 100 nm hole a byproduct cannot block the precursor in the gas phase; it acts only on the
surface. TMA/CH₄: identical to the single-species solver (methane is inert; 99.95 % escapes even at a hypothetical
s_B = 1e-3). TiCl₄/HCl at AR 100: no effect for k_B ≥ 100 /s; +11 to +74 % pulse at k_B = 1 /s depending on s_B and
stoichiometry; failure of saturation at k_B ≤ 0.1 /s or with etching. The HCl rate at 573 K on TiN is not in the
literature; every number is labelled as a range or bound.

### 7.13 Device cells (night 4)
GAA nanosheet stack (IBM 2017: 3 × 5 nm sheets at 10 nm spacing, W 40 nm, Lg 12 nm, voxel 1.25 nm, 50,000
particles per step): everywhere-t90 0.558 s TMA, 0.649 s TDMAH, 0.056 s TiCl₄ at 0.3 Pa; last point the bottom-slot
centre; wall-area ratio 5.7 per wafer at 44 nm CPP, so the chamber-scale loading is small. NAND word-line slab
through a pillar row (5 µm cell): slit floor 9.04 s, slowest recess 11.30 s (recess 86 of 100 at 87 % depth), t99
12.13 s; 2XX-layer cell (12.7 µm, 264 recesses): slowest 33.1 s at 0.3 Pa, wall-area ratio 448 per wafer. These are
applications of the solvers, not validations: no measured per-sheet or per-level ALD thickness exists in the open
literature.

### 7.14 Slit transport correction and corrected refits (night 4)
Transmission Monte Carlo of the infinite slit (4e5 particles per point) confirms the (3/8) ln(L/w) factor; Machball
Trench benchmarks at AR 10 to 100; PillarHall refit with the correction gives s = 1.96e-2 (×/÷ 1.46 sheet to sheet)
at 32.3 Pa; word-line recess times restated 5 to 30 % lower.

---

## 8. Analysis and tooling

- `analysis/doe_rsm.jsl` (JMP), `analysis/doe_gp.m` (MATLAB fitrgp, grid optimum, Monte Carlo), `analysis/doe_gp.py`
  (scikit-learn GP, exported to `figures/doe_gp_python.json` for the browser calculator), Pareto front.
- CLI: `python -m aldconf rule|surrogate|pareto`.
- `coupling/` post-processing scripts, one per study, each writing the CSV its figure is drawn from.
- `paper/wp/wp1 … wp13` are the per-item deliverable notes of nights 3 and 4; `paper/assumptions.md` is the run log
  with every setting and wall time; `CLAIMS.md` traces every quantitative statement to a file and a command.

---

## 9. The live page (docs/)

Static site on GitHub Pages from `docs/`. `index.html` with 13 sections: hero (scroll-driven 16 s pulse), 01 model,
02 rule and pulse-time map, 02b layouts, 03 validation, 04 hole, 05 knobs (DOE surface), 06 calculator (rule plus the
GP surrogate in the browser), 07 water/plasma/spatial, 08 robustness, 09 reproduce, 11 reviewer's questions (13 scenes:
sticking, step coverage, cycle, WPH, exposure, Gakis, 3D wafer, azimuth, taper, recess, chemistry, GAA cell, NAND
cell), 10 layout-to-recipe tool, 12 WebGPU voxel solver. 26 three.js scenes (`scene.js`, `body.js`, `body2.js`), data in
`assets/data.js` (585 KB, generated by `docs/build_assets.py --data` from the figure CSVs), audit figures copied as
`assets/fig_*.png`, dark and light screenshots `assets/viz_*`. QA by `docs/tools/page_qa.py` (Playwright): console
errors, failed requests, images, mobile layout, screenshots in `docs/qa/`; Lighthouse 100 / 100 / 100 / 100 (full
Chromium, `--disable-gpu`). Social preview `assets/social-preview.png` generated by `build_assets.py --preview` and
set on the repository.

---

## 10. The paper (paper/)

`main.tex`, built with Tectonic to `main.pdf`: draft 0.9.5, 42 pages, 25 figures, 17 tables, 54 references,
19,000 words of LaTeX source. Structure: Introduction; Model (chamber; wafer boundary condition, concurrent multiscale;
verification and validation; validation against a measured high-surface-area run; chamber-scale checks; the 3D voxel
Monte Carlo solver); Results (blanket wafer; one-way coupling and shielding; two-way coupling and loading; the
pulse-time map and the rule; H₂O half-cycle; temporal versus spatial ALD; plasma recombination ceiling; operating
optimum and uncertainty; non-uniform layouts; step coverage, purge and throughput; tapered holes, word-line recesses
and other chemistries; a full 3D chamber with a one-sided pump; byproducts; gate-all-around stacks and a NAND cell);
Discussion; Conclusions; code and data availability. `draft.md` is the Markdown draft it
was typeset from; `assumptions.md` the run log; `README.md` the build notes. Title: "A pulse sized on a blanket wafer
is 40× too short for a 3D-NAND wafer: a chamber-to-feature model of supply-limited atomic layer deposition".

Honest framing carried in the paper: the Fluent UDF was never run in Fluent; the GAA and NAND cells are applications,
not validations; the HCl constants are ranges; the TDMAH sticking is a model calibration from the literature, not a
measurement; the chamber-scale non-uniformity of Gakis et al. is not reproduced; isothermal gas, frozen carrier flow,
irreversible Langmuir in the main runs, no film narrowing in the voxel solver, s and N0 are literature-range inputs.

---

## 11. Repository, tests, CI, Docker

- `src/aldconf/` package (`pip install -e '.[test,crosscheck]'`), 129 pytest tests (feature solvers, scaling laws,
  Machball cross-check, slit transport, recesses, device geometries, library bindings, OpenFOAM regression against a
  stored reference to 1e-6, page assets).
- `lib/aldkin` (16 ctest), `lib/aldmc` (10 ctest), `adapters/fluent`.
- `.github/workflows/ci.yml`, four jobs on every push: Python solver and tests; both C++ libraries with their Python
  bindings; the OpenFOAM chain in the OpenCFD container (library build, wmake of the adapters, steady flow, a 2 s
  regression pulse, the regression test, the full suite); live-page asset check. Python 3.12 in the container comes from
  uv (python-build-standalone) and cmake from PyPI, so the jobs do not depend on the Ubuntu package mirrors.
- `Dockerfile` (FROM `opencfd/openfoam-default:2512`) and `reproduce.sh --build | --regress | --quick | --full`.
- Git history was rewritten once (git-filter-repo, force-push approved) to drop 54,364 OpenFOAM output files that had
  been committed by mistake: `.git` went from 1.1 GB to 8.4 MB. Output directories are now ignored; case inputs
  (`0/`, `constant/`, `system/`) and `runtime.txt` are tracked.

---

## 12. Compute and effort

| Item | Value |
|---|---|
| OpenFOAM cases on disk | 145 (57 night 1, 86 night 2 and 3, 2 on 2026-11-11) |
| Solver output on disk | 25,757 time-step folders, 29 GB |
| Typical wall time per case | 84 to 96 s blanket, 104 to 206 s patterned wedge, 1,939 s full 3D chamber, 77 s Fiji F200 steady flow |
| Voxel Monte Carlo | cross-check 24 runs, 6.1e8 particles, 93 min on 7 to 8 threads; device cells 35 runs, 5.9e8 particles, 100 min on 6 threads; longest single run 834 s |
| Other | PillarHall corrected refit 25 min single core; slit transmission Monte Carlo 2 min; Machball runs 1 to 20 s each |
| Total solver time | a few CPU hours in all, on one Mac |
| Code | about 40,000 lines of Python, C++, C, JavaScript and WGSL (excluding the generated data file and vendored libraries) |
| Outputs | 29 CSV tables, 36 figures, 26 page scenes, the paper |
| Commits | 52 over 41 hours of wall clock |

---

## 13. Corrections made along the way (so nobody repeats them)

- Fixed-concentration showerhead inlet over-injected 11×; replaced by the Danckwerts total-flux condition.
- Machball reference values first quoted (1.29 / 4.46 / 15.6 s) were unconverged (Nz = 2 AR); converged values are
  1.158 / 4.476 / 16.19 s, and the 1D model is 3 to 11 % fast, not 0.3 to 7 %. Corrected everywhere, including the
  profile repository.
- The 1D slit Knudsen diffusivity (hydraulic diameter) was 24 % slow at AR 50; the parallel-plate factor fixes it and
  the PillarHall s becomes a range (8.5e-3 to 2.0e-2) depending on the transport assumed.
- The pulse-map caption said carrier flow does not enter the rule; it enters through U c0 (DOE: Q^−0.8).
- "Every blanket-wafer ALD recipe hides a 40× under-dose" overclaimed; the factor is for AR 100 at this density.
- 1/s extrapolation of the blanket pulse at low sticking was wrong by 3.5× (the gas density rises as uptake slows).
- CI: three runs failed on runner DNS to the Ubuntu mirrors; the workflow now avoids apt.

---

## 14. Reproduce

```bash
git clone https://github.com/Soham-ChemE/ald-conformality.git && cd ald-conformality
python3.12 -m venv .venv && .venv/bin/pip install -e '.[test,crosscheck]' && PYTHONPATH=src .venv/bin/python -m pytest -q tests
```

Chamber chain (OpenFOAM v2606 on the Mac, or the Docker image anywhere):

```bash
bash ./reproduce.sh --build && bash ./reproduce.sh --regress     # libraries, adapters, steady flow, regression pulse
docker build -t ald-conformality . && docker run --rm ald-conformality --quick
```

Individual studies: `reactor/Allrun`, `reactor/sweep_pattern_map.py`, `reactor/doe_patterned.py`,
`reactor/sweep_sticking_sensitivity.py`, `reactor/make_layout_cases.py`, `reactor/chamber3d/make_and_run.py`,
`reactor/validation_gakis/`, `reactor/chemistry/make_tdmah_case.py`, `coupling/*.py`, `analysis/doe_gp.{m,py}`;
page assets `python docs/build_assets.py --data --figs --preview --qa --lighthouse`; paper `cd paper && tectonic main.tex`.

---

## 15. Repository layout

```
src/aldconf/        feature solvers (solver, geometry, transport, recess, multispecies, voxel, aldkin binding,
                    surrogate, fit, cli)
lib/aldkin/         C-ABI wall-kinetics library (CMake, 16 ctest, examples, README)
lib/aldmc/          C-ABI voxel Monte Carlo feature solver (CMake, 10 ctest)
adapters/fluent/    Fluent UDF against the kinetics library (compiled against a stub, not run in Fluent)
reactor/            OpenFOAM: lib/aldWallKinetics (adapters), showerhead_* cases, map/, doe/, sticking/,
                    edgering_*, holesize_*, temperature, chamber3d/, chemistry/, validation_gakis/, Allrun,
                    sweep_pattern_map.py, doe_patterned.py, sweep_sticking_sensitivity.py, make_*_cases.py
coupling/           one post-processing or study script per result (reactor_to_feature, validate_supply_rule,
                    fit_pillarhall(_corrected), fab_units, tapered_hole, wordline_recess(_corrected),
                    multispecies, voxel_mc_crosscheck, device_geometries, postprocess_* ...)
analysis/           doe_rsm.jsl (JMP), doe_gp.m (MATLAB), doe_gp.py (scikit-learn), literature compilation
figures/            29 CSV tables and 36 PNG figures, every one written by a script above
tests/              129 pytest tests (19 files), stored regression reference
docs/               live page: index.html, assets/ (scenes, data.js, GPU solver, kit), build_assets.py,
                    tools/page_qa.py, qa/ (screenshots, Lighthouse), partials/, DESIGN.md
paper/              main.tex, main.pdf, draft.md, assumptions.md (run log), wp/ (deliverable notes), README.md
CLAIMS.md           every quantitative claim traced to a file and a command
CITATION.cff, Dockerfile, reproduce.sh, .github/workflows/ci.yml
```
