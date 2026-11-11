# Every blanket-wafer ALD recipe hides a 40x under-dose: a chamber-to-feature model of supply-limited atomic layer deposition

Soham Kavathekar, Department of Chemical and Biomolecular Engineering, University of Pennsylvania

*Draft 0.9.4 (public), 2026-11-11. Every number is from a run documented in `paper/assumptions.md`.*

## Abstract

ALD recipes are tuned on blanket wafers. Product wafers are not blank: a 3D-NAND array of 100 nm holes at aspect
ratio 100 carries 56 times the surface area of the wafer under it, and all of it has to saturate in one pulse. We
built a model of a 300 mm showerhead ALD chamber in which every point on the wafer runs its own feature-scale model
of a hole, coupled both ways to the chamber gas, and asked how long the pulse must be. The answer is a single rule:
the saturating pulse equals 1.2 times (1 + patterned area ratio) times the time the showerhead needs to deliver one
monolayer of precursor. It held within 25 % across aspect ratios 20 to 100 and a twenty-fold range of precursor feed;
on non-uniform layouts the local wall-area ratio of the densest region sets the pulse, not the wafer average.
At a typical operating point the blanket wafer saturates in 0.35 s and the product wafer in 14 s. Chamber pressure,
showerhead gap and carrier flow do not move that number; precursor delivery rate does. Nor does the sticking coefficient across its published range; only the blanket pulse
does, threefold. The water half-cycle, fifty
times less reactive, needs the same pulse on the product wafer and coats it more evenly. The feature model sits within
3 to 11 % of a converged ballistic reference (Argonne's Machball and this work's own 3D Monte Carlo solver, which
agree within 2 %) and fits published thickness profiles from lateral test structures,
and the rule itself reproduces the measured saturation time of a fluidized-bed ALD run on 940 m2 of silica gel
to 7 %. The rule
turns a day of recipe transfer into arithmetic, and says where a tool's throughput actually comes from.

## What this means for a tool (one paragraph for the process engineer)

If your recipe was saturated on a monitor wafer, multiply the pulse by one plus the product wafer's wall-area ratio
before you trust it on product. Size the pulse on the densest or deepest region's own wall-area ratio; the wafer
average under-doses it. Pressure, gap and carrier flow are uniformity knobs; they buy no time. Time comes
from precursor molecules per second: vapour pressure, bubbler flow, injector design, or an architecture that exposes
the wafer to a denser zone (spatial ALD wins by the same factor the rule predicts). Low-sticking reactants are not
the bottleneck on patterned wafers; their cost is on blanket wafers only. And when the deepest feature saturates,
the chamber shows it: the precursor stops being consumed and appears downstream. On a 3D-NAND wafer the pulses are
three quarters of the cycle and the purges a quarter, the reverse of the blanket wafer, so the throughput lever is
precursor delivery and the purge lever is water desorption, not chamber flow.

## 1. Introduction

Three facts set up the problem. First, 3D NAND and DRAM are built from features with aspect ratios of 50 to 100, and
ALD is the only way to coat their insides, so the precursor must reach the bottom of every feature on the wafer
inside a single pulse. Second, the patterned surface area of such a wafer is tens of times the wafer area, so the
wafer consumes precursor tens of times faster than the blank wafer the recipe was developed on. Third, nobody models
these two things together: feature-scale conformality models [Gordon 2003; Yanguas-Gil and Elam; Ylilammi 2018;
Cremers 2019] assume the precursor density at the feature mouth is given, and chamber models assume the wafer's
uptake is given. The quantity that matters, the gas density at the wafer surface while the wafer is eating it,
falls between the two.

This paper closes that gap with a model small enough to run on a laptop and checked hard enough to trust. A chamber
flow model carries, on every wafer face, a feature-scale model of a hole; the two exchange precursor every time
step. The feature model is verified against an independent ballistic code and fitted to published measurements. The
result is not a simulation to be rerun for every recipe but a rule: pulse time is set by precursor supply and
patterned area, and by nothing else that a chamber recipe normally adjusts.

Equipment makers are building exactly this class of model under the name of equipment-and-process digital twins
[Applied Materials, Computational Products and Solutions; AppliedTwin]. This work is an open, reproducible instance
for ALD, with every case file, boundary condition and figure in the repository.

Contributions: (i) a two-way coupled chamber-feature wall condition for OpenFOAM, released; (ii) three-method
verification (continuum, this work's voxel Monte Carlo, Argonne's Machball) and public-data validation of the feature
model; (iii) the supply-limit rule and its range of validity;
(iv) the blanket-to-product recipe-transfer factor, the shielding mechanism and the edge-lag inversion; (v) the
temporal-versus-spatial design chart and the plasma-ALD recombination ceiling.

## 2. Model

### 2.1 Chamber
Axisymmetric showerhead chamber, 300 mm wafer, 15 mm gap, 200 mm chamber radius, annular outlet. N2 carrier at
300 Pa, 573 K (rho 1.76e-3 kg/m3, nu 1.65e-2 m2/s), 300 sccm (0.050 m/s at the showerhead face), Re 0.045,
residence time 0.53 s. Steady laminar flow (simpleFoam), frozen for the species step. TMA as a passive scalar
c/c0, D = 8.4e-3 m2/s at 300 Pa (Chapman-Enskog scaling), gap Peclet number 0.09.
Inlet: Danckwerts total-flux condition U(c_face - c0) = D (c_face - c_cell)/delta. We show that a fixed-concentration
inlet over-injects by 11x at this Peclet number (Section 3.1).

### 2.2 Wafer boundary condition (concurrent multiscale)
Flat surface: irreversible Langmuir, dtheta/dt = (vbar/4) s (1 - theta) c_s/N0.
Holes: cylinder of diameter d and aspect ratio AR, Knudsen diffusivity d vbar/3, quasi-steady gas (Ylilammi 2018)
with sink (vbar/d)[s(1 - theta(x)) + r], coverage advanced per cell; mouth density equals the local wafer-surface
gas density c_s. Net uptake per unit wafer area K_eff c_s = (1 - phi) flat uptake + phi x hole-mouth Knudsen flux,
phi the open-area fraction. Implemented as an OpenFOAM mixed condition with value fraction K_eff/(K_eff + D/delta).
Parameters: s = 5e-3 (TMA) and 1e-4 (H2O) [Arts 2019], N0 = 4e18 m-2, vbar from kinetic theory at 573 K,
d = 100 nm, hole density 1.78e13 m-2 (180 nm hexagonal pitch, 50 % of the wafer patterned): open area 14 %,
A_wall/A_wafer = 11, 28, 56 for AR 20, 50, 100.

### 2.3 Verification and validation
- Scaling tests on the feature solver: penetration depth ~ sqrt(dose); dose for 95 % step coverage ~ AR^2 in the
  diffusion-limited regime (ratio 4.0 for AR 20 -> 40 at s = 0.5; drops to 3.2 at s = 0.05 where the front becomes
  reaction-limited); sticking coefficient recovered from a synthetic profile within 5 % (log).
- Cross-code check: Machball (ballistic, cosine re-emission) vs this Knudsen-continuum solver, 100 nm holes at
  0.3 Pa TMA, 573 K, with Machball's wall discretisation converged (Nz = 8 to 16 AR): 90 % bottom coverage at
  1.035 vs 1.158 s (AR 20), 4.230 vs 4.476 s (AR 50), 15.93 vs 16.19 s (AR 100), 95 % at 1.20 vs 1.33, 4.45 vs
  4.70 and 16.2 vs 16.5 s; the continuum solver is 11, 5.5 and 1.6 % fast against converged Machball and 11, 7 and
  3 % fast against the voxel Monte Carlo solver of Section 2.6, which agrees with converged Machball within 2 %. An
  earlier draft of this paper quoted Machball at Nz = 2 AR (95 % times 1.29 / 4.46 / 15.6 s), which Section 2.6 shows
  is not converged; the values here supersede them.
- Embedded-BC check: with negligible hole density the boundary condition's hole-bottom coverage matches the
  standalone solver driven by the same surface density to 1 % at all times.
- Validation against the Zenodo 3974438 PillarHall Al2O3 profiles (TMA/H2O, 300 C, 500 cycles x 0.1 s, 500 nm
  lateral channel, eight profiles averaged): the slit model with channel narrowing fits the saturation front to an
  rms of 0.038 of the plateau (Figure, pillarhall_fit.png). The front position (AR 215 at half thickness) fixes the
  exposure and its width fixes s, but both depend on the slit's Knudsen diffusivity. With the hydraulic-diameter
  value D_K = 2 w vbar/3 the fit gives s_TMA = 8.5e-3 (profile-to-profile factor 1.5) and 7.5 Pa s per pulse (75 Pa
  TMA for 0.1 s). That value of D_K is exact for a cylinder but not for parallel plates, where grazing flights make
  the free-molecular diffusivity grow as (w vbar/4) ln(L/w) (Clausing; Berman 1965): the 1D slit model is 24 % slow
  against Machball's ballistic trench and the voxel solver at AR 50 and 61 % at AR 100, and with the factor
  (3/8) ln(L/w) on D_K it is within 5 % of Machball from AR 20 to 100 and reproduces the ballistic front width to
  about 10 % at AR 500. With that correction (a factor 2.33 at AR 500) the same data give s_TMA = 1.96e-2
  (sheet-to-sheet x/÷ 1.46) and 3.23 Pa s per pulse (32.3 Pa), at the same rms (Figure, pillarhall_fit_corrected.png):
  the profile cannot distinguish the two transport models, it only re-labels the parameters, s by the factor D rose
  and the exposure by its inverse. Published values are 0.5-2e-3 (Arts 2019, 275 C set-point, extracted with a
  continuum slit model of the hydraulic-diffusivity family), 4e-3 (Vandalon and Kessels 2016, SFG, 300 C) and 5.72e-3
  (Ylilammi 2018, 300 C). The corrected value lies farther from all of them than the uncorrected one; the measured
  front is sharper than a ballistic slit with the literature s would make it. The ideal-plate factor is an upper bound
  for the real cavity, which is 140-200 gaps wide and carries support pillars that interrupt the grazing flights
  (finite width alone removes 8-13 % of the factor), so s from this profile lies between 8.5e-3 and 1.96e-2 depending
  on the cut-off length of the cavity; the chamber runs use 5e-3 and the patterned-wafer conclusions are
  supply-limited and insensitive to it (Section 4).

### 2.4 Validation against a measured high-surface-area ALD run
A patterned wafer with tens of wafer-areas of holes is, to the chamber, a powder bed. The cleanest published
measurement of supply-limited ALD is Argonne's fluidized-bed Al2O3 run on 940 m2 of silica gel [Lu et al. 2022]:
in-situ mass spectrometry put TMA breakthrough at 24 min, and the authors' own accounting gives 0.0115 +/- 0.0020 mol
of TMA supplied in that time against 0.0103 mol needed. In the normalisation used throughout this paper, saturation
arrived at 1.12 (0.92 to 1.31) times the supply-limit time. Our twelve patterned-wafer cases, with nothing tuned to
that experiment, saturate at 1.05 to 1.27 times theirs (Figure, validation_supply_rule.png); the blanket wafers, at
a Damkoehler number of 10 instead of hundreds, sit at 1.5 to 4.3 and are the kinetics-limited exception. Argonne's
plug-flow theory [Lu 2022; Yanguas-Gil and Elam 2024] predicts the same approach to unity from above. A second
measured set, TMA utilization in a rotary-drum reactor at five silica-gel loadings from 40 to 1500 m2 [Coile et al.
2020, Table 2], rises from 25 % to 71 % with area, i.e. saturation moves from 4.0 to 1.4 supply-limit times as the
reactor fills with surface; our blanket and patterned wafers (0.07 to 4 m2) and the fluidized bed (940 m2, 92 %)
fall on the same trend across five orders of magnitude of surface area (Figure, right panel).

### 2.5 Chamber-scale checks
The blanket-wafer pulse of Section 3.1 can be compared with published saturation data once the exposure is defined
with care. The recipe sets a feed exposure, 3 Pa of TMA for 0.345 s, 1.03 Pa s or 7,800 Langmuir; the surface sees
the time integral of its own TMA partial pressure, 0.146 Pa s or 1,100 Langmuir, because the wafer draws the gas
down to 0.14 of the feed value on average while it saturates, and that surface exposure equals the Langmuir value
ln(10) N0/(s F1) at s = 5e-3, as it must. Seven independent sources give either the exposure they applied or a
sticking coefficient or reaction cross-section from which the planar saturating exposure follows (Figure,
blanket_saturation_literature.png; table below). The kinetic values span 0.05 to 1.3 Pa s: 0.27 Pa s from the
sum-frequency cross-section of Vandalon and Kessels, 0.21 from the reactor-fitted sticking coefficient of Gakis et
al., 0.16 from the PillarHall fit of Ylilammi et al., 0.33 to 1.3 from Arts et al. at 220 C, 0.05 at room
temperature in vacuum (Soto and Tysoe, as summarised by Puurunen 2005), 0.08 from our own PillarHall fit. Our
0.146 Pa s sits in the middle of the 300 C cluster, within a factor 1.8 of each of the three sources that state a
sticking coefficient at that temperature; the applied exposures (0.56 to 160 Pa s) are upper bounds that recipes
choose to over-dose and say nothing sharper. Exposure scales as 1/s in this regime, so the model would reproduce the
individual sources at s between 2.7e-3 and 8.7e-3.

| Source | T (C) | Kind of exposure | Exposure (Pa s) | Langmuir | s this model would need |
|---|---|---|---|---|---|
| This work, chamber model, centre: feed p x t90 | 300 | applied | 1.03 | 7,760 | 5e-3 (as run) |
| This work, centre: surface integral of p_s dt to theta = 0.9 | 300 | kinetic | 0.146 | 1,090 | 5e-3 (as run) |
| Ylilammi et al. 2018, Picosun R-150, 300 Pa, 0.1 s pulse: fitted entrance pressure 147 Pa x 0.1 s (flat region saturated) | 300 | applied | 14.7 | 110,000 | upper bound only |
| same, implied by their fitted gamma = 5.72e-3 (q = 5e18 /m2) | 300 | kinetic E90 | 0.155 | 1,160 | 4.7e-3 |
| Arts et al. 2019, Oxford OpAL, set-point 275 C: 190 mTorr s TMA dose, planar surface saturated | ~220 | applied | 25.3 | 190,000 | upper bound only |
| same, 1200 mTorr s (H2O-limited profile) | ~220 | applied | 160 | 1,200,000 | upper bound only |
| same, implied by s0 = 2e-3 / 5e-4 (Table I range; N0 = 4e18 assumed) | ~220 | kinetic E90 | 0.33 / 1.32 | 2,470 / 9,890 | 2.2e-3 / 5.5e-4 |
| Vandalon and Kessels 2019, SFG, sigma_TMA = 6.5 +/- 0.6 e-18 cm2; temperature-independent 100-300 C (2017) | 300 (100) | kinetic E90 from sigma, no site density needed | 0.273 (0.221) | 2,050 (1,650) | 2.7e-3 (3.3e-3) |
| Gakis et al. 2018, Fiji F200, 9.6 Pa, 25 ms valve: wafer impingement 0.0118-0.0130 mol/m2 (+/- 0.0002) x 46.45 Pa s per mol/m2; GPC 0.975-0.998 A/cycle, i.e. >= 97 % saturated | 300 | applied (their CFD) | 0.55-0.60 | 4,100-4,500 | upper bound only |
| Gakis et al. 2019, implied by fitted s0 = 0.004 x 86 % reaction probability, N0 = 4.1e18 from GPC 1.0 A/cycle at 3500 kg/m3 | 300 | kinetic E90 | 0.212 | 1,590 | 3.4e-3 |
| PillarHall record (Zenodo 3974438), this work's fit: mouth exposure 75 Pa x 0.1 s | 300 | applied | 7.5 | 56,000 | upper bound only |
| same, implied by the fitted s = 8.5e-3 | 300 | kinetic E90 | 0.084 | 630 | 8.7e-3 |
| Soto and Tysoe 1991, UHV, as summarised by Puurunen 2005: reaction probability 0.01 | 25 | kinetic E90 (N0 = 4e18 assumed) | 0.051 | 380 | 1.4e-2 |
| Dillon et al. 1995, porous alumina FTIR, as summarised by Puurunen 2005: k_a = 0.0076 /s at 0.01 Torr | 227 | first-order rate on a porous membrane: transport, not planar kinetics | 404 | 3,000,000 | (1.8e-6, not meaningful) |

A second check uses the one published case that combines reactor-scale CFD with a measured thickness map in a
commercial tool: the Ultratech Fiji F200 of Gakis et al. (200 mm wafer, 9.6 Pa Ar, 180 sccm, 300 C, 0.63 mg of TMA
per 25 ms valve opening, measured growth per cycle 0.975 to 0.998 A/cycle along the diameter against 1.00 at 60 ms,
non-uniformity 2.4 % falling to 0.58 %). We rebuilt the chamber as an axisymmetric wedge with dimensions read off their
cross-section (the papers give none), a single axial inlet standing in for their three, the carrier flow frozen,
their computed pulse amount and shape, and their kinetic model reduced to an irreversible sticking coefficient
s = 0.004 x 0.86 = 3.4e-3 and N0 = 4.1e18 /m2 from their growth per cycle (Figure, validation_gakis.png; table
below). The model delivers 0.67 Pa s of TMA impingement to the wafer against 0.55 to 0.60 Pa s in their
three-dimensional CFD and saturates the wafer with the 25 ms pulse (coverage 0.9994, 90 % at 70 ms) and completely
with the 60 ms pulse, as the measurement does to within 2.5 %. It does not reproduce the 2.4 % non-uniformity: their
exposure field is 10.7 % non-uniform because the pulse momentum through an off-axis inlet drives a recirculation over
the wafer and the loading-door purge sweeps one side, ours is 0.3 % because an axisymmetric inlet into a diffusive
gap smooths everything; and with irreversible Langmuir kinetics even their exposure variation would move coverage by
under 0.5 %, so the remainder is their reversible, softer kinetics and the 21 C temperature gradient and cold loading
door they document. The chamber-scale transport that this paper's conclusions rest on, how much of the injected
precursor the wafer takes before the pump does, is reproduced to the accuracy of the dimensions; the within-wafer
uniformity of a particular injector is not a claim of this paper.

| Case | theta after the cycle, centre / min / max (r <= 95 mm) | non-uniformity (max-min)/min | t90 after valve opening, centre / worst | wafer impingement exposure (Pa s) | exposure non-uniformity | TMA captured by the wafer | peak surface TMA (Pa) |
|---|---|---|---|---|---|---|---|
| 25 ms valve, 0.63 mg, s = 3.44e-3 | 0.9994 / 0.9994 / 0.9994 | 0.0007 % | 0.07 / 0.08 s | 0.671-0.673 | 0.28 % | 2.4 % of 5.26e18 (41 monolayer-equivalents delivered) | 5.6 |
| 60 ms, 1.5 x amount | 1.0000 / 1.0000 / 1.0000 | 0.00001 % | 0.06 / 0.06 s | 1.017-1.019 | 0.18 % | 1.6 % | 8.5 |
| 60 ms, 2.4 x amount | 1.0000 / 1.0000 / 1.0000 | 0.0000 % | 0.05 / 0.05 s | 1.639-1.640 | 0.11 % | 1.0 % | 13.7 |
| 25 ms, s = 5e-3 | 1.0000 / 1.0000 / 1.0000 | 0.00003 % | 0.06 / 0.07 s | 0.671-0.673 | 0.28 % | 2.4 % | 5.6 |
| Measured, 25 ms (2019 Fig. 10, digitised, / 1.00 A per cycle at 60 ms) | 0.996 / 0.975 / 0.998 | 2.37 % (2.40 % stated) | | | | | |
| Measured, 60 ms (stated) | 1.00 | 0.58 % | | | | | |
| Their CFD + kinetics, 25 ms / 60 ms | | 2.45 % / 0.83 % | | 0.55-0.60 (2018 Fig. 12) | 10.7 % | | 0.7-3.7 at 25 ms (2019 Fig. 7b) |

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
aspect ratio grows (Figure, voxel_mc_crosscheck.png; table below). The particle budget (25 000 to 400 000 per step),
the voxel size (16 to 32 across the diameter) and the seed change t90 by at most 1 %; the time step enters linearly
through the coverage frozen within a step (0.6 % at dt = t90/200). The frustum with a 60 % bottom lengthens t90 by
the factor 1.19 in the voxel solver and 1.20 in the 1D model. Driven by the chamber's centre-face density history of
the long AR-100 case, the voxel hole reaches 90 % bottom coverage at 14.32 s (14.55 s with the history subsampled to
0.25 s) against 14.03 s in the 1D solver and 14.05 s in the boundary condition, with the same saturated and empty
regions and a somewhat sharper front at 5, 10 and 13 s. For slits the two ballistic codes again agree (AR 50:
1.08-1.11 s voxel, 1.12 s Machball), but the 1D slit model with D_K = 2 w vbar/3 is 24 % slow, the known logarithmic
correction of the parallel-plate Knudsen diffusivity; the word-line recess structure therefore saturates 27 % sooner
in the voxel solver (10.7 s against 14.5 s at 0.3 Pa) with the same ordering of slit bottom and recess ends. The
slowest point of such a stack is not the floor: the per-recess saturation time peaks at 87 to 90 % of the depth
(Section 3.14), because the saturated floor reflects the downward stream back into the lowest recesses, so the
minimum over the recesses (12.06 s) is that recess and the deepest recess (10.88 s) is 10 % ahead of it. The
chamber results rest on the hole model, whose Knudsen diffusivity d vbar/3 is exact; the slit-based numbers of
Section 3.11 are conservative by this factor (the parallel-plate correction of Section 2.3 removes it). The 24 voxel runs behind these numbers used 6.1e8 particles and 93 min
of wall time in all on 7 to 8 threads, no single run above 10 min, against under a second per 1D solve and 1 to
131 s per converged Machball run.

| AR | voxel t90 / t95 (s) | 1D Knudsen t90 / t95 (s) | Machball converged t90 / t95 (s) | Machball Nz = 2 AR t90 / t95 (s) | voxel vs Machball | voxel vs 1D |
|---|---|---|---|---|---|---|
| 20 | 1.148 / 1.316 | 1.035 / 1.204 | 1.158 / 1.334 | 1.123 / 1.289 | -0.9 % / -1.3 % | +10.9 % / +9.3 % |
| 50 | 4.513 / 4.735 | 4.230 / 4.451 | 4.476 / 4.704 | 4.247 / 4.465 | +0.8 % / +0.7 % | +6.7 % / +6.4 % |
| 100 | 16.454 / 16.784 | 15.931 / 16.249 | 16.191 / 16.511 | 15.279 / 15.590 | +1.6 % / +1.7 % | +3.3 % / +3.3 % |

Conditions: s = 5e-3, N0 = 4e18 /m2, 573 K, 0.3 Pa TMA, 100 nm holes; voxel runs 16 voxels across the diameter,
dt 0.005 / 0.0125 / 0.025 s, seed 1; Machball converged at Nz = 8 AR (16 AR at AR 100).

## 3. Results

### 3.1 Blanket wafer
Fixed-concentration inlet absorbs 1179 % of the injected TMA; total-flux inlet 96.5 %. Centre flux plateaus at the
supply rate U c0: at Pe 0.09 the chamber is supply-limited. Dose uniform to 1 % out to r = 100 mm, 7 % dip at
r = 142 mm (radial sweep), edge spike from the un-patterned floor ring. Gas phase clears 0.05 s after the pulse.
Sticking sweep: capture 96 -> 65 %, within-wafer spread 7 -> 0.2 %, purge tail 0.15 -> 0.58 s for s = 1 -> 1e-3.
Self-limiting wafer, s = 5e-3, 1 % TMA: 90 % coverage at 0.35 s (centre), 0.38 s (edge). Operating sweep: pulse
time follows the supply time N0/(U c0) at 1.6 to 2.3x; pressure does not change it at fixed sccm but raises the
edge lag from 5 % (100 Pa) to 14 % (1000 Pa); a 10 mm gap beats 25 mm on time and uniformity.

### 3.2 One-way coupling and the shielding effect
Driving a 100 nm hole with the blanket wafer's surface density: AR 50 bottom saturates at 0.92 s, against 4.45 s at
constant density, because the flat surface stops consuming TMA at 0.35 s and the surface density rises tenfold.
AR 100: 3 % coverage at the bottom after a 1 s pulse. Radius shifts the AR 50 time by only 10 %.

### 3.3 Two-way coupling: loading
3D-NAND-like wafer, AR 100, 1 % TMA: surface density collapses to 0.019 c0 (0.24 c0 blanket), the wafer captures
94 % of the injected TMA during a 1 s pulse, and nothing saturates; a 1 s pulse delivers 8 % of the demand.
Long pulses: hole bottoms reach 90 % at 14.05 s (centre), 13.3 s (edge); flat regions between holes at 1.95 s
(0.35 s blanket). 3 % TMA: 4.75 s.

### 3.4 The pulse-time map and the design rule
Table (worst point on the wafer): see assumptions.md part 7. Patterned wafers sit at 1.05 to 1.25x the supply
limit (1 + A_wall/A_wafer) N0/(U c0) over the whole map; blanket wafers at 1.4 to 3.5x their own limit, rising with
feed fraction as kinetics takes over. Design rule: t_pulse ~ 1.2 (1 + A_wall/A_wafer) N0/(U c0).

### 3.5 H2O half-cycle
Blanket: 1.74 s (5x TMA; s vbar/4 = 0.02 m/s < U). Patterned AR 100: 13.4 s (TMA 14.05 s), flat-to-bottom lag
2.6 s against 12 s for TMA. Cycle time ~30 s at 1 % precursors vs ~2.5 s blanket.

### 3.6 Temporal versus spatial ALD
The design rule gives the molecules required per unit wafer area per half-cycle, Phi_req = 1.2 (1 + A_wall/A_wafer) N0.
A temporal showerhead tool delivers them at U c0 over the whole wafer and pays two purges per cycle; a
separated-chemistry (spatial) tool delivers them in passes through a precursor zone at a flux k times the showerhead
flux and pays rotation time instead of purges. With k = 30, a 1 s rotation and a 25 % zone, the two architectures tie
on blanket wafers (2 s per cycle) and diverge linearly with patterned area: 6.8 vs 2 s at A_wall/A_wafer = 11
(AR 20), 29.9 vs 4 s at 56 (AR 100), 55.9 vs 8 s at 112. The crossover lies at A_wall/A_wafer of 5 to 10 for
k = 30; k is the one quantity a tool maker knows and we do not, so Figure X is a design chart in k. This is the
quantitative form of the argument behind separated-chemistry ALD architectures for 3D NAND [AMAT Olympia release;
AMAT gap-fill patent].

### 3.7 Plasma ALD: the recombination ceiling
With a wall loss probability r for radicals, the time to 90 % bottom coverage in a 100 nm AR-100 hole rises from
8.0 s (r = 0) to 8.6 s (3e-5), 10.2 s (1e-4), 16 s (3e-4) and 64 s (1e-3), and the bottom is not reached at
r >= 3e-3; the largest aspect ratio saturating within a 5 s exposure falls from 79 (r = 1e-5) to 56 (1e-3) and 29
(1e-2). For O atoms on SiO2 and TiO2 (r ~ 6e-5 [Arts 2019]) recombination is a 10 to 25 % penalty at AR 100; on
high-recombination surfaces it is a hard ceiling that no pulse length removes, and radical flux (plasma source
design) or surface chemistry become the levers. The radical sticking coefficient (0.05) is assumed.

### 3.8 Operating optimum and uncertainty
A face-centred central composite design (25 two-way runs; TMA fraction 0.5-3 %, carrier 150-600 sccm, pressure
100-1000 Pa, gap 10-25 mm) on the AR-100 wafer, fitted as a response surface in JMP and as a Gaussian-process
surrogate in MATLAB and scikit-learn (leave-one-out R2 0.999 for log pulse time, 0.988 for utilization), gives the second-order
corrections to the rule. Pulse time scales as x^-1.0, Q^-0.8 and p^-0.12; the gap is not significant. Precursor
utilization up to saturation rises from 68 % at 100 Pa to 88 % at 1000 Pa and falls from 88 % at 150 sccm to 68 %
at 600 sccm: carrier gas that is not needed to deliver precursor carries it past the wafer. The fastest point with
edge lag under 5 % (a constraint never active: lag stays within -3.5 to +1.1 %) is 3 %, 600 sccm, 1000 Pa, 10 mm at
2.38 s. At fixed precursor flow, trading carrier gas for pressure (2.5 %, 150 sccm, 630 Pa, 10 mm) cuts the pulse
from 11.6 to 9.8 s and raises utilization from 85 to 94 %. A Monte Carlo over the site density (lognormal, factor
1.25) puts a 90 % band of 8.2 to 17.2 s on the 11.9 s rule prediction at the design centre; the absolute pulse
inherits the site-density uncertainty, the ratios in this paper do not. The scikit-learn GP reproduces the MATLAB
leave-one-out R2: 0.999 / 0.746 / 0.986 against 0.999 / 0.780 / 0.988 for time / edge lag / utilization. The
cost-of-ownership front of pulse time against TMA per saturating pulse, over the factor ranges with edge lag under
5 % (Figure, doe_pareto.png): at the centre point's pulse time (11.6 s) the front uses 1.044 times the stoichiometric
minimum of TMA against 1.181 at the centre, a 12 % saving per wafer pass at no cost in time; the last 1 % of saving
(to 1.030 at 48.9 s) costs 4.2 times the pulse time.

### 3.9 Non-uniform layouts: the local wall-area ratio sets the pulse
The patterned-surface boundary condition takes hole density and aspect ratio as functions of radius (ARProfile and
holeDensityProfile tables, interpolated at each face; a face with zero hole density is a flat Langmuir surface); the
uniform case rerun with the rebuilt library is identical to the original to the 8 written digits at every face and
write time. Two layouts at 1 % TMA, 300 sccm, 300 Pa (Figure, pattern_layout.png): (a) a dense AR-100 centre, full
hole density for r < 75 mm ramping to bare wafer at 100 mm, wafer-average wall-area ratio 19.2; (b) full density
everywhere with AR 100 at the centre ramping to 60 at the edge, wafer-average ratio 41.0. The rule with the
wafer-average ratio under-predicts both: 5.1 s against the 11.85 s the dense centre needs, 10.6 s against 12.5 s for
the AR ramp. With the local ratio of the most demanding region (55.9 at the centre of both) it predicts 14.4 s and is
conservative by 15 to 21 %, against 2 % on the uniform wafer; the surplus is what the less-demanding neighbourhood
feeds back laterally. A pattern boundary produces a saturation gradient inside the array: 11.85 s at the centre and
8.75 s at r = 74.5 mm, the last face at full density, a 26 % gradient where the uniform wafer varies by 0.4 %. The
bare edge of layout (a) saturates 15 to 30 % slower than a blanket wafer (flat t90 0.40 to 0.50 s against 0.35 to
0.38 s) because the array draws down the gas: surface TMA 0.163 c0 over the bare region at 0.3 s against 0.24 c0 on
a blanket wafer. On the AR ramp the bottom time tracks the local aspect ratio (12.5 s at AR 100 to 7.6 s at AR 60)
while the flat regions are unaffected (1.95 s): a 40 % depth variation gives a 39 % spread in saturation time, so
the recipe is set by the deepest tier and shallower regions are over-dosed by their depth ratio.

### 3.10 Step coverage, purge and throughput
In fab units the rule reads as follows. Step coverage, the hole-bottom coverage over the flat top, follows the
supply integral: at 1 % TMA the bottoms of AR 20, 50 and 100 holes reach 90 % at 2.7, 6.7 and 14.1 s (worst face),
95 % at 2.8, 6.8 and 14.2 s and 99 % at 2.9, 7.0 and 14.3 s (Figure, step_coverage_vs_pulse.png; tables below).
The front is sharp, t99 - t90 is 0.25 to 0.3 s at every aspect ratio, so the rule's pulse (3.1, 7.3, 14.4 s) covers
99 % step coverage as well as 90 % (t99/rule 0.95 to 1.00). The flat top between the holes saturates at 1.8 to 2.0 s
regardless of aspect ratio. Converted to film at 0.103 nm per cycle (the PillarHall record's 51.5 nm after 500
cycles; the review value is about 0.1 nm per cycle [Puurunen 2005]), 50 cycles of the blanket recipe (0.345 s)
leave 1.4 nm of a 5.15 nm target on the top of an AR-100 array and nothing at the bottom, and a 20 % step coverage
on an AR-20 array; 50 cycles of the rule's pulse leave 5.15 nm top and bottom on all three.

| wafer (run) | wall-area ratio | rule 1.2 (1 + A_wall/A_wafer) N0/(U c0) | t90 centre / worst | t95 centre / worst | t99 centre / worst | flat top t90 / t95 / t99 | t99 / rule |
|---|---|---|---|---|---|---|---|
| AR 20 (map, 4.61 s pulse) | 11.2 | 3.08 s | 2.67 / 2.71 s | 2.78 / 2.81 s | 2.92 / 2.94 s | 1.84 / 2.14 / 2.60 s | 0.96 |
| AR 50 (map, 10.26 s) | 28.0 | 7.32 s | 6.73 / 6.73 s | 6.83 / 6.83 s | 6.94 / 6.99 s | 1.94 / 2.35 / 3.06 s | 0.95 |
| AR 100 (map, 19.68 s) | 55.9 | 14.38 s | 14.01 / 14.11 s | 14.21 / 14.21 s | 14.31 / 14.31 s | 1.96 / 2.35 / 3.04 s | 1.00 |
| AR 100 (20 s run, 0.05 s writes) | 55.9 | 14.38 s | 14.05 / 14.10 s | 14.15 / 14.20 s | 14.25 / 14.30 s | 1.95 / 2.30 / 3.05 s | 0.99 |

| pulse | AR 20: step coverage, bottom / top | AR 50: step coverage, bottom / top | AR 100: step coverage, bottom / top |
|---|---|---|---|
| 0.345 s (blanket t90) | 20 %, 0.26 / 1.29 nm | 0.5 %, 0.01 / 1.32 nm | 0 %, 0.00 / 1.44 nm |
| 1 s | 29 %, 0.95 / 3.23 nm | 0.8 %, 0.03 / 3.25 nm | 0 %, 0.00 / 3.16 nm |
| 2 s | 58 %, 2.76 / 4.79 nm | 1.7 %, 0.08 / 4.69 nm | 0 %, 0.00 / 4.67 nm |
| 5 s | 100 %, 5.15 / 5.15 nm | 25 %, 1.27 / 5.15 nm | 0.1 %, 0.00 / 5.15 nm |
| 10 s | 100 %, 5.15 / 5.15 nm | 100 %, 5.15 / 5.15 nm | 4.5 %, 0.23 / 5.15 nm |
| rule (3.08 / 7.32 / 14.38 s) | 100 % | 100 % | 100 % |

The purge has three parts and only one of them is set by the chamber. The hole itself empties by Knudsen
diffusion in microseconds: for 100 nm x 10 um and TMA at 573 K, D_K = 1.4e-5 m2/s, L^2/D_K = 7 us, and a transient
solution with the mouth density set to zero puts the bottom density at 1 % after 14 us. The chamber outlet signal
falls with a time constant of 0.36 to 0.39 s for TMA and 0.46 s for water, between the gap transit (0.3 s) and the
chamber volume over the flow (0.53 s): a decade costs about 1 s, 1 % of the pulse-end level takes 1.9 s (TMA) and
2.2 s (water), 0.1 % about 2.7 and 3.2 s. Desorbing water is the long part: at kdes = 0.1 1/s the exhaust is held
at 1.6 % of the feed level by the surface 3 s after the pulse and reaches 1 % only after about 8 s; at 1 1/s,
4.3 s; at 10 1/s the surface empties within 0.3 s and the purge is gas clearing again (2.4 s). kdes is a parameter,
not a measurement, and the water purge in the cycle below is the kdes = 0.1 1/s value by choice.

(a) Hole purge, 100 nm x 10 um, TMA at 573 K (vbar 410 m/s, D_K = d vbar/3 = 1.37e-5 m2/s), walls saturated, mouth density zero:

| quantity | value |
|---|---|
| L^2/D_K | 7.3 us |
| slowest eigenmode 4 L^2/(pi^2 D_K) (Dirichlet mouth, Neumann bottom) | 3.0 us |
| bottom density to 1 %, first mode (4/pi) exp(-t/tau1) | 14.4 us |
| numeric: bottom density to 1/e | 3.7 us |
| numeric: bottom density to 1 % | 14.4 us |
| numeric: hole-mean density to 1 % | 13.1 us |

(b) Chamber outlet signal after the pulse ends (fraction of the pulse-end outlet value; "extr." = extrapolated with the exponential tail fitted over the last 40 % of the written purge window):

| run | pulse | outlet at pulse end | to 10 % | to 1 % | to 0.1 % | tail decay rate (time constant) |
|---|---|---|---|---|---|---|
| TMA, patterned AR 100, 20 s pulse, purge to 22 s | 20 s | 1.000 c0 | 1.10 s | 1.90 s | 2.73 s (extr.) | 2.75 1/s (0.36 s) |
| TMA, blanket, 0.5 s pulse, purge to 1.5 s | 0.5 s | 0.297 c0 (peak 0.367 at +0.11 s) | 1.25 s (extr.) | 2.14 s (extr.) | 3.04 s (extr.) | 2.57 1/s (0.39 s) |
| TMA, blanket, 1.26 s pulse, purge to 2.26 s (map) | 1.26 s | 0.888 c0 | 1.07 s (extr.) | 1.95 s (extr.) | 2.83 s (extr.) | 2.63 1/s (0.38 s) |
| H2O, blanket, kdes 0, 5 s pulse, purge to 8 s | 5 s | 0.999 c0 | 1.14 s | 2.20 s | 3.25 s (extr.) | 2.18 1/s (0.46 s) |
| reference: gap transit 0.015 m / 0.050 m/s | | | | 0.30 s | | |
| reference: chamber volume / flow, pi 0.2^2 x 0.015 m3 / 3.54e-3 m3/s | | | | 0.53 s | | |

(c) Desorption, water on the blanket wafer, 5 s pulse at 1 %, outlet to 1 % of the pulse-end value:

| kdes | coverage at pulse end / 3 s after | outlet 3 s after the pulse | purge to 1 % | basis |
|---|---|---|---|---|
| 0 | 1.000 / 1.000 | 0.0017 c0 | 2.20 s | run (gas clearing) |
| 0.1 1/s | 0.951 / 0.755 | 0.0167 c0 (tail model 0.0159) | 8.26 s | extrapolated: tail kdes theta N0 A/(Q c0), theta decaying at 0.088 1/s |
| 1 1/s | 0.660 / 0.113 | 0.0254 c0 (tail model 0.0238) | 4.33 s | extrapolated, theta decaying at 0.654 1/s |
| 10 1/s | 0.163 / 0.001 | 0.0030 c0 | 2.40 s | run (surface empties within 0.3 s; gas clearing again) |

The full Al2O3 cycle, TMA pulse, purge to 1 %, water pulse, purge, is 12.2 s on the blanket wafer (0.35 + 1.9 +
1.7 + 8.3 s), 37.6 s on the AR-100 wafer (14.05 + 1.9 + 13.4 + 8.3 s) and 29.4 s at the same-supply optimum of
Section 3.8 (9.8 + 1.9 + 9.4 + 8.3 s, the water pulse scaled by the TMA ratio since the water half-cycle was not in
the design, and the purges held at the baseline values; the optimum's lower carrier flow raises the chamber
residence time 2.8-fold, so its TMA purge would scale to about 5 s and its cycle to about 33 s, a scaling estimate
rather than a run). At 0.103 nm per cycle, 5 nm is 49 cycles: 10 min per wafer and 5.5 wafers per hour with 60 s of
handling on the blanket wafer, 31 min and 1.9 wafers per hour on the AR-100 wafer, 24 min and 2.4 at the
optimum; 10 nm doubles the minutes (Figure, cycle_time_wph.png). With gas clearing only for water (kdes = 0, 2.2 s)
the three cycles are 6.2, 31.6 and 23.3 s. The rule turns this into a throughput formula, 3600 / (N [2 x 1.2
(1 + A_wall/A_wafer) N0/(U c0) + purges] + handling): below a wall-area ratio of about 10 the purges set the cycle
and precursor fraction buys little; above it throughput falls as 1/(1 + A_wall/A_wafer), and tripling the TMA
fraction doubles the AR-100 throughput (1.8 to 3.5 wafers per hour) while leaving the blanket wafer where it was
(5.9 to 6.3). The chamber run at AR 100 sits on the 1 % curve within 2 %.

Al2O3 cycle time and throughput (GPC 0.103 nm per cycle; 49 cycles for 5 nm, 98 for 10 nm; purges TMA 1.90 s, H2O 8.26 s at kdes = 0.1 1/s; in brackets the same cycle with gas clearing only for water, kdes = 0, 2.20 s):

| recipe | TMA + purge + H2O + purge | cycle | 5 nm: min per wafer, WPH, WPH with 60 s handling | 10 nm: min per wafer, WPH, WPH with 60 s handling |
|---|---|---|---|---|
| blanket wafer, baseline (1 %, 300 sccm, 300 Pa, 15 mm) | 0.345 + 1.90 + 1.74 + 8.26 s | 12.2 s (6.2 s) | 10.0 min, 6.0, 5.5 (5.1 min, 11.9, 9.9) | 20.0 min, 3.0, 2.9 (10.1 min, 5.9, 5.4) |
| patterned AR 100, baseline | 14.05 + 1.90 + 13.4 + 8.26 s | 37.6 s (31.6 s) | 30.7 min, 1.95, 1.9 (25.8 min, 2.3, 2.2) | 61.4 min, 0.98, 0.96 (51.5 min, 1.2, 1.1) |
| patterned AR 100, same-supply optimum (2.51 %, 150 sccm, 631 Pa, 10 mm) | 9.83 + 1.90 + 9.38 + 8.26 s | 29.4 s (23.3 s) | 24.0 min, 2.5, 2.4 (19.0 min, 3.2, 3.0) | 48.0 min, 1.25, 1.2 (38.1 min, 1.6, 1.5) |

The rule as throughput: WPH for 5 nm (49 cycles) against wall-area ratio, both half-cycles at the same supply flux, purges 1.90 + 8.26 s:

| A_wall/A_wafer | 1 % TMA: cycle, WPH, WPH with 60 s handling (kdes 0 purge) | 3 % TMA: cycle, WPH, WPH with 60 s handling (kdes 0 purge) |
|---|---|---|
| 1 | 11.2 s, 6.6, 5.9 (11.6) | 10.5 s, 7.0, 6.3 (13.0) |
| 11.2 (AR 20) | 16.3 s, 4.5, 4.2 (6.4) | 12.2 s, 6.0, 5.5 (10.0) |
| 28.0 (AR 50) | 25.1 s, 2.9, 2.8 (3.6) | 15.1 s, 4.9, 4.5 (7.1) |
| 55.9 (AR 100) | 39.3 s, 1.9, 1.8 (2.1) | 19.9 s, 3.7, 3.5 (4.9) |
| 112 | 67.6 s, 1.1, 1.1 (1.2) | 29.3 s, 2.5, 2.4 (3.0) |

### 3.11 Tapered holes, word-line recesses and other chemistries
Real channel holes are not cylinders. The feature solver now takes a frustum (TaperedHole, diameter falling linearly
from the mouth to a fraction of it at the bottom) and solves the conservative form (1/A) d/dx (A D dn/dx) of the
diffusion equation on the varying cross-section; with the sink switched off the density in a frustum is uniform to
1e-9, and bottom_ratio = 1 reproduces the cylinder exactly. Published numbers for the taper of production holes are
scarce: Lam reports less than 0.1 % CD deviation over 10 um with cryogenic etching, about 10 nm on a 100 nm hole,
against about 30 nm for room-temperature etches [Lam 2024; EE Times Japan 2024], and TCAD studies use 0.3 degree
tapers [Lee and Shin 2022]; we therefore take 90 % as the anchored bottom diameter and 80 and 60 % as a sweep (no
peer-reviewed vendor paper with a measured top/bottom CD pair in nanometres was found in open access). Driven by the
chamber's centre-face TMA density of the AR-100 wafer, the bottoms of 100 / 90 / 80 / 60 nm-bottom holes reach 90 %
at 14.03 / 14.21 / 14.36 / 14.62 s (Figure, tapered_hole.png; table below): the taper costs 1 to 4 %, because on a
supply-limited wafer the hole's own conductance is not the bottleneck. At a constant mouth density the same taper
costs 4 to 20 % (15.9 to 19.1 s): the lost conductance (D_K ~ d, A ~ d^2) outweighs the lost wall area. The rule with
the frustum wall area (0.95 / 0.90 / 0.80 of the cylinder's) moves the other way, to 13.7 / 13.0 / 11.6 s, and the
two bracket a tapered-hole wafer: the area reduction buys at most the rule's 5 to 20 %, the transport penalty takes
back 1 to 4 % of it, and for the anchored 90 % taper both are within 5 % of the straight hole.

| bottom diameter | wall area / straight | chamber-driven bottom t90 (one-way) | rule with frustum area, 1.2 (1 + 55.9 A_rel) N0/(U c0) | driven / rule | constant 0.3 Pa mouth, t90 |
|---|---|---|---|---|---|
| 100 nm (straight) | 1.000 | 14.03 s (chamber two-way: 14.05 s first write, 14.01 s interpolated) | 14.39 s | 0.97 | 15.93 s |
| 90 nm (literature-anchored, cryogenic etch) | 0.950 | 14.21 s | 13.68 s | 1.04 | 16.52 s |
| 80 nm | 0.900 | 14.36 s | 12.98 s | 1.11 | 17.23 s |
| 60 nm | 0.800 | 14.62 s | 11.56 s | 1.26 | 19.12 s |

Word lines of a replacement-gate 3D-NAND stack are filled through the gate-line slit into lateral cavities a few tens
of nanometres high and hundreds of nanometres to a micrometre deep. We model this as a nested feature (aldconf.recess):
a slit of width w_s and depth L_s whose walls carry recesses of height g_r, length L_r and pitch p_r, each recess a
Knudsen-Langmuir slit sub-problem fed by the local slit density, with the slit wall's uptake equal to the flat uptake
on the un-recessed fraction plus the recess mouth flux on the recessed fraction, solved quasi-steadily; it is the
patternedSurface construction one level down, and it reduces to the plain slit without recesses (to 1e-9) and to a
slit with the same total wall area when the recesses are short (within 3 % at aspect ratio 1). With literature
dimensions (29 nm word lines at 51 nm pitch [Lee and Shin 2022], 98 pairs in 5 um [Reiter et al. 2024], 250 nm slits
with recesses of 365 nm [TechInsights 2014] or 1 um [Onto Innovation 2023]) and the TMA drive of the AR-100 wafer,
the slit bottom saturates at 5.4 and 8.3 s and the last recess bottom 0.09 and 0.72 s later (2 and 9 %); a 100 nm
slit takes 13.5 and 16.5 s with lags of 0.5 and 1 % (Figure, wordline_recess.png; table below). In every case the
structure saturates within 4 % of the plain slit carrying the same total wall area (14.7 and 39.6 times the slit
wall), so the rule's area-only treatment of recesses holds up to recess aspect ratios of 35. The consequence for the
recipe is in the area: a wafer tiled with these slits at a 1 to 2.3 um pitch carries 150 to 189 wafer-areas of wall,
2.7 to 3.4 x the AR-100 hole wafer, and the rule puts its TMA pulse at 38 to 48 s (fully patterned) at 1 %. Closing
this with a recessed-slit boundary condition in the chamber is the natural next step; the nested model says the
correction it would find inside the structure is below 10 %. The 1D slit transport behind these numbers (the
hydraulic-diameter diffusivity, as in the table) was 24 % slow against the ballistic solvers (Section 2.6); with the
parallel-plate correction of Section 2.3 the chamber-driven times become 5.1 / 8.0 s (250 nm slits) and 11.0 / 15.4 s
(100 nm slits), the recess lags 0.1-0.5 s (1-7 %), the area-equivalent slit holds within 3 %, and the
constant-density w100_Lr365 structure saturates in 10.0 / 10.2 s against the voxel solver's 10.7 / 10.9 s (6 % fast,
the continuum-versus-ballistic offset of an AR-50 slit) instead of 14.5 / 14.7 s (36 % slow); the wall-area ratios
and the rule are geometry only and unchanged.

| slit width w_s | recess length L_r (AR_r) | area factor A_f | chamber-driven t90: plain slit bottom | slit bottom (flat wall) | recess bottom at mouth / mid / stack bottom | delay of the last recess bottom behind the slit bottom | area-equivalent slit bottom | constant 0.3 Pa: slit bottom / last recess bottom; exposure | ratio to the effusion-limit exposure |
|---|---|---|---|---|---|---|---|---|---|
| 250 nm | 365 nm (12.6) | 14.7 | 2.30 s | 5.37 s | 2.09 / 4.46 / 5.46 s | 0.09 s (1.8 %) | 5.39 s | 2.70 / 2.79 s; 0.84 Pa s | 4.6 |
| 250 nm | 1000 nm (34.5) | 39.6 | 2.30 s | 8.33 s | 2.94 / 6.94 / 9.05 s | 0.72 s (8.6 %) | 8.70 s | 6.02 / 7.03 s; 2.11 Pa s | 4.3 |
| 100 nm | 365 nm (12.6) | 14.7 | 3.76 s | 13.46 s | 2.11 / 8.55 / 13.53 s | 0.07 s (0.5 %) | 13.49 s | 14.53 / 14.68 s; 4.40 Pa s | 9.7 |
| 100 nm | 1000 nm (34.5) | 39.6 | 3.76 s | 16.52 s | 2.97 / 13.82 / 16.73 s | 0.20 s (1.2 %) | 16.66 s | 37.2 / 39.2 s; 11.7 Pa s | 9.6 |

TMA (s 5e-3, N0 4e18), 5 um stack, 29 nm recesses at 51 nm pitch. The effusion-limit exposure is the stoichiometric
minimum at the mouth, (A_wall/A_mouth) N0/(n0 vbar/4).

The rule carries no chemistry beyond N0 and c0, and a second chamber chemistry tests that. HfO2 from TDMAH
(M 354.7, vbar 185 m/s at 573 K) was run on the AR-100 wafer at 1 %, 300 Pa, 573 K, 300 sccm with s = 7.0e-3 [Reiter
et al. 2024, calibrated to 3D-NAND HfO2 step coverage; no direct TDMAH measurement exists], N0 = 2.77e18 Hf/m2 from a
GPC of 0.10 nm [Blaschke et al. 2020; Gieraltowska et al. 2023] and the bulk density 9.68 g/cm3 [Delabie et al.
2005], and a gas diffusivity of 3.6e-3 m2/s (hard-sphere Chapman-Enskog relative to TMA). The hole bottoms reach
90 % at 11.6 s against 14.05 s for TMA; the rule gives 9.96 s against 14.38 s (table below). TDMAH sits 17 % above
its rule where TMA sits 2 % below: 0.5 s of the 1.6 s excess is the 2.3 x lower diffusivity (gap Peclet number 0.21
instead of 0.09; a run with the TMA diffusivity gives 11.1 s), the rest the shorter supply time against the fixed
chamber overheads, the trend the DOE found for the rule constant. 573 K is at the upper edge of the TDMAH
self-limiting window [Blaschke et al. 2020; Gieraltowska et al. 2023; Hausmann et al. 2002]; it was kept so that the
flow field, U, c0 and the inlet condition are those of every other run, and the rule at 250 C is 9.3 to 10.0 s
because U c0 at fixed sccm and fraction is temperature-independent. Rule-only predictions for TiN from TiCl4 (s 0.10
[Ylilammi et al. 2018], GPC 0.30 A, density 4.23 g/cm3 [Kinnunen et al. 2017], N0 1.24e18) and W from WF6 (s 0.40
[Elam et al. 2001], GPC 2.5 A [Klaus et al. 2000], 19.25 g/cm3, N0 1.58e19) give 4.4 and 57 s per cycle on the same
wafer at 1 %: the sticking coefficients span 5e-3 to 0.4 and do not enter, since s vbar/4 exceeds the inlet speed by
6 to 400 x for all four precursors; what differs is the number of metal atoms a cycle has to deliver. Per nanometre
of film the pulse scales with the film's metal-atom number density (HfO2 100 s/nm, Al2O3 144, TiN 148, W 227), which
is the number a tool's throughput model needs. The inputs carry their caveats: the TDMAH sticking value is a model
calibration for an undisclosed Hf amide, not a measurement; the TiCl4 value is for TiCl4/H2O on TiO2 at 110 C, not
TiCl4/NH3; the WF6 value is on the Si-terminated surface, not the B-terminated one; the HfO2 and W film densities are
bulk values, not measured ALD densities; N0 therefore carries the GPC x density uncertainty of each chemistry (factor
1.1 to 1.3 across the cited values), the same dominant term the Monte Carlo on the TMA rule found. None of these
moves the patterned-wafer result, which is supply-limited, but each would move a blanket-wafer one.

Chamber runs, AR-100 wafer (wall-area ratio 55.9), 1 % precursor, 300 Pa, 573 K, 300 sccm, 20 s pulse:

| | TMA (reference) | TDMAH | TDMAH with the TMA diffusivity (sensitivity) |
|---|---|---|---|
| M, vbar at 573 K | 72.09, 410 m/s | 354.7, 185 m/s | 354.7, 185 m/s |
| s, N0 | 5e-3, 4.0e18 | 7.0e-3, 2.77e18 | 7.0e-3, 2.77e18 |
| D in N2 at 300 Pa; inlet fraction f | 8.4e-3 m2/s; 1.49e-3 | 3.62e-3 m2/s; 3.45e-3 | 8.4e-3; 1.49e-3 |
| gap Peclet number U g / D | 0.09 | 0.21 | 0.09 |
| flat t90 (centre) | 1.95 s | 1.30 s | 1.20 s |
| hole-bottom t90 centre / worst / r = 149 mm | 14.05 / 14.10 / 13.50 s | 11.60 / 11.60 / 11.10 s | 11.10 / 11.15 / 10.95 s |
| rule 1.2 (1 + 55.9) N0/(U c0) | 14.38 s | 9.96 s | 9.96 s |
| bottom t90 / rule | 0.98 | 1.17 | 1.12 |
| surface density at 5 s (centre) | 0.084 c0 | 0.197 c0 | 0.213 c0 |

Rule-only predictions, AR-100 wafer (55.9), 1 % precursor, baseline chamber (U c0 = 1.90e19 molecules/m2/s):

| film / precursor | T | s (source) | GPC (source) | film density (source) | N0 = GPC rho N_A / M | vbar | s vbar/4 vs U = 0.05 m/s | supply time N0/(U c0) | rule, AR-100 wafer | rule, blanket (x 1.2) |
|---|---|---|---|---|---|---|---|---|---|---|
| Al2O3 / TMA | 300 C | 5e-3 (chamber runs; Arts 2019 0.5-2e-3, PillarHall fit 8.5e-3) | 0.10 nm | (N0 taken as 4e18, one saturated monolayer) | 4.0e18 | 410 m/s | 0.51 | 0.211 s | 14.4 s | 0.25 s |
| HfO2 / TDMAH | 300 C | 7.0e-3 (Reiter 2024, Hf amide, calibrated to 3D-NAND HfO2 step coverage) | 0.10 nm (Blaschke 2020: 1.08 A at 200 C, 0.93 A at 250 C; Gieraltowska 2023: 1.2 A at 180-240 C) | 9.68 g/cm3 bulk monoclinic (Delabie 2005) | 2.77e18 | 185 m/s | 0.32 | 0.146 s | 9.96 s | 0.17 s |
| HfO2 / TDMAH | 250 C | 7.0e-3 | 0.093 nm (Blaschke 2020) | 9.68 g/cm3 | 2.58e18 | 177 m/s | 0.31 | 0.136 s | 9.26 s | 0.16 s |
| TiN / TiCl4 | 325 C | 0.10 (Ylilammi 2018, TiCl4/H2O on TiO2 at 110 C; no TiCl4/NH3 value published) | 0.030 nm (Kinnunen 2017, 325 C; Satta 2002: 0.024 nm at 400 C) | 4.23 g/cm3 XRR (Kinnunen 2017; bulk 5.21) | 1.24e18 | 253 m/s | 6.3 | 0.065 s | 4.44 s | 0.08 s |
| W / WF6 | 300 C | 0.40 (Elam 2001, WF6 on the Si-terminated surface at 573 K; no B-terminated value) | 0.25 nm (Klaus 2000, WF6/Si2H6) | 19.25 g/cm3 (bcc, a = 316.5 pm) | 1.58e19 | 202 m/s | 20 | 0.830 s | 56.7 s | 1.00 s |

### 3.12 A full 3D chamber with a one-sided pump
The wedge is axisymmetric by construction. To test what a real pumping arrangement does to the saturation map, the
same chamber was meshed in full 3D (blockMesh O-grid, 106,000 hexahedra, non-orthogonality below 31 deg, 12 cells
across the 15 mm gap) with the outlet reduced from the full annulus to a 60 degree sector of the side wall at
r = 200 mm, the rest of the side wall a no-slip wall; showerhead, wafer condition, operating point and TMA fraction
are those of Section 3.3 (uniform AR-100 wafer, 1 % TMA, 16 s pulse). The steady carrier flow (Re 0.045) is strongly
one-sided: the gap-averaged outward velocity at r = 140 mm is 0.85 m/s in the port-side quadrant and 0.01 m/s on the
far side, against 0.27 m/s in the axisymmetric chamber, so gas over the far half of the wafer crosses the wafer toward
the port instead of leaving radially (Figure, chamber3d_t90_map.png).

The gas phase sees the asymmetry, the saturation time hardly does. At r = 140 mm the surface TMA density is 9 %
higher on the port side than on the far side at 0.3 s and 8 % higher at 14 s, and the hole-bottom coverage at 14 s is
0.90 against 0.84; but the hole bottoms fill in a steep late wave fed by the integrated supply, so this 6 % coverage
lead is a 0.13 s lead in saturation time. The azimuthal peak-to-peak of the hole-bottom t90 is 0.09 % at r = 50 mm,
0.37 % at 100 mm and 0.96 % at 140 mm, in every case a single smooth lobe earliest at the port and latest opposite it
(Figure, chamber3d_azimuthal.png; table below). The slowest point on the wafer moves from the edge ring of the wedge
to the far side at r = 133 mm, 160 deg from the port, at 14.17 s against 14.05 s in the wedge, and 0.3 % of that
0.8 % is the coarser 3D gap resolution, which shifts the centre value by the same amount. The flat surface between
the holes, which saturates at 2 s while the gas phase is still drawn down to 0.02 c0, is the more sensitive
indicator: 1.92 s on the port side against 2.01 s on the far side at r = 140 mm, a 4.8 % spread.

| radius (mm) | faces in band | mean (s) | min (s) | max (s) | peak-to-peak, raw band (% of mean) | peak-to-peak, 15 deg bins at the exact radius (%) | port side - far side (%) | wedge t90 at r (s) | 3D mean / wedge |
|---|---|---|---|---|---|---|---|---|---|
| 50 | 140 | 14.048 | 14.042 | 14.055 | 0.09 | 0.09 | -0.08 | 14.009 | 1.0028 |
| 100 | 308 | 14.075 | 14.041 | 14.098 | 0.40 | 0.37 | -0.36 | 14.021 | 1.0038 |
| 140 | 360 | 14.089 | 13.975 | 14.156 | 1.29 | 0.96 | -0.90 | 14.006 | 1.0060 |
| 148.9 (outermost row) | 248 | 13.687 | 13.544 | 13.931 | 2.83 | - | -0.70 | 13.280 (r = 149.5 mm) | 1.031 |

The reason is the one behind the whole paper. At gap Peclet 0.09 the precursor reaches the wafer by diffusion across
the gap from the showerhead directly above each point, the patterned wafer captures almost all of what arrives, and
the saturation time is the time to deliver (1 + A_wall/A_wafer) N0 per unit area at the showerhead's own rate. The
cross-flow toward the port can redistribute only the few per cent that escape capture, and it does so over a lateral
diffusion time (R^2/D = 2.7 s) short against the 14 s pulse. A one-sided pump therefore costs about 1 % of pulse time
on this wafer; the recipe margin of the design rule (the factor 1.2) covers it, and the azimuthal lag would need a
convective gap (Peclet of order 1, i.e. ten times the flow or pressure) before it mattered. The same ordering holds in
the wedge's radial direction, where the edge lag is 0.4 %.

### 3.13 Byproducts: what a second species does and does not do
At the pressures of this work the mean free path is 12 to 123 um, so inside a 100 nm feature the Knudsen number is
120 to 1,200 and molecules collide with the walls, not with each other: a reaction byproduct cannot block the
precursor in the gas phase, and each species diffuses independently with its own Knudsen diffusivity (methane 2.1 x
faster than TMA, HCl 2.3 x faster than TiCl4). What a byproduct can do is surface-side: re-adsorb on a free site and
block it until it desorbs, or etch the film. We extended the feature solver to two species (aldconf.multispecies and
the library option): the byproduct is released at nu molecules per adsorbed precursor, solved quasi-steadily with a
source on the right-hand side, re-adsorbs on free sites with sticking s_B, desorbs at k_B and removes a coated site
with probability k_etch per collision; the precursor sees the free fraction 1 - theta_P - theta_B. With s_B = k_etch = 0,
or k_B -> infinity, the single-species solver is reproduced to round-off, and the byproduct balance (produced +
desorbed = escaped + re-adsorbed + etched) closes to 1e-14 in every run. For TMA/H2O the byproduct is methane
(3 per TMA over a cycle [Puurunen 2005]), which does not re-adsorb at these temperatures; the two-species result is
the single-species one, and even a hypothetical re-adsorption at s_B = 1e-3 with 1 ms residence changes the AR-100
pulse by 0.01 % while 99.95 % of the methane leaves the hole (Figure, multispecies.png, panel a; table (a)). Every
result of this paper therefore stands as computed.

For TiCl4/NH3 the byproduct is HCl, 4 per TiCl4 over a cycle [Ritala et al. 1995], about 0.5 of them in the TiCl4
half-step [Tiznado and Zaera 2006], and its re-adsorption is documented on TiN, silica and TiO2 [Tiznado and Zaera
2006; Snyder et al. 2006; Leem et al. 2014]; no sticking coefficient or desorption rate exists for HCl on TiN at
573 K, so we ran the range that the HCl/alumina surface science spans (s_B 1e-3 to 0.4 [Elam et al. 2000; Lostier et
al. 2026]; k_B 0.1 to 1,000 1/s, the 300-650 K TPD span of [Elam et al. 2000] read at 573 K). The desorption rate
decides (panel b, table (b)): at k_B >= 100 1/s HCl is a spectator whatever its sticking; at k_B = 1 1/s the anchored
half-step stoichiometry lengthens the 100 nm x 10 um hole's saturating pulse by 11 to 55 % (7.75 s to 8.6-12.0 s at a
constant 0.3 Pa mouth) and by 11 to 36 % under the chamber's density history (table (c)); at k_B = 0.1 1/s the
blocked sites outlive the pulse and the hole fills three times slower. Releasing the whole cycle's four HCl in the
half-step (a bound) doubles these costs, and an etch probability of 1e-3 per HCl collision (a bound; HCl does not
etch TiO2 [Leem et al. 2014], and the TiN kinetic Monte Carlo of [Kim et al. 2022] lists Ti removal by surface Cl
without a rate) costs 70 to 80 % at nu 0.5 and prevents saturation at nu 4: the hole settles into an etch-limited
steady state. HCl accumulated in the chamber at 10 % of the precursor density adds 4 to 5 % at k_B = 1 1/s and
nothing at 100 1/s, because the HCl produced inside the feature outweighs what the chamber brings in (panel d; table
(d)). The profiles at the rule pulse (panel c) show the mechanism: blocked sites near the mouth let the precursor
pass, so the front is broader and deeper at a given time, but the sites behind it must be released and refilled, and
that is where the time goes.

What this means for a chloride recipe is a measurement, not a number: the supply rule (Section 3.4) is untouched
when k_B >= 100 1/s, carries a 10 to 50 % chemistry margin at k_B ~ 1 1/s, and fails at k_B <= 0.1 1/s or with
etching, so a TiN or HfO2 tool model needs the HCl desorption rate on the growing film at the deposition
temperature (one TPD of HCl from ALD TiN would fix it) before the wall-area rule can be applied with the same
confidence as for Al2O3. The library carries the option behind defaults that reproduce the single-species boundary
conditions bit for bit; a transported HCl field in the chamber is the next step.

(a) TMA/H2O with methane (methane does not re-adsorb; AR-100 hole, bottom t90):

| drive | single-species | two-species, CH4 nu 3, s_B 0 | bound: s_B 1e-3, k_B 1e3 1/s | methane escaping the hole (bound) |
|---|---|---|---|---|
| constant 0.3 Pa mouth | 15.600 s | 15.600 s, max coverage difference 0 | 15.601 s (+0.01 %) | 99.95 % |
| chamber centre-face history | 14.005 s | 14.005 s, max coverage difference 0 | 14.006 s (+0.01 %) | 99.95 % |

(b) TiCl4/NH3 with HCl at a constant 0.3 Pa mouth density; bottom t90 without HCl action 7.754 s. Blocking only (nu, s_B, k_B):

| nu | s_B | k_B 0.1 1/s | k_B 1 | k_B 10 | k_B 100 | k_B 1,000 |
|---|---|---|---|---|---|---|
| 0.5 (anchored) | 1e-3 | 23.51 s (+203 %) | 8.59 (+11 %) | 7.75 | 7.75 | 7.75 |
| 0.5 | 1e-2 | 24.34 (+214 %) | 8.65 (+12 %) | 7.75 | 7.75 | 7.75 |
| 0.5 | 0.1 | 24.51 (+216 %) | 10.06 (+30 %) | 7.93 (+2 %) | 7.75 | 7.75 |
| 0.5 | 0.4 | 27.67 (+257 %) | 12.00 (+55 %) | 8.51 (+10 %) | 7.79 (+0.5 %) | 7.75 |
| 4 (bound) | 1e-3 | 24.73 | 8.41 (+8 %) | 7.75 | 7.75 | 7.75 |
| 4 | 1e-2 | 27.32 | 8.97 (+16 %) | 7.78 | 7.75 | 7.75 |
| 4 | 0.1 | not within 30 s | 13.36 (+72 %) | 8.42 (+9 %) | 7.77 | 7.75 |
| 4 | 0.4 | not within 30 s | 20.39 (+163 %) | 10.16 (+31 %) | 7.97 (+3 %) | 7.75 |

With etching (bound; k_etch per HCl collision on a Ti site), bottom t90 and the state at 30 s:

| nu | s_B | k_B 1/s | k_etch 1e-4 | k_etch 1e-3 | blocking only |
|---|---|---|---|---|---|
| 0.5 | 1e-2 | 1 | 10.44 s | 13.98 s | 8.65 s |
| 0.5 | 1e-2 | 100 | 9.32 s | 12.46 s | 7.75 s |
| 0.5 | 0.1 | 1 | 11.44 s | 14.17 s | 10.06 s |
| 0.5 | 0.1 | 100 | 9.74 s | 13.11 s | 7.75 s |
| 4 | 1e-2 | 1 | never: bottom theta_P 0.045, blocked 0.94 at 30 s | never: 0.001 / 0.98 | 8.97 s |
| 4 | 1e-2 | 100 | never: 0.000 / 0.37 | never: 0.000 / 0.37 | 7.75 s |
| 4 | 0.1 | 1 | never: 0.161 / 0.84 | never: 0.011 / 0.99 | 13.36 s |
| 4 | 0.1 | 100 | never: 0.000 / 0.85 | never: 0.000 / 0.85 | 7.77 s |

(c) The saturating pulse chamber-driven (the TMA run's centre-face history as a stand-in; no TiCl4 chamber run
exists, so only the changes are meaningful) and at a constant 0.3 Pa mouth; bottom t90 without HCl action 9.622 s
driven, 7.754 s constant:

| nu | s_B | k_B 1/s | driven bottom t90 | change | bottom theta_P at the TiCl4 rule pulse 4.44 s | at the TMA rule pulse 14.4 s | constant-mouth t90 | change |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 1e-2 | 1 | 10.72 s | +11 % | 0.000 | 0.997 | 8.65 s | +12 % |
| 0.5 | 1e-2 | 100 | 9.62 s | 0 | 0.000 | 1.000 | 7.75 s | 0 |
| 0.5 | 0.1 | 1 | 11.76 s | +22 % | 0.000 | 0.991 | 10.06 s | +30 % |
| 0.5 | 0.1 | 100 | 9.62 s | 0 | 0.000 | 1.000 | 7.75 s | 0 |
| 0.5 | 0.4 | 1 | 13.11 s | +36 % | 0.000 | 0.959 | 12.00 s | +55 % |
| 0.5 | 0.4 | 100 | 9.67 s | +0.5 % | 0.000 | 1.000 | 7.79 s | +0.5 % |
| 4 | 1e-2 | 1 | 10.89 s | +13 % | 0.012 | 0.997 | 8.97 s | +16 % |
| 4 | 1e-2 | 100 | 9.62 s | 0 | 0.000 | 1.000 | 7.78 s | +0.4 % |
| 4 | 0.1 | 1 | 14.06 s | +46 % | 0.058 | 0.922 | 13.36 s | +72 % |
| 4 | 0.1 | 100 | 9.65 s | +0.3 % | 0.000 | 1.000 | 7.77 s | +0.3 % |
| 4 | 0.4 | 1 | 16.72 s | +74 % | 0.118 | 0.698 | 20.39 s | +163 % |
| 4 | 0.4 | 100 | 9.84 s | +2.3 % | 0.000 | 1.000 | 7.97 s | +3 % |
| 0.5, k_etch 1e-4 | 0.1 | 1 | 12.51 s | +30 % | - | - | 11.44 s | +48 % |
| 0.5, k_etch 1e-3 | 0.1 | 1 | 13.74 s | +43 % | - | - | 14.17 s | +83 % |

(d) Byproduct at the mouth (chamber accumulation, a parameter): bottom t90 for a byproduct mouth density of 0 / 1 / 3 / 10 % of the precursor mouth density:

| case | 0 % | 1 % | 3 % | 10 % |
|---|---|---|---|---|
| TiCl4 constant 0.3 Pa, nu 0.5, s_B 0.1, k_B 1 | 10.059 s | 10.111 (+0.5 %) | 10.215 (+1.6 %) | 10.568 (+5.1 %) |
| TiCl4 constant, nu 0.5, s_B 0.1, k_B 100 | 7.749 s | 7.749 | 7.749 | 7.749 |
| TiCl4 chamber-driven, nu 0.5, s_B 0.1, k_B 1 | 11.759 s | 11.808 (+0.4 %) | 11.905 (+1.2 %) | 12.234 (+4.0 %) |
| TMA constant, CH4 bound s_B 1e-3, k_B 1e3 | 15.601 s | 15.601 | 15.601 | 15.601 |

### 3.14 Gate-all-around stacks and a current NAND cell
The voxel solver of Section 2.6 takes any shape, so two device geometries that the one-dimensional models cannot
describe were run with it: the replacement-gate step of a stacked-nanosheet transistor and a 3D-NAND word-line cell
with the channel-hole pillars standing in the recesses (Figure, device_geometries.png). The nanosheet cell follows
the IBM device of Loubet et al. [2017] (three 5 nm sheets at 10 nm spacing, 12 nm gate, 40 nm sheets [Semiconductor
Engineering 2018]) under a 40 nm gate trench, periodic across the sheet width at a 60 nm pitch, with a four-sheet
7 nm-spacing variant [Wong and Kakushima 2022] and the IBM cell after 3 nm of oxide on every surface (4 nm slots,
6 nm trench) for the TiN step; the NAND cells are the word-line slit of Section 3.11 with two rows of 100 nm pillars
per side, and a 2XX-layer cell built from a teardown table [TechInsights 2022] (12.7 um slit depth, 48 nm word-line
pitch) and a metrology note [Onto Innovation 2023] (250 nm slit, 1 um recesses, 100 nm holes). Each was run with the
TMA, TDMAH and TiCl4 parameters at a constant 0.3 Pa and under the centre-face history of the AR-100 wafer, with 8
voxels across every slot or recess, particle doublings and grid changes moving the saturation times by at most 1 %
(nanosheets) and 3 % (NAND).

The nanosheet stack is, for ALD, a flat surface with a short slot: the last point, the centre of the slot under the
lowest sheet, reaches 90 % at 0.56 s of TMA at 0.3 Pa against 0.47 s for a flat surface, and 99 % at 1.04 against
0.95 s; the slots fill from their two open edges inward with a centre-to-edge lag of 1 % (TMA) to 8 % (TiCl4), 28 % in
the 4 nm slots (Figure, gaa_coverage.png; table below). A wafer tiled with such devices carries 5.7 wafer areas of
surface (the gate trench walls, not the sheets, dominate), so the rule gives a 1.7 s TMA pulse at 1 %, by which time
the cell is at 99.5 % everywhere; the exposure for 90 % everywhere is 0.17 Pa s (TMA), 0.20 (TDMAH) and 0.017 Pa s
(TiCl4), 30 x and 6 x the stoichiometric minimum, i.e. reaction-limited like a blanket wafer. The word-line cells are
the opposite case (second table below). The 5 um cell with pillars saturates its last recess at 11.3 s at 0.3 Pa
(3.4 Pa s), 6 % sooner than without pillars because a 100 nm pillar in a 30 nm recess removes more plate area than it
adds in wall; the 2XX cell needs 33 s (9.9 Pa s of TMA, 15 Pa s of TDMAH, 9 Pa s of TiCl4), eight times its effusion
limit, and a wafer of it carries 448 wafer areas, so the rule puts its TMA pulse at 113 s at 1 %, 3.4 x the feature's
own need; for TiN the margin is only 17 %. The nested 1D model of Section 3.11 with the parallel-plate slit
diffusivity reproduces the deepest recess within 2-6 % for TMA and TDMAH (10.2 vs 10.4 s; 29.9 vs 28.2 s) and is
18-29 % off for TiCl4, where the recess rather than the slit limits. What the 1D model cannot give is where the last
point sits: in every voxel run it is at 87-90 % of the depth, about five slit widths above the floor, and the deepest
recess is 9-22 % ahead of it, because the saturated floor returns the downward stream into the lowest recesses; the
same turn-up is present at the bottom of the AR-100 hole of Section 2.6 and absent from the continuum profile.

These are applications of a solver validated in Section 2.6, not validations: no published measurement of film
thickness per sheet or per word line was found for either structure (the one step-coverage datum on a word-line
stack, Reiter et al.'s 85 % HfO2, fixed the TDMAH sticking coefficient in Section 3.11 and cannot also test the
transport), and at saturating pulses both solvers predict uniform films, which is what the device papers'
cross-sections show and what any model would. The predictions that would test the transport are the 3-8 %
(TMA-like) to 30 % (TiCl4, 4 nm slot) deficit at the slot centre of the lowest sheet after a pulse one third of
saturation, and the 10-20 % dose-equivalent lead of the lowest word lines over those five slit widths above them
after a two-thirds pulse; a TEM or EELS thickness at those positions after a deliberately short pulse is the
experiment.

GAA cells (573 K; constant 0.3 Pa mouth density or the AR-100 wafer's centre-face history; M 50 000 particles per
step, dt = t_flat/100; the last point is the centre of the slot under the lowest sheet in every run; times in s to
90 / 99 % area-weighted coverage; the rule's pulse is 1.2 (1 + A_wall/A_wafer) N0/(U c0) at 1 % precursor):

| cell | chemistry, drive | top of top sheet t90 / t99 | everywhere (last point) t90 / t99 | exposure 90 / 99 % (Pa s) | wall-area ratio cell / wafer | rule pulse |
|---|---|---|---|---|---|---|
| IBM 2017 (3 x 5 nm sheets, 10 nm slots) | TMA, 0.3 Pa | 0.490 / 0.966 | 0.558 / 1.038 | 0.167 / 0.311 | 17.5 / 5.69 | 1.69 s |
| | TMA, chamber | 1.974 / 3.050 | 2.152 / 3.184 | - | 17.5 / 5.69 | 1.69 s |
| | TDMAH, 0.3 Pa | 0.546 / 1.068 | 0.649 / 1.179 | 0.195 / 0.354 | 17.5 / 5.69 | 1.17 s |
| | TDMAH, chamber | 2.121 / 3.238 | 2.376 / 3.430 | - | 17.5 / 5.69 | 1.17 s |
| | TiCl4, 0.3 Pa | 0.018 / 0.032 | 0.056 / 0.070 | 0.017 / 0.021 | 17.5 / 5.69 | 0.52 s |
| | TiCl4, chamber | 0.130 / 0.208 | 0.341 / 0.413 | - | 17.5 / 5.69 | 0.52 s |
| four sheets, 7 nm slots | TMA, 0.3 Pa | 0.494 / 0.970 | 0.569 / 1.052 | 0.171 / 0.316 | 19.5 / 6.15 | 1.81 s |
| | TMA, chamber | 1.982 / 3.056 | 2.180 / 3.204 | - | 19.5 / 6.15 | 1.81 s |
| | TDMAH, 0.3 Pa | 0.550 / 1.074 | 0.665 / 1.199 | 0.200 / 0.360 | 19.5 / 6.15 | 1.25 s |
| | TDMAH, chamber | 2.133 / 3.247 | 2.419 / 3.462 | - | 19.5 / 6.15 | 1.25 s |
| | TiCl4, 0.3 Pa | 0.018 / 0.033 | 0.062 / 0.076 | 0.019 / 0.023 | 19.5 / 6.15 | 0.56 s |
| | TiCl4, chamber | 0.133 / 0.213 | 0.374 / 0.447 | - | 19.5 / 6.15 | 0.56 s |
| IBM 2017 after 3 nm oxide (4 nm slots, 6 nm trench) | TiCl4, 0.3 Pa | 0.024 / 0.038 | 0.142 / 0.157 | 0.043 / 0.047 | 25.6 / 4.35 | 0.42 s |
| | TiCl4, chamber | 0.165 / 0.246 | 0.761 / 0.827 | - | 25.6 / 4.35 | 0.42 s |
| | TDMAH, 0.3 Pa | 0.576 / 1.105 | 0.898 / 1.450 | 0.269 / 0.435 | 25.6 / 4.35 | 0.94 s |
| | TDMAH, chamber | 2.198 / 3.301 | 2.913 / 3.852 | - | 25.6 / 4.35 | 0.94 s |

Word-line cells with channel-hole pillars (regions: slit floor; deepest recess, deep half; slowest recess = the recess
whose deep half is the last to pass 90 %, which is "everywhere"; for the 2XX cell its index of 264 and depth fraction
are given; times in s to 90 / 99 %):

| cell | chemistry, drive | slit floor | deepest recess | slowest recess (everywhere) | exposure 90 / 99 % (Pa s) | wall-area ratio cell / wafer; rule pulse at 1 % |
|---|---|---|---|---|---|---|
| 5 um, 100 nm slit, 29 nm recesses at 51 nm, 365 nm long, two pillar rows | TMA, 0.3 Pa | 9.04 / 11.39 | 10.33 / 11.78 | 11.30 / 12.13 | 3.39 / 3.64 | 1422 / 171; 43.5 s |
| | TMA, chamber | 10.45 / 11.88 | 11.25 / 12.10 | 11.85 / 12.32 | - | same |
| | TDMAH, 0.3 Pa | 12.42 / 16.59 | 15.29 / 17.49 | 17.10 / 18.08 | 5.13 / 5.42 | 1422 / 171; 30.1 s |
| | TDMAH, chamber | 12.47 / 14.34 | 13.87 / 14.51 | 14.45 / 14.63 | - | same |
| | TiCl4, 0.3 Pa | 0.82 / 1.49 | 4.38 / 4.95 | 6.80 / 6.94 | 2.04 / 2.08 | 1422 / 171; 13.5 s |
| | TiCl4, chamber | 2.62 / 3.77 | 6.95 / 7.45 | 8.88 / 8.99 | - | same |
| 2XX-layer: 250 nm slit, 12.7 um, 48 nm pitch, 1 um recesses, four pillar rows | TMA, 0.3 Pa | 14.60 / 23.63 | 28.18 / 32.26 | 236 (0.90): 33.09 / 34.20 | 9.93 / 10.26 | 4009 / 448; 113 s |
| | TMA, chamber | 13.55 / 15.30 | 15.78 / 16.47 | 228 (0.87): 16.89 / 17.19 | - | same |
| | TDMAH, 0.3 Pa | 17.78 / 30.22 | 41.31 / 47.20 | 233 (0.88): 50.40 / 51.60 | 15.12 / 15.50 | 4009 / 448; 78.5 s |
| | TDMAH, chamber | 14.55 / 16.08 | 17.82 / 18.57 | 228 (0.87): 19.04 / 19.24 | - | same |
| | TiCl4, 0.3 Pa | 1.25 / 2.35 | 23.15 / 24.68 | 231 (0.88): 29.88 / 30.37 | 8.96 / 9.11 | 4009 / 448; 35.1 s |
| | TiCl4, chamber | 2.99 / 4.43 | 20.25 / - | not reached by the end of the 20 s pulse: recesses from index 172 (65 % depth) down are below 90 % except the two next to the floor; 13.8 % of the area below 90 % at 22 s | - | same |

## 4. Discussion
- Throughput on patterned wafers is bought with precursor delivery rate (vapour pressure, bubbler carrier flow,
  injector design), not with chamber pressure, gap or carrier flow; those knobs move uniformity.
- Recipe transfer: the blanket-to-product factor is the wall-area ratio, 11 to 56 here; monitor-wafer recipes
  under-dose product wafers by that factor at the same supply.
- Precursor utilization: 60 % (TMA) and 50 % (H2O) over the saturating pulse on the patterned wafer; the remainder
  is pumped after saturation. Shorter pulses with feedback (end-point on surface density) would recover it [TODO
  estimate].
- The TMA sticking coefficient is the least certain input of the chamber model: the published values at or near this
  temperature span 0.5-2e-3 (Arts et al. 2019, lateral-high-aspect-ratio front slope at a 275 C set-point), 4 +/- 1e-3
  (Vandalon and Kessels 2016, sum-frequency generation at 300 C), 5.72e-3 (Ylilammi et al. 2018, fit to a lateral
  channel profile at 300 C) and 8.5e-3 to 1.96e-2 (this work, fit to the Zenodo 3974438 profiles at 300 C, factor
  1.5 sheet-to-sheet; the lower value with the hydraulic-diameter slit diffusivity of the literature extractions, the
  upper with the parallel-plate diffusivity that matches the ballistic solvers, Section 2.3). We reran the baseline operating point (1 % TMA, 300 sccm, 300 Pa, 573 K) at s = 5e-4, 1e-3, 2e-3,
  5e-3 and 8.5e-3 for the blanket wafer and the AR-100 patterned wafer (Figure, sticking_sensitivity.png). On the
  patterned wafer the hole bottoms reach 90 % coverage at 14.00-14.15 s (centre) and 14.05-14.20 s (worst face) for
  every value, 1.17-1.18 times the supply-limit time of 11.99 s, with 84-85 % precursor utilization: a seventeen-fold
  change in s moves the saturating pulse by three write intervals. Low sticking changes the order of events inside
  the hole (at 10 s the bottom coverage is 0.38 at s = 5e-4 and 0.02 at 8.5e-3) but not the finish, because the last
  tenth of the bottom coverage waits for molecules that the showerhead has not yet delivered. The design rule, the
  blanket-to-product transfer factor and the pulse-time map are therefore independent of which published value is
  right.
- The blanket wafer is not: low sticking is a conformality asset and a throughput cost only on blanket wafers. Its
  saturation time runs from 0.33 s at s = 8.5e-3 to 1.01 s at 5e-4, about s^-0.4 rather than 1/s, because the TMA
  density at the wafer rises as the uptake slows (0.065 c0 against 0.23 c0 at the centre 0.1 s into the pulse) and
  partly compensates; the flat surface between the holes of the patterned wafer behaves the same way (1.55 to
  6.00 s). Every blanket-wafer time in this paper therefore carries a factor-3 uncertainty from the sticking
  coefficient, whereas the patterned-wafer times carry only the site-density uncertainty of the supply limit. The
  8.5e-3 of our PillarHall fit is a slope-constrained value: the front slope fixes s and the penetration depth fixes
  the exposure per site, so the unknown TMA partial pressure trades against the assumed site density, not against s.
  It sits a factor 1.5 to 2 above the two other 300 C values and a factor 4 above the upper bound of the 275 C
  set-point data, a spread that also reflects three different kinetic definitions (initial sticking probability,
  lumped Langmuir coefficient with an equilibrium constant, irreversible first-order Langmuir on 4e18 sites).
  Correcting the slit transport does not close that gap: the 1D slit diffusivity was 24 % low at AR 50 and a factor
  2.3 low at the channel's AR 500, and with the exact ballistic front width the same measured profile needs a larger
  s, not a smaller one. A profile grown at a known TMA partial pressure at 300 C, with the site density taken from the measured growth per
  cycle, would remove the trade-off and fix s to the sheet-to-sheet scatter.
- Byproducts: at Knudsen numbers of 100 to 1,000 a byproduct cannot block the precursor in the gas, only on the
  surface. For TMA/H2O (methane) the single-species model is exact to the purposes of this paper; for chlorides the
  HCl desorption rate at the deposition temperature decides between no effect (k_B >= 100 1/s), a 10 to 50 % longer
  pulse (k_B ~ 1 1/s) and failure of the supply rule (k_B <= 0.1 1/s or etching), and that rate is not in the
  literature for TiN or HfO2 at 573 K; the two-species solver and the library option (aldkin 1.1) are ready for it.
- Robustness checks: replacing the bare floor by a 25 mm edge ring with pumping below the wafer plane changes no
  saturation time by more than one write interval and leaves the within-wafer spread within 3 % (blanket 9-11 %,
  patterned 4 %); the pulse-time results do not depend on the chamber outboard of the wafer at this Peclet number.
- Robustness checks continued: 200, 300 and 400 C at fixed sccm and feed fraction give the same pulse times within one
  write interval (the supply U c0 is temperature-independent; temperature acts through the chemistry, held fixed
  here); 50, 100 and 200 nm holes at fixed open area and aspect ratio saturate within 1 % of each other, so the
  wall-area ratio is the rule's only geometric input.
- Limitations: isothermal gas (temperature enters only through properties); frozen carrier flow; irreversible Langmuir in the main runs; a first-order desorption term was added for water and shows the
  purge-tail mechanism (desorbing water holds the exhaust above 1 % for seconds at kdes = 0.1 to 1 1/s) and the
  equilibrium-coverage limit (above ~0.3 1/s at 1 % H2O the surface never reaches 90 %), but kdes is a parameter,
  not a measurement, and a two-site chemisorbed/physisorbed model is the next step; single pattern density (14 % open area, half the wafer patterned); Knudsen D_K = d vbar/3 without
  Clausing correction; s and N0 are literature-range inputs (the sensitivity to s is quantified above; the
  patterned-wafer times carry only the site-density uncertainty).

## 5. Conclusions

1. On a patterned wafer, the ALD pulse that saturates the deepest feature is set by precursor supply:
   t_pulse = 1.2 (1 + A_wall/A_wafer) N0/(U c0), within 25 % over aspect ratios 20-100 and a twenty-fold range of
   feed fraction, and independent of hole diameter at fixed wall area, of temperature at fixed sccm and fraction,
   of the chamber geometry outboard of the wafer, of the TMA sticking coefficient across its published range
   (5e-4 to 8.5e-3: 1.17 to 1.18 supply-limit times, three write intervals of movement) and, in a full 3D chamber
   with a one-sided pump, of the pumping arrangement (azimuthal spread of the saturation time under 1 %, within the
   rule's factor 1.2).
2. The same normalisation puts a measured fluidized-bed ALD run at 1.12 supply-limit times and orders the measured
   utilization of a rotary-drum reactor with loading; wafers and powders sit on one curve of utilization against
   surface area spanning five orders of magnitude.
3. The chamber model reproduces the published Fiji F200 case where that case can be represented: the TMA exposure
   delivered to the wafer per pulse agrees with the authors' three-dimensional CFD to 12 to 22 % and the saturation
   state of the 25 and 60 ms pulses with the measured growth per cycle; the 2.4 % measured non-uniformity, which
   comes from an off-axis pulse recirculation, a cold loading door and softer-than-Langmuir kinetics, is not
   reproduced by an axisymmetric frozen-flow model. The blanket saturation exposure, 0.146 Pa s at the surface, sits
   inside the 0.05 to 1.3 Pa s spread of kinetic saturation exposures from seven independent sources.
4. A recipe saturated on a blanket wafer under-doses a product wafer by the wall-area ratio, 11 to 56 here. Blanket
   wafers are the kinetics-limited exception to the rule, not the norm.
5. On a non-uniform layout the pulse is set by the local wall-area ratio of the most demanding region, not by the
   wafer average: the average under-predicts by 2.3x (dense centre) and 1.2x (aspect-ratio ramp), the local ratio
   is conservative by 15 to 21 %, and a pattern boundary puts a 26 % saturation gradient inside the array.
6. Real feature shapes and other chemistries stay inside the rule: tapered channel holes cost 1 to 4 % on the
   supply-limited wafer (the anchored 90 % taper is within 5 % of the straight hole either way); word-line recesses
   fed through a slit saturate within 4 % of a plain slit with the same wall area up to recess aspect ratio 35, so
   their cost is in the wall area alone (150 to 189 wafer-areas, 38 to 48 s at 1 %); and the rule is
   chemistry-agnostic through s, N0 and c0, with the HfO2/TDMAH chamber run at 1.17 times its rule and TiN and W
   predicted at 4.4 and 57 s per cycle on the same wafer.
7. A 3D voxel Monte Carlo feature solver agrees with converged Machball within 2 % and places the 1D hole model 3 to
   11 % fast between AR 100 and AR 20; the 1D slit model with the hydraulic-diameter diffusivity is 24 % slow, and
   with the parallel-plate correction the word-line recess numbers come within 6 % of the ballistic solver. A second species cannot block the precursor in the gas at Knudsen numbers of 100 to
   1,000: for TMA/H2O the methane byproduct leaves every result as computed, while for chlorides the HCl desorption
   rate decides between no effect (k_B >= 100 1/s), a 10 to 50 % longer pulse (k_B ~ 1 1/s) and failure of the rule
   (k_B <= 0.1 1/s or etching).
8. Two device geometries run in the voxel solver: gate-all-around nanosheet cells are reaction-limited (the last
   point, the slot centre under the lowest sheet, 18 % behind a flat surface at 90 % for TMA and 5x for TiCl4) and a
   wafer of them carries 4 to 6 times the flat area, so the rule's pulse is 1.7 s at 1 % TMA; word-line stacks with
   channel-hole pillars carry 171 (5 um cell) to 448 (2XX cell) wafer areas, their slowest recess sits at 87 to 90 %
   of the depth, about five slit widths above the floor, and on the 2XX cell at 1 % the TiN rule pulse exceeds the
   feature's own need by only 17 %.
9. Second-order knobs: pressure up and carrier flow down raise precursor utilization from about 55 % to 95 % at
   fixed precursor flow; gap and pressure set within-wafer edge lag (5 to 18 %) but not time.
10. The water half-cycle needs the same pulse as TMA on patterned wafers despite a 50x lower sticking coefficient,
   and coats deep features more evenly; reversible adsorption, not chamber flushing, governs its purge.
11. In fab units the rule's pulse gives 99 % step coverage as well as 90 %, because the saturation front is 0.25 to
   0.3 s wide at every aspect ratio; 50 cycles of the blanket recipe leave 1.4 nm of a 5.15 nm film on top of an
   AR-100 array and nothing at the bottom. The Al2O3 cycle on the AR-100 wafer is 3.1 times the blanket cycle (37.6
   against 12.2 s; 1.9 against 5.5 wafers per hour for 5 nm), its pulses are three quarters of the cycle, and above
   a wall-area ratio of about 10 throughput falls as 1/(1 + A_wall/A_wafer); the purge is set by water desorption,
   not by the chamber, which clears to 1 % in 2 s while a 10 um hole empties in microseconds.
12. Plasma ALD adds a hard ceiling: above a wall recombination probability of about 3e-3, aspect ratio 100 is
   unreachable by dose. Separated-chemistry (spatial) architectures beat temporal pulsing by the same wall-area
   factor the rule predicts.
13. All of this runs on a laptop in open-source tools, with the feature model verified against two independent
   ballistic descriptions and fitted to public measurements, and every number traceable to a script in the repository.

## Acknowledgements
Validation data by R. L. Puurunen and O. M. E. Ylivaara (Aalto University, Zenodo 3974438, CC BY 4.0). Machball by
A. Yanguas-Gil and J. W. Elam (Argonne). OpenFOAM v2606 (OpenCFD) via the gerlero/openfoam-app build.

## Code and data availability
OpenFOAM v2606 cases, the `aldWallKinetics` boundary-condition library, the Python `aldconf` solver and its command
line (`python -m aldconf` evaluates the rule, the surrogate and the Pareto front), tests and all
figures: https://github.com/Soham-ChemE/ald-conformality. The repository also ships a layout-to-recipe tool on its
live page (pattern zones in, per-zone local rule and recipe pulse out, calibrated on the two non-uniform runs of
Section 3.9) and continuous integration that rebuilds the boundary conditions and reruns a regression case in an
OpenFOAM container (129 tests). The kinetics behind the OpenFOAM boundary conditions live in `lib/aldkin`, a C-ABI
library with Fluent and plain-C hooks, and the voxel Monte Carlo feature solver of Section 2.6 in `lib/aldmc`.
Validation data: Zenodo 3974438 (CC BY 4.0).

## References (to format)
Gordon, Hausmann, Kim, Shepard, Chem. Vap. Deposition 9, 73 (2003). Yanguas-Gil and Elam, Chem. Vap. Deposition
18, 46 (2012); Machball, machball.readthedocs.io. Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018).
Arts et al., J. Vac. Sci. Technol. A 37, 030908 (2019). Arts et al., J. Phys. Chem. C (2019); (2021). Cremers,
Puurunen, Knoops, Appl. Phys. Rev. 6, 021302 (2019). Yim et al., PCCP (2020), doi 10.1039/D0CP03358H; Zenodo
3974438. Jarvilehto et al., PCCP (2023). Applied Materials, "Semiconductor equipment and processes need digital
twins" (blog); Applied Materials, Olympia ALD release (2015); Centris Spectral SiN ALD release. OpenFOAM v2606;
ViennaPS 4.7.0 (TU Wien) [used for feature topography, TODO].
