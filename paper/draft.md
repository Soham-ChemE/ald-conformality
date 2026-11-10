# Every blanket-wafer ALD recipe hides a 40x under-dose: a chamber-to-feature model of supply-limited atomic layer deposition

Soham Kavathekar, Department of Chemical and Biomolecular Engineering, University of Pennsylvania

*Draft 0.9.2 (public), 2026-11-10. Every number is from a run documented in `paper/assumptions.md`.*

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
times less reactive, needs the same pulse on the product wafer and coats it more evenly. The feature model agrees
with an independent Argonne code to 7 % and fits published thickness profiles from lateral test structures,
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

Contributions: (i) a two-way coupled chamber-feature wall condition for OpenFOAM, released; (ii) two-code
verification and public-data validation of the feature model; (iii) the supply-limit rule and its range of validity;
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
- Two-code check: Machball (ballistic, cosine re-emission) vs this Knudsen-continuum solver, 95 % bottom coverage
  time at 0.3 Pa TMA, 573 K: 1.20 vs 1.29 s (AR 20), 4.45 vs 4.46 s (AR 50), 16.2 vs 15.6 s (AR 100).
- Embedded-BC check: with negligible hole density the boundary condition's hole-bottom coverage matches the
  standalone solver driven by the same surface density to 1 % at all times.
- Validation against the Zenodo 3974438 PillarHall Al2O3 profiles (TMA/H2O, 300 C, 500 cycles x 0.1 s, 500 nm
  lateral channel, eight profiles averaged): the slit model with channel narrowing fits the saturation front to an rms
  of 0.038 of the plateau with s_TMA = 8.5e-3 (profile-to-profile factor 1.5) and a mouth exposure of 7.5 Pa s per
  pulse (75 Pa TMA for 0.1 s). The front position (AR 215 at half thickness) fixes the exposure; its width fixes s.
  Published values: 0.5-2e-3 (Arts 2019, 275 C set-point); our value is higher, with temperature, the site-density
  assumption and lumped-vs-explicit kinetics as candidate reasons; s = 1e-3 at the same exposure is visibly excluded
  by the data. The chamber runs use 5e-3; the patterned-wafer conclusions are supply-limited and insensitive to it.

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
  channel profile at 300 C) and 8.5e-3 (this work, fit to the Zenodo 3974438 profiles at 300 C, factor 1.5
  sheet-to-sheet). We reran the baseline operating point (1 % TMA, 300 sccm, 300 Pa, 573 K) at s = 5e-4, 1e-3, 2e-3,
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
  lumped Langmuir coefficient with an equilibrium constant, irreversible first-order Langmuir on 4e18 sites). A
  profile grown at a known TMA partial pressure at 300 C, with the site density taken from the measured growth per
  cycle, would remove the trade-off and fix s to the sheet-to-sheet scatter.
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
6. Second-order knobs: pressure up and carrier flow down raise precursor utilization from about 55 % to 95 % at
   fixed precursor flow; gap and pressure set within-wafer edge lag (5 to 18 %) but not time.
7. The water half-cycle needs the same pulse as TMA on patterned wafers despite a 50x lower sticking coefficient,
   and coats deep features more evenly; reversible adsorption, not chamber flushing, governs its purge.
8. In fab units the rule's pulse gives 99 % step coverage as well as 90 %, because the saturation front is 0.25 to
   0.3 s wide at every aspect ratio; 50 cycles of the blanket recipe leave 1.4 nm of a 5.15 nm film on top of an
   AR-100 array and nothing at the bottom. The Al2O3 cycle on the AR-100 wafer is 3.1 times the blanket cycle (37.6
   against 12.2 s; 1.9 against 5.5 wafers per hour for 5 nm), its pulses are three quarters of the cycle, and above
   a wall-area ratio of about 10 throughput falls as 1/(1 + A_wall/A_wafer); the purge is set by water desorption,
   not by the chamber, which clears to 1 % in 2 s while a 10 um hole empties in microseconds.
9. Plasma ALD adds a hard ceiling: above a wall recombination probability of about 3e-3, aspect ratio 100 is
   unreachable by dose. Separated-chemistry (spatial) architectures beat temporal pulsing by the same wall-area
   factor the rule predicts.
10. All of this runs on a laptop in open-source tools, with the feature model verified against an independent code
   and fitted to public measurements, and every number traceable to a script in the repository.

## Acknowledgements
Validation data by R. L. Puurunen and O. M. E. Ylivaara (Aalto University, Zenodo 3974438, CC BY 4.0). Machball by
A. Yanguas-Gil and J. W. Elam (Argonne). OpenFOAM v2606 (OpenCFD) via the gerlero/openfoam-app build.

## Code and data availability
OpenFOAM v2606 cases, the `aldWallKinetics` boundary-condition library, the Python `aldconf` solver and its command
line (`python -m aldconf` evaluates the rule, the surrogate and the Pareto front), tests and all
figures: https://github.com/Soham-ChemE/ald-conformality. The repository also ships a layout-to-recipe tool on its
live page (pattern zones in, per-zone local rule and recipe pulse out, calibrated on the two non-uniform runs of
Section 3.9) and continuous integration that rebuilds the boundary conditions and reruns a regression case in an
OpenFOAM container (77 tests). Validation data: Zenodo 3974438 (CC BY 4.0).

## References (to format)
Gordon, Hausmann, Kim, Shepard, Chem. Vap. Deposition 9, 73 (2003). Yanguas-Gil and Elam, Chem. Vap. Deposition
18, 46 (2012); Machball, machball.readthedocs.io. Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018).
Arts et al., J. Vac. Sci. Technol. A 37, 030908 (2019). Arts et al., J. Phys. Chem. C (2019); (2021). Cremers,
Puurunen, Knoops, Appl. Phys. Rev. 6, 021302 (2019). Yim et al., PCCP (2020), doi 10.1039/D0CP03358H; Zenodo
3974438. Jarvilehto et al., PCCP (2023). Applied Materials, "Semiconductor equipment and processes need digital
twins" (blog); Applied Materials, Olympia ALD release (2015); Centris Spectral SiN ALD release. OpenFOAM v2606;
ViennaPS 4.7.0 (TU Wien) [used for feature topography, TODO].
