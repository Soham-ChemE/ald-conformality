# Every blanket-wafer ALD recipe hides a 40x under-dose: a chamber-to-feature model of supply-limited atomic layer deposition

Soham Kavathekar, Department of Chemical and Biomolecular Engineering, University of Pennsylvania

*Draft 0.9.1 (public), 2026-11-10. Every number is from a run documented in `paper/assumptions.md`.*

## Abstract

ALD recipes are tuned on blanket wafers. Product wafers are not blank: a 3D-NAND array of 100 nm holes at aspect
ratio 100 carries 56 times the surface area of the wafer under it, and all of it has to saturate in one pulse. We
built a model of a 300 mm showerhead ALD chamber in which every point on the wafer runs its own feature-scale model
of a hole, coupled both ways to the chamber gas, and asked how long the pulse must be. The answer is a single rule:
the saturating pulse equals 1.2 times (1 + patterned area ratio) times the time the showerhead needs to deliver one
monolayer of precursor. It held within 25 % across aspect ratios 20 to 100 and a twenty-fold range of precursor feed;
on non-uniform layouts the local wall-area ratio of the densest region sets the pulse, not the wafer average.
At a typical operating point the blanket wafer saturates in 0.35 s and the product wafer in 14 s. Chamber pressure,
showerhead gap and carrier flow do not move that number; precursor delivery rate does. The water half-cycle, fifty
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
the chamber shows it: the precursor stops being consumed and appears downstream.

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

## 4. Discussion
- Throughput on patterned wafers is bought with precursor delivery rate (vapour pressure, bubbler carrier flow,
  injector design), not with chamber pressure, gap or carrier flow; those knobs move uniformity.
- Recipe transfer: the blanket-to-product factor is the wall-area ratio, 11 to 56 here; monitor-wafer recipes
  under-dose product wafers by that factor at the same supply.
- Precursor utilization: 60 % (TMA) and 50 % (H2O) over the saturating pulse on the patterned wafer; the remainder
  is pumped after saturation. Shorter pulses with feedback (end-point on surface density) would recover it [TODO
  estimate].
- Low sticking is a conformality asset and a throughput cost only on blanket wafers.
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
  Clausing correction; s and N0 are literature-range inputs (PillarHall fit 8.5e-3; Arts 2019 0.5-2e-3).

## 5. Conclusions

1. On a patterned wafer, the ALD pulse that saturates the deepest feature is set by precursor supply:
   t_pulse = 1.2 (1 + A_wall/A_wafer) N0/(U c0), within 25 % over aspect ratios 20-100 and a twenty-fold range of
   feed fraction, and independent of hole diameter at fixed wall area, of temperature at fixed sccm and fraction,
   and of the chamber geometry outboard of the wafer.
2. The same normalisation puts a measured fluidized-bed ALD run at 1.12 supply-limit times and orders the measured
   utilization of a rotary-drum reactor with loading; wafers and powders sit on one curve of utilization against
   surface area spanning five orders of magnitude.
3. A recipe saturated on a blanket wafer under-doses a product wafer by the wall-area ratio, 11 to 56 here. Blanket
   wafers are the kinetics-limited exception to the rule, not the norm.
4. On a non-uniform layout the pulse is set by the local wall-area ratio of the most demanding region, not by the
   wafer average: the average under-predicts by 2.3x (dense centre) and 1.2x (aspect-ratio ramp), the local ratio
   is conservative by 15 to 21 %, and a pattern boundary puts a 26 % saturation gradient inside the array.
5. Second-order knobs: pressure up and carrier flow down raise precursor utilization from about 55 % to 95 % at
   fixed precursor flow; gap and pressure set within-wafer edge lag (5 to 18 %) but not time.
6. The water half-cycle needs the same pulse as TMA on patterned wafers despite a 50x lower sticking coefficient,
   and coats deep features more evenly; reversible adsorption, not chamber flushing, governs its purge.
7. Plasma ALD adds a hard ceiling: above a wall recombination probability of about 3e-3, aspect ratio 100 is
   unreachable by dose. Separated-chemistry (spatial) architectures beat temporal pulsing by the same wall-area
   factor the rule predicts.
8. All of this runs on a laptop in open-source tools, with the feature model verified against an independent code
   and fitted to public measurements, and every number traceable to a script in the repository.

## Acknowledgements
Validation data by R. L. Puurunen and O. M. E. Ylivaara (Aalto University, Zenodo 3974438, CC BY 4.0). Machball by
A. Yanguas-Gil and J. W. Elam (Argonne). OpenFOAM v2606 (OpenCFD) via the gerlero/openfoam-app build.

## Code and data availability
OpenFOAM v2606 cases, the `aldWallKinetics` boundary-condition library, the Python `aldconf` solver and its command
line (`python -m aldconf` evaluates the rule, the surrogate and the Pareto front), tests and all
figures: https://github.com/Soham-ChemE/ald-conformality. Validation data: Zenodo 3974438 (CC BY 4.0).

## References (to format)
Gordon, Hausmann, Kim, Shepard, Chem. Vap. Deposition 9, 73 (2003). Yanguas-Gil and Elam, Chem. Vap. Deposition
18, 46 (2012); Machball, machball.readthedocs.io. Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018).
Arts et al., J. Vac. Sci. Technol. A 37, 030908 (2019). Arts et al., J. Phys. Chem. C (2019); (2021). Cremers,
Puurunen, Knoops, Appl. Phys. Rev. 6, 021302 (2019). Yim et al., PCCP (2020), doi 10.1039/D0CP03358H; Zenodo
3974438. Jarvilehto et al., PCCP (2023). Applied Materials, "Semiconductor equipment and processes need digital
twins" (blog); Applied Materials, Olympia ALD release (2015); Centris Spectral SiN ALD release. OpenFOAM v2606;
ViennaPS 4.7.0 (TU Wien) [used for feature topography, TODO].
