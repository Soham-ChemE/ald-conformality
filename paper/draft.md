# Supply-limited atomic layer deposition on high-aspect-ratio wafers: a chamber-to-feature digital twin

Soham Kavathekar, Department of Chemical and Biomolecular Engineering, University of Pennsylvania

*Draft 0.2, 2026-11-10. Numbers below are from the runs in `paper/assumptions.md`; sections marked [TODO] await runs.*

## Abstract

Atomic layer deposition (ALD) recipes are developed on blanket wafers and then transferred to product wafers whose
patterned surface area can exceed the wafer area by one to two orders of magnitude. This paper presents a
concurrent multiscale model of a single-wafer showerhead ALD chamber in which every point on the wafer carries a
feature-scale model of a high-aspect-ratio hole, two-way coupled to the chamber gas phase. The chamber is solved
with OpenFOAM (laminar carrier flow, pulsed precursor transport with a total-flux inlet); the wafer boundary
condition integrates Langmuir coverage on the flat surface and Knudsen diffusion with Langmuir uptake inside the
holes. The feature model is verified against an independent ballistic view-factor code (Machball) to within 0.3 to
7 % and [TODO] validated against published PillarHall thickness profiles. Across 16 operating points spanning TMA
feed fractions of 0.5 to 10 % and hole aspect ratios of 20 to 100 on a 3D-NAND-like array, the pulse time needed to
saturate the hole bottoms everywhere on a 300 mm wafer lies within 5 to 25 % of a pure supply limit,
t_pulse = 1.2 (1 + A_wall/A_wafer) N0/(U c0), independent of aspect ratio except through the wall-area ratio, while
blanket wafers are kinetics-limited. A recipe saturating a blanket wafer in 0.35 s under-doses an aspect-ratio-100
product wafer by a factor of 40. The H2O half-cycle, with a 50x lower sticking coefficient, needs the same pulse on
the patterned wafer but coats it far more evenly. The result is a closed-form design rule for pulse time and
precursor utilization that depends on precursor delivery rate and patterned area, not on chamber pressure, gap or
carrier flow, with direct consequences for recipe transfer and for the temporal-versus-spatial ALD trade-off [TODO].

## 1. Introduction

- ALD in 3D NAND and DRAM: channel holes and word-line recesses at aspect ratios of 50 to 100 and beyond; the film
  must reach the bottom of every feature on the wafer within one pulse.
- The gap in current practice: recipes tuned on blanket monitor wafers; feature-scale models (Gordon 2003,
  Yanguas-Gil and Elam, Ylilammi 2018, Cremers 2019 review) assume a prescribed precursor density at the feature
  mouth; chamber CFD models assume a prescribed wafer uptake. Neither sees the loading effect, in which the patterned
  wafer depresses the gas-phase density that the features depend on.
- Equipment makers are moving to equipment-and-process digital twins (Applied Materials' AppliedTwin and its
  Computational Products and Solutions group; GPU-accelerated topography simulation) [cite blog posts]. This work is
  a small, open, fully reproducible example of that genre for ALD.
- Contributions: (i) a two-way coupled chamber-feature boundary condition for OpenFOAM; (ii) two-code verification
  and public-data validation of the feature model; (iii) the supply-limit design rule and its range of validity;
  (iv) the blanket-to-product recipe-transfer factor; (v) [TODO] the temporal-vs-spatial ALD crossover and the
  plasma-ALD recombination ceiling.

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
- [TODO] Validation against the Zenodo 3974438 PillarHall Al2O3 profiles (TMA/H2O 300 C, 500 cycles): fit s_TMA,
  compare with Arts 2019 and Ylilammi 2018.

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

## 4. Discussion
- Throughput on patterned wafers is bought with precursor delivery rate (vapour pressure, bubbler carrier flow,
  injector design), not with chamber pressure, gap or carrier flow; those knobs move uniformity.
- Recipe transfer: the blanket-to-product factor is the wall-area ratio, 11 to 56 here; monitor-wafer recipes
  under-dose product wafers by that factor at the same supply.
- Precursor utilization: 60 % (TMA) and 50 % (H2O) over the saturating pulse on the patterned wafer; the remainder
  is pumped after saturation. Shorter pulses with feedback (end-point on surface density) would recover it [TODO
  estimate].
- Low sticking is a conformality asset and a throughput cost only on blanket wafers.
- Limitations: isothermal; frozen carrier flow; irreversible Langmuir (no desorption, which governs the real water
  purge); single hole size and density; Knudsen D_K = d vbar/3 without Clausing correction; s and N0 from the
  literature, not fitted (until the PillarHall validation is in); idealised bare floor at the wafer edge.

## 5. Conclusions
[Write after 3.6, 3.7 and the validation.]

## Code and data availability
OpenFOAM v2606 cases, the `aldWallKinetics` boundary-condition library, the Python `aldconf` solver, tests and all
figures: GitHub (to be published). Validation data: Zenodo 3974438 (CC BY 4.0).

## References (to format)
Gordon, Hausmann, Kim, Shepard, Chem. Vap. Deposition 9, 73 (2003). Yanguas-Gil and Elam, Chem. Vap. Deposition
18, 46 (2012); Machball, machball.readthedocs.io. Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018).
Arts et al., J. Vac. Sci. Technol. A 37, 030908 (2019). Arts et al., J. Phys. Chem. C (2019); (2021). Cremers,
Puurunen, Knoops, Appl. Phys. Rev. 6, 021302 (2019). Yim et al., PCCP (2020), doi 10.1039/D0CP03358H; Zenodo
3974438. Jarvilehto et al., PCCP (2023). Applied Materials, "Semiconductor equipment and processes need digital
twins" (blog); Applied Materials, Olympia ALD release (2015); Centris Spectral SiN ALD release. OpenFOAM v2606;
ViennaPS 4.7.0 (TU Wien) [used for feature topography, TODO].
