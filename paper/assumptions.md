# Reactor-scale model: assumptions and first results (2026-11-09, night 1)

## Geometry and operating point
Axisymmetric single-wafer showerhead chamber, 5 degree wedge (OpenFOAM v2606, blockMesh, 6000 cells, orthogonal).
r = 0..200 mm, gap 15 mm. Wafer r < 150 mm on the floor; an un-absorbing floor from 150 to 200 mm; annular outlet
at r = 200 mm over the full gap. Showerhead face over r < 150 mm as a uniform downward inflow.
N2 carrier, 300 Pa, 573 K: rho = 1.76e-3 kg/m3, mu = 2.9e-5 Pa s, nu = 1.65e-2 m2/s. 300 sccm -> Q = 3.54e-3 m3/s
at chamber conditions, inlet speed 0.050 m/s. Re = U g / nu = 0.045 (creeping flow), residence time V/Q = 0.53 s.
TMA diffusivity in N2 at these conditions: D = 8.4e-3 m2/s (Chapman-Enskog scaling of ~8e-6 m2/s at 1 atm, 300 K;
a parameter to vary, not a measured value). Gap Peclet number U g / D = 0.09: transport across the gap is diffusive.

## Equations and boundary conditions
1. Steady laminar incompressible N2 flow, simpleFoam, frozen for the species step (TMA fraction small; its injection does
   not change the carrier flow). Mass balance closes to machine precision: showerhead in = outlet out = 4.90e-5 m3/s
   per 5 degree wedge = 3.53e-3 m3/s per 360 degrees.
2. TMA as a passive scalar, scalarTransportFoam, dimensionless concentration c/c0, c0 = TMA concentration in the
   showerhead feed during the pulse. 0.1 s pulse, 4 s purge, dt = 1 ms, written every 5 ms.
   - Showerhead: Danckwerts total-flux inlet, U (c_face - c0(t)) = D (c_face - c_cell)/delta, implemented as
     uniformMixed with value fraction f = U/(U + D/delta) = 1.49e-3. A fixed-concentration inlet is WRONG here: with
     D/L = 0.56 m/s versus U = 0.05 m/s it feeds 11x more TMA by diffusion than the gas actually carries (the first
     run absorbed 1179 % of the injected TMA). With the Danckwerts inlet the wafer captures 96.5 %, the rest leaves
     through the outlet.
   - Wafer: c = 0, a perfectly absorbing surface, the sticking -> 1 limit. The surface-kinetics (finite s, Langmuir
     saturation) version follows; it will lower the capture and lengthen the tail.
   - Floor, top wall: zero flux. Outlet: inletOutlet. Isothermal.
   - Wafer flux = D c_cell / delta from the wafer-adjacent cell row (delta = 0.25 mm), which is the finite-volume
     wall gradient on this orthogonal mesh. (OpenFOAM's patch sampler with nearCellValue returned the patch value, so
     the fields are read directly.)

## First results (figures/reactor_dose_vs_r.png, figures/reactor_dose_vs_r.csv)
- Centre flux plateaus at 0.048 m/s x c0, i.e. at the supply rate U c0: at Pe = 0.09 the wafer takes everything the
  showerhead delivers; the process is supply-limited, not transport-limited.
- Dose per pulse at the centre: 4.79e-3 m x c0 (multiply by c0 in mol/m3 for mol/m2).
- Within-wafer dose: flat to within 1 % out to r = 100 mm, a 7 % dip centred at r = 142 mm where the carrier
  accelerates radially toward the outlet and sweeps TMA past the surface, then a rise in the last 2 mm (to 1.3x at
  r = 149 mm) fed by lateral diffusion from the un-absorbing floor region. Over r < 140 mm the spread is 7.1 %.
  The edge spike is an artefact of the idealised bare floor; a real edge ring / focus ring geometry is the first
  geometry refinement.
- Purge: centre flux falls below 1 % of its peak 0.05 s after the pulse ends (0.15 s). At this gap and pressure a 4 s
  purge is ~80x longer than the gas-phase clearing time; purge length is set by surface desorption and dead volumes,
  not by the chamber flow. This is a throughput statement the feature-scale model will have to confirm.

## Not yet modelled (ordered)
finite sticking and Langmuir saturation on the wafer; pressure and gap sweep; edge-ring geometry; temperature;
coupling the dose(r) curve into the feature-scale model (aldconf solver and ViennaPS) for step coverage vs radius;
H2O half-cycle; DOE and GP surrogate in MATLAB.

## Night 1, part 2: finite sticking and the self-limiting wafer
Sticking sweep (mixed BC, j = s (vbar/4) c_face, vbar = 410 m/s for TMA at 573 K; figures/reactor_sticking_sweep.png):
| wafer | captured % | centre dose/c0 (m) | spread r<140 mm | gas clears (flux < 1 % peak) |
| absorbing | 96.5 | 4.79e-3 | 7.1 % | 0.15 s |
| s = 1 | 96.5 | 4.79e-3 | 6.9 % | 0.15 s |
| s = 0.1 | 96.0 | 4.77e-3 | 5.6 % | 0.16 s |
| s = 0.01 | 91.9 | 4.57e-3 | 1.6 % | 0.21 s |
| s = 0.001 | 65.4 | 3.26e-3 | 0.2 % | 0.58 s |
Reading: the lower the sticking, the more uniform the dose (the gas phase homogenises when uptake is slow) and the
longer the purge tail (unreacted TMA sits in the chamber). Realistic TMA sticking is 1e-3 to 1e-2 (Arts 2019,
Ylilammi 2018), so expect ~65-90 % capture, sub-2 % within-wafer dose spread, and a 0.2-0.6 s gas-phase purge.

Self-limiting wafer: custom OpenFOAM boundary condition `langmuirSticking` (reactor/lib/aldWallKinetics, compiled with
wmake against OpenFOAM v2606): per-face coverage theta integrated in time, uptake (vbar/4) s (1 - theta) c_face,
optional recombination term r for plasma radicals. Parameters used: s = 5e-3, N0 = 4e18 sites/m2 (about one
TMA-saturated monolayer on Al2O3), c0 = 3.79e20 /m3 (3 Pa TMA at 573 K). Both N0 and the TMA partial pressure are
assumptions to be swept; the sticking value is literature-range, not fitted yet.
Results (figures/reactor_langmuir_coverage.png):
- The showerhead delivers U c0 = 1.9e19 molecules/m2/s: one monolayer's worth every 0.21 s at 100 % capture. A 0.1 s
  pulse therefore cannot saturate the wafer at this flow and partial pressure: final coverage 0.40 everywhere, of
  which 0.25 arrived during the pulse and the rest during the purge from TMA still in the chamber.
- A 0.5 s pulse saturates: 90 % coverage at 0.345 s in the centre and 0.385 s at r = 140 mm; the edge saturates
  last, 12 % later than the centre. The minimum pulse for a saturated wafer at this operating point is set by the
  wafer edge, about 0.4 s, not by the 0.1 s recipe value from the tube-reactor dataset.
- Interpretation: at 300 Pa, 300 sccm and 3 Pa TMA the process is supply-limited; pulse time scales inversely with
  TMA partial pressure x flow. This is the knob map the DOE will quantify (pressure, flow, TMA fraction, gap).

## Night 1, part 3: operating-point sweep with the self-limiting wafer (figures/reactor_sweep_operating.png, .csv)
One factor at a time around 300 Pa / 300 sccm / 15 mm / 1 % TMA, s = 5e-3, N0 = 4e18, 1.0 s pulse + 1.0 s purge.
| case | Pe | monolayer supply time N0/(U c0) | t90 centre | t90 worst on wafer | edge lag | theta after 0.1 s |
| base 300 Pa, 300 sccm, 15 mm, 1 % | 0.090 | 0.21 s | 0.35 s | 0.38 s | 8.6 % | 0.25 |
| 100 Pa | 0.090 | 0.21 | 0.41 | 0.43 | 4.9 % | 0.24 |
| 1000 Pa | 0.090 | 0.21 | 0.35 | 0.40 | 14.3 % | 0.18 |
| 150 sccm | 0.045 | 0.42 | 0.56 | 0.64 | 14.3 % | 0.13 |
| 600 sccm | 0.179 | 0.105 | 0.23 | 0.24 | 4.3 % | 0.43 |
| 10 mm gap | 0.060 | 0.21 | 0.31 | 0.34 | 9.7 % | 0.30 |
| 25 mm gap | 0.149 | 0.21 | 0.40 | 0.47 | 17.5 % | 0.16 |
| 0.5 % TMA | 0.090 | 0.42 | 0.60 | 0.65 | 8.3 % | 0.13 |
| 3 % TMA | 0.090 | 0.07 | 0.16 | 0.19 | 18.8 % | 0.63 |
Reading:
- Pulse time is set by supply. Doubling TMA fraction or N2 flow halves the monolayer supply time and t90 follows it
  (t90 = 1.6 to 2.3 x the supply time; the excess is the finite sticking plus the Langmuir (1 - theta) tail).
- Pressure barely moves t90 at fixed sccm (U and D both scale as 1/p, so the gap Peclet number is unchanged), but it
  moves uniformity: the edge lags the centre by 5 % at 100 Pa and 14 % at 1000 Pa, because lateral mixing slows.
- Gap: 10 mm beats 25 mm on both counts, 0.31 vs 0.40 s and 10 % vs 18 % edge lag. Residence volume and diffusion
  distance both shrink.
- The within-wafer saturation-time spread (edge lag) is the uniformity figure of merit for the feature-scale step: the
  last point to saturate sets the recipe pulse time, so every knob that raises the edge lag costs throughput twice.
- Caveats: one-factor-at-a-time, s and N0 from literature ranges, bare-floor edge, isothermal, frozen carrier flow,
  no H2O half-cycle. A central composite design over (Q, x, g, p) with the GP surrogate is the next layer.

## Night 1, part 4: reactor -> feature coupling (coupling/reactor_to_feature.py, figures/coupling_hole_coverage.png)
Method: the TMA number density at the wafer surface from the chamber model (wall value of the langmuirSticking patch,
times c0) at four radii drives the feature-scale Knudsen + Langmuir model (aldconf.simulate_driven) of a 100 nm hole
at AR 20, 50, 100, s = 5e-3, N0 = 4e18. One-way coupling: holes do not feed back on the chamber (no loading effect yet).
Benchmarks of the feature solver alone at constant 0.3 Pa TMA (the surface density of the unsaturated base case):
95 % step coverage needs 1.20 s (AR 20), 4.45 s (AR 50), 16.2 s (AR 100): ratio 3.7 per doubling, the AR^2 law.
Results, base case (300 sccm, 1 % TMA, 1 s pulse):
- Flat wafer saturates at 0.35 s. AR 20 hole bottom: 0.49 s. AR 50: 0.92 s. AR 100: 3 % coverage at the bottom after
  the full 1 s pulse plus 1 s purge; the saturation front stops about 60 diameters deep.
- The hole fills far faster than the constant-density benchmark predicts (AR 50: 0.92 s vs 4.45 s). Mechanism: once
  the flat surface saturates at 0.35 s, the wafer stops consuming TMA, the gas density at the surface rises from 0.1 c0
  to c0, and the holes see ten times the precursor. The flat wafer "shields" its own holes until it is full. Neither
  model alone shows this; it is a reactor-feature interaction.
- Radius matters little for the holes: centre vs r = 149 mm shifts the AR 50 bottom t90 from 0.92 to 1.01 s (10 %),
  the same edge lag as the flat wafer. At this operating point the aspect ratio, not the wafer position, is the
  bottleneck.
- 3 % TMA fills AR 100 in 0.89 s (centre) to 0.99 s (edge); 600 sccm at 1 % does not (6 % at the bottom). Supply of
  precursor molecules, not carrier flow, is what deep features need: tripling the TMA fraction did what doubling the
  flow could not.
Caveats: single hole, no loading (a real 3D NAND wafer with ~1e9 holes of AR 100 has ~400 x the flat area and would
depress the surface density far more, which makes the "shielding" effect stronger and the pulse longer); H2O half-cycle
not modelled; hole-scale D_K uses d vbar/3 without the Clausing correction; Machball cross-check still to be wired
(its saturation() needs a Structure built from view factors, see machball.structure / vf).

## Night 1, part 5: concurrent multiscale wall (reactor/lib/aldWallKinetics, patternedSurface BC)
Every wafer face now carries a flat Langmuir surface plus a 200-cell Knudsen/Langmuir model of a 100 nm hole; the
face's net uptake K_eff c_face = (1 - phi) flat uptake + phi x hole-mouth flux, two-way coupled to the chamber.
Verification: with a negligible hole density (one-way limit) the BC's hole-bottom coverage matches the Python
aldconf.simulate_driven model driven by the same face density to within 1 % at every time (t90 0.92 s in both).
3D-NAND-like array: 100 nm holes, 180 nm hex pitch, half the wafer patterned -> 1.78e13 holes/m2, 14 % open area,
wall area 56 x wafer area at AR 100 (28 x at AR 50).
Loading result, 1 s pulse at 1 % TMA: the surface density collapses to 0.019 c0 (0.24 c0 on a blanket wafer), the
wafer captures 94 % of everything injected, and nothing saturates: flat regions and hole bottoms both stay far from
90 % after 1 s pulse + 5 s purge. Supply arithmetic: the pulse delivers 1.9e19 molecules/m2 against a demand of
57 x 4e18 = 2.3e20, i.e. 8 % of a monolayer-equivalent. A blanket-wafer recipe (0.35 s) under-doses a 3D NAND wafer
by more than an order of magnitude at this precursor supply. Long-pulse runs (20 s at 1 %, 8 s at 3 %) are next.

Long-pulse results for the 3D-NAND-like AR-100 wafer (figures/coupling_patterned_long.png):
| supply | blanket wafer t90 | patterned wafer, flat regions t90 | hole bottom t90 centre | hole bottom t90 r = 149 mm | TMA captured over the pulse |
| 1 % TMA, 300 sccm, 20 s pulse | 0.35 s | 1.95 s | 14.05 s | 13.3 s | 60 % |
| 3 % TMA, 300 sccm, 8 s pulse | 0.16 s | 0.70 s | 4.75 s | 4.55 s | 50 % |
Reading:
- The pulse that saturates the bottom of AR-100 holes across a 3D-NAND-like wafer is 14 s at 1 % TMA: 40 x the
  blanket-wafer value of 0.35 s. Tripling the TMA fraction cuts it to 4.75 s, almost exactly 1/3: the process is
  supply-limited end to end, so pulse time x precursor supply rate is roughly constant (14.05 x 1 = 14.0, 4.75 x 3 =
  14.3 "percent-seconds"). The knob is precursor molecules per second, not carrier flow, pressure or gap.
- The flat regions between holes take 1.95 s instead of 0.35 s: the holes starve the flat surface (loading).
- The wafer edge saturates slightly earlier than the centre on the patterned wafer (13.3 vs 14.05 s), the opposite of
  the blanket wafer (edge 8 % later): with the surface density collapsed to 0.02 c0 everywhere, the un-patterned floor
  ring outside the wafer feeds the edge laterally. An edge ring that is also a sink would remove this.
- Captured fraction over the pulse: 60 % and 50 %, the rest pumped away mostly after the surfaces saturate.
Caveats unchanged (isothermal, frozen carrier, no H2O half-cycle, Knudsen D_K = d vbar/3, s and N0 literature-range,
hole density and pattern coverage assumed). Next: a 2D map of saturating pulse time versus TMA fraction and AR,
the H2O half-cycle (s ~ 1e-4, which will dominate cycle time), the edge-ring geometry, and DOE/GP over the knob set.

## Cross-check of the hole model against Machball (Argonne, ballistic view-factor model)
Same inputs (s = 5e-3, N0 = 4e18, 573 K, 0.3 Pa TMA, 100 nm cylinder), no fitting. Time for 95 % coverage at the
bottom: AR 20: aldconf 1.20 s vs Machball 1.29 s (7 %); AR 50: 4.45 s vs 4.46 s (0.3 %); AR 100: 16.2 s vs 15.6 s (4 %). Two different transport
descriptions (Knudsen continuum with Langmuir sink vs ballistic transport with cosine re-emission) agree, which is the
same kind of two-method check as the Monte Carlo comparison in Jarvilehto et al. 2023. Test: tests/test_machball_crosscheck.py.

## Night 1, part 6: the H2O half-cycle (reactor/make_h2o_cases.py; figures/coupling_patterned_h2o_vs_tma.png)
H2O at 573 K: vbar = 821 m/s, D(H2O-N2) = 2.5e-2 m2/s at 300 Pa (Chapman-Enskog scaling of 2.4e-5 at 1 atm, 300 K,
an estimate), s = 1e-4 on the TMA-terminated surface (Arts et al. 2019 JVSTA give 0.8-2e-4), N0 = 4e18 (one H2O
reaction per TMA site, a simplification), 1 % H2O in the feed. Irreversible Langmuir only: reversible water
adsorption and its slow desorption, which is what really sets the water purge in practice, is not modelled.
| | blanket wafer t90 (centre / r = 140) | patterned AR-100 wafer: flat regions t90 | hole bottom t90 (centre / r = 140 / r = 149) |
| TMA, 1 % | 0.35 / 0.38 s | 1.95 s | 14.05 / 14.0 / 13.3 s |
| H2O, 1 % | 1.74 / 1.82 s | 10.8 s | 13.4 / 13.6 / 13.4 s |
Reading:
- Blanket wafer: water is 5 x slower than TMA (1.74 vs 0.35 s) because with s = 1e-4 the kinetic uptake coefficient
  s vbar/4 = 0.02 m/s is below the inlet speed 0.05 m/s: the surface, not the supply, is the bottleneck (29 % capture
  at zero coverage, falling as the surface fills).
- Patterned wafer: water and TMA need almost the same pulse (13.4 vs 14.05 s). With 57 monolayer-equivalents of area
  to fill, both half-cycles are supply-limited and the sticking coefficient stops mattering for the total.
- Conformality inverts between the two: TMA's hole bottoms lag its flat regions by 12 s (1.95 -> 14.05 s), water's by
  only 2.6 s (10.8 -> 13.4 s). Low sticking penetrates deep features more evenly (the classic result), which here shows
  up as a small flat-to-bottom lag.
- Cycle time at this operating point (1 % precursors, 300 sccm, 300 Pa): about 14 + 13.4 s of pulses plus purges,
  i.e. ~30 s per cycle for an AR-100 3D-NAND-like wafer, against ~2.5 s for a blanket wafer. At 3 % precursors the TMA
  half drops to 4.75 s; the water half should scale the same way (to be run).

## Night 1, part 7: the pulse-time map (reactor/sweep_pattern_map.py, figures/pulse_time_map.png, .csv)
16 two-way multiscale runs: blanket wafer and 3D-NAND-like wafers with AR 20 / 50 / 100 holes (wall area 11 / 28 / 56
x wafer area) at 0.5, 1, 3, 10 % TMA; baseline chamber (300 Pa, 300 sccm, 15 mm), s = 5e-3, N0 = 4e18.
| TMA | blanket t90 | AR 20 bottoms | AR 50 bottoms | AR 100 bottoms | (worst point on the wafer) |
| 0.5 % | 0.65 s | 5.37 s | 13.4 s | 28.1 s |
| 1 % | 0.38 s | 2.71 s | 6.73 s | 14.1 s |
| 3 % | 0.18 s | 0.95 s | 2.32 s | 4.79 s |
| 10 % | 0.09 s | 0.33 s | 0.75 s | 1.51 s |
Design rule that falls out: the patterned wafers sit at 1.05-1.25 x the pure supply limit everywhere on the map,
    t_pulse ~ 1.2 x (1 + A_wall/A_wafer) x N0 / (U c0),
with U c0 the precursor molecules delivered per unit wafer area per second. Every patterned point scales as 1/x within
a few percent, and the aspect ratio enters only through the wall-area ratio: the Knudsen transport inside a 100 nm
hole is not the bottleneck even at AR 100 once the chamber is supply-limited. The blanket wafer is the exception: it
sits at 1.4-3.5 x its (much smaller) supply limit because with so little area to fill the finite sticking (5e-3) and
the Langmuir tail dominate; at 10 % TMA it is purely kinetics-limited.
Why this matters to a tool: throughput on a 3D-NAND wafer is bought with precursor delivery rate (vapour pressure,
carrier flow through the bubbler, injector design), not with chamber pressure, gap or carrier flow; and a recipe tuned
on blanket monitor wafers under-doses the product wafer by the wall-area ratio (11 to 56 x here).
Caveats: single hole size and density (14 % open area, 50 % patterned); s, N0 literature-range; isothermal; frozen
carrier flow; irreversible Langmuir; no desorption; H2O map not yet run (expect the same 1/x law from part 6).

## Night 1, part 8: plasma-ALD recombination ceiling (coupling/peald_ceiling.py, figures/peald_ceiling.png)
Feature solver, 100 nm hole, radical mass 16, s = 0.05 (assumed: radicals react readily on the methyl-terminated
surface), 0.3 Pa-equivalent radical density at the mouth, loss probability r on every wall collision.
Time to 90 % bottom coverage (s):
| r | 0 | 1e-5 | 3e-5 | 1e-4 | 3e-4 | 1e-3 | 3e-3 | 1e-2 |
| AR 20 | 0.31 | 0.31 | 0.31 | 0.31 | 0.32 | 0.35 | 0.46 | 1.07 |
| AR 50 | 1.90 | 1.91 | 1.94 | 2.03 | 2.31 | 3.54 | 10.6 | 206 |
| AR 100 | 7.98 | 8.18 | 8.59 | 10.2 | 16.2 | 64.5 | not reached | not reached |
Largest AR whose bottom reaches 90 % within a 5 s exposure: 79 (r = 1e-5), 75 (1e-4), 56 (1e-3), 29 (1e-2).
Reading: for O atoms on SiO2/TiO2 (r ~ 6e-5, Arts 2019) recombination costs AR-100 holes only ~10-25 % extra time;
the ceiling bites between r = 3e-4 and 3e-3, i.e. on surfaces with high recombination (metals, some nitrides), where
AR 100 becomes unreachable by dose alone. This is the regime where plasma source design (radical flux) and surface
choice, not pulse time, decide conformality. Caveat: s for radicals assumed; results scale with s vbar/4 vs D_K.

## Night 1, part 9: temporal vs spatial ALD from the design rule (coupling/temporal_vs_spatial.py)
Molecules required per unit wafer area per half-cycle Phi_req = 1.2 (1 + A_wall/A_wafer) N0. Temporal: whole wafer
under the showerhead flux J_T = U c0 (1.9e19 /m2/s at 1 %), two half-cycles plus two 0.5 s purges. Spatial
(separated-chemistry, Olympia-type): precursor zone flux J_S = k J_T for a dwell t_d = f_zone tau_rot per pass,
passes n = ceil(Phi_req/(J_S t_d)), no purge steps; k = 10, 30, 100, f_zone = 0.25, tau_rot = 1 s (assumptions to
sweep, not measurements).
| A_wall/A_wafer | temporal 1 % | temporal 3 % | spatial k = 30 |
| 1 (blanket) | 2.0 s | 1.3 s | 2 s |
| 11 (AR 20) | 6.8 s | 2.9 s | 2 s |
| 28 (AR 50) | 16.3 s | 6.1 s | 4 s |
| 56 (AR 100) | 29.9 s | 10.6 s | 4 s |
| 112 | 55.9 s | 19.3 s | 8 s |
Reading: on blanket wafers the two architectures tie (purge overhead vs rotation overhead); the spatial advantage
grows linearly with patterned area because the zone delivers precursor at a flux the showerhead cannot match without
raising the feed fraction. This is the quantitative form of the argument in Applied's 3D-NAND deposition patent
literature and the premise of Olympia's separated-chemistry design; the crossover sits at A_wall/A_wafer ~ 5-10 for
k = 30. The zone flux k is the number a tool maker knows and we do not; the figure is a design chart in k.

## Night 1, part 10: PillarHall validation (coupling/fit_pillarhall.py, figures/pillarhall_fit.png)
Zenodo 3974438, 8 long-channel profiles averaged (L = 0.5, 1, 5 mm), plateau 51.5 nm after 500 cycles (103 pm/cycle),
half-thickness depth 107 um = AR 215 in the 500 nm channel. Two-parameter fit of the aldconf slit model with channel
narrowing (s and the mouth TMA density; N0 = 4e18 fixed, so the pressure is really an exposure/N0 estimate):
s_TMA = 8.5e-3 (sheet-to-sheet factor 1.5), p_TMA = 75 Pa (factor 1.02), exposure 7.5 Pa s per 0.1 s pulse, rms
0.038 of the plateau. Comparison: Arts et al. 2019 (275 C set-point, lumped model) give 0.5-2e-3 for TMA; a thesis
citing Ylilammi 2018 estimates that study's TMA partial pressure at about 43 Pa. Our s is 4x above Arts' upper
bound: candidates are the higher temperature (300 C), the N0 assumption (which trades against p but also sets the
front sharpness through the narrowing term), and the lumped vs explicit kinetics. The curves for s = 1e-3 and 1e-2
at the same exposure bracket the data visibly (figure): the data exclude 1e-3 at this exposure. The value used in
the chamber runs (5e-3) lies between the two literature-anchored estimates; the patterned-wafer results are
supply-limited and insensitive to it, the blanket-wafer numbers are not (t90 ~ 1.6-2.3x the supply time at 5e-3).
Data points beyond AR 300 rise again (membrane remnants near the far pillars, as the record warns) and were excluded.

## Night 2, part 1 (2026-11-10 ~02:10): validation of the supply-limit rule against a measurement
Source: Lu, Yanguas-Gil, Kang, Darapaneni, Mane, Marshall, Elam, "Scalable synthesis of supported catalysts using
fluidized bed ALD", J. Vac. Sci. Technol. A 40 (2022) (OSTI 1902316). Al2O3 ALD on 10 g silica gel, 940 m2 total
area, Forge Nano fluidized bed, 50 sccm N2, 14.5 Torr, TMA at 1.6 +/- 0.3 Torr, 200 C. Measured by in-situ QMS: TMA
breakthrough (saturation) at 24 min, H2O at 13 min. Their numbers: TMA supplied during the 24 min dose 0.0115 +/-
0.0020 mol; TMA needed for saturation 0.0103 mol (from the deposited Al2O3); utilization 91.8 +/- 14.6 %.
Normalised: t0 = sites / supply rate = 0.0103 / (0.0115 / 24 min) = 21.5 min, so tau_sat = 24 / 21.5 = 1.12
(0.92 to 1.31 over their stated supply uncertainty). Their Damkoehler number ~1e4; their own plug-flow model (Eqs. 3-9 of
the report; Yanguas-Gil and Elam, arXiv 2408.13116) predicts tau_sat -> 1 from above as Da grows.
Our two-way chamber model, same normalisation: tau_sat = 1.05 to 1.27 for the twelve patterned-wafer cases
(Da = 124 to 581 by their definition with S the total reactive area and V the chamber volume), and 1.5 to 4.3 for
the blanket wafers (Da = 10, kinetics-limited). Figure: figures/validation_supply_rule.png.
Reading: a wafer carrying tens of wafer-areas of holes behaves like a powder bed: saturation arrives when roughly
1.1 to 1.2 site-equivalents of precursor have been delivered, in the measurement and in the model alike. This is
the measured anchor for the "1.2" in the design rule; the model was not tuned to it (N0, s and the chamber were fixed
before the comparison). One point is one point: the Strempel et al. 2018 fixed-bed data (Davisil silica, 506 m2/g,
atmospheric, 50 mL/min) has the breakthrough curves in its reference 15 and would add a second, lower-Da point.

## Night 2, part 1b: second measured set, Argonne rotary drum (Coile et al., JVST A 38, 052403, 2020; OSTI 1660751)
Hot-wall tube reactor, 200 sccm carrier at 0.8 Torr, 200 C, silica gel in a rotating mesh drum at 120 rpm; TMA
saturation by in-situ QMS (95 % of the integrated CH4 signal). Table 2, TMA utilization vs loading: 2 g 24.7 %,
5 g 47.9 %, 15 g 46.7 %, 30 g 55.3 %, 75 g 70.7 % (about 20 m2/g: 50 g ~ 1000 m2, 75 g ~ 1500 m2 per the paper).
Utilization = reacted/delivered up to saturation = 1/tau_sat, so tau_sat falls from 4.0 (40 m2) to 1.4 (1500 m2):
the more area the reactor holds, the closer saturation sits to the supply limit. Same trend as our wafers (blanket
1.5-4.3 at 0.07 m2 of area; patterned 1.05-1.27 at 0.8-4 m2) and the fluidized bed (1.09 at 940 m2). The drum does
not reach the fluidized bed's utilization because precursor must diffuse through the mesh into a 5 mm bed (the
authors' explanation), i.e. transport to the surface, not supply, caps it; its Damkoehler number by the flow
definition is not meaningful, so the drum points appear on the utilization-vs-area panel only.
Figure: validation_supply_rule.png now has two panels: tau_sat vs Da (chamber, theory, fluidized bed) and precursor
utilization vs total reactive area across five orders of magnitude (0.07 m2 blanket wafer to 1500 m2 of powder).

## Night 2, part 2: edge-ring chamber geometry (reactor/edgering_*, figures/edgering_comparison.png)
Geometry change: a 25 mm non-absorbing edge ring at the wafer plane (r 150-175 mm) and an annular pumping slot whose
outlet sits 25 mm below the wafer plane (r 175-200 mm), replacing the bare floor with the outlet at gap height.
Mass balance exact (4.90e-5 m3/s per wedge in and out). Same baseline conditions, 1 % TMA.
| | blanket t90 centre / r=140 / r=149 (s) | spread | patterned AR-100 bottom t90 centre / r=140 / r=149 (s) | spread |
| bare floor | 0.350 / 0.380 / 0.380 | 8.6 % | 14.01 / 14.01 / 13.33 | 4.2 % |
| edge ring | 0.350 / 0.380 / 0.390 | 11.4 % | 14.05 / 14.05 / 13.35 | 3.9 % |
Reading: the outboard geometry moves the saturation times by at most one write interval (0.01 s blanket, 0.05 s
patterned). The blanket edge lag stays 9-11 %, the patterned edge still saturates slightly early (gas over the
un-absorbing ring feeds it laterally), and the within-wafer spread is unchanged within 3 %. The dose-profile edge
spike seen in part 1 was a feature of the dose integral, not of the saturation time. Conclusion for the paper: the
pulse-time results are insensitive to the chamber geometry outboard of the wafer at this Peclet number; the lever
on edge lag is gap and pressure (part 3), not the ring design.

## Night 2, part 3: reversible water adsorption (kdes parameter of langmuirSticking, developed in a side library and since merged into reactor/lib/aldWallKinetics; figures/h2o_desorption.png)
The flat-surface BC now carries a first-order desorption rate kdes (1/s): dtheta/dt = (vbar/4) s (1-theta) c_s/N0 -
kdes theta, with the desorbed flux returned to the gas through the mixed-condition reference value. Blanket wafer,
H2O at 1 % (s = 1e-4, vbar 821 m/s), 5 s pulse, 3 s purge:
| kdes (1/s) | t90 | coverage at pulse end | 3 s after | outlet H2O below 1 % of the pulse level |
| 0 (irreversible) | 1.74 s | 1.000 | 1.000 | 2.2 s after the pulse (pure gas-phase flushing, 4.6 residence times) |
| 0.1 | 2.04 s | 0.951 | 0.755 | not within 3 s: desorbing water holds the exhaust above 1 % |
| 1 | never | 0.660 | 0.113 | not within 3 s |
| 10 | never | 0.163 | 0.001 | 2.4 s (surface empties within 0.3 s; purge is gas-phase again) |
Reading: with reversible adsorption the surface has an equilibrium coverage K c/(1 + K c) set by the water partial
pressure, and at 1 % H2O a desorption rate above ~0.3 1/s stops the surface reaching 90 %: water would then need a
higher partial pressure or a lower temperature, which is the familiar ALD trade-off. Between 0.1 and 1 1/s the
desorbing layer feeds the exhaust for seconds after the pulse, which is the mechanism behind long water purges that
the irreversible model cannot show. A faithful water model needs two site types, an irreversible chemisorbed
hydroxyl fraction plus a reversible physisorbed layer (Yanguas-Gil's soft-saturation form); kdes here stands in for
the physisorbed part and is a parameter, not a measured value. The TMA results in this paper are unaffected (TMA
chemisorption on hydroxyls is irreversible at 300 C).

## Night 2, part 4: temperature and hole-size sweeps (reactor/temperature, reactor/holesize_*; figures/temperature_sweep.csv)
Temperature 200 / 300 / 400 C at fixed 300 Pa, 300 sccm, 1 % TMA, 15 mm gap. Temperature enters through carrier
density and viscosity, inlet speed (U ~ T at fixed sccm), diffusivity (~T^1.75), precursor density (c0 ~ 1/T) and
thermal speed (~sqrt T); the product U c0, the supply, is the same at every temperature. Sticking kept at 5e-3.
| T | blanket t90 centre / worst | patterned AR-100 bottom t90 centre / worst | flat regions |
| 200 C | 0.350 / 0.390 s | 14.02 / 14.02 s | 1.97 s |
| 300 C | 0.350 / 0.380 s | 14.02 / 14.15 s | 1.97 s |
| 400 C | 0.340 / 0.380 s | 14.15 / 14.15 s | 1.97 s |
Reading: at fixed sccm and fraction the pulse time is temperature-independent to within one write interval, for
the blanket wafer as well as the patterned one, because the supply is; temperature acts on ALD through the chemistry
(s, N0, desorption), which this sweep holds fixed by construction. A recipe transferred between temperatures at
constant sccm and fraction keeps its dose; one transferred at constant partial pressure does not (c0 ~ 1/T).
Hole size 50 / 100 / 200 nm at fixed open area (14 %) and AR 100, so the wall-area ratio is 55.9 in all three
(density 7.1e13, 1.8e13, 4.5e12 /m2):
| d | bottom t90 centre / worst | flat |
| 50 nm | 14.10 / 14.10 s | 2.00 s |
| 100 nm | 14.01 / 14.11 s | 1.96 s |
| 200 nm | 14.10 / 14.10 s | 2.00 s |
Reading: identical within 1 %. The rule's only geometric input is the wall-area ratio; intra-hole Knudsen transport
(D_K ~ d, depth ~ d) is not the bottleneck in the supply-limited regime even at 50 nm and AR 100.

## Night 2, part 5: DOE, response surface, GP surrogate and uncertainty (items 2 and 4)
Face-centred central composite design, 25 runs of the two-way multiscale model on the AR-100 3D-NAND-like wafer:
TMA fraction 0.5-3 % (log), N2 flow 150-600 sccm (log), pressure 100-1000 Pa (log), gap 10-25 mm. Driver
reactor/doe_patterned.py; results figures/doe_results.csv and .png; JMP fits figures/doe_rsm_report.txt
(analysis/doe_rsm.jsl, generated by analysis/make_doe_jsl.py); MATLAB GP and Monte Carlo figures/doe_gp_report.txt
(analysis/doe_gp.m). Three corner cases (0.5 %, 600 sccm, 100 Pa; 3 %, 600 sccm, 100 Pa, 25 mm) needed a 3.2x
supply-time pulse to saturate and were rerun.
JMP response surface (coded units, RS model): log10 t90: R2 0.993, c_x -0.383 (vs -0.389 for exact 1/x), c_Q -0.240
(vs -0.301 for exact 1/Q: t ~ Q^-0.8, the extra carrier flow partly bypasses the wafer), c_p -0.062 (t ~ p^-0.12),
gap not significant (p = 0.16). Utilization: R2 0.88, +0.104 per half-decade of pressure, -0.101 per coded flow
unit, fraction and gap not significant. Edge lag: R2 0.54, range -3.5 to +1.1 %, never near the 5 % constraint.
MATLAB GP (ARD squared-exponential, linear basis), leave-one-out: R2 0.999 (log t90), 0.988 (utilization), 0.78
(edge lag). Grid optimum: fastest point with edge lag <= 5 % is 3 %, 600 sccm, 1000 Pa, 10 mm: 2.38 s at 85 %
utilization (centre point 11.6 s, 85 %). At the SAME precursor supply as the centre (x Q fixed): 2.5 %, 150 sccm,
631 Pa, 10 mm gives 9.83 s and 94 % utilization against 11.6 s and 85 %: 15 % less time and 9 points less waste
for the same precursor flow, by trading carrier gas for pressure.
Marginal means: utilization 68 % at 100 Pa, 84 % at 316 Pa, 88 % at 1000 Pa; 88 % at 150 sccm, 83 % at 300, 68 %
at 600; tau_sat 1.55 at 100 Pa vs 1.14 at 1000 Pa.
Uncertainty: the rule's constant is 1.03-1.25 across the 22 points at p >= 316 Pa and up to 1.96 at 100 Pa with
600 sccm (bypass). Monte Carlo on the rule at the centre point with N0 lognormal (factor 1.25): 11.9 s, 90 % band
8.2 to 17.2 s. The site density is the dominant uncertainty in any absolute pulse time; the ratios (recipe
transfer factor, 1/x scaling, utilization ordering) do not depend on it.
Reading for a tool: the supply-limit rule is the first-order answer; the second-order knobs are pressure (up) and
carrier flow (down), which together move utilization from ~55 % to ~95 % without changing precursor flow. Gap is a
uniformity knob only. These are the two sentences the DOE adds to the paper.

## Night 3, part 1: spatially varying pattern (reactor/make_layout_cases.py, coupling/postprocess_layout.py; figures/pattern_layout.png, .csv)
BC change: patternedSurface keeps the scalar AR and holeDensity as wafer-wide defaults and accepts optional
ARProfile and holeDensityProfile entries, each a table of (radius_m value) pairs, e.g.
`ARProfile table ((0 100) (0.075 100) (0.1 60) (0.15 60));`, interpolated linearly at each face centre's radius
sqrt(y^2 + z^2) (wedge axis x) and clamped outside the table. The per-face values enter the hole length AR d, the
open-area fraction holeDensity pi d^2/4, the blending of flat and hole fluxes, and the wall-area ratio
holeDensity pi d AR d. A face with holeDensity 0 (or AR 0) is a flat Langmuir surface: no hole is solved there, its
hole coverages stay 0 and the open-area fraction is 0. The per-face values are written as ARField and
holeDensityField beside thetaFlat etc.; the scalar entries are written unchanged, so coupling/patterned_lib.py
reads old and new cases alike (it now also returns ARField and holeDensityField). Face-centre radii are the area
centroids of the wedge faces (0.67 mm for the innermost face, r + 0.003 mm elsewhere, times cos 2.5 deg), so a
table value at radius r lands within 0.1 % of r.
Regression: reactor/showerhead_pulse_layout_regress is a copy of showerhead_pulse_patterned_check (AR 50,
holeDensity 1, 1 s pulse, 2 s run) rerun with the rebuilt library. thetaBottom, thetaFlat, thetaMid, Keff and the
face value are identical to the original case to the 8 written digits at all 201 write times and all 150 faces
(maximum absolute difference 0). The uniform path is unchanged. 32 s wall time.
Two layouts, both copies of showerhead_pulse_patterned_AR100_long (1 % TMA, 300 sccm, 300 Pa, 15 mm gap, same
frozen flow in 0/U and 0/phi, s 5e-3, N0 4e18, c0 3.79e20, d 100 nm, 200 hole cells) with the pulse shortened to
16 s and the run to 17 s, writes every 0.05 s, 133 and 139 s wall time each:
(a) dense centre: AR 100 everywhere, 1.78e13 holes/m2 for r < 75 mm ramping linearly to 0 between 75 and 100 mm,
    bare wafer beyond (wafer-average open area 4.8 %, wafer-average wall-area ratio 19.2);
(b) AR ramp: 1.78e13 holes/m2 everywhere, AR 100 at the centre ramping linearly to 60 at r = 150 mm (area-averaged
    wall-area ratio 41.0).
Rule t = 1.2 (1 + A_wall/A_wafer) N0/(U c0) with U = 0.0501 m/s, c0 = 3.79e20 /m3 (monolayer supply time 0.211 s),
evaluated with the wafer-average ratio and with the local ratio of the face in question:
| layout | wafer-average ratio, rule | hole-bottom t90 centre (deepest holes) | worst t90 on wafer | rule (wafer avg) / worst | local rule at the worst face, ratio to actual | hole-bottom t90 at r = 149.5 mm vs local rule | flat t90 centre / edge |
| uniform AR 100 (reference run, 20 s pulse) | 55.9, 14.39 s | 14.05 s | 14.10 s (r = 128.5 mm) | 1.02 | 14.39 s, 1.02 | 13.30 s vs 14.39 s | 1.95 / 1.55 s |
| (a) dense centre | 19.2, 5.10 s | 11.85 s | 11.85 s (centre) | 0.43 | 14.39 s, 1.21 | bare (no holes); last holes at 99.5 mm: 5.10 s | 1.95 / 0.45 s |
| (b) AR ramp | 41.0, 10.62 s | 12.50 s | 12.50 s (centre) | 0.85 | 14.36 s, 1.15 | 7.60 s vs 8.76 s (AR 60.2, local ratio 33.7) | 1.95 / 1.55 s |
Within-layout profiles (figures/pattern_layout.csv): (a) hole-bottom t90 11.85 s at the centre, 11.60 at 25 mm,
10.75 at 50 mm, 8.75 at 74.5 mm (the last face at full density), 5.10 s at 99.5 mm (density 4.2e11 /m2); flat
t90 1.95 s at the centre, 1.35 s at 74.5 mm, 0.50 s at 100.5 mm, 0.40-0.45 s over the bare wafer beyond 110 mm.
(b) hole-bottom t90 12.50 s at the centre, 11.40 at 50 mm, 10.45 at 80 mm, 9.05 at 125 mm, 8.50 at 140 mm, 7.60 s
at 149.5 mm; flat t90 1.95 s everywhere except the outermost faces (1.55 s), as on the uniform wafer.
Reading:
- The rule with the wafer-average wall-area ratio fails for a non-uniform layout: it predicts 5.1 s for the
  dense-centre wafer whose array needs 11.85 s, and 10.6 s for the AR ramp whose centre needs 12.5 s. The showerhead
  supplies precursor per unit area, so the pulse time is set by the LOCAL wall-area ratio at the most demanding
  location, not by the wafer total. With the local ratio the rule is conservative by 15-21 % on both layouts
  (14.4 s predicted vs 11.85 and 12.5 s measured) against 2 % on the uniform wafer: the surplus is what the
  less-demanding neighbourhood feeds back laterally. For a tool recipe: size the pulse on the densest/deepest
  region's own wall-area ratio; the wafer average under-doses it.
- Open-field question, layout (a): the bare edge does not saturate at the blanket time. Its flat t90 is 0.50 s
  next to the array (r = 100.5 mm) and 0.40-0.45 s beyond 110 mm, against 0.35 s (centre) / 0.38 s (edge) on a
  blanket wafer: 15-30 % slower. The surface TMA density at 0.3 s is 0.163 c0 over the bare region (0.108 c0 at
  100.5 mm, 0.019 c0 over the array) against 0.24 c0 on a blanket wafer, i.e. the dense centre draws down the gas
  that would have reached the edge. The reverse transfer is larger in absolute terms: the array's hole bottoms
  saturate in 11.85 s at the centre and 8.75 s at its outer edge instead of 14.05 s on the uniform wafer, because
  the bare region leaves precursor that diffuses back inward (gap Peclet 0.09; lateral diffusion over 50 mm takes
  about L^2/D = 0.3 s, short against 14 s). A pattern boundary therefore produces a saturation gradient inside the
  array: 26 % between the centre and r = 74.5 mm, where the uniform wafer varies by 0.4 %.
- Layout (b): the bottom t90 tracks the local AR (12.5 s at AR 100 down to 7.6 s at AR 60); the flat regions are
  unaffected (1.95 s), because the open area and therefore the loading of the gas phase are the same everywhere.
  A tier-to-tier depth variation of 40 % across the wafer gives a 39 % spread in saturation time, so the recipe must
  be set by the deepest tier; shallower regions are over-dosed by exactly their depth ratio.
Caveats: 16 s pulse (all faces saturate before 12.6 s, so the shortened pulse does not truncate anything); t90
resolution 0.05 s; the face radius uses the wedge face centroid; the chamber floor outside the wafer is a
non-absorbing zeroGradient wall in both layouts as in all earlier runs.

## Night 3, part 2: sticking-coefficient sensitivity (reactor/sweep_sticking_sensitivity.py; figures/sticking_sensitivity.png, .csv)
Nine new runs at the baseline point (300 Pa, 573 K, 300 sccm, 1 % TMA, 15 mm gap, N0 4e18, c0 3.79e20, frozen flow)
with s = 5e-4, 1e-3, 2e-3, 8.5e-3 on the AR-100 patterned wafer (24 s pulse, endTime 26 s, writes 0.05 s; the 20 s
long case is the 5e-3 point) and s = 5e-4, 1e-3, 2e-3, 5e-3, 8.5e-3 on the blanket wafer (8 s pulse, endTime 9 s,
writes 0.01 s). Cases are copies of showerhead_pulse_patterned_AR100_long and showerhead_pulse_langmuir_p0.5 with the
wafer BC line, pulse table and endTime/writeInterval substituted (each regex asserted to hit once), under
reactor/sticking/. Wall time 84-96 s per blanket case, 197-206 s per patterned case, 21 min in all.
| s | blanket t90 centre / worst | patterned flat t90 | patterned bottom t90 centre / worst | bottom t90 / 11.99 s (worst) | utilization patterned / blanket |
| 5e-4 | 0.92 / 1.01 s | 6.00 s | 14.15 / 14.20 s | 1.184 | 84.4 % / 20.9 % |
| 1e-3 | 0.64 / 0.71 | 4.25 | 14.05 / 14.10 | 1.176 | 85.0 / 29.7 |
| 2e-3 | 0.47 / 0.53 | 3.00 | 14.00 / 14.05 | 1.172 | 85.3 / 39.7 |
| 5e-3 | 0.35 / 0.39 | 1.95 | 14.05 / 14.10 | 1.176 | 85.0 / 54.0 |
| 8.5e-3 | 0.30 / 0.33 | 1.55 | 14.05 / 14.10 | 1.176 | 85.0 / 63.8 |
Reading:
- The patterned-wafer pulse is independent of s across the whole published range: 0.15 s (three write intervals)
  of movement for a seventeen-fold change, tau_sat 1.17-1.18 throughout, utilization 84-85 %. The surface density over
  the wafer at 10 s is 0.154-0.162 c0 for all five values: once the wafer takes essentially the whole feed, the gas
  phase does not know the sticking coefficient. Inside the holes the bottoms fill earlier at low s (coverage at 10 s:
  0.384 at 5e-4, 0.017 at 8.5e-3) and the flat-to-bottom lag shrinks from 12.5 to 8.15 s, but the last tenth waits for
  the supply in every case.
- The blanket wafer scales as about s^-0.4 (worst 1.01 to 0.33 s) and the flat surface of the patterned wafer as
  about s^-0.5 (6.00 to 1.55 s), weaker than 1/s because the TMA density at the wafer rises as the uptake slows
  (centre, 0.1 s: 0.229 c0 at 5e-4, 0.065 c0 at 8.5e-3). At s = 5e-4 the uptake coefficient s vbar/4 = 0.051 m/s
  equals the inlet speed, so the surface sees about half the feed and the first-order time ln(10) N0/(s vbar/4 c0)
  = 0.47 s becomes 0.92 s. The 1/s extrapolation (3.5 s at 5e-4) was wrong by a factor 3.5; the 8 s pulse was
  longer than needed by a factor 8.
- Blanket 5e-3 with 0.01 s writes gives 0.35 / 0.39 s against 0.345 / 0.385 s with 0.005 s writes: one write interval.
- Literature check (full texts of Arts 2019 JVST A and Ylilammi 2018 JAP): Arts give TMA s0 = (0.5-2)e-3 at a 275 C
  set-point (about 220 C substrate) from the front slope, and tabulate the 300 C values 2-7e-3 (Ylilammi data, slope
  relation), 5.72e-3 (Ylilammi fit; their fitted partial pressure 147 Pa, K = 219 1/Pa) and 4 +/- 1e-3 (Vandalon and
  Kessels 2016, SFG). Our PillarHall fit (8.5e-3 x/÷ 1.5 at 300 C) is a factor 1.5-2 above the 300 C values, not
  4x above the field; the front slope fixes s and the penetration depth fixes exposure per site, so the unknown TMA
  pressure trades against N0, not against s. Arts' J. Phys. Chem. C papers (2019, 2021) are on O-atom recombination
  and give no TMA value.
Caveats: single operating point (the DOE showed tau_sat 1.03-1.25 at p >= 316 Pa and up to 1.96 at 100 Pa with
600 sccm, where carrier bypass, not sticking, is the cause); s held constant with coverage (first-order Langmuir);
t90 resolution 0.05 s patterned, 0.01 s blanket; the Cremers 2019 review and Yim 2020 full texts were not readable
through the tools available, so no numbers were taken from them.

## Night 3, part 3: fab units, purge and throughput (coupling/fab_units.py; figures/step_coverage_vs_pulse.png, .csv, step_coverage_table.csv, purge_times.csv, cycle_time_wph.png, .csv, wph_vs_wall_ratio.csv)
No new chamber runs. Inputs: the 1 % map runs at AR 20 / 50 / 100 (pulses 4.61 / 10.26 / 19.68 s, writes 0.023 /
0.051 / 0.098 s), the 20 s AR-100 run (0.05 s writes, purge to 22 s), the blanket TMA runs (0.5 s pulse to 1.5 s;
1.26 s pulse to 2.26 s), the water kdes 0 / 0.1 / 1 / 10 runs (5 s pulse to 8 s), `figures/pulse_time_map.csv`
(t_supply column, already (1 + A_wall/A_wafer) N0/(U c0)) and `figures/doe_gp_report.txt` (same-supply optimum
9.83 s). GPC 0.103 nm per cycle from the PillarHall plateau (night 1 part 10).
Step coverage SC = thetaBottom/thetaFlat at a face during the pulse, read as the coverage a pulse of that length
leaves (irreversible Langmuir; the chamber holds 0.019 c0 over the patterned wafer at the valve close, so the
post-pulse addition is negligible; the long run's 14.05 s and the map's 14.01 s agree to one write). Worst face over
faces 0-148 as in the map.
| AR | rule | t90 centre / worst | t95 | t99 | flat t90 / t95 / t99 | t99 / rule |
| 20 | 3.08 s | 2.67 / 2.71 | 2.78 / 2.81 | 2.92 / 2.94 | 1.84 / 2.14 / 2.60 | 0.96 |
| 50 | 7.32 s | 6.73 / 6.73 | 6.83 / 6.83 | 6.94 / 6.99 | 1.94 / 2.35 / 3.06 | 0.95 |
| 100 (map) | 14.38 s | 14.01 / 14.11 | 14.21 / 14.21 | 14.31 / 14.31 | 1.96 / 2.35 / 3.04 | 1.00 |
| 100 (20 s run) | 14.38 s | 14.05 / 14.10 | 14.15 / 14.20 | 14.25 / 14.30 | 1.95 / 2.30 / 3.05 | 0.99 |
Thickness after 50 cycles at the blanket pulse 0.345 s, bottom / top: AR 20 0.26 / 1.29 nm (SC 20 %), AR 50
0.01 / 1.32 nm, AR 100 0.00 / 1.44 nm; saturated 5.15 nm.
Purge (a), hole 100 nm x 10 um, TMA 573 K, vbar 410 m/s, D_K = d vbar/3 = 1.368e-5 m2/s: L^2/D_K 7.31 us; slowest
eigenmode 4 L^2/(pi^2 D_K) 2.96 us; bottom to 1 % by the first mode (4/pi) exp(-t/tau1): 14.36 us. Implicit 1D
transient (200 cells, 20 ns steps, n = 1 at t = 0, n(0) = 0, dn/dx(L) = 0): bottom to 1/e 3.70 us, to 1 % 14.42 us,
hole mean to 1 % 13.08 us. Hole clearing is 1.3e5 times faster than the chamber's 1 % time.
Purge (b), outlet patch mean after the pulse, fraction of the pulse-end value; tails fitted log-linearly over the
last 40 % of the written purge, extrapolated where the window ends first:
| run | outlet at pulse end | 10 % | 1 % | 0.1 % | tail decay (time constant) |
| TMA patterned AR 100, 20 s | 1.000 c0 | 1.10 s | 1.90 s | 2.73 s extr. | 2.75 1/s (0.36 s) |
| TMA blanket 0.5 s (1 s window; outlet still rising at pulse end, peak 0.367 c0 at +0.11 s) | 0.297 c0 | 1.25 extr. | 2.14 extr. | 3.04 extr. | 2.57 1/s (0.39 s) |
| TMA blanket 1.26 s (map, 1 s window) | 0.888 c0 | 1.07 extr. | 1.95 extr. | 2.83 extr. | 2.63 1/s (0.38 s) |
| H2O blanket kdes 0, 5 s | 0.999 c0 | 1.14 s | 2.20 s | 3.25 extr. | 2.18 1/s (0.46 s) |
References: gap transit 0.015 / 0.050 = 0.30 s; V/Q = pi 0.2^2 x 0.015 / 3.54e-3 = 0.53 s. The outlet time constant
sits between the two; the wafer flux itself clears within 0.05 s (part 1).
Purge (c), water outlet to 1 % of the pulse-end level (metric of coupling/postprocess_desorption.py): kdes 0 2.20 s
(run); 0.1 1/s 8.26 s and 1 1/s 4.33 s, extrapolated beyond the 3 s window as the desorption-fed level kdes theta
N0 A_wafer/(Q c0) (checked against the simulated outlet at 3 s: 0.0159 vs 0.0167 c0 and 0.0238 vs 0.0254 c0) with
theta decaying at the simulated rate over the last 2 s (0.088 and 0.654 1/s, below kdes because of re-adsorption);
10 1/s 2.40 s (run).
Cycle (TMA, purge 1.90 s, H2O, purge 8.26 s at kdes 0.1; in brackets kdes 0, 2.20 s): blanket 0.345 + 1.74 s pulses,
12.2 s (6.2); patterned AR 100 14.05 + 13.4 s, 37.6 s (31.6); same-supply optimum 9.83 + 9.38 s (H2O scaled by
9.83/14.05, not in the DOE; purges kept at baseline), 29.4 s (23.3). Caveat on the optimum: at 150 sccm, 631 Pa,
10 mm the volumetric flow is 8.4e-4 m3/s, the showerhead speed 0.012 m/s, gap transit 0.84 s and V/Q 1.5 s, 2.8
times the baseline; a TMA purge scaling with those would be about 5.3 s and the cycle about 32.8 s (estimate, not
run); the desorption-fed water purge scales with the precursor supply, held within 25 % of the baseline.
5 nm = 49 cycles, 10 nm = 98 at 0.103 nm per cycle. Minutes per wafer, 5 / 10 nm: 10.0 / 20.0, 30.7 / 61.4,
24.0 / 48.0. Wafers per hour at 5 nm, no overhead / 60 s handling: 6.0 / 5.5, 1.95 / 1.9, 2.5 / 2.4; at 10 nm
3.0 / 2.9, 0.98 / 0.96, 1.25 / 1.2.
Rule as throughput, 5 nm, both half-cycles at the same supply flux (constants of coupling/temporal_vs_spatial.py),
purges 1.90 + 8.26 s, 60 s handling: 1 % TMA 5.9 WPH at wall ratio 1, 4.2 at 11.2, 2.8 at 28, 1.8 at 55.9, 1.1 at
112; 3 % TMA 6.3, 5.5, 4.5, 3.5, 2.4. The AR-100 chamber run (37.6 s, 1.9 WPH) lies on the 1 % curve (39.3 s,
1.8 WPH).
Caveats: the purge criterion is 1 % of the pulse-end outlet level, a choice; kdes is a parameter; the blanket TMA
purge numbers are extrapolations from 1 s windows; the DOE optimum's purges are taken from the baseline chamber;
the time resolution of the step-coverage table is one write interval; no gas-phase reaction, no dead volumes, no
downstream plumbing, which in a real tool lengthen the purge.

## Night 3, part 4: chamber-scale validation (reactor/validation_gakis, analysis/blanket_saturation_literature.py; figures/blanket_saturation_literature.png, .csv, figures/validation_gakis.png, .csv)
(a) Blanket saturation against published data. No new chamber run. From reactor/showerhead_pulse_langmuir_p0.5:
feed exposure p_feed t90 = 3.00 Pa x 0.345 s = 1.03 Pa s (7,760 L) at the centre, 1.16 Pa s at r = 149 mm; surface
exposure (integral of the face TMA partial pressure to theta = 0.9) 0.146 Pa s (1,090 L) at the centre and 0.146 at
the edge, equal to ln(10) N0/(s F1) = 0.142 Pa s (F1 = 1.30e22 /m2/s/Pa at 573 K). Literature exposures compiled with
their DOIs in the csv: applied exposures 0.56 (Gakis, impingement from their CFD), 7.5 (our PillarHall fit), 14.7
(Ylilammi, fitted entrance pressure), 25 and 160 Pa s (Arts, stated doses); kinetic E90 implied by stated s or sigma:
0.05 (Soto and Tysoe, 25 C, via Puurunen 2005), 0.08 (PillarHall fit, s = 8.5e-3), 0.16 (Ylilammi, gamma = 5.72e-3,
q = 5e18), 0.21 (Gakis, s0 = 0.004 x 0.86), 0.27 (Vandalon and Kessels, sigma = 6.5e-18 cm2, 300 C; 0.22 at 100 C),
0.33 to 1.3 (Arts, s0 = 2e-3 to 5e-4, 220 C); Dillon 1995 via Puurunen (k = 0.0076 /s at 0.01 Torr, 227 C) gives
404 Pa s on a porous membrane and is transport-limited, listed but excluded. Values read from figures: Gakis 2018
Fig. 12 colour scale (0.0118-0.0130 mol/m2, +/- 0.0002), Gakis 2019 Fig. 10 squares (digitised by marker centroids,
+/- 0.003 A/cycle; digitised non-uniformity 2.37 % vs 2.40 % stated), Gakis 2018 Fig. 8 pulse shape (+/- 20 % per
ordinate, integral renormalised). Not obtained: Ott 1997, Groner 2004, Elam 2002 full texts (not open; secondary
sources give pulse times without partial pressures). Reading: our surface exposure sits in the middle of the 300 C
cluster (0.08-0.27 Pa s), within a factor 1.8 of each explicit-sticking source; the feed exposure is 7x larger because
the blanket wafer draws the surface density down to 0.14 c0 while it saturates (the supply-limited chamber), and
literature feed exposures overstate surface exposures for the same reason. s needed to match each source: 2.7e-3
(Vandalon), 3.4e-3 (Gakis), 4.7e-3 (Ylilammi), 2.2e-3 to 5.5e-4 (Arts), 8.7e-3 (PillarHall fit); s = 5e-3 stays.
(b) Gakis et al. Fiji F200. Geometry measured on the 2018 Fig. 5a cross-section (1.65 px/mm on the 200 mm chuck top,
+/- 3 mm): inlet tube r = 20 mm, cone from r = 20 mm at +199 mm to r = 128 mm at +35 mm, wide chamber r = 128 mm to
-50 mm, chuck cylinder r = 100 mm to -50 mm (real dome 124 mm deep), duct r = 106 mm to -205 mm, outlet on its floor.
One axial inlet for their three (top 100, door 50, side 30 sccm), no loading-door box, carrier flow frozen, TMA a
passive scalar (its inlet partial pressure peaks at 11.5 Pa against 9.6 Pa of Ar: the real pulse is a slug of
nearly pure TMA; their outlet pressure rises 18 %), isothermal 573 K, D = 0.263 m2/s (N2 value scaled 1/p),
lambda = 1.4 mm with 2 mm first cells. Q = 0.0603 m3/s (180 sccm at 24 C, 1 bar -> 9.6 Pa, 573 K), U_in = 48 m/s,
Re 4.3, nu = 0.447 m2/s, residence 0.23 s. Mesh 13,020 cells, max non-orthogonality 56 deg, max skewness 2.5
(checkMesh OK). simpleFoam: the showerhead template's consistent SIMPLE at relaxation 0.9 diverged at iteration 98
(11 s); standard SIMPLE with one non-orthogonal corrector, p 0.3 / U 0.7, converged in 2000 iterations (Ux residual
3.8e-7, p 2.5e-8, ClockTime 77 s); inlet and outlet flux 0.0602 m3/s. Pulses: s = 3.44e-3 (their s0 = 0.004 times
the 86 % reaction probability of adsorbed TMA at 300 C), N0 = 4.1e18 /m2 (their 1.0 A/cycle at 3500 kg/m3), vbar
410 m/s; 25 ms amount 8.74e-6 mol (their 0.63 mg) with their Fig. 8 shape (peak at 12 ms, 5 % by 50 ms, tail to
0.3 s), 60 ms at 1.5x and 2.4x (not quantified in the paper), one case at s = 5e-3; Danckwerts inlet fraction
f = 48/(48 + 128) = 0.27; dt 0.4 ms to 1.0 s, writes every 10 ms; 104 to 117 s wall time each (1.7 to 1.9 min,
sharing the CPU with another solver; a first attempt at dt 0.2 ms to 1.5 s was stopped at 0.25 s after 3 min).
| case | theta after the cycle (centre / min / max over |r| <= 95 mm) | non-uniformity | t90 centre / worst | wafer impingement exposure | exposure non-uniformity | captured | peak surface TMA |
| 25 ms, s = 3.44e-3 | 0.9994 / 0.9994 / 0.9994 | 0.0007 % | 0.07 / 0.08 s | 0.671-0.673 Pa s | 0.28 % | 2.4 % | 5.6 Pa at 50 ms |
| 60 ms x1.5 | 1.0000 | 0.00001 % | 0.06 / 0.06 s | 1.017-1.019 | 0.18 % | 1.6 % | 8.5 Pa |
| 60 ms x2.4 | 1.0000 | 0 % | 0.05 / 0.05 s | 1.639-1.640 | 0.11 % | 1.0 % | 13.7 Pa |
| 25 ms, s = 5e-3 | 1.0000 | 0.00003 % | 0.06 / 0.07 s | 0.671-0.673 | 0.28 % | 2.4 % | 5.6 Pa |
| measured (2019 Fig. 10 / stated) | 25 ms: 0.996 / 0.975 / 0.998 of the 60 ms value; 60 ms: 1.00 | 2.40 % / 0.58 % | | their CFD 0.55-0.60 Pa s | 10.7 % | | 0.7-3.7 Pa at 25 ms (2019 Fig. 7b) |
Reading: exposure to the wafer per pulse agrees with their 3D CFD to 12-22 %; saturation state agrees (25 ms at the
knee, 60 ms past it); non-uniformity not reproduced (0.0007 % vs 2.4 %): their recirculating off-axis injection and
cold loading door give a 10.7 % exposure spread we cannot carry axisymmetrically, and irreversible Langmuir turns even
that into < 0.5 % at 3.2 E90, so their reversible kinetics, the 21 C wafer gradient and the H2O purge residue (the
0.58 % at 60 ms is a saturated-TMA value) carry the rest. By 1/s, s_eff = 1.7e-3 to 3.0e-3 would put the 25 ms pulse
at the measured 0.975 to 0.998. Caveats: dimensions from a figure; frozen flow during a pulse that is 50 % of the
gas inventory; passive-scalar TMA above the carrier pressure; isothermal; chuck dome truncated; 60 ms amount assumed.

## Night 3, part 6: 3D chamber (reactor/chamber3d/make_and_run.py, coupling/postprocess_chamber3d.py; figures/chamber3d_t90_map.png, chamber3d_azimuthal.png, chamber3d_t90.csv, chamber3d_summary.txt)
Full 3D cylinder, r 0-200 mm, gap 15 mm, axis along x as in the wedge so the patternedSurface face radius
sqrt(y^2 + z^2) and the wedge conventions carry over. blockMesh O-grid, nine blocks (central square of half side
75 mm with its edges bulged to r = 82 mm, four curved blocks to r = 150 mm, four ring blocks to 200 mm, arc edges on
both circles): 48 cells per quarter arc, 20 + 14 radial cells, 12 across the gap; 105,984 hexahedra, max
non-orthogonality 30.7 deg (average 5.1), max skewness 0.794, max aspect ratio 5.04, checkMesh OK. Patches:
showerhead and wafer (r < 150 mm, 6144 faces each, the wafer as one patch), topWall and floor (2688 each), sideWall
(1920), and outlet = the 60 deg sector of the side wall centred on +y, 384 faces, carved with topoSet (patchToFace
sideWall, subset boxToFace |z| < 100 mm, y > 170 mm; the sector boundary falls on face edges with 1.875 deg cells) and
createPatch -overwrite. Steady flow: simpleFoam with the wedge's settings and BCs (uniform -0.050 m/s inflow, noSlip
walls, totalPressure 0 at the port), converged by the residual controls in 93 iterations, 18 s wall (last initial
residuals Ux 9.6e-8, Uy 1.0e-9, Uz 8.6e-9, p 4.1e-8); showerhead in 3.53366094e-3 m3/s, port out 3.53366097e-3 m3/s;
mean port speed 1.125 m/s. Gap-averaged outward velocity at r = 140 mm: port side +0.854 m/s, 90 deg +0.073, far
side +0.013 (axisymmetric 0.268); at r = 100 mm the far side flows inward, -0.081 m/s. Pulse: scalarTransportFoam on
the frozen flow, 0/T as showerhead_pulse_patterned_AR100_long (uniformMixed Danckwerts inlet with f recomputed for
delta 0.625 mm: 3.706e-3; patternedSurface wafer, s 5e-3, N0 4e18, c0 3.79e20, d 100 nm, AR 100, 1.78e13 /m2, 200
hole cells; inletOutlet port; zeroGradient elsewhere), 16 s pulse, 17 s run, dt 1 ms, writes every 0.1 s (25 MB
each, 3.4 GB; the wedge wrote every 0.05 s), 1939 s wall serial (a 0.2 s timing copy took 26 s). Face centres from
postProcess -func writeCellCentres -time 0 (boundary values of C). t90 = 0.9 crossing interpolated between writes
(first-write values kept in the csv); the wedge reread the same way gives 14.007 s centre (14.05 first write),
14.051 s worst at r = 130.5 mm (14.10), 13.28 s outermost face, flat 1.923 s.
Hole-bottom t90 by radius (faces within +-2.5 mm; 15 deg bins detrended to the exact radius):
| r (mm) | mean (s) | min / max (s) | azimuthal p2p (%) | port - far (%) | wedge (s) | 3D / wedge |
| 50 | 14.048 | 14.042 / 14.055 | 0.09 | -0.08 | 14.009 | 1.0028 |
| 100 | 14.075 | 14.041 / 14.098 | 0.37 | -0.36 | 14.021 | 1.0038 |
| 140 | 14.089 | 13.975 / 14.156 | 0.96 (raw band 1.29, of which 0.58 is radial) | -0.90 | 14.006 | 1.0060 |
| 148.9, outermost row | 13.687 | 13.544 / 13.931 | 2.83 | -0.70 | 13.280 (149.5 mm) | 1.031 |
Innermost face 14.044 s (first write 14.1) vs wedge 14.007 (+0.26 %, the 12-cell gap); worst face 14.168 s (14.2) at
r = 132.7 mm, azimuth -160 deg, vs wedge worst 14.051 (+0.83 %); wafer-wide spread 4.44 % of 14.059 s; face-by-face
3D minus wedge at the same radius +0.44 % mean, -0.09 to +1.37 %; faint O-grid block imprint of about 0.1 % in the
deviation map. Flat t90: 1.929 s at r = 50 mm (p2p 0.02 %), 1.934 at 100 mm (0.46 %), 1.983 at 140 mm (4.75 %:
1.917 port side, 2.011 far side; wedge 1.969). Gas side at r = 140 mm, port vs far: value/c0 0.0186 vs 0.0172 at
0.3 s (p2p 9.2 %), 0.128 vs 0.126 at 8 s, 0.253 vs 0.237 at 14 s; flat coverage at 0.3 s 0.213 vs 0.191; bottom
coverage 0.0084 vs 0.0075 at 8 s, 0.899 vs 0.835 at 14 s (p2p 10.8 %; 0.7 % at r = 50 mm). Reading: the one-sided
pump moves the gas (three times the axisymmetric outward speed on the port side, almost none on the far side) and
the instantaneous surface density and coverage follow it at the 10 % level, but the hole-bottom saturation time is
the integral of a locally delivered, locally captured supply across a Pe 0.09 gap, so the azimuthal lag is 1 % at
r = 140 mm and the worst point moves to the far side (r 133 mm, 160 deg from the port) at +0.8 % of the wedge value,
0.3 % of it mesh. The flat surface, saturating at 2 s under a strongly drawn-down gas phase, is the more sensitive
indicator (4.8 % spread at r = 140 mm). Caveats: 12 gap cells; writes 0.1 s apart (interpolated); frozen flow; a plain
sector port at uniform pressure; non-absorbing floor and side wall; Stokes flow, so only the Peclet numbers change
the result with operating point.

## Night 3, part 5: geometry and chemistry (reactor/chemistry/make_tdmah_case.py, coupling/tapered_hole.py, coupling/wordline_recess.py, src/aldconf/recess.py)

Tapered holes: TaperedHole(diameter, length, bottom_ratio) in aldconf.geometry; the solver uses the conservative
(1/A) d/dx (A D dn/dx) form on the frustum cross-section when the geometry is tapered and the previous uniform-area
form otherwise (film narrowing included), so no earlier number changes. Checks: zero sink gives a uniform density
to 1e-9 in a 0.5 frustum; bottom_ratio 1 is bitwise the Hole; the area form at bottom_ratio 0.999 agrees with the
cylinder within 0.01 in coverage; frustum area matches a numerical integral to 1e-4. AR-100 holes (100 nm mouth,
10 um) with 100 / 90 / 80 / 60 nm bottoms driven by the centre-face value of showerhead_pulse_patterned_AR100_long
times c0 (aldconf.simulate_driven, 5 substeps per 0.05 s write), s 5e-3, N0 4e18:
| bottom | wall area / straight | driven bottom t90 | constant 0.3 Pa t90 | rule 1.2 (1 + 55.9 A_rel) N0/(U c0) |
| 100 nm | 1.000 | 14.03 s | 15.93 s | 14.39 s |
| 90 nm | 0.950 | 14.21 s | 16.52 s | 13.68 s |
| 80 nm | 0.900 | 14.36 s | 17.23 s | 12.98 s |
| 60 nm | 0.800 | 14.62 s | 19.12 s | 11.56 s |
(Chamber two-way bottom t90 at the centre: 14.05 s at the first write above 0.9, 14.01 s interpolated.) 8 s wall
time for the script. Reading: supply-limited wafer, taper nearly invisible (1-4 %); feature-only penalty 4-20 %; the
rule's area reduction goes the other way (5-20 % shorter) and the two bracket a tapered-hole wafer.
Word-line recesses: aldconf.recess, RecessedSlit(w_s, L_s, g_r, L_r, p_r), 200 slit cells x 50 recess cells,
vectorised Thomas solves per substep; limits tested (plain slit to 1e-9; AR-1 recesses vs area-equivalent slit
within 3 % on the front). Dimensions: g_r 29 nm, p_r 51 nm (Lee and Shin 2022), L_s 5 um (98 pairs, Reiter 2024),
w_s 250 nm (Onto 2023) and 100 nm, L_r 365 nm (TechInsights 2014: 730 nm slit-to-slit) and 1 um (Onto 2023: slit
pitch > 2 um). Drive A = centre-face TMA density of the AR-100 long case (one-way); drive B = constant 0.3 Pa:
| w_s | L_r | A_f | A: plain slit | slit bottom | recess bottoms mouth / mid / deep | lag | area-eq slit | B: slit bottom / last recess; exposure | effusion limit |
| 250 nm | 365 nm | 14.7 | 2.30 s | 5.37 s | 2.09 / 4.46 / 5.46 s | 0.09 s | 5.39 s | 2.70 / 2.79 s; 0.84 Pa s | 0.61 s |
| 250 nm | 1 um | 39.6 | 2.30 s | 8.33 s | 2.94 / 6.94 / 9.05 s | 0.72 s | 8.70 s | 6.02 / 7.03 s; 2.11 Pa s | 1.63 s |
| 100 nm | 365 nm | 14.7 | 3.76 s | 13.46 s | 2.11 / 8.55 / 13.53 s | 0.07 s | 13.49 s | 14.53 / 14.68 s; 4.40 Pa s | 1.52 s |
| 100 nm | 1 um | 39.6 | 3.76 s | 16.52 s | 2.97 / 13.82 / 16.73 s | 0.20 s | 16.66 s | 37.2 / 39.2 s; 11.7 Pa s | 4.08 s |
Wall-area ratio of a wafer tiled at pitch w_s + 2 L_r: 150 / 176 / 178 / 189 (fully patterned); rule 38.3 / 44.8 /
45.2 / 48.0 s at 1 % TMA (19.3-24.1 s at 50 % patterned). 16 s wall time for the script. Reading: the recess bottoms
lag the slit bottom by 0.5-9 % and the structure is an area-equivalent slit within 4 %; the rule's area-only treatment
holds to recess AR 35; the chamber loading of such a wafer (not run; no recessed-slit BC) would be 2.7-3.4 x the hole
wafer's.
TDMAH chamber run (reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah, copy of the AR-100 long case with
vbar 185, s 7e-3, N0 2.77e18, DT 3.62e-3, f 3.445e-3; 20 s pulse, 22 s, writes every 0.05 s): 132 s wall time;
bottom t90 11.60 s centre / 11.60 worst / 11.10 at r = 149 mm, flat 1.30 s, surface density at 5 s 0.197 c0; rule
9.96 s, ratio 1.165. Sensitivity ..._tdmah_Dtma (DT 8.4e-3, f 1.489e-3, TDMAH kinetics): 185 s; bottom t90 11.10 /
11.15 / 10.95 s, flat 1.20 s, ratio 1.115. TMA reference from the same post-processing: 14.05 / 14.10 / 13.50 s, flat
1.95 s, rule 14.38 s, ratio 0.977. Reading: the chemistry enters through N0 (and c0); the lower TDMAH diffusivity
(gap Peclet 0.21) adds 0.5 s, the shorter supply time against the chamber's fixed overheads the remaining 1.1 s.
573 K kept for comparability (upper edge of the TDMAH window); rule at 523 K 9.26-9.96 s (U c0 temperature-
independent at fixed sccm and fraction). Rule-only: TiN/TiCl4 (s 0.10, GPC 0.30 A, 4.23 g/cm3, N0 1.24e18) 4.44 s;
W/WF6 (s 0.40, GPC 2.5 A, 19.25 g/cm3, N0 1.58e19) 56.7 s; per nanometre of film HfO2 100, Al2O3 144, TiN 148,
W 227 s/nm (proportional to the film's metal-atom number density). Caveats: TDMAH sticking is a model calibration for an undisclosed Hf amide (Reiter 2024), not a
measurement; TiCl4 sticking is for TiCl4/H2O on TiO2 at 110 C, not TiCl4/NH3; WF6 sticking is on the Si-terminated
surface, not the B-terminated one; film densities for ALD HfO2 (bulk used) and ALD W (bulk used) are not measured
values; N0 therefore carries the GPC x density uncertainty (factor 1.1-1.3). Tests: 86 pass (9 new) in 14 s.

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

## Night 4, part 3: two-species kinetics (src/aldconf/multispecies.py, coupling/multispecies.py, lib/aldkin 1.1; figures/multispecies.png, .csv)
Physics first: at 100-1,000 Pa the mean free path is 123-12 um (N2, 3.8 A, 573 K), Kn = 1,233-123 in a 100 nm hole,
so a byproduct cannot block the precursor in the gas phase; each species diffuses with its own Knudsen diffusivity
(CH4 2.12 x TMA, HCl 2.28 x TiCl4) and a byproduct acts only at the surface (re-adsorption, etching). The feature
solver gained a two-species variant: precursor P (irreversible, sticking on the free fraction 1 - theta_P - theta_B),
byproduct B released at nu per adsorbed P, quasi-steady gas for both (two tridiagonal solves per substep, the second
with the wall source nu x uptake + k_B N0 theta_B on the right-hand side; _solve_gas gained an optional source
argument, None leaves the earlier arithmetic unchanged), re-adsorption s_B, desorption k_B, etch probability k_etch
per collision on a coated site, byproduct mouth density as a fraction of the precursor's; the coverage pair advanced
per substep by the exact closed-form solution of its linear 2 x 2 system at frozen gas (eigenvalue of smaller
magnitude from det/l2, no cancellation at k_B -> infinity). Limits tested: s_B = 0 or nu = 0 reproduces simulate
(difference 0), k_B = 1e12 to 1e-9, the driven variant equals simulate_driven to 1e-12; the byproduct balance
produced + desorbed = escaped + re-adsorbed + etched closes to 1e-16 to 1e-14 (gas side, frozen coverages); the
surface inventory agrees with it to first order in the substep, and the HCl rates (a_B up to 50 1/s) set the substep
at 1 ms (t90 10.111 / 10.074 / 10.059 / 10.052 s at 4 / 2 / 1 / 0.5 ms). Stoichiometries: 3 CH4 per TMA per cycle
(Puurunen 2005); 4 HCl per TiCl4 per cycle (Ritala 1995), 0.5 in the TiCl4 half-step (Cl:Ti 4 -> 3.5, Tiznado and
Zaera 2006); nu 4 run as the bound. HCl kinetics: no value on TiN or TiO2 at 573 K; s_B 1e-3 (reactive sticking on
alpha-Al2O3(0001) at 298 K, Elam 2000) to 0.4 (Knudsen-cell initial uptake 0.38 on alumina at 296 K, Lostier 2026);
k_B 0.1-1,000 1/s from the 300-650 K TPD span of Elam 2000 (E_des 0.8-1.8 eV, k 1e-3-1e6 1/s at 573 K); k_etch 0
with bounds 1e-4, 1e-3 (HCl does not etch TiO2, Leem 2014; Ti removal by surface Cl in the TiN kMC of Kim 2022, no
rate). Results, AR-100 hole: (a) TMA/CH4 at s_B 0 identical to the single-species solver (15.600 s constant 0.3 Pa,
14.005 s chamber-driven), the hypothetical bound s_B 1e-3, k_B 1e3 +0.01 % with 99.95 % of the methane escaping;
(b) TiCl4/HCl constant 0.3 Pa: 7.754 s without HCl action; k_B >= 100: 0 to +0.5 %; k_B 1: 8.59 / 8.65 / 10.06 /
12.00 s for s_B 1e-3 / 1e-2 / 0.1 / 0.4 (nu 0.5), 8.41-20.39 s (nu 4); k_B 0.1: 23.5-27.7 s (nu 0.5), no saturation
in 30 s at nu 4 and s_B >= 0.1; etch bounds 1e-4 / 1e-3: +20-35 % / +70-80 % at nu 0.5, no saturation at nu 4
(etch-limited steady state); (c) under the chamber's centre-face history (the TMA run as a stand-in; only changes
are meaningful): 9.622 s without HCl action, +11 / +22 / +36 % (nu 0.5, k_B 1, s_B 1e-2 / 0.1 / 0.4), +46 / +74 %
(nu 4, s_B 0.1 / 0.4), 0 to +2.3 % at k_B 100, etch bounds +30 / +43 %; the relative cost is 20-40 % below the
constant-mouth one (the rising late-pulse density refills released sites); profiles at the rule pulse show a
broader, deeper front with blocking and the blocked coverage peaking ahead of it; (d) HCl at the mouth at 1 / 3 /
10 % of the precursor density: +0.5 / +1.6 / +5.1 % constant and +0.4 / +1.2 / +4.0 % driven at k_B 1, nothing at
k_B 100, nothing for the TMA bound. Library: aldkin 1.1.0 with nu, vbarB, sB, kB, kEtch, cB0 appended to
aldkin_params (defaults 0 = the single-species path), blocked-coverage selectors, aldkin_set_byproduct_mouth,
aldkin_byproduct_flux, aldkin_byproduct_balance, aldkin_hole_byproduct_gas; exponential scheme uses the same
closed-form 2 x 2 advance, Euler one forward step of both; 5 new ctest cases (16 passed): inactive fields ==
bit-identical for both models and schemes, k_B 1e14 within 1e-9, balance to 1e-9 per step, blocking and etching
ordering with Euler/exponential within 0.01, error paths; the binding exposes the fields and 3 new tests compare the
library with aldconf.multispecies within 1e-9 over 2,000 steps (a free-fraction inconsistency in the library's
precursor sink was caught by this test and fixed; the single-species path was not affected). pytest: 106 passed
(95 + 8 + 3; tests/test_voxel_mc.py excluded while the other package builds). ./reproduce.sh --regress through the
launcher with the 1.1 library: ctest 16 passed (aldkin) and 10 passed (aldmc), wmake, steady flow, regression pulse,
"3 passed" (the boundary conditions run the single-species path unchanged). The boundary conditions do not read
the new fields yet (next step: dictionary keywords and a transported HCl field). Caveats: every HCl kinetic number
is a range or a bound, labelled in the tables; the TiCl4 sticking is the TiCl4/H2O value; the chamber-driven TiCl4
runs use the TMA run's density history; the etch product is not tracked; byproduct produced on the flat wafer
surface is not fed back to the feature unless mouth_fraction > 0.
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

## Night 4, part 4: slit transport correction (src/aldconf/transport.py slit_knudsen_factor, geometry.py Slit.knudsen_correction, recess.py; coupling/fit_pillarhall_corrected.py, coupling/wordline_recess_corrected.py; figures/pillarhall_fit_corrected.png, .csv, wordline_recess_corrected.csv, .png; tests/test_slit_transport.py)
The 24 % slit error of WP9 (f) is the parallel-plate Knudsen diffusivity: D_K = d_h vbar/3 with d_h = 2 w is exact for
a cylinder and not for plates, where grazing flights give D_slit = (w vbar/4)[ln(L/w) + c], i.e. (3/8)[ln(L/w) + c]
times D_K (Clausing slit transmission; Berman 1965; Steckelmacher 1986). A transmission Monte Carlo of the infinite
slit (4e5 particles per point) gives W = 0.6844 / 0.5425 / 0.3572 / 0.2399 / 0.1545 / 0.0802 / 0.0472 / 0.0270 / 0.0124
at L/w = 1 / 2 / 5 / 10 / 20 / 50 / 100 / 200 / 500, so W L/w - ln(L/w) = +0.10 from L/w 10 to 200. Machball
Trench(AR, 8 AR) at s 5e-3, N0 4e18, 573 K, 0.3 Pa, 100 nm gap: t90 0.5335 / 0.6231 / 1.1170 / 2.6305 s at AR 10 / 20 /
50 / 100 (Nz convergence 0.2-0.3 %) against the 1D slit model 0.5087 / 0.6149 / 1.3842 / 4.2305 s (bottom coverage
0.9, 200 cells, 200 substeps; 400/400 moves AR 50 by 0.2 %); the knudsen_factor that matches Machball is 0.59 / 0.95 /
1.405 / 1.721 (AR 10-20 reaction-limited, kinetic floor 0.47 s, factor ill-conditioned) against (3/8) ln(AR) = 0.863 /
1.123 / 1.467 / 1.727; the 4 % at AR 50 is the continuum-vs-ballistic offset of WP9 (a). Implemented as
Slit(knudsen_correction="plates") (keyword-only, default "none" reproduces every earlier number), applied in simulate,
simulate_driven, multispecies and simulate_recessed (slit and recesses separately), factor (3/8) ln(L/w) without
the constant: corrected 1D t90 0.5143 / 0.5993 / 1.0891 / 2.6234 s, -3.6 / -3.8 / -2.5 / -0.3 % against Machball (with
c = 0.10: -3.9 / -4.5 / -3.9 / -2.0 %). Front shape (Machball Trench vs 1D at equal times): AR 500, s 8.5e-3, 0.8 -> 0.2
width 40-48 gaps ballistic over z50 60-187, 42.1 corrected (6-13 % narrow), 27.6 uncorrected (31-43 % narrow); front
position corrected 14-30 % deep, uncorrected 15-25 % shallow. PillarHall refit (same data, model, nine starts, per-sheet
refits; factor 2.330 at AR 500): s = 1.96e-2, 32.3 Pa, 3.23 Pa s per pulse, rms 0.0381 (unchanged), sheet-to-sheet
s = 1.99e-2 x/÷ 1.46 (geometric mean of the eight sheets; range 7.9e-3 to 2.7e-2), p = 32.4 Pa x/÷ 1.02; the nine starts agree within 1 %. s rises by the factor D rose and the pressure falls by it; against the 300 C
literature (5.72e-3, 4e-3) the corrected value is 3.4-4.9 x high (was 1.5-2.1 x), i.e. the correction moves the fit
away from the literature; the Arts slope relation itself is a hydraulic-diffusivity continuum result, so 8.5e-3 is the
like-for-like number against it; a fully ballistic fit would give about 2.6e-2; the ideal-plate factor is an upper bound
for the 140-200-gap-wide, pillar-supported cavity (duct Monte Carlo with side walls: 1.94 / 1.89 / 1.70 / 1.49 at b/w
200 / 140 / 50 / 20 against 2.17 infinite at L/w 250; 2.21 / 1.90 / 1.88 at b/w 200 / 140 / 50 vs 2.30 at L/w 500). Word-line recesses
restated (factors 1.123 / 1.467 slit, 0.950 / 1.328 recess): chamber-driven slit bottom 5.09 / 7.96 / 10.99 / 15.36 s
(was 5.37 / 8.33 / 13.46 / 16.52), last recess +0.10 / +0.53 / +0.08 / +0.13 s, plain slit 2.26 / 3.27 s, area-equivalent
within 0.5-3 %; constant 0.3 Pa w100_Lr365 10.03 / 10.17 s vs voxel 10.70 / 10.88 s (6 % fast; was 36 % slow); rule
and wall ratios unchanged. Everything else that used Slit: voxel_mc_crosscheck part slit/recess references (stay as
the uncorrected values they measure), docs/build_assets.py PillarHall panel (not touched), test_scaling and
test_recess (transport-independent). pytest 124 passed in 32 s (118 + 6). Wall time: Machball runs 1-20 s each,
transmission Monte Carlo 2 min, recess script 10 s, refit about 25 min single core. Caveats: the correction is for an
ideal, infinitely wide slit closed at one end with the cut-off at the channel length; finite width, pillars or an
open far end change the cut-off; it is a few per cent fast against ballistic solvers at AR 20-50 (the diffusion model's
own offset) and 14-30 % deep in front position at AR 500; the factor falls below 1 for AR < 14.4, where the slit is
reaction-limited and the factor hardly matters; no narrowing in the Machball comparisons.

## Night 4, part 5: device geometries (coupling/device_geometries.py; figures/device_geometries.png, gaa_coverage.png, nand_cell_coverage.png, device_geometries.csv, device_geometries_recess_profile.csv)
Two builders on top of aldconf.voxel.from_array: gaa_stack (N sheets t_s at spacing g_s, width W, under a gate trench
of length Lg closed by the spacers, H_top above the top sheet, g_bot to the substrate, periodic along the width at
W + gap_x; plain voxel faces, exact for the axis-aligned cell, inventory 1.000) and nand_cell (voxel.recessed_slit's
slit and recesses plus channel-hole pillars as solid cylinders of round(pi r^2) voxels per slice in staggered rows at
the hole pitch along the slit; area correction with a 2-voxel normal kernel so the 5-6-voxel recesses are not blended
across, inventory 0.991 / 0.994). Cells: gaa_ibm2017 (3 x 5 nm at 10 nm, W 40, Lg 12, gap 20, H_top 40, g_bot 10;
h 1.25 nm; Loubet 2017 via Wong and Kakushima 2022 for N, t_s, g_s, Lg; W from an IBM figure caption; gap, H_top,
g_bot assumed), gaa_scaled4 (4 x 5 nm at 7 nm, W 25, h 0.875), gaa_ibm2017_after_hk (3 nm removed from every
surface: 4 nm slots, 6 nm trench, h 0.5), nand_wp4_pillars (wp4's w100_Lr365 at h 5 nm with 2 rows of 100 nm holes at
200 nm pitch), nand_2xx (250 nm slit, 12.7 um, 48 nm pitch split 30/18, 1 um recesses, 4 rows of 100 nm holes, h 6 nm;
TechInsights YMTC 232L table for depth and pitch, Onto 2023 for slit, slit pitch and hole; split, hole pitch and rows
assumed). Regions: sheet tops, slot centres (central 20 % of W) and edges (outer 20 %), floors; slit floor, deep halves
and end walls of the first, middle and deepest recess, and the minimum over all recesses' deep halves. Chemistries TMA
5e-3 / 4e18 / 410 m/s, TDMAH 7e-3 / 2.77e18 / 185, TiCl4 0.1 / 1.24e18 / 253; drives constant 0.3 Pa and the long
AR-100 case's centre-face history; GAA dt = t_flat/100, M 50 000; NAND dt 0.1 / 0.25 s, M 100 000 / 60 000; stop at
99.5 % everywhere or the end of the history; wall-clock cap 18 min never reached.
GAA (everywhere t90 / t99 at 0.3 Pa; flat Langmuir 0.474 / 0.947 s TMA): IBM cell TMA 0.558 / 1.038 s, TDMAH 0.649 /
1.179, TiCl4 0.056 / 0.070 (exposure 0.17 / 0.20 / 0.017 Pa s); order top sheet < slots from the top down < bottom slot
centre; centre-edge lag 1 % TMA, 8 % TiCl4; four-sheet 7 nm cell +2-3 %; narrowed cell TiCl4 0.142 s (28 % centre lag),
TDMAH 0.898 s. Wall-area ratio 17.5 / 19.5 / 25.6 per mouth, 5.7 / 6.2 / 4.4 per wafer at 44 nm CPP; rule 1.69 / 1.17 /
0.52 s; all cells at 99.5 % at the rule's pulse. Chamber drive: 2.15 / 3.18 s TMA, 2.38 / 3.43 TDMAH, 0.34 / 0.41 TiCl4
(flat top of that wafer 1.95 s). Convergence M x2 0.558, h 0.625 nm 0.563, dt/2 0.556 s.
NAND 5 um cell (TMA, 0.3 Pa): slit floor 9.04 s, deepest recess end 10.37, slowest recess 11.30 (recess 86 of 100 at
4.3 um; the deepest is 9 % ahead because the saturated floor reflects the downward stream), t99 12.13; wp9 without
pillars 10.70 / 10.88 / 12.06 (the 12.06 was called noise there; it is this effect); nested model nominal 14.68 s, with
knudsen_correction="plates" 10.17 s; TDMAH 12.42 / 15.38 / 17.10 (nested 22.34 / 15.42); TiCl4 0.82 / 4.55 / 6.80 (8.00 /
5.57); chamber-driven TMA 10.45 / 11.27 / 11.85 (nested 13.53 / 11.07), TDMAH 12.47 / 13.92 / 14.45, TiCl4 2.62 / 7.09 /
8.88. Pillars take 8 % of the area (two discs out of the plates per pillar and recess, one wall in). Wall-area ratio
1422 per mouth, 171 per wafer; rule 43.5 / 30.1 / 13.5 s; exposure 3.39 / 5.13 / 2.04 Pa s, 7.6-9.2 x the effusion
limit. M x2 11.28 s; h 10 nm 10.99 s.
NAND 2XX cell: TMA 0.3 Pa slit floor 14.60, deepest recess 28.18 (end 28.76), slowest (recess 236 of 264, 90 % depth)
33.09 / 34.20 s, 9.93 Pa s; TDMAH 17.78 / 41.31 / 50.40 s (15.1 Pa s); TiCl4 1.25 / 23.15 / 29.88 s (9.0 Pa s); nested
deepest recess 43.6 / 66.2 / 26.9 s (wp4 diffusivity), 29.9 / 45.5 / 19.0 s (plates); chamber-driven 16.90 / 19.04 s,
TiCl4 not reached by the end of the 20 s pulse (recesses from 172 on below 90 % except the two at the floor; slowest
at 0.09, deepest at 0.97; 13.8 % of the area below 90 %); nested chamber 17.0 / 19.3 / 15.0 and 15.7 / 17.2 / 13.5 s. Wall-area ratio 4009 per mouth,
448 per wafer at the 2.25 um slit pitch; rule 113 / 78.5 / 35.1 s at 1 % (3.4 / 1.6 / 1.2 x the 0.3 Pa need); effusion
limit 4.1 / 6.3 / 2.1 s. M x2 32.81 s; h 8 nm 33.90 s. 35 runs, 5.9e8 particles, 100 min on 6 threads, longest 834 s;
every particle bookkept. Literature: no measured ALD thickness per sheet on a GAA stack or per level on a 2XX word-line
stack found (patents and a thesis describe TiN pinch-off between sheets qualitatively); Reiter 2024's 85 % HfO2 step
coverage fixed s_TDMAH in wp4 and is not an independent test. Everything here is application; the testable
predictions are the slot-centre deficit after a one-third pulse (3-8 % TMA-like, 30 % TiCl4 in 4 nm slots) and the
10-20 % dose lead of the lowest word lines over those five slit widths above them after a two-thirds pulse.
Tests: tests/test_device_geometries.py, 5 new; 129 pass in 13 s.
```
