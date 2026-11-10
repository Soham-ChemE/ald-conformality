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

## Night 2, part 3: reversible water adsorption (reactor/lib/aldWallKineticsDev, kdes; figures/h2o_desorption.png)
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
