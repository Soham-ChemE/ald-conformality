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
