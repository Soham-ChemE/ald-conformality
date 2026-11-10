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
