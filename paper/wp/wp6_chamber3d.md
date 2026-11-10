# Work package 6: a full 3D chamber with a one-sided pumping port

Question: the axisymmetric wedge cannot show azimuthal non-uniformity. A real single-wafer chamber pumps through one
port on the side, so the carrier gas crosses the wafer from the far side toward the port. Does that cross-flow produce
a measurable azimuthal lag in the saturation time of the patterned (AR 100) wafer at the baseline operating point
(N2 300 Pa, 573 K, 300 sccm, 1 % TMA, s 5e-3, N0 4e18 /m2, d 100 nm, holeDensity 1.78e13 /m2)?

Answer in one paragraph: yes, but it is small. With the pump on one side the carrier gas over the far half of the
wafer barely moves outward (0.01 m/s at r = 140 mm against 0.27 m/s axisymmetric) while the port-side quadrant drains
at 0.85 m/s, and the instantaneous gas-side quantities are visibly asymmetric (surface TMA density 9 % and hole-bottom
coverage 11 % higher on the port side at r = 140 mm). The saturation time of the hole bottoms integrates that supply
and barely notices: its azimuthal peak-to-peak is 0.09 % at r = 50 mm, 0.37 % at 100 mm and 0.96 % at 140 mm (port
side earlier), the worst face on the wafer is 14.17 s at r = 133 mm, 160 deg from the port, against 14.05 s on the
wedge, and of that +0.8 % about +0.3 % is the coarser 3D gap resolution seen at the centre as well. The reason is the
one the wedge already gave: at gap Peclet 0.09 the TMA reaches the wafer by diffusion across the gap from the
showerhead directly above, the wafer captures almost all of it locally, and the cross-flow can only redistribute the
few per cent that escape. The flat surface between the holes, which saturates at 2 s while the gas phase is still being
drawn down hard, shows a larger 4.8 % azimuthal spread at r = 140 mm (1.92 s port side, 2.01 s far side).

## What was done

### Mesh (reactor/chamber3d/steady/system/blockMeshDict, topoSetDict, createPatchDict)
Full 3D cylinder, radius 200 mm, gap 15 mm, axis along x as in the wedge (wafer at x = 0, showerhead at x = 15 mm,
face radius sqrt(y^2 + z^2)), so the patternedSurface boundary condition and the wedge conventions carry over
unchanged. blockMesh O-grid of nine blocks: one central square block (half side 75 mm, its edges bulged to r = 82 mm
at their mid points to round the O-grid corners), four curved blocks out to r = 150 mm (wafer / showerhead radius) and
four ring blocks from 150 to 200 mm, with arc edges on both circles. 48 cells along each quarter arc and across the
square (1.875 deg per cell at the side wall), 20 radial cells from the square to r = 150 mm, 14 from 150 to 200 mm,
12 cells across the gap (1.25 mm, first cell centre 0.625 mm from the showerhead; the wedge had 30 cells, 0.25 mm).
Patches: showerhead (x = 15 mm, r < 150 mm, 6144 faces), topWall (x = 15 mm, 150 < r < 200 mm, 2688 faces), wafer
(x = 0, r < 150 mm, one patch of 6144 faces, so one patternedSurface condition covers it), floor (x = 0,
150 < r < 200 mm, 2688 faces), sideWall (r = 200 mm, 1920 faces) and outlet, the pumping port: the 60 degree sector
of the side wall centred on +y (azimuth 0), 384 faces (32 x 12). The port was carved after meshing with topoSet
(faceSet from the sideWall patch, subset to the box |z| < 100 mm, y > 170 mm; with 1.875 deg cells the sector
boundary at +-30 deg falls on face edges, face centres sit at |z| <= 97.1 mm inside and >= 102.8 mm outside, so the
box selects exactly the sector) and createPatch -overwrite, which moved those faces into a new patch of type patch.
Method chosen because a 60 degree sector does not fit the 4-fold block structure of an O-grid; the alternative, 19
blocks with the square and the +-y ring blocks split at +-30 deg, gives the same faces.
checkMesh: 105,984 hexahedra, 116,077 points, 327,936 faces; max aspect ratio 5.04; non-orthogonality max 30.7 deg,
average 5.1 deg; max skewness 0.794; all checks OK. Total volume 1.8846e-3 m3 (pi 0.2^2 0.015 = 1.8850e-3).
Port area 3.1415e-3 m2 (2 pi 0.2 0.015 / 6 = 3.1416e-3), showerhead area 7.0673e-2 m2 (pi 0.15^2 = 7.0686e-2).

### Steady carrier flow (reactor/chamber3d/steady)
simpleFoam, laminar, nu 1.65e-2 m2/s, the wedge's fvSchemes, fvSolution (SIMPLEC, relaxation 0.9, residualControl
p 1e-6, U 1e-7) and boundary conditions: uniform inflow (-0.050 0 0) m/s on the showerhead, noSlip on wafer, floor,
topWall and sideWall, totalPressure p0 = 0 with pressureInletOutletVelocity on the port. Converged by the residual
controls after 93 iterations (18 s wall): initial residuals at the last iteration Ux 9.6e-8, Uy 1.0e-9, Uz 8.6e-9,
p 4.1e-8 (postProcessing/residuals/0/solverInfo.dat). Mass balance: showerhead in 3.53366094e-3 m3/s, port out
3.53366097e-3 m3/s (flowRatePatch), i.e. closed to 1e-8 relative; mean port speed 1.125 m/s. U, p and phi are in
steady/93 and were copied to the pulse case as its frozen flow.
The flow is strongly one-sided (gap-averaged radial velocity from the frozen field, positive outward; the
axisymmetric value Q/(2 pi r g) in brackets): at r = 140 mm, port side +0.854 m/s, 90 deg +0.073 m/s, far side
+0.013 m/s (0.268 m/s); at r = 100 mm the far side is slightly inward, -0.081 m/s (0.375 m/s); at r = 199 mm the port
sector carries +1.12 m/s and the rest of the side wall 0; |U| in the outer 15 mm of the chamber is 1.13 m/s in the
port sector and 0.019 m/s in the far sector. Gas over the far half of the wafer therefore moves slowly across the
wafer toward the port; gas over the port-side quadrant leaves at about three times the axisymmetric speed.

### Pulse (reactor/chamber3d/pulse)
scalarTransportFoam on the frozen flow with the 0/T of reactor/showerhead_pulse_patterned_AR100_long: Danckwerts
uniformMixed inlet with the value fraction recomputed for this mesh, f = U/(U + D/delta) = 0.050/(0.050 + 13.44) =
3.706e-3 (delta 0.625 mm; the wedge's 1.49e-3 at delta 0.25 mm); patternedSurface wafer with s 5e-3, vbar 410 m/s,
N0 4e18 /m2, c0 3.79e20 /m3 (1 % TMA), d 100 nm, AR 100, holeDensity 1.78e13 /m2, 200 hole cells; zeroGradient on
topWall, floor and sideWall; inletOutlet on the port. DT 0.0084 m2/s. 16 s pulse, run to 17 s, dt 1 ms, written every
0.1 s (the wedge wrote every 0.05 s; each 3D write is 25 MB because the BC writes its 200 hole cells for each of
the 6144 faces, 3.4 GB for the run, so the interval was doubled; t90 is interpolated between writes, see below).
Face centres of every patch were written with `postProcess -func writeCellCentres -time 0` (the boundary values of
C are the patch face centres, in the same face order as every other patch field).
Timing: a 0.2 s copy of the case (reactor/chamber3d/pulse_timing) took 26 s wall, predicting 0.6 h serial for 17 s;
the full run took 1939 s (32.3 min, a steady 1.9 min per simulated second, slightly under the 37 min the timing copy
predicted), serial on one core. Open MPI 5.0.11 is in the app, but reconstructPar
would zero the per-face hole-coverage array of the BC (its mapping constructor does not carry thetaHole), so a
parallel run would have needed a processor-wise reader; not needed at this wall time. Nothing was coarsened.

### Post-processing (coupling/postprocess_chamber3d.py)
Reads the wafer patch's thetaBottom, thetaFlat and value lists from every time directory (same face order as the
wafer block of 0/C; cached in pulse/wafer_fields.npz), computes the 0.9 crossing per face by linear interpolation
between write times (and keeps the first write time at or above 0.9, the definition used for the wedge numbers in
earlier parts), azimuthal statistics of faces within +-2.5 mm of r = 50, 100 and 140 mm in 15 degree bins (raw band
statistics, bin means, and bin values detrended to the exact radius by a linear fit in r inside each bin, since the
O-grid rows are not circles), the outermost row, the worst face on the wafer, the deviation of every face from the
wedge profile at the same radius (the wedge case reread with coupling/patterned_lib.py and the same interpolation),
the surface TMA density and coverages by azimuth at 0.3, 8 and 14 s, and the steady radial velocity by azimuth from
the frozen U field and the cell centres. Writes figures/chamber3d_t90.csv (6144 rows: face centre, radius, azimuth,
both t90 definitions for bottoms and flat, deviation from the wedge, value/c0 at 0.3 and 8 s),
figures/chamber3d_t90_map.png, figures/chamber3d_azimuthal.png and figures/chamber3d_summary.txt (every number below).

## Results

Hole-bottom t90 (0.9 crossing interpolated between 0.1 s writes; wedge values from the same interpolation on
reactor/showerhead_pulse_patterned_AR100_long, whose first-write values are the 14.05 s centre / 14.10 s worst quoted
in earlier parts; port at azimuth 0):

| radius (mm) | faces in band | mean (s) | min (s) | max (s) | peak-to-peak, raw band (% of mean) | peak-to-peak, 15 deg bins at the exact radius (%) | port side - far side (%) | wedge t90 at r (s) | 3D mean / wedge |
|---|---|---|---|---|---|---|---|---|---|
| 50 | 140 | 14.048 | 14.042 | 14.055 | 0.09 | 0.09 | -0.08 | 14.009 | 1.0028 |
| 100 | 308 | 14.075 | 14.041 | 14.098 | 0.40 | 0.37 | -0.36 | 14.021 | 1.0038 |
| 140 | 360 | 14.089 | 13.975 | 14.156 | 1.29 | 0.96 | -0.90 | 14.006 | 1.0060 |
| 148.9 (outermost row) | 248 | 13.687 | 13.544 | 13.931 | 2.83 | - | -0.70 | 13.280 (r = 149.5 mm) | 1.031 |

The raw band spread at r = 140 mm includes the radial gradient inside the +-2.5 mm band (0.58 % on the wedge over the
same band); the detrended bin values remove it: 14.00 s in the port bins, 14.13 s in the far bins, minimum 13.996 s
at 8 deg and maximum 14.131 s at 158 deg. The azimuthal profile is a smooth single-lobed curve at every radius
(figures/chamber3d_azimuthal.png), earliest at the port and latest opposite it, with amplitude growing from 0.012 s at
r = 50 mm to 0.05 s at 100 mm and 0.13 s at 140 mm.

Whole wafer: innermost face (r = 2.4 mm) 14.044 s (first write 14.1 s) against the wedge centre 14.007 s (14.05 s),
+0.26 %: the systematic offset of the 3D mesh (12 cells across the gap against 30, first cell centre 0.625 mm against
0.25 mm). Worst face 14.168 s (first write 14.2 s) at r = 132.7 mm, azimuth -160 deg (far side), against the wedge
worst 14.051 s (14.10 s) at r = 130.5 mm: +0.83 %, of which about 0.3 % is the mesh offset and 0.5 % the one-sided
pump. Earliest face 13.544 s on the outermost row (r = 148.9 mm; wedge outermost face 13.28 s). Wafer-wide spread
4.44 % of the mean 14.059 s, almost all of it the radial edge drop that the wedge also has. Face by face, 3D minus
wedge at the same radius: mean +0.44 %, range -0.09 to +1.37 %; the deviation map shows a faint imprint of the
O-grid block boundaries (about 0.1 %, a mesh artefact) under the smooth port-to-far-side gradient. All 6144 faces
reach 90 % at the bottom and on the flat surface within the 16 s pulse.

Flat-surface t90 (same bands): r = 50 mm 1.929 s (peak-to-peak 0.02 %), 100 mm 1.934 s (0.46 %, port side 0.38 %
earlier), 140 mm 1.983 s (4.75 %; 1.917 s minimum on the port side, 2.011 s maximum on the far side; wedge 1.969 s).
The flat surface saturates at 2 s, while the gas phase is still drawn down to 0.02 c0 at the surface and lateral
transport still matters; the port-side cross-flow brings it TMA that the far side does not get.

Gas-side asymmetry during the pulse, r = 140 mm band, port side (|azimuth| <= 30 deg) against far side
(>= 150 deg): surface TMA density value/c0 at 0.3 s 0.0186 against 0.0172 (peak-to-peak 9.2 % of the mean; 0.01 % at
r = 50 mm), flat coverage at 0.3 s 0.213 against 0.191 (12.8 %); at 8 s surface density 0.128 against 0.126 (2.5 %),
hole-bottom coverage 0.0084 against 0.0075 (15 %); at 14 s surface density 0.253 against 0.237 (8.4 %), hole-bottom
coverage 0.899 against 0.835 (10.8 %; 0.7 % at r = 50 mm). The bottom coverage climbs from 0.008 at 8 s to 0.85 at
14 s, so a 6 % coverage lead at 14 s is a 0.13 s lead in saturation time: the instantaneous asymmetry of order 10 %
collapses to 1 % in t90 because the hole bottoms fill in a steep late wave fed by the integrated supply.

Transport numbers for the reading: gap Peclet number U g / D = 0.089; radial Peclet at the wafer edge
u_r g / D = 0.45 axisymmetric (1.5 on the port side, 0.02 on the far side of the 3D chamber); lateral diffusion time
over the wafer radius R^2/D = 2.7 s against the 14 s saturation.

Caveats: 3D gap resolution 12 cells (offset +0.26 % at the centre against the wedge; the azimuthal numbers are
differences on the same mesh and do not carry it); t90 interpolated between 0.1 s writes (first-write values are in
the csv); the frozen carrier flow; the port is a plain 60 degree sector of the side wall with a uniform-pressure
condition, no pumping plenum or baffle; the floor beyond the wafer and the side wall do not absorb; Stokes flow
(Re 0.045), so the flow pattern is the same at any flow rate and the only operating-point dependence of the
azimuthal lag is through the Peclet numbers.

## Figures
- figures/chamber3d_t90_map.png: polar map of the hole-bottom t90 over the wafer at the 6144 face centres (left,
  absolute; right, deviation from the wedge value at the same radius in %), the 60 degree port marked in red at
  r = 200 mm, the wafer edge dashed.
- figures/chamber3d_azimuthal.png: hole-bottom (left) and flat-surface (right) t90 against azimuth from the port at
  r = 50, 100 and 140 mm (faces within +-2.5 mm as dots, 15 degree bin values at the exact radius as lines, the
  wedge value at that radius dashed).
- figures/chamber3d_t90.csv: per-face table. figures/chamber3d_summary.txt: every number in this note.

## Proposed paper text

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
(Figure, chamber3d_azimuthal.png). The slowest point on the wafer moves from the edge ring of the wedge to the far
side at r = 133 mm, 160 deg from the port, at 14.17 s against 14.05 s in the wedge, and 0.3 % of that 0.8 % is the
coarser 3D gap resolution, which shifts the centre value by the same amount. The flat surface between the holes,
which saturates at 2 s while the gas phase is still drawn down to 0.02 c0, is the more sensitive indicator: 1.92 s on
the port side against 2.01 s on the far side at r = 140 mm, a 4.8 % spread.

The reason is the one behind the whole paper. At gap Peclet 0.09 the precursor reaches the wafer by diffusion across
the gap from the showerhead directly above each point, the patterned wafer captures almost all of what arrives, and
the saturation time is the time to deliver (1 + A_wall/A_wafer) N0 per unit area at the showerhead's own rate. The
cross-flow toward the port can redistribute only the few per cent that escape capture, and it does so over a lateral
diffusion time (R^2/D = 2.7 s) short against the 14 s pulse. A one-sided pump therefore costs about 1 % of pulse time
on this wafer; the recipe margin of the design rule (the factor 1.2) covers it, and the azimuthal lag would need a
convective gap (Peclet of order 1, i.e. ten times the flow or pressure) before it mattered. The same ordering holds in
the wedge's radial direction, where the edge lag is 0.4 %.

## CLAIMS rows

| Full 3D chamber, one-sided 60 degree pump port: azimuthal non-uniformity of the AR-100 wafer's hole-bottom saturation time | peak-to-peak 0.09 % at r = 50 mm, 0.37 % at 100 mm, 0.96 % at 140 mm (port side earlier); worst face 14.17 s at r = 133 mm, 160 deg from the port, vs 14.05 s on the wedge (+0.8 %, of which +0.3 % is the 3D gap resolution, seen at the centre too) | figures/chamber3d_t90.csv, figures/chamber3d_azimuthal.png, figures/chamber3d_t90_map.png, figures/chamber3d_summary.txt | reactor/chamber3d/make_and_run.py (mesh, simpleFoam, scalarTransportFoam, 32 min serial); coupling/postprocess_chamber3d.py |
| Same chamber, gas-side asymmetry of the one-sided pump | outward velocity at r = 140 mm 0.85 m/s port side vs 0.01 m/s far side (0.27 m/s axisymmetric); surface TMA density 9 % and hole-bottom coverage at 14 s 11 % higher on the port side at r = 140 mm; flat-surface t90 at r = 140 mm 1.92 s port side vs 2.01 s far side (4.8 %) | figures/chamber3d_summary.txt | same |
| 3D mesh (12 gap cells) vs wedge (30) at the same radius | hole-bottom t90 +0.26 % at the centre, +0.44 % mean over the wafer; mass balance of the 3D steady flow closed to 1e-8 | figures/chamber3d_summary.txt; reactor/chamber3d/steady/log.checkMesh, postProcessing | same |

## Run log (for paper/assumptions.md)

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
