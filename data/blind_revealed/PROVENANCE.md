# Blind test: revealed curves and provenance

Released by the custodian on 2026-11-26 after verifying tag `blind-predictions-v1` (commit 6279d07, committed 2026-11-26T15:55:06-05:00). The CSVs are byte-identical copies of the sealed set; `SHA256SUMS` holds the hashes computed at sealing time (before the tag), check with `shasum -a 256 -c SHA256SUMS`.

Columns in every CSV: `x_um`, the y value as plotted (name in the header), `x_over_H` (x divided by the nominal gap, or hole/pore diameter for Rose and Elam), `thickness_over_plateau` (plateau as reported in docs/blindtest/inputs.json; background-subtracted for Gonsalves), plus flag columns where present (`occluded_by_other_curve`, `uncertain_marker`, `band_low_nm`/`band_high_nm`). Script names below refer to the custodian folder (~/Documents/MATLAB/ald-blind-custody/scripts). The sections were written while the curves were sealed and are reproduced unchanged.

---

# PROVENANCE 1: Gao, Arpiainen, Puurunen, JVST A 33, 010601 (2015), Fig. 5(b)  [SEALED]

Source file: papers_inbox/"Microscopic silicon-based lateral high-aspect-ratio structures for thin film conformality analysis.pdf", page index 4.
Fig. 5(a) was NOT digitised (already used by the predictor).

## Method
Vector extraction with pymupdf `page.get_drawings()` (scripts/gao_explore.py, gao_extract.py, gao_finalise.py).
- 500 nm series: filled black 2.9 pt squares (`re` items), 50 points.
- 1000 nm series: white-filled 2.9 pt squares, 79 points.
- 200 nm series: crosses, each drawn as two 4.1 pt diagonal strokes; centre = bbox centre, 23 points.
Legend markers (x > 215 pt, y < 570 pt) excluded.

## Axis calibration (page points)
- x (aspect ratio, bottom axis of panel b): major ticks at 112.80, 150.30, 187.81, 225.31, 262.82 pt = AR 0, 50, 100, 150, 200. Spacing 37.50-37.51 pt (residual 0.01 pt = 0.013 AR).
- y (thickness nm): major ticks at 658.01, 635.13, 612.16, 589.27, 566.38, 543.41 pt = 0, 20, ..., 100 nm. Spacing 22.88-22.97 pt, i.e. residual up to 0.05 pt = 0.04 nm.
- x_um = AR x nominal gap (0.2, 0.5, 1.0 um). The scan step comes out as 2.0 um in all three series (10, 4, 2 AR units), which confirms the conversion.

## Errors
x: +/-0.1 AR unit (0.02 / 0.05 / 0.1 um). y: +/-0.1 nm (0.1 % of plateau). These are vector coordinates; the real uncertainty is the reflectometer's.

## Plateau and x50
The curves have no flat part: thickness falls from the first scan point. Plateau was therefore defined as the mean of the scan points at 0-4 um (three points, about one spot diameter of 5 um): 95.7 (sd 1.8), 102.2 (sd 0.9), 97.7 (sd 0.5) nm for 200, 500, 1000 nm.
x50 = first downward crossing of half that value, linear interpolation: 17.8, 70.4, 118.1 um (AR 89, 141, 118).

## Oddities
- Many points sit at exactly 0 nm (reflectometer returned nothing) interleaved with non-zero points near the film end; a late non-zero point appears at AR ~200 (200 nm) and ~196 (500 nm) at ~17 nm. All kept in the CSVs.
- 1000 nm series: one zero-valued point at AR 116 sits inside the front, between 58 nm (AR 114) and 50 nm (AR 118); then the next point (AR 120) is missing and AR 122 is 15 nm. The zero was ignored for x50 (otherwise crossing would be at AR 114.3). x50 here rests on the 118 -> 122 interval, hence the larger +/-2 um.
- The paper text says the 500 nm curve is at 50 % of the original value around AR 130; my 141 uses the 0-4 um plateau (102 nm).
- Channel length for these chips is not stated.

---

# PROVENANCE 2: Arts et al., JVST A 37, 030908 (2019), Fig. 3(a)  [SEALED]

Source: papers_inbox/arts2019.pdf, page index 4, embedded image xref 282 (1000 x 1339 px RGB), saved as renders/arts_fig3.png. Fig. 2 was NOT used.

## Method
Raster colour selection on the native embedded image (scripts/arts_extract.py, arts_finalise.py, overlay_check.py).
Masks: 150 C navy (B>90, R<90, G<90), 275 C orange (R>170, 80<G<185, B<90), 400 C red (R>150, G<80, B<80). Legend box (x 230-340, y 340-500 px) excluded. For every pixel column the largest vertical cluster of mask pixels (gaps <= 12 px, to bridge the white core of the open-circle markers) was taken and its median row used. One point per pixel column (655 / 647 / 704 points). A flag column marks columns where another curve's colour lies within 3 px of the cluster (the curve is then partly hidden: draw order is blue, orange, red on top).

## Axis calibration (pixels)
- x: tick marks every 100 units of distance/cavity height; x = 0 at 209.5 px, 1600 at 951.5 px (ticks found at 256.5, 302.5, ... 766.5, 905.5, 951.5; spacing 46.375 px, residual <= 0.5 px = 1.1 x/H = 0.5 um). x_um = (x/H) x 0.5 um; the top axis (um, 0-800) confirms H = 500 nm.
- y: ticks every 0.1; 0.0 at 571.5 px, 1.2 at 121.5 px (37.5 px per 0.1, residual <= 0.5 px = 0.0013).

## Errors
x: 1 px = 1.08 um; total +/-1.5 um. y: +/-0.01 (line + marker width ~8 px, median taken); +/-0.02 in flagged (occluded) columns, where the visible part of the lower curve is biased away from the curve on top.
Overlay check: renders/arts_overlay.png (digitised points drawn over the faded figure) follows all three traces.

## Plateau and x50
The plotted quantity is already normalised by the authors (1.0 at the cavity opening). None of the three curves is flat: all decline from the entrance (the authors call it soft saturation), most strongly at 150 C. Plateau was taken as the authors' 1.0. My reading over 0-10 um: 0.990 (150 C), 0.968 (275 C, mostly hidden under the red curve there), 0.996 (400 C).
x50 (crossing of 0.5): 533.8, 546.1, 599.3 um (x/H 1068, 1092, 1199). Single crossing in each case, also after a 7-column running median.

## Oddities
- Points at x < 0 (grey band = open area) are kept.
- Grey straight lines in the figure are the authors' linear fits at theta = 1/2; they are grey and were not picked up.
- 150 C curve has two small upward steps near x/H ~350 and ~520 (visible in the figure; kept).
- Near the film end (x/H > 1300) the three curves overlap heavily; values there are unreliable (flagged).
- Cycles/thickness (400 cycles, ~46 nm) are stated once for "the profiles reported in this work"; doses are given only for Fig. 2.

---

# PROVENANCE 3: Gonsalves et al., arXiv:2609.12460v1 (11 Sep 2026), Fig. 2(d-f)  [SEALED]

Source: downloaded 2026-11-26 from https://arxiv.org/pdf/2609.12460 into downloads_arxiv_2609.12460/2609.12460.pdf (921,376 bytes, 26 pages). Page index 7.
The ESI (Fig. S10 normalised profiles, Tables S11/S13/S14) is not part of the arXiv PDF and was not obtained.

## Method
Vector extraction (scripts/gons_explore.py, gons_extract.py, gons_finalise.py). Every marker is a small rectangle:
- DEZ/H2O: black filled squares (1.33 pt); some neighbouring squares are merged into one even-odd path, each sub-rectangle was taken as a marker.
- MeOH/DEZ/H2O: blue stroked squares (1.2 pt in panel d, 1.0 pt in e, f).
Legend markers (top right of each panel) excluded. Each panel is a mirrored structure: grey band = open area (Zone I), x = 0 at each edge of the band, distance increasing outwards. Each scan was split into left channel, right channel and Zone I: 12 half-profiles, 411-449 points each.

## Axis calibration (page points)
Bottom axis (x/H), major ticks every 200:
- panel d (150 C): left zero 161.10, left 800 at 112.06; right zero 173.36, right 800 at 222.40
- panel e (175 C): 310.30 / 261.27; 322.57 / 371.60
- panel f (200 C): 457.31 / 408.27; 469.56 / 518.60
(12.26 pt per 200 x/H = 100 um; residual of intermediate ticks <= 0.01 pt.) Top axis in um confirms H = 500 nm.
y (X-ray counts): 0 at 322.55 pt, 2000 at 256.48 pt, ticks every 500 at 306.00, 289.57, 273.02 (residual <= 0.03 pt = 1 count).

## Errors
Marker coordinates are quantised on a 0.11 pt grid: 0.9 um in x (= the ~1 um scan step) and 3.3 counts in y. So x +/-0.5 um, y +/-3 counts (0.3 % of net plateau without inhibitor, 0.6 % with). A few black markers have clipped widths (0.55-1.2 pt instead of 1.33) from the path merging; their centres may be off by up to 0.3 pt (2.7 um). Black series have 5-8 % fewer markers than blue (lost in merges).
The data themselves scatter: sd of plateau counts = 2.6-3.4 % of net plateau (DEZ/H2O), 4.1-4.8 % (MeOH).

## Plateau, background, x50
Followed the authors' normalisation (Methods 9.5.3): plateau = mean counts over 20-60 um in the channel, except DEZ/H2O at 150 C where they (and I) use 10-20 um because the signal already declines beyond that. Background = mean over 350-450 um (authors use 500-700 um, outside the plotted range). thickness_over_plateau = (counts - bg)/(plateau - bg).
x50 = first crossing of 0.5 on a 9-point running mean of the background-subtracted counts.

| curve | plateau | bg | Zone I | x50 um | x50/H |
|---|---|---|---|---|---|
| DEZ 150 left | 1482 | 91 | 1465 | 168.9 | 338 |
| DEZ 150 right | 1463 | 93 | 1465 | 129.1 | 258 |
| MeOH 150 left | 631 | 96 | 535 | 220.4 | 441 |
| MeOH 150 right | 644 | 98 | 535 | 214.4 | 429 |
| DEZ 175 left | 1413 | 92 | 1380 | 150.1 | 300 |
| DEZ 175 right | 1428 | 96 | 1380 | 147.9 | 296 |
| MeOH 175 left | 696 | 88 | 568 | 214.6 | 429 |
| MeOH 175 right | 701 | 92 | 568 | 214.9 | 430 |
| DEZ 200 left | 1294 | 86 | 1248 | 145.2 | 290 |
| DEZ 200 right | 1301 | 90 | 1248 | 144.2 | 288 |
| MeOH 200 left | 662 | 88 | 578 | 219.2 | 438 |
| MeOH 200 right | 665 | 89 | 578 | 239.7 | 479 |

Check against Table 1 (authors' mean of six scans, x/H): 297/430, 298/428, 288/454. My left-right means: 298/435, 298/430, 289/459. Agreement within a few units.

## Oddities
- DEZ/H2O at 150 C is strongly asymmetric left/right (x50 169 vs 129 um) and is not flat in the channel.
- MeOH/200 C is asymmetric too (219 vs 240 um).
- With inhibitor, the in-channel plateau is HIGHER than Zone I (ratio net plateau / net Zone I = 1.17-1.28); without, 1.00-1.05. The first ~20 um inside the channel are lower than further in for the MeOH scans.
- Background is not zero (86-98 counts) and drifts down by ~10 counts from 250 to 450 um.
- The date on the arXiv listing (Sep 2026) is as found.

---

# PROVENANCE 4: Dendooven et al., JES 156, P63 (2009) Fig. 5; JES 157, G111 (2010) Fig. 2(a)  [SEALED]

## 2009, Fig. 5 (papers_inbox/dendooven2009.pdf, page index 4)
Vector extraction (scripts/dend09_extract.py, dend09_finalise.py). Squares (AR 100) = `re` items, circles (AR 50) = 4-Bezier paths, triangles (AR 33) = 3-line paths; 20 points each plus one legend marker each (x ~405 pt, y 574/583/592), excluded. Centres = bbox centre (squares, circles), vertex centroid (triangles). All three first points coincide at y = 1.000, which confirms the triangle anchor.
Calibration (page points): x = 0 mm at 392.44, 20 mm at 547.39 (minor ticks every 1 mm, 7.75 pt, residual <= 0.06 pt = 0.008 mm). y = 0 at 599.07, 1 at 557.18 (ticks at 588.58, 578.16, 567.67).
Errors: x +/-0.02 mm, y +/-0.005.
Plateau = authors' relative scale, 1.0 at the first point (0.5 mm). None of the curves is flat: all decline from the entrance (to ~0.85-0.90 before the front for AR 100 and 50; AR 33 declines to 0.85 at 19.5 mm without a front).
x50 (linear interpolation, points 1 mm apart): AR 100: 10.68 mm; AR 50: 18.78 mm (only the last two points are below 0.76, the last is 0.42, so the front is cut by the slot end); AR 33: never reaches 0.5.
Gap heights are not stated; inferred from a = L p/(4A): 100, 200, 300 um.
All other figures in the paper are simulations or schematics.

## 2010, Fig. 2(a) (papers_inbox/dendooven2010.pdf, page index 1)
Vector extraction (scripts/dend10_extract.py). Panel (a) only (x < 198 pt). AR ~22 = green up-triangles (26 points), AR ~50 = purple down-triangles (27 points); each marker appears as fill + stroke copies, de-duplicated. Centre = vertex centroid. Legend markers excluded.
Calibration: x = -5 mm at 99.57 pt, 0 at 115.76, 25 at 196.70 (16.19 pt per 5 mm, residual <= 0.04 pt). y = 0 nm at 733.96 pt, 35 nm at 659.10 (10.69 pt per 5 nm).
Errors: x +/-0.03 mm; y +/-0.05 nm from coordinates, but whether the plotting program anchors triangles at the centroid or the bbox centre is unknown (difference 0.4 pt = 0.19 nm, opposite sign for up and down triangles), so +/-0.2 nm.
Plateau = mean of points at 1-3 mm: AR 22: 17.88 nm (sd 0.20); AR 50: 18.02 nm (sd 0.34). Points at negative depth (outside the slot) are 17.2-17.5 (AR 22) and 18.2-18.5 (AR 50).
x50: AR 50: 13.04 mm. AR 22: no crossing inside the slot; thickness is 15.4 nm at 19 mm, 11.6 at 20 mm (the slot end) and ~1-2 nm beyond.
Spacer thickness for AR ~50 is not given in the 2010 paper (0.2 mm inferred). Fig. 2(b) and Figs. 3-10 are plasma processes (not requested).

---

# PROVENANCE 5: Haimi et al., JVST A 44, 022410 (2026), Figs. 5 and 6  [SEALED]
(The predictor has read this paper's text.)

Source: papers_inbox/"Atomic layer deposition of zinc oxide films on lateral high-aspect-ratio test structures using diethylzinc and water as precursors.pdf" (accepted manuscript). Fig. 5 = page index 13, image xref 249; Fig. 6 = page index 15, xref 254; both 668 x 512 px RGB.
Listed but not digitised: Fig. 3 (raw Zn/O/Si X-ray intensity for sample B), Fig. 8(a,b) (samples A and B again under simulated curves), Figs. 2 and 7 (maps).

## Method
Raster hue selection on the native images (scripts/haimi_extract.py, haimi_finalise.py, overlay_check.py). Hue windows (saturation >= 0.45, 0.30 for yellow): A red 350-12 deg, B orange 20-33, C yellow 38-52; D blue 225-255, E cyan 170-195. A first attempt with RGB distance mixed anti-aliased orange with yellow; hue does not. Per pixel column: largest vertical cluster (gaps <= 6 px), median row; the cluster's upper/lower rows are stored as band_low/high. Flag column = another curve's colour within 2 px.
Legend area (rows < 170, cols > 340) excluded.

## Axis calibration (pixels; identical layout in both figures)
x ticks at 119.0, 184.5, 249.0, 314.5, 379.5, 444.5, 509.5, 574.5 px = -100 ... 600 um. Linear fit 1.5371 um/px, residuals <= 0.44 um.
y ticks at 426.5, 374.0, ..., 59.0 px = 0 ... 35 nm; 0.09524 nm/px, residual 0.

## Errors
x +/-2 um (1 px = 1.54 um, the data were measured at 1 um steps so each column mixes 1-2 points). y +/-0.15 nm (line 2-3 px thick plus noise band; the column median sits in the middle of the noise band). In flagged columns up to +/-0.3 nm.
Overlay: renders/haimi_fig5_overlay.png follows all three traces.

## Plateau and x50
Plateau = mean over 20-60 um inside the cavity (uniform rule for this source, chosen to skip the entrance step at 0-20 um). Open area = mean over -80 to -10 um.
| sample | plateau nm (sd) | open area nm | x50 um |
|---|---|---|---|
| A 125 C | 24.7 (0.9) | 28.9 | 369 |
| B 200 C | 22.1 (0.3) | 22.1 | 343 |
| C 250 C | 18.3 (1.3) | 13.2 | 361 |
| D 125 C 1 s | 8.4 (0.35) | 9.4 | 371 (raw crossings 365-371) |
| E 125 C 0.5 s | 7.7 (0.4) | 8.6 | 283 (raw 283-293) |
Paper values: PD50% of B "appr. 340 um" (mine 343); D and E 355 and 255 um (mine 371 and 283: the authors evidently normalise to a higher reference level).

## Oddities
- Only B has a flat part. A, D, E drop at the entrance and then decline all the way to the front; the 20-60 um mean is still falling inside the window (sd 0.9 nm for A).
- C is irregular in the cavity (several local maxima), and is thicker in the cavity than in the open area.
- Curves cross each other; flagged columns are 30-40 % of the points.
- A residual level of 0.3-0.6 nm persists beyond the front in all curves.

---

# PROVENANCE 6: Rose & Bartha, Appl. Surf. Sci. 255, 6620 (2009) Fig. 3; Rose, Bartha, Endler, Appl. Surf. Sci. 256, 3778 (2010) Figs. 1, 3  [SEALED]

Cylindrical holes (130 nm diameter, 2 um deep), NOT slits. Each experimental profile has only four TEM points.

## Method
Raster: embedded images rose2009 xref 22 (682 x 508), rose2010 xref 17 (669 x 1688, Fig. 1 a-c) and xref 34 (667 x 540, Fig. 3). Filled markers found by 4x4 erosion of dark pixels (sum RGB < 330), reconstruction by dilation, centroid (scripts/rose_explore.py). Legend markers and arrowheads rejected by position/size. Coordinates hard-coded in scripts/rose_finalise.py.

## Calibration (pixels)
- 2009 Fig. 3: x ticks 69.0 (0 nm) ... 665.0 (2000 nm), every 250 nm = 74.5 px (residual <= 0.5 px = 1.7 nm). y ticks 8.0 (21 nm) ... 421.0 (13 nm), 51.6 px/nm.
- 2010 Fig. 1: x ticks 68.5 (250 nm) ... 637.5 (1250 nm), every 125 nm = 71.1 px. y: panel a 29.5 (20 nm) - 470.0 (0); b 602.5 - 1044.0; c 1178.5 - 1619.5 (22.05 px/nm).
- 2010 Fig. 3: x as Fig. 1 (637.0 at 1250); y 30.0 (15 nm) - 471.0 (11 nm), 110.25 px/nm.

## Points (depth nm, thickness nm)
- 2009 TiO2/TDMAT 180 C (squares): (740, 19.58) (1315, 17.46) (1649, 15.37) (1782, 13.63). First depth matches the 740 nm quoted in the text.
- 2010 HfO2 180 C: (387, 18.89) (707, 17.76) (1022, 13.40) (1220, 9.90)
- 2010 HfO2 230 C: (405, 18.53) (616, 15.77) (793, 14.31) (971, 7.03)
- 2010 HfO2 270 C: (389, 17.89) (611, 16.45) (775, 6.02) (923, 3.11)
- 2010 TiO2/Cp*Ti(OMe)3 270 C: (470, 13.56) (741, 13.45) (930, 13.01) (1250, 12.14)

## Errors
x +/-5 nm digitising (authors: depth +/-10 nm). y: triangles' anchor (centroid vs bbox centre) uncertain by ~1.5 px: +/-0.07 nm (Fig. 1), +/-0.02 nm (Fig. 3); squares +/-0.03 nm.

## Plateau and x50
Plateau = shallowest point (no way to verify a flat part from four points; for 2009 the text says 20 nm at 740 nm equals the planar thickness). x50 by linear interpolation exists only for HfO2 230 C (0.916 um) and 270 C (0.729 um); the other three never fall to half within the measured depths. With four points these x50 values are crude (+/-0.05 um).
The simulated curves in the figures (grey/black/coloured lines) were not digitised.

---

# PROVENANCE 7: Elam, Routkevitch, Mardilovich, George, Chem. Mater. 15, 3507 (2003), Figs. 5(a), 6(a)  [SEALED]

Anodic alumina pores (cylinders open at both ends, L = 50 um), NOT slits. The figures show 0-25 um from one face.
No measured Al2O3 thickness-vs-depth profile exists in the paper (Fig. 3 is SEM at top/middle/bottom).

## Method
The figures are 1-bit image masks (page index 4, xrefs 46 and 47). Page regions were rendered at zoom 8 (renders/elam_fig5_z8.png, elam_fig6_z8.png; ~40 px per um) and markers detected with scipy.ndimage (scripts/elam_markers.py, elam_finalise.py):
- filled markers (circles, triangles): 9x9 erosion removes the connecting lines; blobs reconstructed by dilation; centroid.
- open markers (squares, diamonds): enclosed white regions 5-22 px; squares have fill ratio > 0.85 of their bbox, diamonds 0.5-0.7.
Series assignment: Fig. 5(a): 1 s = filled circles at integer um; 10 s = filled triangles at half-integer um; 3 s = open squares; 30 s = open diamonds. Fig. 6(a): 46 nm = open squares at half-integer um; 19 nm (circles) and 65 nm (triangles) both at integer um, separated by value: lower = 19 nm for x >= 3 um; for x = 1, 2 um the circle is the upper one and at x = 0 the lower one (assigned by eye; flagged uncertain in the CSV).

## Calibration (pixels at zoom 8)
- Fig. 5(a): x = 0 um at 226.0 px, 25 um at 1227.0 (major ticks every 5 um: 226, 426, 625.5, 826, 1026, 1227; residual <= 0.5 px = 0.012 um). y = 1.4 at 74.5 px, 0 at 1022.5 (major ticks every 0.2 at 74.5, 209.5, 346, 481.5, 617.5, 753, 888; residual <= 1 px = 0.0015).
- Fig. 6(a): x = 0 at 228.0, 25 um at 1266.0 (228, 435, 642, 850, 1058, 1266). y = 1.4 at 82.5, 0 at 1064.25 (82.5, 223, 364, 505, 645).

## Errors
x +/-0.05 um; y +/-0.01 for clean markers, +/-0.03 where markers merge with lines or each other (flagged). Scan points are 1 um apart and the EPMA resolution is 3-5 um, which dominates.
Missing points: Fig. 5(a) 1 s at x = 0; 10 s at x = 0.5 (merged blobs at the entrance); Fig. 6(a) 46 nm at x = -0.5 (if any).

## Plateau and x50
Curves are normalised by the authors to the maximum measured Zn concentration (one common maximum), so short exposures peak below 1. Rule: plateau = mean of the three consecutive points (within 1-10 um) with the highest mean; if that window has sd > 0.05 (no flat part), the peak value is used.
| curve | plateau | window um | x50 um |
|---|---|---|---|
| Fig 5a 1 s | 0.875 (peak) | 1 | 5.83 |
| Fig 5a 3 s | 0.896 | 2-4 | 9.72 |
| Fig 5a 10 s | 0.970 | 3.5-5.5 | 18.44 |
| Fig 5a 30 s | 1.004 | 4-6 | none (stays >= 0.80) |
| Fig 6a 19 nm | 0.999 (peak) | 2 | 3.98 |
| Fig 6a 46 nm | 0.988 | 2.5-4.5 | 10.01 |
| Fig 6a 65 nm | 0.980 | 4-6 | 16.65 |

## Oddities
- All curves rise over the first 1-2 um (membrane face + probe size): an entrance artefact, excluded from the plateau.
- The 1 s curve has no flat part at all (peak at 1 um, then falls); the 19 nm curve has only two points near 1.0.
- Baselines differ: ~0.03 (1 s), ~0.10 (3 s), ~0.02-0.04 (Fig. 6a); not subtracted.
- Fig. 6(a) x axis starts slightly below 0.

