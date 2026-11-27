# Second blind test: revealed curves (released by the custodian after tag `blind-predictions-v2`, 2026-11-26T16:43:49-05:00)

Columns: x_um, x_over_H, thickness_over_plateau, raw value as plotted. Byte-identical copies of the sealed CSVs (SHA256SUMS). Sealed originals, scripts, overlays: ~/Documents/MATLAB/ald-blind-custody-v2/.

Plateau = median of the digitised points between 0.1 and 0.4 of x50 (x >= 5 um); x50 = first downward crossing of half the plateau (3-point running median, linear interpolation). Modes: grid = marker centres read on the 2 um scan grid; ridge = markers not resolved along x, values in 2 um bins; line = line trace, one value per pixel column.

## Entering curves (scored)

| file | figure | role | gap nm | cycles | plateau (nm or figure units) | x50 um | mode |
|---|---|---|---|---|---|---|---|
| yim2020_s01.csv | Yim 2020 ESI Fig. S4 | predict | 500 | 500 | 52.52 | 140.4 | grid |
| yim2020_s02.csv | Yim 2020 Fig. 8(b) | predict | 500 | 500 | 54.03 | 141.3 | grid |
| yim2020_s03.csv | Yim 2020 ESI Fig. S4 | predict | 500 | 500 | 53.6 | 110.6 | grid |
| yim2020_s04.csv | Yim 2020 Fig. 7 | predict | 500 | 500 | 47.19 | 105.1 | grid |
| yim2020_s09.csv | Yim 2020 Fig. 9(a) | calibration | 2000 | 500 | 48.53 | 376.6 | ridge |
| yim2020_s11.csv | Yim 2020 Fig. 9(d) | predict | 500 | 1000 | 99.95 | 96.1 | grid |
| yim2020_s12.csv | Yim 2020 Fig. 11(b) | predict | 500 | 500 | 50.15 | 97.3 | grid |
| yim2020_s13.csv | Yim 2020 Fig. 11(a) | predict | 500 | 500 | 49.35 | 133.3 | ridge |
| yim2020_s15.csv | Yim 2020 Fig. 11(b) | predict | 500 | 500 | 51.9 | 99.8 | grid |
| yim2020_s16.csv | Yim 2020 Fig. 11(b) | predict | 500 | 500 | 48.75 | 110.8 | grid |
| haimi2021_A_opt.csv | Haimi 2021 Fig. 8 (opt) | predict | 500 | 1080 | 0.9813 | 97.8 | line |

## Not entered (reported, not scored)

| file | figure | reason |
|---|---|---|
| not_entered/yim2020_s05.csv | Yim 2020 Fig. 7 | criterion 2 (plateau not flat: value 0.059, limit 0.05) |
| not_entered/yim2020_s06.csv | Yim 2020 Fig. 7 | criterion 2 (plateau not flat: value 0.155, limit 0.05) |
| not_entered/yim2020_s07.csv | Yim 2020 Fig. 9(c) | criterion 2 (plateau not flat: value 0.309, limit 0.05); criterion 3 (fewer than six points between 0.25 and 0.75: 5 recoverable) |
| not_entered/yim2020_s10.csv | Yim 2020 Fig. 9(d) | criterion 3 (fewer than six points between 0.25 and 0.75: 4 recoverable) |
| not_entered/yim2020_s14.csv | Yim 2020 Fig. 11(a) | criterion 2 (plateau not flat: value 0.131, limit 0.05) |
| not_entered/haimi2021_A_xray.csv | Haimi 2021 Fig. 8 (x-ray) | criterion 5 (spot or step size not stated in the source) |
| not_entered/utriainen2022_Al2O3_R2.csv | Utriainen 2022 Fig. 4 (Al2O3) | criterion 5 (spot or step size not stated in the source); criterion 6 (no second entering curve of the same process, temperature and reactor) |
| not_entered/utriainen2022_TiO2_R2.csv | Utriainen 2022 Fig. 4 (TiO2) | criterion 2 (plateau not flat: value 0.103, limit 0.05); criterion 5 (spot or step size not stated in the source); criterion 6 (no second entering curve of the same process, temperature and reactor) |

Notes: yim2020_s10 fails criterion 3 because its markers below about 0.4 of the plateau are hidden under other series in the figure. yim2020_s05 is borderline on criterion 2 (0.059; about 0.04 with the authors' tabulated plateau); the rule was not changed. yim2020_s07 uses the authors' normalisation of Fig. 9(c).

---

# PROVENANCE 1: Yim et al., Phys. Chem. Chem. Phys. 22, 23107 (2020) and ESI  [released]

Source PDFs: ald-conformality/data/papers_inbox/"Saturation profile based conformality analysis ... .pdf" and d0cp03358h1_suppl.pdf (CC BY 3.0). All figures are embedded raster images; extracted at native resolution (renders/yim2020_p*_x*.png, yim2020esi_p*_x*.png).
Al2O3, TMA/H2O, 300 C, Picosun R-150, ca. 3 hPa, N2 150 sccm; PillarHall-3; reflectometer FilmTek 2000M, 2 um step, spot 5-6 um (50x), 25 um (10x), 50 um (5x). Fig. 9 curves are averages over repeated scans and lateral lengths; series B is the L = 5 mm scan.

## Method (scripts/yim_extract.py)
Colour mask per series -> 2-D correlation with the marker shape taken from the figure legend -> per-column ridge of the correlation.
'grid' mode (>= 4 px per 2 um step): values read at the 2 um scan grid, the grid phase chosen to maximise the summed correlation; a point is kept if the correlation is >= 0.55 (0.35 for s01, 0.40 for s06, both partly hidden).
'ridge' mode (markers not resolved along x: Fig. 9(a) at 0.92 um/px, Fig. 11(a) at 0.64 um/px): ridge averaged in 2 um bins centred on odd micrometres.
Sample 7 is read from Fig. 9(c) (authors' normalised thickness against x/H, markers 20 units apart and individually resolved).
Axis calibration: tick marks detected on the frame, linear fit; residuals in the table. Overlays: renders/overlay_<panel>.png.

| curve | panel | mask | mode | points | x unit/px | y unit/px | tick residual x | tick residual y |
|---|---|---|---|---|---|---|---|---|
| s01_S4_v1a | S4 | magenta | grid | 55 | 0.375 | 0.1556 | 0.06 | 0.064 |
| s02_fig8b | fig8b | lgreen | grid | 84 | 0.358 | 0.1623 | 0.16 | 0.128 |
| s03_S4_v2a | S4 | black | grid | 100 | 0.375 | 0.1556 | 0.06 | 0.064 |
| s04_fig7_x50 | fig7 | dblue | grid | 91 | 0.322 | 0.1191 | 0.11 | 0.045 |
| s05_fig7_x10 | fig7 | black | grid | 61 | 0.322 | 0.1191 | 0.11 | 0.045 |
| s06_fig7_x5 | fig7 | gray | grid | 53 | 0.322 | 0.1191 | 0.11 | 0.045 |
| s07_fig9a | fig9a | lblue9 | ridge | 62 | 0.923 | 0.1384 | 0.28 | 0.046 |
| s07_fig9c | fig9c | lblue9 | grid | 22 | 1.031 | 0.0028 | 0.00 | 0.001 |
| s09_fig9a | fig9a | green9 | ridge | 200 | 0.923 | 0.1384 | 0.28 | 0.046 |
| s10_fig9d | fig9d | orange | grid | 66 | 0.460 | 0.2985 | 0.07 | 0.000 |
| s11_fig9d | fig9d | dblue | grid | 97 | 0.460 | 0.2985 | 0.07 | 0.000 |
| s12_fig11b | fig11b | red | grid | 67 | 0.378 | 0.1099 | 0.26 | 0.000 |
| s13_fig11a | fig11a | teal | ridge | 62 | 0.642 | 0.1095 | 0.38 | 0.044 |
| s14_fig11a | fig11a | black | ridge | 200 | 0.642 | 0.1095 | 0.38 | 0.044 |
| s15_fig11b | fig11b | navy | grid | 81 | 0.378 | 0.1099 | 0.26 | 0.000 |
| s16_fig11b | fig11b | black | grid | 94 | 0.378 | 0.1099 | 0.26 | 0.000 |

(s07_fig9a is an alternative reading of sample 7 from Fig. 9(a), kept in renders/ only.)

## Errors
x: 1.5 px + 0.3 um tick residual (+1 um in ridge mode, where the scan phase is not recovered). y: 1.5 px. Plateau: about 1 % (differences from Table 1 S_IIb are within 3 %).
Points hidden under other series are missing, not interpolated.

## Not digitised
Sample 8 everywhere (Zenodo 3974438, seen). ESI S7 (individual repeats behind Fig. 9), S15 (sample 11 again), S16 (sample 8), S8 (SEM-EDS of samples 10, 2, 11: entrance not marked), S17 (V0005: nine overlapping line traces).

# PROVENANCE 2: Haimi et al., Appl. Surf. Sci. Adv. 5, 100102 (2021), Fig. 8  [released]
Source: downloads/haimi_2021_apsadv/ (Aalto open copy, CC BY-NC-ND). Fig. 8 = image xref 262, 800 x 555 px.
Specimen A: Al2O3, TMA/H2O, 300 C, ca. 300 Pa, Picosun R-150, 0.1-4.0-0.1-4.0 s, 1080 cycles, gap 500 nm, L 1000 um; thickness 103 nm (optical), 100 nm (ED-EPMA).
Traces: 'opt' (blue; reflectometer, about 5 um spot, 100 points at 2 um) and 'x-ray' (orange; ED-EPMA at 5 kV, average of 10 scans; step and probe size not stated -> criterion 5).
Method (scripts/lines_extract.py): colour mask, per pixel column the mid-row of the largest cluster. x axis is x/H: ticks 122.5 ... 775 px = 0 ... 500 (0.767 per px = 0.38 um); y ticks 379 ... 11.5 px = 0 ... 1.25 (0.0034 per px). Tick residual 0.28 x-units.
The blue trace is hidden where the orange one runs on top of it (372 of about 486 columns recovered). Both traces are smoothed spreadsheet lines through the measured points.
Errors: x 0.9 um, y 0.006 (line width 3 px), plateau 1 %.
Fig. 7 (x-ray trace alone with the slope fit) and Fig. 3 (raw intensities) not digitised.

# PROVENANCE 3: Utriainen et al., Proc. SPIE 12008, 120080D (2022), Fig. 4  [released]
Source: downloads/utriainen_2022_spie/. Fig. 4 = image xref 169, 724 x 548 px. Relative thickness from greyscale microscope images (Olympus MX51), PillarHall LHAR4, gap 500 nm, Beneq TFS-200: Al2O3 500 cycles (64 nm, TMA/H2O 300 C, 0.15 s / 1.0 s, 2.1-2.2 mbar) and TiO2 500 cycles (32 nm, TiCl4, 120 C). Pixel size of the analysis not stated (criterion 5); one curve per process (criterion 6).
Calibration: x ticks 152.5 ... 675.5 px = 0 ... 250 um (0.477 um/px); y ticks 457.5 ... 28.5 px = 0 ... 1.0. The arrow annotations cross the traces near half height.

