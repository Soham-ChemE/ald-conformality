# Work package 10: gate-all-around nanosheet stacks and a current 3D-NAND cell in the voxel solver

Repository state: no commits made; everything is in the working tree for the integrator. README.md, CLAIMS.md,
paper/assumptions.md, paper/draft.md, main.tex and docs/ are untouched (their text is proposed at the end). The
files of the slit-transport refit (src/aldconf/transport.py, geometry.py, solver.py, coupling/fit_pillarhall*.py,
coupling/wordline_recess*.py) were not edited; aldconf.recess is only called.

## What was done

1. Literature geometry (Section "Sources"): published dimensions of stacked-nanosheet gate-all-around (GAA)
   devices (IBM 2017, imec 2016/2018, Intel 2024, Samsung 2018, the Wong and Kakushima 2022 scaling review, the
   Mohapatra 2021 design study) and of 3D-NAND word-line stacks beyond 200 layers (TechInsights teardown tables of
   YMTC 232L and of the 128/176-layer products; imec 2021-2025 research stacks; the Onto Innovation metrology
   notes), with what each source measured and how. No source was found that reports the thickness of an ALD film
   (interfacial oxide, HfO2 or TiN) measured sheet by sheet on a GAA stack or inside an inner spacer; the GAA part
   of this package is therefore an application of the solver, not a validation. For NAND the only published
   step-coverage number on a word-line stack is the 85 % HfO2 step coverage of Reiter et al. 2024 on Lam's 76/98-layer
   test structure, already used in wp4 to calibrate the TDMAH sticking value.
2. Two geometry builders in `coupling/device_geometries.py` (occupancy arrays through `voxel.from_array`):
   `gaa_stack` (N sheets of thickness t_s at spacing g_s and width W under a gate trench of length Lg closed by the
   spacers, periodic along the sheet width with period W + gap_x, mouth at the trench top) and `nand_cell` (the
   word-line slit with recesses of `voxel.recessed_slit`, plus the channel-hole pillars standing in the recesses as
   solid cylinders in staggered rows at the hole pitch along the slit). Face-region classifiers (sheet tops, slot
   centres and edges per slot, floor; slit floor, recess deep halves and end walls per recess, the minimum over all
   recesses) give area-weighted coverage histories. `figures/device_geometries.png` pictures every cell.
3. Runs with lib/aldmc at the chamber's TMA conditions (s 5e-3, N0 4e18, 573 K), HfO2/TDMAH (7e-3, 2.77e18, mass
   354.7) and TiN/TiCl4 (0.1, 1.24e18, 189.7), at constant 0.3 Pa mouth density and under the chamber's centre-face
   density history of the long AR-100 case (the drive of coupling/voxel_mc_crosscheck.py part `driven`): coverage
   against time of every region, the time (and exposure) for 90 and 99 % everywhere, the last region to saturate,
   the wall-area ratio of a wafer tiled with the cell and the rule's pulse for it, the step coverage at that pulse;
   one M doubling and one resolution change per geometry. The nested 1D recess model (aldconf.recess) was run on the
   two NAND cells without pillars as the reference, with wp4's hydraulic-diameter slit diffusivity and with the
   parallel-plate correction the refit package added (`knudsen_correction="plates"`).
4. Tests `tests/test_device_geometries.py` (5): builder volumes, periodicity, pillar discs and area formulas; a
   tiny-budget GAA run and a tiny NAND run complete with every particle bookkept.

## Sources

What each source gives, how it was measured, and the uncertainty as read here. "Teardown" means a cross-section
measurement by a reverse-engineering firm; "device paper" means the authors' own TEM and process description;
"review" and "TCAD" give assumed or compiled values.

### Gate-all-around nanosheet stacks

| source (DOI) | type | what was read | uncertainty / remark |
|---|---|---|---|
| Loubet et al. (IBM), Symposium on VLSI Technology 2017, doi 10.23919/VLSIT.2017.7998183 | device paper | first stacked-nanosheet GAA: 3 sheets, sheet thickness 5 nm, vertical sheet-to-sheet spacing 10 nm, gate length 12 nm, CPP 44/48 nm; HfO2/TiN replacement gate with multiple work functions, "aggressive sheet to sheet spacing for wide stacked sheets". The abstract was read from IEEE Xplore and the IBM Research record; the 5 nm / 10 nm / 3-sheet values were read from two secondary texts that cite this paper (Wong and Kakushima 2022, Section 2.3, "the sheet height is reduced to 5 nm ... the intersheet spacing is 10 nm, the same value reported by the IBM group"; Tirumaladass 2022, Section 2.4) | the paper itself is paywalled; the values are nominal process targets as quoted, no TEM scale bar was read here (+/- 1 nm assumed) |
| Wong, Kakushima, Nanomaterials 12, 1739 (2022), doi 10.3390/nano12101739 | review (scaling analysis) | gate-oxide physical thickness 3 nm and metal gate 2 nm as the thinnest practical gate stack, hence "the intersheet/interwire spacing should not be less than 8 nm" and "the minimum intersheet spacing should be in the range of 7 to 8 nm"; 5-nm sheets in an 80 nm stack headroom; 8 nm sheets at 12 / 8 / 5.3 nm spacing for 4 / 5 / 6 sheets; the seven-level Leti stack (150 nm Si0.7Ge0.3/Si multilayer, 8 nm sheets) | analysis, not measurement; used for the scaled cell's 7 nm spacing and for the gate-stack budget |
| Mohapatra, Dash, Jena, Das, Maiti, SN Applied Sciences 3, 540 (2021), doi 10.1007/s42452-021-04539-y | TCAD design study | gate stack SiO2 0.5 nm + HfO2 1.5 nm, TiN gate; sheet thickness 4-6 nm, widths 3-12 nm (simulation) ; calibrated to the IBM device (Lg 12 nm, thickness 5 nm, width 50 nm) | assumed values; used only for the film budget (interfacial oxide about 1 nm, HfO2 1.5-2 nm) |
| Mertens et al. (imec), Symposium on VLSI Technology 2016, doi 10.1109/VLSIT.2016.7573416 | device paper | 8 nm diameter vertically stacked horizontal Si nanowires, replacement metal gate on bulk Si, Lg 24 nm (SS 65 mV/dec, DIBL 42 mV/V); imec's release of the paper gives 45 nm lateral pitch and 20 nm vertical separation | the pitch and separation are from imec's press text, not read from a TEM here |
| Ritzenthaler et al. (imec), IEDM 2018, doi 10.1109/IEDM.2018.8614528 | device paper | "vertically stacked gate-all-around Si nanowire CMOS transistors with reduced vertical nanowires separation, new work function metal gate solutions"; the separation value is in the paywalled text | title and metadata only; cited for the trend (the imec separation went down from 20 nm) |
| Bae et al. (Samsung), IEDM 2018, doi 10.1109/IEDM.2018.8614629 | device paper | 3 nm MBCFET: stacked nanosheets, variable sheet widths, "multiple threshold voltages ... despite challengingly limited vertical spacing between channels", SS 65 mV/dec | no dimensions in the abstract; the SemiAnalysis teardown of the SF2 product (2026, no DOI) reads sheet widths of 19-30 nm (NPU rows) and 37-50 nm (CU cell) and three sheets |
| Agrawal et al. (Intel), IEDM 2024, doi 10.1109/IEDM50854.2024.10873367 | device paper | silicon RibbonFET at 6 nm gate length and 45 nm CPP on a single-nanoribbon vehicle; cross-section TEMs at Tsi 5.5, 3.1 and 1.7 nm (1.5 nm in Intel's own summary); mobility holds to about 3 nm Tsi | ribbon spacing not stated; Intel 18A products carry four ribbons (Intel 18A technology brief; SemiAnalysis teardown) |
| Yeap et al. (TSMC), IEDM 2024, doi 10.1109/IEDM50854.2024.10873475 | device paper | N2 nanosheet platform; sheet thickness "tuned", no value | no dimensions published |
| Barraud et al. (Leti), IEDM 2017, doi 10.1109/IEDM.2017.8268473 | device paper | stacked nanowire/nanosheet FETs up to seven levels (the 150 nm Si/SiGe multilayer quoted by Wong and Kakushima) | text not accessible here (HAL blocked); cited through the review |
| Semiconductor Engineering, "New transistor structures at 3nm/2nm" (2018, no DOI) | trade press | IBM figure caption "Wsheet = 40 nm"; TEL: the inner spacer is "a critical ~5 nm recess etch with no line of sight"; Intel SiGe nanosheet pMOS with 5 nm sheets at 25 nm gate length | secondary; W = 40 nm is the value used for the IBM cell |
| Tirumaladass, KTH MSc thesis 2022 (DiVA diva2:1710391, no DOI) | TCAD thesis | restates Loubet's 3 x 5 nm sheets at 10 nm spacing, Lg 12 nm, CPP 44/48 nm; its own device: sheets 7 nm, Lg 13.5 nm, spacing 8 nm (Table 5.1); shows tungsten pinch-off between the sheets as the TiN thickness grows | secondary; the pinch-off remark is the only published statement found on film fill between sheets |

ALD film thickness measured on a GAA stack: none found. Searches for TEM-measured HfO2 or TiN thickness per sheet
(top / middle / bottom), for inner-spacer ALD conformality and for sheet-to-sheet work-function-metal thickness
returned patents (US 10084055, US 10297671: TiN "pinches off" the space between sheets so that no further work-function
metal is deposited there), a variability paper that names "metal thickness variation" as a Vt source without TEM
numbers (Harsha Vardhan et al. 2019, cited in the thesis) and process descriptions; no peer-reviewed measurement of
film thickness at several positions on a nanosheet stack was found. This is stated again under "Validation versus
application".

### 3D-NAND word-line stacks beyond 200 layers

| source (DOI) | type | what was read | uncertainty / remark |
|---|---|---|---|
| TechInsights, "YMTC 232L TLC 3D NAND" product brief (November 2022, no DOI; PDF table "Device comparison") | teardown | YMTC 232L (die from a HikSemi CC700 SSD): word-line pitch 48 nm, bit-line pitch 39 nm, vertical-channel height 12.4 um, slit depth 12.7 um, 253 total gates (232 active), 2 decks; the 128L generation: 58 nm, 39 nm, 8.5 um, 8.8 um | the table gives no channel-hole diameter; pitch to 1 nm, heights to 0.1 um as printed |
| TechInsights, "Comparison: latest 3D NAND products from YMTC, Samsung, SK hynix and Micron" (January 2023, no DOI), Table 1 | teardown | word-line pitch 44 / 43 nm (Samsung 128L / 176L), 52 / 45 nm (SK hynix 128L / 176L), 56 / 56 nm (Micron 128L / 176L), 58 nm (YMTC 128L); slit depths 6.7-11 um; bit-line pitch 38-39 nm | no hole diameter; products of 2020-2022 |
| imec, "Unlocking z-pitch scaling for next-generation 3D NAND flash" (7 November 2025, no DOI), citing Rachidi et al., IMW 2025, doi 10.1109/IMW61990.2025.11026936 | research article | "the pitch between adjacent word-lines ... today is about 40 nm"; products with over 300 layers and a stack about 30 um thick; imec test device at 30 nm pitch (15 nm gate, 15 nm oxide) with an 80 nm memory hole | research test device, not a product |
| imec (Rachidi et al.), IMW 2022, doi 10.1109/IMW52921.2022.9779303 | research device | 25 nm z-pitch with 10 nm word lines (already in wp4) | scaled research stack |
| imec, Ru and Mo word lines in 40 nm pitch 3D-NAND, Symposium on VLSI Technology 2021 (IEEE Xplore record 9508691; the DOI did not resolve here) | research device | 40 nm word-line pitch demonstrator with TEM of a 5-tier Mo stack | pitch only |
| Keller et al. (Onto Innovation), Semiconductor Digest, July 2023 (no DOI; already in wp4) | metrology note | 2XX-layer product: slit width about 250 nm, slit pitch above 2 um, channel holes about 100 nm, W recess about 50 nm | order-of-magnitude values from a vendor note |
| Antonelli (Onto Innovation), Semiconductor Digest, March 2021 (no DOI; SPIE counterpart doi 10.1117/12.2583786) | metrology note | channel holes of 100 nm diameter and more than 6 um depth as the metrology target | typical, not a specific product |
| Reiter et al., Solid-State Electronics 211, 108816 (2024), doi 10.1016/j.sse.2023.108816 (already in wp4) | device data + model | Lam's 76/98 ON-layer test structure (4-5 um, AR about 50): HfO2 step coverage 85 %, from which the Hf-precursor sticking 7.0e-3 was calibrated | the only published step-coverage number on a word-line stack found; it fixes s, so it cannot also validate the transport |
| Lee and Shin 2022 (doi 10.3390/mi13071139), Kihara et al. 2023 (doi 10.23919/VLSITechnologyandCir57934.2023.10185160), Shen et al. 2023 (doi 10.35848/1347-4065/accbc7), Oehrlein et al. 2024 (doi 10.1116/6.0003579) | TCAD / etch papers (in wp4) | word line 29 / oxide 22 nm; 10 um holes of about 100 nm for 400+ layers | as in wp4 |

The 2XX-layer cell built below takes the slit depth (12.7 um) and word-line pitch (48 nm) from the YMTC 232L
teardown table, the slit width (250 nm), slit pitch (above 2 um, so recesses of 1 um) and hole diameter (100 nm)
from the Onto note, and splits the 48 nm pitch into a 30 nm word line and an 18 nm oxide (wp4's 29/22 ratio at the
new pitch, rounded to 6 nm voxels); the hole pitch along the slit (200 nm) and the number of hole rows between slits
(four per side) are this work's assumptions consistent with those numbers. No published TEM of a 2XX-layer product
was available to read the recess length or the hole pitch from.

## Geometry table (as built; `meta` of each VoxelGeometry records the rounded values)

GAA cells (periodic along x, the sheet-width direction; the mouth is the top of the gate trench; the y ends are the
gate spacers / inner spacers; plain voxel faces, which are exact for this axis-aligned geometry, inventory 1.000):

| cell | source of the numbers | sheets N x t_s | slot g_s | sheet width W | gate trench Lg | stack-to-stack gap, period | trench above the top sheet H_top; below the lowest sheet g_bot | voxel h, grid, faces |
|---|---|---|---|---|---|---|---|---|
| gaa_ibm2017 | Loubet 2017 (N, t_s, g_s, Lg), W from the IBM figure caption; gap, H_top, g_bot assumed | 3 x 5 nm | 10 nm (8 voxels) | 40 nm | 12 nm (built 12.5) | 20 nm, period 60 nm | 40 nm; 10 nm | 1.25 nm, 48 x 12 x 69, 8 400 faces |
| gaa_scaled4 | four ribbons (Intel 18A), 7 nm spacing (Wong and Kakushima minimum), W 25 nm (inside the teardown range) | 4 x 5 nm (built 5.25) | 7 nm (8 voxels) | 25 nm (built 25.4) | 12 nm (built 12.25) | 15 nm, period 40.25 nm | 40 nm; 7 nm | 0.875 nm, 46 x 16 x 103, 12 556 faces |
| gaa_ibm2017_after_hk | the IBM cell with 1 nm interfacial oxide + 2 nm HfO2 already on every surface (the TiN step sees this; the solver has no narrowing) | 3 x 11 nm | 4 nm (8 voxels) | 46 nm | 6 nm | 14 nm, period 60 nm | 37 nm; 4 nm | 0.5 nm, 120 x 14 x 165, 36 864 faces |

NAND cells (periodic along the slit at the hole pitch; mouth = slit top; area correction with a 2-voxel normal
kernel so that the 5-6-voxel recesses are not blended across, site inventory 0.991 / 0.994 of the smooth structure;
the pillar walls alone carry 7 % too many sites at this kernel, the plates are exact):

| cell | source of the numbers | slit w_s x depth L_s | word line g_r / pitch p_r (recesses per wall) | recess length L_r | holes: diameter, pitch along the slit, rows per side | voxel h, grid, faces | reference without pillars |
|---|---|---|---|---|---|---|---|
| nand_wp4_pillars | wp4's w100_Lr365 cell (Lee and Shin 29/51 nm; Reiter 5 um; TechInsights 2014 730 nm slit-to-slit) + 100 nm holes (Onto 2023) at 200 nm pitch, two staggered rows | 100 nm x 5.0 um (AR 50) | 30 / 50 nm (100) | 365 nm | 100 nm, 200 nm, 2 (4 pillars per cell) | 5 nm, 40 x 168 x 1001, 1 189 600 faces | wp9 voxel recess run (10.70 / 10.88 s) and wp4 nested model (14.53 / 14.68 s) |
| nand_2xx | TechInsights YMTC 232L: slit depth 12.7 um, WL pitch 48 nm; Onto 2023: slit 250 nm, slit pitch > 2 um (L_r 1 um), holes 100 nm; WL/oxide split 30/18 assumed | 250 nm (built 252) x 12.7 um (AR 50) | 30 / 48 nm (264) | 1.0 um | 100 nm, 200 nm (built 198), 4 (8 pillars per cell) | 6 nm, 33 x 378 x 2118, 5 742 132 faces | nested model run here (none, plates) |

## Results

Conditions: 573 K; constant mouth density 0.3 Pa (n = 3.79e19 m^-3) or the chamber's centre-face density history of
the long AR-100 case (1 % TMA, 300 Pa, 20 s pulse; n0 = T_face x 3.79e20 m^-3, the same series as wp9 part `driven`);
TMA s 5e-3 / N0 4e18 / vbar 410 m/s, TDMAH 7e-3 / 2.77e18 / 185 m/s, TiCl4 0.1 / 1.24e18 / 253 m/s. Times are to 90 and
99 % area-weighted coverage of a face region; "everywhere" is the pointwise minimum over all regions (for NAND also
over every recess's deep half). Exposure = 0.3 Pa x time. The rule's pulse is 1.2 (1 + A_wall/A_wafer) N0/(U c0) at
1 % precursor (U c0 = 1.90e19 m^-2 s^-1) with the wafer wall-area ratio of the cell (GAA: the cell's exposed area plus
the flat ILD over the rest of a 44 nm CPP, per footprint; NAND: per slit-pitch footprint w_s + 2 L_r). "Step coverage
at the rule" is the slowest region's coverage over the top region's coverage at that pulse (GAA: top of the top sheet;
NAND: the first recess's deep half). The effusion limit is (A_wall/A_mouth) N0 / (n vbar/4) at 0.3 Pa.

### Table 1. GAA cells (figures/device_geometries.csv; every run M 50 000 particles per step, dt = t_flat/100)

t_flat is the flat-surface Langmuir time to 90 %, 2.3 N0/(s n vbar/4): 0.474 s TMA, 0.520 s TDMAH, 0.012 s TiCl4 at
0.3 Pa. Columns: time to 90 / 99 % in seconds.

| cell, chemistry, drive | top of top sheet | first slot, centre | bottom slot, centre (last everywhere) | bottom slot, edge | floor between stacks | everywhere t90 / t99 | exposure 90 / 99 % (Pa s) | wall-area ratio cell / wafer | rule pulse | step coverage at the rule (min / top) | wall time |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gaa_ibm2017, TMA, 0.3 Pa | 0.490 / 0.966 | 0.535 / 1.015 | 0.558 / 1.038 | 0.553 / 1.033 | 0.547 / 1.027 | 0.558 / 1.038 | 0.167 / 0.311 | 17.5 / 5.69 | 1.69 s | 0.999 (0.995 / 0.997) | 13 s |
| gaa_ibm2017, TMA, chamber | 1.974 / 3.050 | 2.094 / 3.138 | 2.152 / 3.184 | 2.137 / 3.172 | 2.123 / 3.161 | 2.152 / 3.184 | - | 17.5 / 5.69 | 1.69 s | 0.938 (0.789 / 0.841) | 35 s |
| gaa_ibm2017, TDMAH, 0.3 Pa | 0.546 / 1.068 | 0.615 / 1.144 | 0.649 / 1.179 | 0.640 / 1.170 | 0.631 / 1.161 | 0.649 / 1.179 | 0.195 / 0.354 | 17.5 / 5.69 | 1.17 s | 0.996 (0.990 / 0.994) | 13 s |
| gaa_ibm2017, TDMAH, chamber | 2.121 / 3.238 | 2.294 / 3.366 | 2.376 / 3.430 | 2.355 / 3.413 | 2.336 / 3.399 | 2.376 / 3.430 | - | 17.5 / 5.69 | 1.17 s | 0.823 (0.519 / 0.631) | 36 s |
| gaa_ibm2017, TiCl4, 0.3 Pa | 0.018 / 0.032 | 0.041 / 0.057 | 0.056 / 0.070 | 0.052 / 0.067 | 0.047 / 0.064 | 0.056 / 0.070 | 0.017 / 0.021 | 17.5 / 5.69 | 0.52 s | 0.995 (0.995 / 1.000) | 17 s |
| gaa_ibm2017, TiCl4, chamber | 0.130 / 0.208 | 0.258 / 0.346 | 0.341 / 0.413 | 0.317 / 0.395 | 0.295 / 0.380 | 0.341 / 0.413 | - | 17.5 / 5.69 | 0.52 s | 0.995 (0.995 / 1.000) | 8 s |
| gaa_scaled4, TMA, 0.3 Pa | 0.494 / 0.970 | 0.534 / 1.014 | 0.569 / 1.052 | 0.564 / 1.047 | 0.560 / 1.042 | 0.569 / 1.052 | 0.171 / 0.316 | 19.5 / 6.15 | 1.81 s | 0.998 (0.995 / 0.997) | 16 s |
| gaa_scaled4, TMA, chamber | 1.982 / 3.056 | 2.095 / 3.140 | 2.180 / 3.204 | 2.170 / 3.197 | 2.157 / 3.188 | 2.180 / 3.204 | - | 19.5 / 6.15 | 1.81 s | 0.941 (0.814 / 0.866) | 47 s |
| gaa_scaled4, TDMAH, 0.3 Pa | 0.550 / 1.074 | 0.613 / 1.142 | 0.665 / 1.199 | 0.659 / 1.192 | 0.652 / 1.184 | 0.665 / 1.199 | 0.200 / 0.360 | 19.5 / 6.15 | 1.25 s | 0.997 (0.992 / 0.995) | 23 s |
| gaa_scaled4, TDMAH, chamber | 2.133 / 3.247 | 2.296 / 3.368 | 2.419 / 3.462 | 2.404 / 3.450 | 2.386 / 3.437 | 2.419 / 3.462 | - | 19.5 / 6.15 | 1.25 s | 0.814 (0.536 / 0.659) | 62 s |
| gaa_scaled4, TiCl4, 0.3 Pa | 0.018 / 0.033 | 0.039 / 0.056 | 0.062 / 0.076 | 0.059 / 0.074 | 0.056 / 0.071 | 0.062 / 0.076 | 0.019 / 0.023 | 19.5 / 6.15 | 0.56 s | 0.995 (0.995 / 1.000) | 13 s |
| gaa_scaled4, TiCl4, chamber | 0.133 / 0.213 | 0.247 / 0.343 | 0.374 / 0.447 | 0.357 / 0.434 | 0.338 / 0.421 | 0.374 / 0.447 | - | 19.5 / 6.15 | 0.56 s | 0.995 (0.995 / 1.000) | 7 s |
| gaa_ibm2017_after_hk (4 nm slots, 6 nm trench), TiCl4, 0.3 Pa | 0.024 / 0.038 | 0.100 / 0.119 | 0.142 / 0.157 | 0.111 / 0.138 | 0.092 / 0.122 | 0.142 / 0.157 | 0.043 / 0.047 | 25.6 / 4.35 | 0.42 s | 0.995 (0.995 / 1.000) | 56 s |
| gaa_ibm2017_after_hk, TiCl4, chamber | 0.165 / 0.246 | 0.561 / 0.656 | 0.761 / 0.827 | 0.618 / 0.741 | 0.522 / 0.666 | 0.761 / 0.827 | - | 25.6 / 4.35 | 0.42 s | 0.131 (0.131 / 1.000) | 25 s |
| gaa_ibm2017_after_hk, TDMAH, 0.3 Pa | 0.576 / 1.105 | 0.796 / 1.343 | 0.898 / 1.450 | 0.846 / 1.396 | 0.812 / 1.362 | 0.898 / 1.450 | 0.269 / 0.435 | 25.6 / 4.35 | 0.94 s | 0.934 (0.914 / 0.979) | 28 s |
| gaa_ibm2017_after_hk, TDMAH, chamber | 2.198 / 3.301 | 2.709 / 3.693 | 2.913 / 3.852 | 2.809 / 3.774 | 2.745 / 3.725 | 2.913 / 3.852 | - | 25.6 / 4.35 | 0.94 s | 0.474 (0.237 / 0.500) | 58 s |

Convergence (gaa_ibm2017, TMA, 0.3 Pa; everywhere t90 / t99): M 50 000: 0.558 / 1.038 s; M 100 000: 0.558 / 1.039 s;
voxel 1.25 -> 0.625 nm (16 voxels per slot, 96 x 21 x 137, 33 072 faces): 0.563 / 1.046 s (+0.9 %; the finer grid builds
Lg as 11.9 instead of 12.5 nm and the wall-area ratio as 18.1 instead of 17.5); dt t_flat/100 -> t_flat/200: 0.556 /
1.037 s (-0.4 %). Every run: escaped + stuck = launched, left = 0; plain-face inventory 1.000 of the smooth cell.

Reading of Table 1.
- The GAA stack at the replacement-gate step is a nearly flat surface for ALD. The slowest point, the centre of the
  slot under the lowest sheet, reaches 90 % at 0.558 s against 0.474 s for a flat surface at the same density (+18 %)
  and 99 % at 1.038 s against 0.947 s (+10 %); the top of the top sheet (0.490 s) is itself 3 % behind the flat value
  because it shares the entering flux with the walls of the trench. The slot centres lag their edges by 1 % (TMA) to 8 %
  (TiCl4): with s = 0.1 the 10 nm slot of aspect ratio 2 (20 nm from each open edge) is mildly transport-limited, with
  s = 5e-3 it is not. The ordering is always top of the top sheet, then the slots from the top one down, then the
  bottom slot centre; the floor between the stacks sits between the first and the last slot.
- Exposure for 90 / 99 % everywhere at 0.3 Pa: TMA 0.17 / 0.31 Pa s, TDMAH 0.20 / 0.35 Pa s, TiCl4 0.017 / 0.021 Pa s
  (the IBM cell); 2 to 3 % more for the four-sheet 7 nm-slot cell and 60 % (TMA-like chemistries) to 2.5 x (TiCl4) more
  for the narrowed cell with 4 nm slots and a 6 nm trench. For comparison the AR-100 hole needs 4.9 Pa s for 90 % (wp9: 16.45 s at
  0.3 Pa). The effusion limit of the cells is 0.018 / 0.028 / 0.009 s (TMA / TDMAH / TiCl4): the cells saturate at
  30 x (TMA) and 6 x (TiCl4) their stoichiometric minimum, i.e. they are reaction-limited like a flat surface, not
  transport-limited like the hole.
- Wall-area ratio: the cell carries 17.5 (IBM), 19.5 (four sheets) or 25.6 (narrowed) times its mouth area, but the
  mouth is only the 12 nm gate trench of a 44 nm CPP; per wafer area a wafer tiled with these devices has 5.7 / 6.2 /
  4.4 times its area in exposed surface (the trench walls dominate: the three sheets of the IBM cell add 1.4 wafer areas,
  the 85 nm deep trench 3.3). The rule's TMA pulse for such a wafer is 1.69 / 1.81 / 1.65 s at 1 %, against 14.4 s for the
  AR-100 hole wafer and 0.25 s for a blanket wafer; HfO2 1.17 / 1.25 s, TiN 0.52 / 0.56 s. At a constant 0.3 Pa mouth
  density the cells are at 99.5 % everywhere well before the rule's pulse (step coverage 0.995 to 0.999), which is the
  rule's claim (the pulse set by supply covers the feature's own transport) carried over to a feature with a 2:1 slot.
- Chamber drive: under the AR-100 wafer's centre-face history (a wafer with ten times the GAA wafer's area, which
  holds the surface density at a few per cent of c0 for the first seconds) the IBM cell reaches 90 % everywhere at 2.15 s
  (TMA), 2.38 s (TDMAH) and 0.34 s (TiCl4), 99 % at 3.18 / 3.43 / 0.41 s, with the same ordering of regions and lags of
  8-9 % (TMA, TDMAH) and 160 % (TiCl4) between the top sheet and the bottom slot centre; the flat top of that wafer
  reached 90 % at 1.95 s in the chamber run, so the stack costs 10 % of pulse time on top of a flat surface under this
  drive. The "step coverage at the rule's pulse" column of the chamber rows (0.94 TMA, 0.82 TDMAH, 0.13 TiCl4 in the
  narrowed cell) compares the GAA wafer's own rule pulse with a history that belongs to a different, much more loaded
  wafer, in which the density at 0.4-1.7 s is still rising; it says that the history is the wrong one for this wafer,
  not that the rule fails. A chamber run with the GAA wafer's wall-area ratio (5.7) would have a surface density close to
  the blanket wafer's and would saturate the cell in about the rule's 1.7 s; it was not run here (the patternedSurface
  boundary condition carries a hole model, not this cell).
- Narrowed cell (the TiN step after 1 nm interfacial oxide and 2 nm HfO2 on every surface, slots 4 nm, trench 6 nm):
  the TiCl4 exposure for 90 % everywhere rises from 0.017 to 0.043 Pa s and the slot centre lags its edge by 28 %
  (0.142 vs 0.111 s); still reaction-limited (3 x the effusion limit) and still saturated at 99.5 % within the rule's
  0.42 s at 0.3 Pa. Beyond this the solver's lack of film narrowing matters: 2-3 nm of TiN on both faces of a 4 nm slot
  pinches it off, which is the state the patents and the KTH thesis describe for the work-function metal, and the model
  then predicts the pinch-off order: slot edges before slot centres, top slot before bottom slot, by the lags above.

### Al2O3 / HfO2 / TiN thickness after N cycles (application; no measurement to compare with)

With the cells saturated to 99.5 % everywhere at the rule's pulse, the film after N cycles is N x GPC everywhere
within 0.5 %: 50 cycles of TMA/H2O give 5.13 nm on the top sheet and 5.12 nm at the bottom slot centre (GPC
0.103 nm); HfO2 at 0.10 nm per cycle 4.97 / 4.95 nm; TiN at 0.030 nm per cycle 1.500 / 1.493 nm
(figures/device_geometries.csv, columns `film_50cyc_at_rule_*`). A pulse one third of the rule's (0.56 s of TMA at
0.3 Pa) leaves 92.8 % on the top sheet and 90.0 % at the bottom slot centre, i.e. 4.78 against 4.64 nm after 50
cycles: a 0.14 nm top-to-bottom difference, below what a TEM scale bar resolves on a 5 nm film. The measurement that would test the prediction is a HfO2 or TiN
thickness read by TEM/EELS at the slot centre and slot edge of the lowest sheet after a deliberately short pulse
(one third of saturation), where the model predicts 3-8 % (TMA-like chemistries) to 30 % (TiCl4 in a 4 nm slot) less
film at the centre; at a saturating pulse it predicts no measurable difference, which is consistent with the uniform
gate stacks shown in the device papers' TEMs but does not test the transport.

### Table 2. The 5 um word-line cell with channel-hole pillars (nand_wp4_pillars; h 5 nm, dt 0.1 s, M 100 000)

Regions: slit floor (the 100 x 200 nm floor of the slit itself); recess "deep half" = every face of a recess more
than L_r/2 from the slit (ceiling, floor, pillar walls, end wall), for the first ("top"), the middle and the deepest
recess; "slowest recess" = the recess whose deep half is the last to pass 90 % (the pointwise minimum over the 100
recesses is what "everywhere" uses). Times in s to 90 / 99 %.

| chemistry, drive | slit floor | top recess | mid recess | deepest recess (end wall) | slowest recess = everywhere | exposure 90 / 99 % (Pa s) | nested model, deepest recess: wp4 diffusivity / parallel-plate (nominal 29/51 nm dims) | wp9 voxel without pillars: deepest recess end / min over recesses | wall-area ratio cell / wafer; rule; effusion limit | step coverage at the rule | wall time |
|---|---|---|---|---|---|---|---|---|---|---|---|
| TMA, 0.3 Pa | 9.04 / 11.39 | 0.75 / 1.30 | 7.09 / 8.77 | 10.33 / 11.78 (10.37) | 11.30 / 12.13 | 3.39 / 3.64 | 14.68 / 10.17 (slit bottom 14.53 / 10.03) | 10.88 / 12.06 | 1422 / 171; 43.5 s; 1.46 s | 0.995 | 310 s |
| TMA, chamber | 10.45 / 11.88 | 2.58 / 3.61 | 9.16 / 10.27 | 11.25 / 12.10 (11.27) | 11.85 / 12.32 | - | 13.53 / 11.07 (13.46 / 10.99) | - | same | 0.998 | 238 s |
| TDMAH, 0.3 Pa | 12.42 / 16.59 | 0.91 / 1.55 | 10.23 / 12.29 | 15.29 / 17.49 (15.38) | 17.10 / 18.08 | 5.13 / 5.42 | 22.34 / 15.42 | - | 1422 / 171; 30.1 s; 2.25 s | 0.996 | 420 s |
| TDMAH, chamber | 12.47 / 14.34 | 2.91 / 3.97 | 11.22 / 12.41 | 13.87 / 14.51 (13.92) | 14.45 / 14.63 | - | 15.02 / 13.87 | - | same | 0.999 | 243 s |
| TiCl4, 0.3 Pa | 0.82 / 1.49 | 0.26 / 0.30 | 3.68 / 3.86 | 4.38 / 4.95 (4.55) | 6.80 / 6.94 | 2.04 / 2.08 | 8.00 / 5.57 | - | 1422 / 171; 13.5 s; 0.74 s | 1.000 | 198 s |
| TiCl4, chamber | 2.62 / 3.77 | 0.75 / 0.89 | 6.13 / 6.30 | 6.95 / 7.45 (7.09) | 8.88 / 8.99 | - | 9.33 / 7.60 | - | same | 0.998 | 127 s |

Convergence (TMA, 0.3 Pa, everywhere t90 / t99): M 100 000: 11.30 / 12.13 s; M 200 000: 11.28 / 12.11 s (-0.2 %);
h 5 -> 10 nm (recess 3 voxels, pillars 10 across, 20 x 84 x 501, 297 800 faces): 10.99 / 11.82 s (-2.7 %; the coarse
grid also changes the inventory to 0.975). Site inventory 0.991 of the smooth cell; every particle bookkept.

Where the last point sits (figures/device_geometries_recess_profile.csv, a dedicated run at M 50 000 recording all
100 recesses; figure nand_cell_coverage.png, top right): the saturation time rises almost linearly with depth from
0.74 s at the first recess to 11.29 s at recesses 85-89 (4.3-4.5 um of the 5 um stack) and then falls again to
10.33 s at the deepest recess next to the floor. The last point of the structure to saturate is therefore not the
bottom but 0.5-0.7 um (five to seven slit widths) above it: once the slit floor and the lowest recess mouths have
saturated, molecules that pass the bottom recesses on the way down are reflected by the floor and get a second
chance on the way up, whereas a recess higher in the front loses the molecules that pass it to the recesses below.
The same profile was present in wp9's run without pillars, where the 11 % gap between the deepest recess (10.88 s)
and the minimum over recesses (12.06 s) was attributed to noise; the per-recess record here shows it is systematic
(neighbouring recesses agree to 0.05 s), so the wp9 sentence should be corrected (proposal below). The 1D nested model
cannot show it: its slit density is monotonic in depth, so its last recess is always the deepest.

Reading of Table 2.
- Pillars versus no pillars (TMA, 0.3 Pa): the deepest recess end saturates at 10.37 s with the pillars against
  10.88 s without (wp9), 5 % sooner, and the slowest recess at 11.30 against 12.06 s. The pillars remove more area
  than they add (a 100 nm pillar in a 30 nm recess takes two 7 854 nm2 discs out of the plates and adds 9 425 nm2 of
  wall), 8 % of the cell's area, and obstruct the lateral path into the recess; the net effect at this recess length
  is the area, with the obstruction worth about 3 % (area alone would give 7-8 %).
- Against the nested 1D model: with wp4's hydraulic-diameter slit diffusivity the model is 30 % slow (14.68 vs
  11.30 s for the last recess; 14.53 vs 10.37 s for slit bottom vs deepest recess end), the 24-27 % of wp9 plus the
  pillars' 5 %. With the parallel-plate correction of the refit package the nested model gives 10.17 s for the deepest
  recess, 2 % under the voxel deepest recess (10.37 s) and 10 % under the voxel's true last point (11.30 s), which the
  nested model cannot represent. For TDMAH the same holds (15.42 vs 15.38 s deepest; 17.10 s last point). For TiCl4
  (s = 0.1, a sharp front and recess filling limited by the recess's own transport) the corrected nested model gives
  5.57 s against the voxel deepest recess 4.55 s and last point 6.80 s: the recess-level transport differs more when
  the recess, not the slit, is the bottleneck.
- Under the chamber drive the ordering and the lags are the same (floor 10.45 s, deepest recess 11.27 s, last point
  11.85 s for TMA; nested 13.46 / 13.53 s with wp4's diffusivity, 10.99 / 11.07 s corrected). The exposure for
  90 % everywhere at 0.3 Pa is 3.39 Pa s (TMA), 5.13 (TDMAH) and 2.04 Pa s (TiCl4), 7.7 / 7.6 / 9.2 times the
  effusion limit of the cell.
- Wall-area ratio and rule: the cell carries 1 422 times its slit-mouth area and a wafer tiled with it at the 830 nm
  slit pitch carries 171 wafer areas (wp4's plain recessed slit: 178), so the rule's pulse is 43.5 s (TMA), 30.1 s
  (HfO2) and 13.5 s (TiN) at 1 %, three to four times the constant-0.3 Pa saturation time; at the rule's pulse the
  structure is at 99.5 % everywhere (step coverage 0.995-1.000).

### Table 3. The 2XX-layer cell (nand_2xx: 250 nm slit, 12.7 um, 48 nm pitch, 1 um recesses with four rows of 100 nm pillars; h 6 nm, dt 0.25 s, M 60 000)

| chemistry, drive | slit floor | top recess | mid recess | deepest recess (end) | slowest recess: index of 264 (depth fraction), t90 / t99 = everywhere | exposure 90 / 99 % (Pa s) | nested model, deepest recess: wp4 diffusivity / parallel-plate (rounded dims) | wall-area ratio cell / wafer; rule; effusion limit | step coverage at the rule | wall time |
|---|---|---|---|---|---|---|---|---|---|---|
| TMA, 0.3 Pa | 14.60 / 23.63 | 1.52 / 2.19 | 19.21 / 21.71 | 28.18 / 32.26 (28.76) | 236 (0.90): 33.09 / 34.20 | 9.93 / 10.26 | 43.55 / 29.94 (slit bottom 41.46 / 28.67) | 4009 / 448; 113 s; 4.12 s | 0.997 | 413 s |
| TMA, chamber | 13.55 / 15.30 | 3.82 / 4.67 | 14.83 / 15.21 | 15.78 / 16.47 (15.91) | 228 (0.87): 16.89 / 17.19 | - | 17.03 / 15.69 (16.84 / 15.57) | same | 1.000 | 155 s |
| TDMAH, 0.3 Pa | 17.78 / 30.22 | 2.04 / 2.80 | 28.42 / 31.19 | 41.31 / 47.20 (42.24) | 233 (0.88): 50.40 / 51.60 | 15.12 / 15.50 | 66.17 / 45.46 | 4009 / 448; 78.5 s; 6.33 s | 0.997 | 607 s |
| TDMAH, chamber | 14.55 / 16.08 | 4.53 / 5.41 | 16.20 / 16.60 | 17.82 / 18.57 (17.96) | 228 (0.87): 19.04 / 19.24 | - | 19.28 / 17.23 | same | 1.000 | 138 s |
| TiCl4, 0.3 Pa | 1.25 / 2.35 | 1.44 / 1.71 | 16.91 / 17.46 | 23.15 / 24.68 (24.16) | 231 (0.88): 29.88 / 30.37 | 8.96 / 9.11 | 26.86 / 18.95 | 4009 / 448; 35.1 s; 2.07 s | 1.000 | 339 s |
| TiCl4, chamber | 2.99 / 4.43 | 2.62 / 2.95 | 16.41 / 16.91 | 20.25 / - (21.25) | not reached: from recess 172 (65 % depth) down the recesses are below 90 % when the pulse ends, except the two next to the floor (deepest 0.97); slowest (index 170-200) at 0.09; 13.8 % of the area below 90 % at 22 s | - | 14.96 / 13.54 | same | 0.09 | 170 s |

Convergence (TMA, 0.3 Pa, everywhere t90 / t99): M 60 000: 33.09 / 34.20 s; M 120 000: 32.81 / 33.90 s (-0.8 %); h 6 -> 8 nm
(word line 32 / oxide 16 nm, 25 x 284 x 1589, 3.0e6 faces): 33.90 / 34.99 s (+2.4 %; inventory 0.994 -> 0.996, wall-area
ratio 448 -> 450). Slowest recess at index 230-236 in all three. Every particle bookkept (left = 0). Longest run 834 s
(M 120 000), 3 300-5 700 particles/s per thread, 6 800-12 000 voxel segments per particle.

Reading of Table 3 and comparison with wp4's nested-model numbers.
- The 2XX cell needs 9.9 Pa s of TMA for 90 % everywhere at 0.3 Pa (33.1 s), 2.9 x the 5 um cell's 3.39 Pa s, for
  2.5 x the depth, 2.5 x the slit width and 2.7 x the recess length. The slit floor saturates at 14.6 s and the deepest
  recess at 28.2 s: with a 250 nm slit and 1 um recesses the recess filling, not the slit, is the bottleneck (the
  TiCl4 case makes it plain: slit floor 1.25 s, mid recess 16.9 s). Against the effusion limit the structure is
  transport-limited by 8.0 x (TMA), 8.0 x (TDMAH) and 14 x (TiCl4).
- Nested model: with wp4's hydraulic-diameter diffusivities the deepest recess is at 43.6 s (TMA), 32 % above the
  voxel's last point and 55 % above the voxel's deepest recess; with the parallel-plate correction 29.9 s, 6 % above
  the voxel deepest recess (28.2 s) and 10 % below the voxel's last point (33.1 s). The corrected nested model is the
  better reference for the deepest recess in all three chemistries (TDMAH 45.5 vs 41.3 s, TiCl4 19.0 vs 23.2 s) and
  its remaining error is of the size of the pillars' effect; the point it cannot give is the 10-15 % between the
  deepest recess and the true last point at 87-90 % depth.
- Chamber drive (the AR-100 wafer's history): the cell reaches 90 % everywhere at 16.9 s (TMA) and 19.0 s (TDMAH), and
  TiCl4 does not get there before the pulse ends at 20 s (174 of 264 recesses saturated: the first 172 and the two
  next to the floor, which the floor-reflection effect fills while 65-98 % depth is at 9-90 %). The nested model under the same drive
  gives 17.0 / 19.3 / 15.0 s with wp4's diffusivity and 15.7 / 17.2 / 13.5 s corrected; the chamber drive's rising
  density after 6 s (above 0.3 Pa from 7 s on) compresses the differences between the models. The TiCl4 difference
  (nested 13.5-15.0 s, voxel not reached by 20 s) is the recess-level transport at s = 0.1 again.
- Wall-area ratio and rule: 4 009 slit-mouth areas per cell, 448 wafer areas for a wafer tiled at the 2.25 um slit
  pitch (wp4's plain 250 nm / 1 um recessed slit: 176; the 48 nm pitch with 30 nm word lines, 264 recesses in 12.7 um
  and the 1 um recesses on both sides make 2.5 x the area per wafer area). The rule's pulse at 1 % is 113 s (TMA),
  78.5 s (HfO2) and 35.1 s (TiN), 3.4 / 1.6 / 1.2 x the constant-0.3 Pa saturation time; at the rule's pulse the cell is
  at 99.7-100 % everywhere (step coverage 0.997-1.000). The TiN rule (35 s) is only 17 % above the 29.9 s the cell
  needs at 0.3 Pa, the smallest margin in this package: for a chloride step on such a wafer the recess transport is
  close to mattering, and wp12's HCl re-adsorption would add to it.
- Where the last point sits: at 87-90 % of the depth in every run of both cells (recess 86 of 100 in the 5 um cell,
  228-236 of 264 here), 1.2-1.5 um above the floor for the 250 nm slit and 0.5-0.7 um for the 100 nm slit, i.e. about
  five slit widths; the deepest recess saturates 9 % (5 um cell), 15 % (TMA, 2XX cell) to 22 % (TiCl4, 2XX cell) sooner
  than that point. The same turn-up exists in the plain AR-100 hole of wp9 (figures/voxel_mc_profiles.csv: at 5 s the
  voxel wall coverage rises from 0.0022 at 87.5 % depth to 0.0041 at the bottom while the 1D continuum profile keeps
  falling, 0.0014 to 0.0006). The reason is the reflecting floor: once it has saturated it returns the whole downward
  stream, so a recess mouth within a few slit widths of the floor is hit from below as well as from above, whereas a
  mouth higher up receives from below only what the absorbing recesses beneath it let through. The 1D model carries
  this only as dn/dz = 0 at the floor and keeps the deepest recess slowest. A TEM of the lowest ten word lines against
  those 1-1.5 um higher, after a deliberately short pulse (two thirds of saturation), would test it: the voxel model
  predicts the lowest word lines ahead of the ones above them by the equivalent of 10-20 % of the dose.

## Validation versus application

Validation in this package is limited to what published data support: (i) the solver itself was validated in wp9
against Machball (ballistic, within 2 %) and against the 1D solver for holes, slits and the plain recessed slit; (ii)
for the 5 um word-line cell the only external anchor is wp4's calibration of s_TDMAH = 7e-3 to the 85 % HfO2 step
coverage of Reiter et al. 2024 on Lam's 76/98-layer structure, and because that datum fixed s it cannot also test the
transport; (iii) no published film-thickness profile on a GAA stack (per sheet, slot centre vs edge, or inside an
inner spacer) was found, and no TEM of a 2XX-layer word-line stack with film thickness per level is public. Everything
in Tables 1-3 is therefore application of a solver validated elsewhere, with the geometry taken from the sources in
the tables and the gaps (GAA trench height, stack-to-stack gap, hole pitch and rows, the 30/18 nm split of the 48 nm
pitch) stated as assumptions. The predictions that a measurement would test are: the 3-8 % (TMA-like) to 30 % (TiCl4
in a 4 nm slot) slot-centre deficit after a pulse one third of saturation on a nanosheet stack, and the 10-20 %
dose-equivalent lead of the lowest word lines over those five slit widths above them after a two-thirds pulse on a
word-line stack. At saturating pulses both models predict uniform films, consistent with the uniform gate stacks and
word lines in the device papers' cross-sections, which is consistency, not validation.

## Figures

- `figures/device_geometries.png`: the five cells as occupancy slices: the three GAA cells across the sheet width at
  mid-gate (black Si sheets suspended in the white gate trench over the substrate; period 60 / 40 / 60 nm), the plan
  view inside a recess of both NAND cells (slit in the middle, recesses either side, pillar discs in staggered rows,
  x stretched), and a vertical section through a pillar row of the 2XX cell for the first six word lines. Caption:
  "Device geometries as voxel occupancy arrays: GAA nanosheet unit cells (3 x 5 nm sheets at 10 nm, 4 x 5 nm at 7 nm,
  and the first after 3 nm of oxide on every surface) and 3D-NAND word-line cells with channel-hole pillars (5 um,
  100 nm slit; 12.7 um, 250 nm slit, 48 nm pitch)."
- `figures/gaa_coverage.png`: top row, coverage maps of the underside of the lowest sheet of the narrowed IBM cell
  under TiCl4 at 0.3 Pa at t90/3, 2 t90/3 and t90 (47, 94, 142 ms): the film grows in from the two open edges toward
  the slot centre; bottom row, region histories of the IBM cell under TMA at 0.3 Pa, the slowest region of all GAA
  cells and chemistries at 0.3 Pa (log time), and the IBM cell under the chamber's centre-face history. Caption:
  "GAA stacks fill from the sheet edges inward; the slot centre under the lowest sheet is the last point, 18 % (TMA)
  to 5 x (TiCl4) behind a flat surface at 90 % and within the rule's pulse in every case."
- `figures/nand_cell_coverage.png`: histories of slit floor, deepest recess and slowest recess for both NAND cells,
  three chemistries, constant 0.3 Pa (with the two nested-model references) and chamber-driven; the per-recess
  saturation time along the 5 um stack (dedicated M 50 000 run) and the 2XX cell; the slice coverage along the depth
  of the 2XX cell at the end of each run (the slice mean alternates between saturated oxide/slit-wall slices and
  unsaturated recess slices where the TiCl4 chamber run stopped, and recovers next to the floor). Caption: "Word-line cells with pillars: the last point to saturate sits at
  87-90 % of the depth, five slit widths above the floor; the nested 1D model with the parallel-plate correction
  matches the deepest recess within 2-6 % (TMA, TDMAH) and misses the last point by 10-15 %."

## Timings

35 voxel runs, 5.9e8 particles, 100 min of wall time on 6 threads (ALDMC_THREADS=6) of an Apple M4; the longest run
834 s (2XX cell, M 120 000); GAA runs 7-62 s (130 000-210 000 particles/s per thread, 130-300 voxel segments per
particle), 5 um NAND cell 127-706 s (5 800-12 400 p/s/thread, 3 300-5 300 segments), 2XX cell 138-834 s (3 300-5 700
p/s/thread, 6 800-12 000 segments; 5.7e6 faces, 1.5 GB resident). The nested-model references took 24 s in all; the
dedicated per-recess profile run 7.5 min on 4 threads while the main runs were going. Geometry creation 1-2 s for the
2XX cell. pytest: 129 passed in 13 s (5 new; `PYTHONPATH=src .venv/bin/python -m pytest -q tests`).

## Limits

- No film narrowing: the TiN step was run on an explicitly narrowed GAA cell (3 nm removed from every surface), and
  the pinch-off of a 4 nm slot by 2-3 nm of TiN is outside the model; for the NAND cells the films (a few nm in 30 nm
  recesses) do not narrow the transport paths appreciably within one half-cycle.
- The chamber drive is the AR-100 wafer's history, not the history a GAA wafer (wall-area ratio 5.7) or a 2XX NAND
  wafer (448) would produce; the GAA cell would saturate in about the rule's 1.7 s under its own loading, the NAND
  wafer's pulse is the rule's 113 s. The patternedSurface boundary condition carries a hole model, so neither wafer
  was run in the chamber.
- Site inventory of the NAND cells 0.991-0.994 with the 2-voxel normal kernel (the plates exact, the pillar walls
  7 % over); GAA cells exact (plain faces).
- Per-face coverages are noisy (a 6 nm face has 144 TMA sites and a TiCl4 particle carries 15 % of a face's sites);
  every number above is an area-weighted mean over at least 330 faces (recess end walls) and usually thousands, and
  the slowest-recess profiles are smooth to 0.05 s between neighbours.
- Geometry assumptions stated in the tables; the 2XX cell's hole pitch, rows and word-line/oxide split are not from a
  measurement.

## Correction to wp9 (proposal)

wp9, part (d) and the recess run: "(The dotted 'slowest recess end' curve in the figure is the pointwise minimum over
98 noisy recess histories, 12.06 s; it is biased late by the noise and the deepest recess is the honest number.)"
should read: "The slowest recess is not the deepest: the per-recess saturation time peaks at 87-90 % of the depth
(wp10, figures/device_geometries_recess_profile.csv) because the saturated floor reflects the downward stream back
into the lowest recesses; the minimum over recesses (12.06 s) is that recess, not noise, and the deepest recess
(10.88 s) is 10 % ahead of it."

## Proposed paper text (paper/draft.md, new subsection after 3.13)

```
### 3.14 Gate-all-around stacks and a current NAND cell

The voxel solver of Section 2.6 takes any shape, so two device geometries that the one-dimensional models cannot
describe were run with it: the replacement-gate step of a stacked-nanosheet transistor and a 3D-NAND word-line cell
with the channel-hole pillars standing in the recesses. The nanosheet cell follows the IBM device of Loubet et al.
(three 5 nm sheets at 10 nm spacing, 12 nm gate, 40 nm sheets) under a 40 nm gate trench, periodic across the sheet
width at a 60 nm pitch, with a four-sheet 7 nm-spacing variant and the IBM cell after 3 nm of oxide on every surface
(4 nm slots, 6 nm trench) for the TiN step; the NAND cells are the word-line slit of Section 3.11 with two rows of 100
nm pillars per side, and a 2XX-layer cell built from a teardown table (12.7 um slit depth, 48 nm word-line pitch) and a
metrology note (250 nm slit, 1 um recesses, 100 nm holes). Each was run with the TMA, TDMAH and TiCl4 parameters at a
constant 0.3 Pa and under the centre-face history of the AR-100 wafer, with 8 voxels across every slot or recess,
particle doublings and grid changes moving the saturation times by at most 1 % (nanosheets) and 3 % (NAND).

The nanosheet stack is, for ALD, a flat surface with a short slot: the last point, the centre of the slot under the
lowest sheet, reaches 90 % at 0.56 s of TMA at 0.3 Pa against 0.47 s for a flat surface, and 99 % at 1.04 against
0.95 s; the slots fill from their two open edges inward with a centre-to-edge lag of 1 % (TMA) to 8 % (TiCl4), 28 % in
the 4 nm slots. A wafer tiled with such devices carries 5.7 wafer areas of surface (the gate trench walls, not the
sheets, dominate), so the rule gives a 1.7 s TMA pulse at 1 %, by which time the cell is at 99.5 % everywhere; the
exposure for 90 % everywhere is 0.17 Pa s (TMA), 0.20 (TDMAH) and 0.017 Pa s (TiCl4), 30 x and 6 x the stoichiometric
minimum, i.e. reaction-limited like a blanket wafer. The word-line cells are the opposite case. The 5 um cell with
pillars saturates its last recess at 11.3 s at 0.3 Pa (3.4 Pa s), 6 % sooner than without pillars because a 100 nm
pillar in a 30 nm recess removes more plate area than it adds in wall; the 2XX cell needs 33 s (9.9 Pa s of TMA,
15 Pa s of TDMAH, 9 Pa s of TiCl4), eight times its effusion limit, and a wafer of it carries 448 wafer areas, so
the rule puts its TMA pulse at 113 s at 1 %, 3.4 x the feature's own need; for TiN the margin is only 17 %. The nested
1D model of Section 3.11 with the parallel-plate slit diffusivity reproduces the deepest recess within 2-6 % for TMA
and TDMAH (10.2 vs 10.4 s; 29.9 vs 28.2 s) and is 18-29 % off for TiCl4, where the recess rather than the slit limits.
What the 1D model cannot give is where the last point sits: in every voxel run it is at 87-90 % of the depth, about
five slit widths above the floor, and the deepest recess is 9-22 % ahead of it, because the saturated floor returns the
downward stream into the lowest recesses; the same turn-up is present at the bottom of the AR-100 hole of Section
2.6 and absent from the continuum profile.

These are applications of a solver validated in Section 2.6, not validations: no published measurement of film
thickness per sheet or per word line was found for either structure (the one step-coverage datum on a word-line
stack, Reiter et al.'s 85 % HfO2, fixed the TDMAH sticking coefficient in Section 3.11 and cannot also test the
transport), and at saturating pulses both solvers predict uniform films, which is what the device papers' cross-
sections show and what any model would. The predictions that would test the transport are the 3-8 % (TMA-like) to
30 % (TiCl4, 4 nm slot) deficit at the slot centre of the lowest sheet after a pulse one third of saturation, and the
10-20 % dose-equivalent lead of the lowest word lines over those five slit widths above them after a two-thirds
pulse; a TEM or EELS thickness at those positions after a deliberately short pulse is the experiment.
```

## README text (one paragraph for the "Night three/four" results list)

```
- **Device geometries in the voxel solver.** A stacked-nanosheet gate (IBM's 3 x 5 nm sheets at 10 nm, a four-sheet
  7 nm variant, and the 4 nm slots left for TiN after the high-k) and 3D-NAND word-line cells with the channel-hole
  pillars in the recesses (the 5 um cell of night three and a 2XX-layer cell: 12.7 um slit, 48 nm pitch, 1 um
  recesses) run with TMA, TDMAH and TiCl4 at 0.3 Pa and under the chamber's centre-face history. The nanosheet stack
  is a flat surface with a 2:1 slot (last point 18 % behind a blanket at 90 %, saturated well inside the rule's 1.7 s
  pulse for a 5.7x wafer); the 2XX NAND cell needs 9.9 Pa s of TMA and a 448x wafer, 113 s by the rule; the last point
  of a word-line stack sits five slit widths above the floor, not at it, which the 1D nested model cannot show
  (`figures/device_geometries.png`, `figures/gaa_coverage.png`, `figures/nand_cell_coverage.png`,
  `paper/wp/wp10_device_geometries.md`). No published per-sheet or per-word-line film thickness exists to compare
  with: application, not validation.
```

## CLAIMS rows (CLAIMS.md format)

```
| GAA nanosheet cell (3 x 5 nm sheets at 10 nm, W 40 nm, Lg 12 nm, 60 nm pitch) in the voxel solver at 0.3 Pa, 573 K | last point (bottom slot centre) t90 / t99 0.558 / 1.038 s TMA (flat 0.474 / 0.947 s), 0.649 / 1.179 s TDMAH, 0.056 / 0.070 s TiCl4; exposure 90 % 0.17 / 0.20 / 0.017 Pa s; four-sheet 7 nm cell +2-3 %; 4 nm slots after 3 nm oxide: TiCl4 0.142 s, centre 28 % behind edge; M x2, h/2, dt/2 within 1 % | figures/device_geometries.csv, figures/gaa_coverage.png | PYTHONPATH=src .venv/bin/python coupling/device_geometries.py --part gaa,maps |
| GAA wafer wall-area ratio and rule | 5.7 (IBM cell, 44 nm CPP), 6.2 (four sheets), 4.4 (narrowed): rule 1.69 / 1.17 / 0.52 s (TMA / HfO2 / TiN) at 1 %; cells at 99.5 % everywhere at the rule's pulse (step coverage 0.995-0.999) | figures/device_geometries.csv | same |
| GAA cell under the AR-100 wafer's centre-face history | 90 / 99 % everywhere at 2.15 / 3.18 s TMA, 2.38 / 3.43 s TDMAH, 0.34 / 0.41 s TiCl4 (flat top of that wafer 1.95 s) | figures/device_geometries.csv | same |
| 5 um word-line cell with channel-hole pillars (100 nm slit, 30/50 nm, 365 nm recesses, 2 rows of 100 nm holes) at 0.3 Pa | slit floor / deepest recess / slowest recess t90 9.04 / 10.37 / 11.30 s TMA (wp9 without pillars 10.70 / 10.88 / 12.06 s; nested model 14.68 s, 10.17 s with the parallel-plate correction); TDMAH 12.42 / 15.38 / 17.10 s; TiCl4 0.82 / 4.55 / 6.80 s; chamber-driven TMA 10.45 / 11.27 / 11.85 s; M x2 -0.2 %, h 10 nm -2.7 % | figures/device_geometries.csv, figures/nand_cell_coverage.png | same, --part nand5,refs |
| Where the last point of a word-line stack sits | recess 86 of 100 (4.3 um of 5 um) at 11.29 s vs the deepest at 10.33 s; 2XX cell recess 228-236 of 264 (87-90 %), deepest 9-22 % ahead; the saturated floor reflects the downward stream into the lowest recesses | figures/device_geometries_recess_profile.csv, figures/device_geometries.csv (t90_slowest_recess, t90_deepest_recess) | same |
| 2XX-layer word-line cell (250 nm slit, 12.7 um, 48 nm pitch, 1 um recesses, 4 rows of 100 nm holes) at 0.3 Pa | 90 / 99 % everywhere 33.1 / 34.2 s TMA (9.9 Pa s), 50.4 / 51.6 s TDMAH, 29.9 / 30.4 s TiCl4; deepest recess 28.2 / 41.3 / 23.2 s vs nested 43.6 / 66.2 / 26.9 s (wp4 diffusivity) and 29.9 / 45.5 / 19.0 s (parallel-plate); chamber-driven 16.9 / 19.0 s, TiCl4 not reached in 20 s; M x2 -0.8 %, h 8 nm +2.4 % | figures/device_geometries.csv | same, --part nand13,refs |
| 2XX NAND wafer wall-area ratio and rule | 448 wafer areas (2.25 um slit pitch); rule 113 / 78.5 / 35.1 s (TMA / HfO2 / TiN) at 1 %, 3.4 / 1.6 / 1.2 x the feature's constant-0.3 Pa need; effusion limit 4.1 / 6.3 / 2.1 s | figures/device_geometries.csv | same |
| Device-geometry literature | no published TEM film thickness per sheet on a GAA stack or per word line on a 2XX stack; dimensions from Loubet 2017, Wong and Kakushima 2022, Agrawal 2024, TechInsights 2022/2023, Onto 2021/2023, imec 2021-2025 | paper/wp/wp10_device_geometries.md (Sources) | literature |
```

## Run log for paper/assumptions.md

```
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
