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
