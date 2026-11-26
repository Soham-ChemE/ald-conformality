# Second blind test: entry log (custodian, 2026-11-26)

Entry rule of `../PROTOCOL_v2.md`, applied by the custodian before release. Criteria: 1 not seen before; 2 flat plateau;
3 at least six measured points between 0.25 and 0.75; 4 half depth reached and between 60 and 1000 gap heights;
5 gap, cycles (or thickness) and spot or step known; 6 member of a series of at least two entering curves.
A curve can fail more than one criterion. Curves sealed in `~/Documents/MATLAB/ald-blind-custody-v2/`.

Totals: 22 sources examined (11 more sought but not obtained), 19 curves digitised, 11 entering, in 1 series (1 calibration + 10 to predict).
By the protocol's own threshold (six entering curves in two series) the test is NOT informative as it stands: only one series.

## Sources with digitised curves

| source | figures digitised | digitised | entered | failed, number per criterion |
|---|---|---|---|---|
| Yim et al., PCCP 22, 23107 (2020) and ESI | Figs. 7, 8(b), 9(a), 9(c), 9(d), 11(a), 11(b); ESI Fig. S4 (one curve per Table 1 sample; sample 8 excluded as already seen) | 15 | 10 | criterion 2: 4; criterion 3: 2 |
| Haimi et al., Appl. Surf. Sci. Adv. 5, 100102 (2021) | Fig. 8 (two traces of specimen A) | 2 | 1 | criterion 5: 1 |
| Utriainen et al., Proc. SPIE 12008, 120080D (2022) | Fig. 4 (two traces, second reactor) | 2 | 0 | criterion 2: 1; criterion 5: 2; criterion 6: 2 |

## Sources examined, nothing digitised (reason is a property of the source, not of a curve shape)

| source | reason |
|---|---|
| Yim 2020 ESI Figs. S7, S15, S16 | repeats of curves already digitised from the main figures, or sample 8 (seen) |
| Yim 2020 ESI Fig. S8 (SEM-EDS line scans) | second measurement of samples already entered; channel entrance not marked on the axis, so no half depth |
| Yim 2020 ESI Fig. S17 (sample V0005) | nine overlapping line traces that cannot be separated in a raster figure |
| J. Yim, doctoral thesis, Aalto (2024) | Figs. 17, 18 reproduce Yim 2020 Fig. 9 |
| K. Arts, doctoral thesis, TU/e (2021) | thermal profiles are those of JVST A 37, 030908 (seen); the rest is plasma ALD |
| Schwille et al., JVST A 35, 01B118 (2017) | radial cavity fed by a central hole (criterion 1), cycle number and film thickness not stated (criterion 5) |
| Schwille et al., JVST A 35, 01B119 (2017) | same structure; one curve per temperature (criteria 1, 5, 6) |
| Schwille et al., Modelling Simul. Mater. Sci. Eng. 25, 035008 (2017) | republished SiO2 cavity data, same structure |
| M. Mattinen, MSc thesis, Helsinki (2016) (data of Langmuir 32, 10559) | half depth below 60 gap heights on the published axis (criterion 4) |
| J. Järvilehto, MSc thesis, Aalto (2023) | XPS with a spot of several hundred micrometres; the thesis states that half-thickness depths could not be extracted |
| Kia et al., Nanomaterials 9, 1035 (2019) | uncalibrated ToF-SIMS intensity, one curve (criteria 5, 6) |
| Emara et al., Appl. Surf. Sci. (2025) | XPS at 10 um spacing, one curve per process (criteria 3, 6) |
| Madadi et al., ACS Appl. Electron. Mater. (2024) | LiPON, one curve per process, conditions only in an SI not obtained (criteria 5, 6) |
| Philip et al., Nanoscale (2026); Philip et al., Chem. Commun. (2025) | ALD/MLD hybrid films, one curve per precursor (criterion 6) |
| Klabbers, MSc thesis, TU/e (2026) | re-plots the Haimi 2021 profile; no new measurement |
| Choolakkal et al., J. Phys. Chem. Lett. (2025); Heikkinen, MSc thesis, Aalto (2023) | images only, no thickness profile |
| Yim 2022, Järvilehto 2023, Gonsalves 2024 (PCCP, local, with ESIs) | simulated profiles only |
| Zenodo search (PillarHall, LHAR, saturation profile) | only record 3974438 (seen) holds measured profiles; the others are simulated |

## Not obtained (no open copy)

Puurunen and Gao, BALD 2016 (IEEE); Mattinen et al., Langmuir 32, 10559 (2016) (thesis used instead); Verkama MSc thesis;
Gonsalves thesis; Korpelainen et al. 2019; Manninen MSc thesis (Aalto, 2022, file not public); Utriainen et al. 2022 workshop paper
(captcha); Fan 2025, Deng 2025, Al Chimali 2026, Kim 2026 (closed).

## What is missing

A second series: a second process, temperature or reactor with at least two entering curves. The candidates found are single
curves or do not meet criteria 2 to 5.
