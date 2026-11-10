# Work package 4: real device geometry and other chemistries

Tapered channel holes, word-line recesses filled through a slit, and a second chemistry (HfO2 from TDMAH) in the
chamber, plus rule-only predictions for TiN and W. Everything below comes from the runs listed in the run log or from
the cited sources; nothing is estimated from memory.

## What was done

1. **Tapered holes.** `aldconf.geometry.TaperedHole(diameter, length, bottom_ratio)` (a `Hole` subclass; `bottom_ratio`
   = 1 is the straight hole and takes the unchanged code path, so every earlier result is reproduced bit for bit; test
   `test_straight_tapered_hole_is_the_hole`). The diameter falls linearly from the mouth to `bottom_ratio x diameter`
   at the bottom. For a tapered geometry the solver uses the conservative form of the 1D diffusion equation,
   (1/A) d/dx (A D dn/dx) = (P/A) (vbar/4) [s (1 - theta) + r] n, with the frustum cross-section A(x) on the cell faces
   (`_solve_gas` gained optional `A_face`/`A_cell` arguments; with a zero wall sink the density in a frustum comes out
   uniform to 1e-9, `test_area_weighted_solve_conserves_flux`). Film narrowing keeps the previous uniform-area form so
   that the PillarHall fit and the chamber runs are unchanged. `wall_area()` returns the frustum lateral area
   pi (r_top + r_bot) sqrt(L^2 + (r_top - r_bot)^2), checked against a numerical integral to 1e-4.
   `coupling/tapered_hole.py` drives AR-100 holes (100 nm mouth, 10 um deep) with bottom diameters 100, 90, 80 and
   60 nm by the chamber's centre-face precursor density from `reactor/showerhead_pulse_patterned_AR100_long` (1 % TMA,
   300 Pa, 573 K, 300 sccm, 20 s pulse), exactly as `docs/build_assets.py` drives the straight hole, and also at a
   constant 0.3 Pa mouth density (feature-only benchmark, as in the Machball check).
2. **Word-line recesses.** New module `src/aldconf/recess.py`: `RecessedSlit(width, depth, recess_gap, recess_length,
   recess_pitch)`, `AreaEquivalentSlit`, `simulate_recessed` (time-varying mouth density), `simulate_recessed_constant`,
   `time_to_coverage`. A vertical slit (gap w_s, depth L_s) whose two walls carry lateral recesses (gap g_r, length
   L_r, vertical pitch p_r); a fraction f_r = g_r/p_r of each wall is recess mouth. Each recess is a 1D slit
   sub-problem (Knudsen diffusivity 2 g_r vbar/3, Langmuir sink on both faces, closed end) fed by the local slit gas
   density; the slit wall's uptake per unit wall area per unit gas density is K_wall(x) = (1 - f_r)(vbar/4)[s(1 -
   theta_s) + r] + f_r j_r(x), with j_r the recess mouth flux per unit mouth density from the quasi-steady sub-problem:
   the patternedSurface construction one level down. All recesses are solved per substep with a vectorised Thomas
   algorithm (200 slit cells x 50 recess cells). Limits verified by tests: no recesses reproduces the plain slit to
   1e-9 (`test_no_recess_is_the_plain_slit`); recesses of aspect ratio 1 (30 nm gap, 30 nm long, 60 nm pitch, area
   factor 1.5) reproduce the area-equivalent slit (wall sink x 1.5) within 3 % on the unsaturated front
   (`test_short_recesses_are_an_area_equivalent_slit`); recess bottoms lag the slit wall and order with depth
   (`test_recess_bottoms_lag_the_slit_wall_and_order_with_depth`). Dimensions from the 3D-NAND literature (Sources):
   word-line height 29 nm and inter-word-line oxide 22 nm (51 nm vertical pitch), 5 um stack (98 pairs), slit width
   250 nm and 100 nm, lateral recess length 365 nm and 1 um. `coupling/wordline_recess.py` runs the four geometries
   under the chamber's centre-face density (one-way, TMA, same drive as the holes) and at a constant 0.3 Pa.
3. **TDMAH chemistry.** `reactor/chemistry/make_tdmah_case.py` builds, runs and post-processes
   `reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah`: a copy of the long AR-100 case (1 %, 300 Pa, 573 K,
   300 sccm, 20 s pulse, 22 s run, frozen carrier flow) with vbar = sqrt(8RT/(pi M)) = 185 m/s for M = 354.7 g/mol,
   s = 7.0e-3, N0 = 2.77e18 Hf/m2, D = 3.62e-3 m2/s (hard-sphere Chapman-Enskog estimate for TDMAH-N2 relative to
   TMA-N2, applied to the chamber's TMA value 8.4e-3: ratio of sqrt(1/mu) / sigma_12^2 with the TDMAH diameter 1008 pm
   from its liquid density by the rule Ylilammi et al. 2018 use for TMA (591 pm) and the N2 Lennard-Jones diameter
   3.798 A; the absolute estimates are 8.14e-3 for TMA and 3.51e-3 for TDMAH) and the Danckwerts inlet fraction
   f = U/(U + D/delta) = 3.445e-3. A sensitivity variant `..._tdmah_Dtma` keeps the TMA diffusivity and inlet
   fraction with the TDMAH kinetics. The same script writes the rule-only table for TiN (TiCl4) and W (WF6)
   (`figures/chemistry_rule_table.csv`).
4. Tests: `tests/test_tapered_hole.py` (5) and `tests/test_recess.py` (4); 86 pass in 14 s.

## Results

### Table 1: tapered AR-100 holes (100 nm mouth, 10 um deep), 1 % TMA, baseline chamber (figures/tapered_hole.csv)

| bottom diameter | wall area / straight | chamber-driven bottom t90 (one-way) | rule with frustum area, 1.2 (1 + 55.9 A_rel) N0/(U c0) | driven / rule | constant 0.3 Pa mouth, t90 |
|---|---|---|---|---|---|
| 100 nm (straight) | 1.000 | 14.03 s (chamber two-way: 14.05 s first write, 14.01 s interpolated) | 14.39 s | 0.97 | 15.93 s |
| 90 nm (literature-anchored, cryogenic etch) | 0.950 | 14.21 s | 13.68 s | 1.04 | 16.52 s |
| 80 nm | 0.900 | 14.36 s | 12.98 s | 1.11 | 17.23 s |
| 60 nm | 0.800 | 14.62 s | 11.56 s | 1.26 | 19.12 s |

Reading. Under the chamber drive the taper costs 1 to 4 % in pulse time (14.03 -> 14.21 / 14.36 / 14.62 s): the
wafer is supply-limited and the hole's own transport is not the bottleneck, so the smaller bottom conductance (D_K ~ d,
A ~ d^2) is almost invisible. At a constant mouth density the same taper costs 4 to 20 % (15.9 -> 19.1 s), which is the
feature-only penalty. The rule with the frustum wall area moves the other way (less area, shorter pulse: 14.4 -> 11.6 s)
and under-predicts the driven time by up to 26 % at 60 % bottom diameter. The one-way drive keeps the straight wafer's
loading (wall-area ratio 55.9); on a wafer of tapered holes the loading itself falls by the area ratio, which the rule
captures and the one-way drive cannot, so the two numbers bracket the tapered-wafer result: the area reduction buys at
most the rule's 5 to 20 % and the transport penalty takes back 1 to 4 % of it. For the literature-anchored taper
(bottom 90 % of the mouth) both effects are within 5 % of the straight hole. Literature anchor: Lam Research reports
less than 0.1 % top-to-bottom CD deviation (defined as (CD_max - CD_min)/depth) in 10 um deep 3D-NAND channel holes with
its cryogenic etch, i.e. about 10 nm of taper on a 100 nm hole, with the EE Times Japan write-up quoting 108 nm top
and 99 nm bottom CD; Lee and Shin 2022 use a 0.3 degree taper (56.5 -> 51.5 nm channel radius over 1 um) in TCAD and
state that the bottom channel radius in mass production is 40-50 nm; room-temperature etches were reported with about
30 nm of CD loss (secondary source). No peer-reviewed vendor paper with a measured top/bottom CD pair in nm was found
in open access, so 90 % is the anchored value and 80 / 60 % remain a parametric sweep (60 % corresponds to the
simulated bottom/top ratio of Huang et al. 2019, 84/156 nm at AR 40, a reactive-ion-etch simulation, not a measurement).

### Table 2: word-line recesses through a slit, TMA (s 5e-3, N0 4e18), 5 um stack, 29 nm recess, 51 nm pitch (figures/wordline_recess.csv)

| slit width w_s | recess length L_r (AR_r) | area factor A_f | chamber-driven t90: plain slit bottom | slit bottom (flat wall) | recess bottom at mouth / mid / stack bottom | delay of the last recess bottom behind the slit bottom | area-equivalent slit bottom | constant 0.3 Pa: slit bottom / last recess bottom; exposure | ratio to the effusion-limit exposure |
|---|---|---|---|---|---|---|---|---|---|
| 250 nm | 365 nm (12.6) | 14.7 | 2.30 s | 5.37 s | 2.09 / 4.46 / 5.46 s | 0.09 s (1.8 %) | 5.39 s | 2.70 / 2.79 s; 0.84 Pa s | 4.6 |
| 250 nm | 1000 nm (34.5) | 39.6 | 2.30 s | 8.33 s | 2.94 / 6.94 / 9.05 s | 0.72 s (8.6 %) | 8.70 s | 6.02 / 7.03 s; 2.11 Pa s | 4.3 |
| 100 nm | 365 nm (12.6) | 14.7 | 3.76 s | 13.46 s | 2.11 / 8.55 / 13.53 s | 0.07 s (0.5 %) | 13.49 s | 14.53 / 14.68 s; 4.40 Pa s | 9.7 |
| 100 nm | 1000 nm (34.5) | 39.6 | 3.76 s | 16.52 s | 2.97 / 13.82 / 16.73 s | 0.20 s (1.2 %) | 16.66 s | 37.2 / 39.2 s; 11.7 Pa s | 9.6 |

Rule with the total wall area, for a wafer tiled with these slits at pitch w_s + 2 L_r (recesses from neighbouring
slits meeting mid-block): wall-area ratio 150 / 176 / 178 / 189 (fully patterned) and 75 / 88 / 89 / 94 (50 % patterned
as the reference wafer); rule pulse 38.3 / 44.8 / 45.2 / 48.0 s and 19.3 / 22.5 / 22.7 / 24.1 s at 1 % TMA. For
comparison the AR-100 hole wafer has 55.9 and 14.4 s. The effusion-limit exposure is the stoichiometric minimum at the
mouth, (A_wall/A_mouth) N0 / (n0 vbar/4), i.e. molecules needed per unit mouth area divided by the impingement rate.

Reading. Under the chamber drive the recess bottoms lag the slit bottom by 0.5 to 9 %, and the whole structure
saturates within 4 % of a plain slit carrying the same wall area (area-equivalent slit): the rule's area-only
treatment of the recesses holds for recess aspect ratios up to 35 (1 um long, 29 nm high). The lag grows with L_r/g_r
and is largest for the wide slit, whose recesses see a higher gas density and fill their mouths before their ends. The
word-line stack costs 2.3 to 4.4 x the plain-slit time (A_f 14.7 to 39.6 at f_r 0.57), and a 100 nm slit with 1 um
recesses takes 16.5 s under the same drive that fills the AR-100 holes in 14.0 s. At a constant mouth density the
recess structure is transport-limited by only 4 to 10 x over the effusion limit, against 39 x for the AR-100 hole
(15.9 s vs 0.41 s), because the slit is wide and the recesses short; the recess bottom lag is 1 to 17 % there.
What the one-way drive does not include: a wafer of such slits carries 2.7 to 3.4 x the wall area of the hole wafer
and would depress the surface density further, so its pulse is set by the rule at 150 to 189 wall areas (38 to 48 s
fully patterned at 1 % TMA), not by the 5 to 17 s above; the nested model's result is that the transport inside the
recess structure adds at most 9 % to that. A chamber boundary condition for recessed slits (the recess model inside
the patternedSurface BC) is the step that would close this loop and is not part of this package.
Geometry sources: word-line 29 nm / oxide 22 nm (Lee and Shin 2022, TCAD of a production-like cell); 25 to 30 nm oxide
(Lee et al. 2021); 20/20 nm (Miyaji et al. 2014); imec's scaled devices at 25 nm z-pitch with 10 nm word lines
(Rachidi et al. 2022); 98 oxide/nitride pairs for a 5 um test stack at AR about 50 (Reiter et al. 2024, Lam data);
slit width about 250 nm, slit pitch above 2 um, channel holes about 100 nm in a 2XX-layer device (Onto Innovation,
Semiconductor Digest 2023, no DOI); 730 nm slit-to-slit space in Samsung's 32-layer V-NAND (TechInsights/3DInCites
2014, no DOI), hence L_r = 365 nm; the lateral channel of the PillarHall structure (Gao et al. 2015; Ylilammi et al.
2018) is the measured analogue of a single recess.

### Table 3: second chemistry in the chamber and rule-only chemistries (reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah/results.json, figures/chemistry_rule_table.csv)

Chamber runs, AR-100 wafer (wall-area ratio 55.9), 1 % precursor, 300 Pa, 573 K, 300 sccm, 20 s pulse:

| | TMA (reference) | TDMAH | TDMAH with the TMA diffusivity (sensitivity) |
|---|---|---|---|
| M, vbar at 573 K | 72.09, 410 m/s | 354.7, 185 m/s | 354.7, 185 m/s |
| s, N0 | 5e-3, 4.0e18 | 7.0e-3, 2.77e18 | 7.0e-3, 2.77e18 |
| D in N2 at 300 Pa; inlet fraction f | 8.4e-3 m2/s; 1.49e-3 | 3.62e-3 m2/s; 3.45e-3 | 8.4e-3; 1.49e-3 |
| gap Peclet number U g / D | 0.09 | 0.21 | 0.09 |
| flat t90 (centre) | 1.95 s | 1.30 s | 1.20 s |
| hole-bottom t90 centre / worst / r = 149 mm | 14.05 / 14.10 / 13.50 s | 11.60 / 11.60 / 11.10 s | 11.10 / 11.15 / 10.95 s |
| rule 1.2 (1 + 55.9) N0/(U c0) | 14.38 s | 9.96 s | 9.96 s |
| bottom t90 / rule | 0.98 | 1.17 | 1.12 |
| surface density at 5 s (centre) | 0.084 c0 | 0.197 c0 | 0.213 c0 |
| wall time | (earlier run) | 132 s | 185 s |

Reading. TDMAH fills the AR-100 wafer in 11.6 s against 14.05 s for TMA at the same 1 % feed, because one HfO2
cycle needs 2.77e18 Hf/m2 against 4e18 Al-site equivalents for TMA; the rule tracks this (9.96 s) but the run sits
17 % above it where TMA sits 2 % below. 0.5 s of the 1.6 s excess is the lower gas diffusivity (gap Peclet 0.21
instead of 0.09: less of the injected precursor reaches the wafer before the outlet; the sensitivity run with the TMA
diffusivity gives 11.10 s, ratio 1.12); the rest is the shorter supply time (0.146 s per monolayer-equivalent against
0.211 s) against the fixed chamber overheads (0.53 s residence time, the Langmuir tail), the same trend the DOE showed
for the rule constant. Both chemistries are inside the 1.05 to 1.27 band of the twelve TMA cases once the
diffusivity is matched. Temperature: 573 K (300 C) is at the upper edge of the TDMAH/H2O self-limiting window.
Blaschke et al. 2020 find the GPC minimum at 300 C (0.88 A/cycle) and attribute the rise at 350 C to thermal
decomposition of the precursor; Gieraltowska et al. 2023 report ligand decomposition raising the GPC above 240 C;
Hausmann et al. 2002 deposited from 50 to 500 C with 100 % step coverage above AR 35. The chamber temperature was
kept at 573 K so that the carrier flow field (0/U, 0/phi), U, c0 and the Danckwerts inlet are those of every other
run in the paper and the comparison isolates the chemistry. For the record, the rule at 523 K (250 C, inside the
window) is 9.26 s with the 250 C GPC of 0.93 A (N0 2.58e18) and 9.96 s with 1.0 A, because U c0 at fixed sccm and
fraction does not depend on temperature (U ~ T, c0 ~ 1/T; night 2 part 4).

Rule-only predictions, AR-100 wafer (55.9), 1 % precursor, baseline chamber (U c0 = 1.90e19 molecules/m2/s):

| film / precursor | T | s (source) | GPC (source) | film density (source) | N0 = GPC rho N_A / M | vbar | s vbar/4 vs U = 0.05 m/s | supply time N0/(U c0) | rule, AR-100 wafer | rule, blanket (x 1.2) |
|---|---|---|---|---|---|---|---|---|---|---|
| Al2O3 / TMA | 300 C | 5e-3 (chamber runs; Arts 2019 0.5-2e-3, PillarHall fit 8.5e-3) | 0.10 nm | (N0 taken as 4e18, one saturated monolayer) | 4.0e18 | 410 m/s | 0.51 | 0.211 s | 14.4 s | 0.25 s |
| HfO2 / TDMAH | 300 C | 7.0e-3 (Reiter 2024, Hf amide, calibrated to 3D-NAND HfO2 step coverage) | 0.10 nm (Blaschke 2020: 1.08 A at 200 C, 0.93 A at 250 C; Gieraltowska 2023: 1.2 A at 180-240 C) | 9.68 g/cm3 bulk monoclinic (Delabie 2005) | 2.77e18 | 185 m/s | 0.32 | 0.146 s | 9.96 s | 0.17 s |
| HfO2 / TDMAH | 250 C | 7.0e-3 | 0.093 nm (Blaschke 2020) | 9.68 g/cm3 | 2.58e18 | 177 m/s | 0.31 | 0.136 s | 9.26 s | 0.16 s |
| TiN / TiCl4 | 325 C | 0.10 (Ylilammi 2018, TiCl4/H2O on TiO2 at 110 C; no TiCl4/NH3 value published) | 0.030 nm (Kinnunen 2017, 325 C; Satta 2002: 0.024 nm at 400 C) | 4.23 g/cm3 XRR (Kinnunen 2017; bulk 5.21) | 1.24e18 | 253 m/s | 6.3 | 0.065 s | 4.44 s | 0.08 s |
| W / WF6 | 300 C | 0.40 (Elam 2001, WF6 on the Si-terminated surface at 573 K; no B-terminated value) | 0.25 nm (Klaus 2000, WF6/Si2H6) | 19.25 g/cm3 (bcc, a = 316.5 pm) | 1.58e19 | 202 m/s | 20 | 0.830 s | 56.7 s | 1.00 s |

Reading. The rule is chemistry-agnostic: s does not enter it at all once the wafer is supply-limited (s vbar/4 is
6 to 400 x the inlet speed for all four precursors, against 0.4 x for H2O), and the chemistry enters only through N0,
the metal atoms one cycle must deliver per unit wall area, and c0. TiN needs the shortest pulse per cycle (4.4 s,
0.3 A/cycle of a light film) and W the longest (57 s, 2.5 A/cycle of a dense metal at 1 % WF6). Per nanometre of film
the pulse is proportional to the film's metal-atom number density over c0: HfO2 100 s/nm (2.8e28 Hf/m3), Al2O3 144
(4.0e28), TiN 148 (4.1e28), W 227 s/nm (6.3e28), which is the throughput number a tool sees. The N0 uncertainty is the GPC x density uncertainty of each chemistry (factor 1.1 to 1.3 across the cited
values), the same dominant term the Monte Carlo on the TMA rule found. TDMAH's sticking value is the weakest input:
no direct measurement exists; the 7.0e-3 is a model calibration to Lam 3D-NAND HfO2 step-coverage data for an
undisclosed Hf amide precursor, bracketed by the TMA-like range 1e-3 to 1e-2; it does not move the patterned-wafer
result (supply-limited) but would move a blanket-wafer one.

## Figures

- `figures/tapered_hole.png`: (left) bottom coverage vs time of AR-100 holes with 100 / 90 / 80 / 60 nm bottoms
  driven by the chamber's centre-face TMA density (1 %, 20 s pulse), with the chamber's own two-way 14.01 s marked;
  (middle) coverage along the hole at 4, 8 and 12 s, straight (solid) vs 60 % bottom (dashed); (right) time to 90 %
  bottom coverage vs bottom/mouth diameter: chamber-driven, the rule with the frustum wall area, and the constant
  0.3 Pa feature-only benchmark. Caption: "Tapered channel holes. Under the supply-limited chamber the taper costs
  1-4 %; at constant mouth density 4-20 %; the rule with the reduced frustum area predicts 5-20 % less."
- `figures/wordline_recess.png`: base geometry (250 nm slit, 5 um deep, 29 nm recesses 1 um long at 51 nm pitch,
  TMA), chamber-driven. Top row: coverage in the recesses (x = depth along the slit, y = distance into the recess) at
  3, 6 and 9 s. Bottom left: coverage of the slit wall (solid) and of the recess bottoms (dashed) along the slit at
  the same times. Bottom middle: slit mouth, slit bottom and recess bottoms at the mouth and at the stack bottom vs
  time under the chamber drive. Bottom right: the same at a constant 0.3 Pa mouth density. Caption: "Word-line
  recesses fed through a slit. The saturation front runs down the slit and into each recess behind it; the last
  recess bottom lags the slit bottom by 0.7 s (9 %) under the chamber drive, and the structure saturates within 4 %
  of a slit with the same total wall area."

## Proposed paper text

### 3.11 Tapered holes, word-line recesses and other chemistries

Real channel holes are not cylinders. The feature solver now takes a frustum (TaperedHole, diameter falling linearly
from the mouth to a fraction of it at the bottom) and solves the conservative form (1/A) d/dx (A D dn/dx) of the
diffusion equation on the varying cross-section; with the sink switched off the density in a frustum is uniform to
1e-9, and bottom_ratio = 1 reproduces the cylinder exactly. Published numbers for the taper of production holes are
scarce: Lam reports less than 0.1 % CD deviation over 10 um with cryogenic etching, about 10 nm on a 100 nm hole,
against about 30 nm for room-temperature etches [Lam 2024; EE Times Japan 2024], and TCAD studies use 0.3 degree
tapers [Lee and Shin 2022]; we therefore take 90 % as the anchored bottom diameter and 80 and 60 % as a sweep.
Driven by the chamber's centre-face TMA density of the AR-100 wafer, the bottoms of 100 / 90 / 80 / 60 nm-bottom holes
reach 90 % at 14.03 / 14.21 / 14.36 / 14.62 s (Figure, tapered_hole.png): the taper costs 1 to 4 %, because on a
supply-limited wafer the hole's own conductance is not the bottleneck. At a constant mouth density the same taper
costs 4 to 20 % (15.9 to 19.1 s): the lost conductance (D_K ~ d, A ~ d^2) outweighs the lost wall area. The rule with
the frustum wall area (0.95 / 0.90 / 0.80 of the cylinder's) moves the other way, to 13.7 / 13.0 / 11.6 s, and the
two bracket a tapered-hole wafer: the area reduction buys at most the rule's 5 to 20 %, the transport penalty takes
back 1 to 4 % of it, and for the anchored 90 % taper both are within 5 % of the straight hole.

Word lines of a replacement-gate 3D-NAND stack are filled through the gate-line slit into lateral cavities a few tens
of nanometres high and hundreds of nanometres to a micrometre deep. We model this as a nested feature (aldconf.recess):
a slit of width w_s and depth L_s whose walls carry recesses of height g_r, length L_r and pitch p_r, each recess a
Knudsen-Langmuir slit sub-problem fed by the local slit density, with the slit wall's uptake equal to the flat uptake
on the un-recessed fraction plus the recess mouth flux on the recessed fraction, solved quasi-steadily; it is the
patternedSurface construction one level down, and it reduces to the plain slit without recesses (to 1e-9) and to a
slit with the same total wall area when the recesses are short (within 3 % at aspect ratio 1). With literature
dimensions (29 nm word lines at 51 nm pitch [Lee and Shin 2022], 98 pairs in 5 um [Reiter et al. 2024], 250 nm slits
with recesses of 365 nm [TechInsights 2014] or 1 um [Onto Innovation 2023]) and the TMA drive of the AR-100 wafer,
the slit bottom saturates at 5.4 and 8.3 s and the last recess bottom 0.09 and 0.72 s later (2 and 9 %); a 100 nm
slit takes 13.5 and 16.5 s with lags of 0.5 and 1 % (Figure, wordline_recess.png). In every case the structure
saturates within 4 % of the plain slit carrying the same total wall area (14.7 and 39.6 times the slit wall), so the
rule's area-only treatment of recesses holds up to recess aspect ratios of 35. The consequence for the recipe is in
the area: a wafer tiled with these slits at a 1 to 2.3 um pitch carries 150 to 189 wafer-areas of wall, 2.7 to 3.4 x
the AR-100 hole wafer, and the rule puts its TMA pulse at 38 to 48 s (fully patterned) at 1 %. Closing this with a
recessed-slit boundary condition in the chamber is the natural next step; the nested model says the correction it
would find inside the structure is below 10 %.

The rule carries no chemistry beyond N0 and c0, and a second chamber chemistry tests that. HfO2 from TDMAH
(M 354.7, vbar 185 m/s at 573 K) was run on the AR-100 wafer at 1 %, 300 Pa, 573 K, 300 sccm with s = 7.0e-3 [Reiter
et al. 2024, calibrated to 3D-NAND HfO2 step coverage; no direct TDMAH measurement exists], N0 = 2.77e18 Hf/m2 from a
GPC of 0.10 nm [Blaschke et al. 2020; Gieraltowska et al. 2023] and the bulk density 9.68 g/cm3 [Delabie et al.
2005], and a gas diffusivity of 3.6e-3 m2/s (hard-sphere Chapman-Enskog relative to TMA). The hole bottoms reach
90 % at 11.6 s against 14.05 s for TMA; the rule gives 9.96 s against 14.38 s. TDMAH sits 17 % above its rule where
TMA sits 2 % below: 0.5 s of the 1.6 s excess is the 2.3 x lower diffusivity (gap Peclet number 0.21 instead of 0.09; a
run with the TMA diffusivity gives 11.1 s), the rest the shorter supply time against the fixed chamber overheads,
the trend the DOE found for the rule constant. 573 K is at the upper edge of the TDMAH self-limiting window [Blaschke
et al. 2020; Gieraltowska et al. 2023; Hausmann et al. 2002]; it was kept so that the flow field, U, c0 and the inlet
condition are those of every other run, and the rule at 250 C is 9.3 to 10.0 s because U c0 at fixed sccm and
fraction is temperature-independent. Rule-only predictions for TiN from TiCl4 (s 0.10 [Ylilammi et al. 2018], GPC
0.30 A, density 4.23 g/cm3 [Kinnunen et al. 2017], N0 1.24e18) and W from WF6 (s 0.40 [Elam et al. 2001], GPC 2.5 A
[Klaus et al. 2000], 19.25 g/cm3, N0 1.58e19) give 4.4 and 57 s per cycle on the same wafer at 1 %: the sticking
coefficients span 5e-3 to 0.4 and do not enter, since s vbar/4 exceeds the inlet speed by 6 to 400 x for all four
precursors; what differs is the number of metal atoms a cycle has to deliver. Per nanometre of film the pulse
scales with the film's metal-atom number density (HfO2 100 s/nm, Al2O3 144, TiN 148, W 227), which is the number a
tool's throughput model needs.

## CLAIMS rows

| Claim | Value | Source in repo | How to reproduce |
|---|---|---|---|
| Tapered AR-100 holes (bottom 90 / 80 / 60 % of the 100 nm mouth), chamber-driven at 1 % TMA | bottom t90 14.21 / 14.36 / 14.62 s vs 14.03 s straight (one-way; chamber two-way 14.05 s); constant 0.3 Pa 16.5 / 17.2 / 19.1 vs 15.9 s; rule with frustum area 13.7 / 13.0 / 11.6 s | figures/tapered_hole.csv, figures/tapered_hole.png | PYTHONPATH=src .venv/bin/python coupling/tapered_hole.py |
| Frustum solver | zero-sink density uniform to 1e-9; bottom_ratio 1 identical to Hole; wall area matches the integral to 1e-4 | tests/test_tapered_hole.py | pytest |
| Word-line recess model limits | no recesses = plain slit to 1e-9; AR-1 recesses = area-equivalent slit within 3 % | tests/test_recess.py | pytest |
| Word-line recesses (29 nm, 51 nm pitch, 5 um stack) under the AR-100 wafer's TMA drive | slit bottom t90 5.37 / 8.33 s (250 nm slit, L_r 365 / 1000 nm), 13.46 / 16.52 s (100 nm slit); last recess bottom 0.09 / 0.72 / 0.07 / 0.20 s later (0.5-9 %); area-equivalent slit within 4 % | figures/wordline_recess.csv, figures/wordline_recess.png | PYTHONPATH=src .venv/bin/python coupling/wordline_recess.py |
| Rule for a word-line-slit wafer | wall-area ratio 150-189 (fully patterned), rule 38-48 s at 1 % TMA | figures/wordline_recess.csv | same |
| HfO2 from TDMAH on the AR-100 wafer, 1 %, 300 Pa, 573 K | bottom t90 11.60 s centre / 11.10 s edge, flat 1.30 s; rule 9.96 s (ratio 1.17); with the TMA diffusivity 11.10 s (1.12) | reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah/results.json | .venv/bin/python reactor/chemistry/make_tdmah_case.py --run --d-sensitivity |
| Rule-only chemistries on the AR-100 wafer at 1 % | TiN/TiCl4 4.44 s, W/WF6 56.7 s, HfO2/TDMAH 9.96 s (9.26 s at 250 C), Al2O3/TMA 14.4 s | figures/chemistry_rule_table.csv | .venv/bin/python reactor/chemistry/make_tdmah_case.py --post |

## Run log: Night 3, part 5: geometry and chemistry (reactor/chemistry/make_tdmah_case.py, coupling/tapered_hole.py, coupling/wordline_recess.py, src/aldconf/recess.py)

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

## Sources

Chemistry
- Reiter, Aguinsky, Rodrigues, Weinbub, Hoessinger, Filipovic, Solid-State Electronics 211, 108816 (2024),
  doi 10.1016/j.sse.2023.108816. Table 3: calibrated initial sticking probability beta0 = 7.0e-3 for the hafnium
  step of HfO2 ALD on Lam's 3D-NAND test structure (76/98 ON layers, 4-5 um, AR about 50, 85 % step coverage);
  precursor assumed TEMAH; hard-sphere diameters d_TMA = 591 pm used with Bosanquet interpolation.
- Blaschke, Munnik, Grenzer, Rebohle, Zahn, Schmidt, Gemming, Applied Surface Science 506, 144188 (2020),
  doi 10.1016/j.apsusc.2019.144188. TDMAH/H2O: 500 cycles give 70.7 / 62.0 / 53.8 / 46.7 / 44.2 / 46.1 nm at
  100-350 C in 50 C steps (1.41 to 0.88 A/cycle; 1.08 A at 200 C, 0.93 A at 250 C); minimum 0.88 A/cycle at 300 C;
  the rise at 350 C attributed to thermal decomposition of the precursor.
- Gieraltowska et al., Materials 16, 4077 (2023), doi 10.3390/ma16114077. TDMAH/H2O GPC 1.2 A/cycle at 180-240 C;
  ligand decomposition raises the GPC above 240 C.
- Hausmann, Kim, Becker, Gordon, Chemistry of Materials 14, 4350 (2002), doi 10.1021/cm020357x. HfO2 and ZrO2 from
  metal amides and water, 50-500 C; 100 % step coverage in holes of aspect ratio above 35.
- Delabie et al., Journal of Applied Physics 97, 064104 (2005), doi 10.1063/1.1856221. Bulk HfO2 density
  9.68 g/cm3; one HfO2 monolayer 9.15 Hf/nm2; HfCl4/H2O steady growth 1.1 Hf/nm2 per cycle.
- Wang et al., Chemistry of Materials 19, 3127 (2007), doi 10.1021/cm061761p. Bulk HfO2 9.69 g/cm3 = 8.316e22
  atoms/cm3; TEMAH/D2O 0.11 nm/cycle on SiO2.
- Ylilammi, Ylivaara, Puurunen, Journal of Applied Physics 123, 205301 (2018), doi 10.1063/1.5028178. Table I:
  lumped sticking coefficient 0.00572 for TMA/H2O on Al2O3 at 300 C and 0.10 for TiCl4/H2O on TiO2 at 110 C;
  hard-sphere diameters from condensed-phase density, d_TMA = 591 pm, d_TiCl4 = 703.9 pm.
- Arts, Vandalon, Puurunen, Utriainen, Gao, Kessels, Knoops, Journal of Vacuum Science and Technology A 37, 030908
  (2019), doi 10.1116/1.5093620. TMA (0.5-2)e-3, H2O (0.8-2)e-4 at a 275 C set-point.
- Kinnunen, Malm, Arstila, Lahtinen, Sajavaara, Nuclear Instruments and Methods B 406, 152 (2017),
  doi 10.1016/j.nimb.2016.12.032. TiCl4/NH3 TiN: 0.30 A/cycle at 325 C, 0.36 A/cycle at 450 C; XRR density
  4.23 g/cm3 (bulk 5.21 g/cm3).
- Satta et al., Journal of Applied Physics 92, 7641 (2002), doi 10.1063/1.1522485. TiCl4/NH3 TiN 24 +/- 3 A after
  about 100 cycles at 400 C (0.24 A/cycle).
- Elam, Nelson, Grubbs, George, Surface Science 479, 121 (2001), doi 10.1016/S0039-6028(01)00969-4. WF6 reactive
  sticking coefficient 0.4 at 573 K on the Si2H6-saturated surface; Si2H6 about 5e-2.
- Klaus, Ferro, George, Thin Solid Films 360, 145 (2000), doi 10.1016/S0040-6090(99)01076-7. W ALD from WF6/Si2H6,
  2.5 A per AB cycle, 425-600 K.
- Poling, Prausnitz, O'Connell, The Properties of Gases and Liquids, 5th ed., McGraw-Hill (2001), Appendix B
  (no DOI): Lennard-Jones diameter of N2 3.798 A.
- TDMAH liquid density 1.098 g/mL at 25 C: vendor property compilation, ChemicalBook entry for CAS 19782-68-4
  (no DOI; used only for the hard-sphere diameter, to which the chamber result is insensitive: the diffusivity
  sensitivity run changes the bottom t90 by 0.5 s for a 2.3 x change in D).
- Tungsten bcc lattice constant 316.5 pm (Kittel, Introduction to Solid State Physics, 8th ed., Wiley 2005, Table 1.4;
  no DOI): 2 atoms per cell gives 6.31e28 W/m3 = 19.25 g/cm3.

Geometry
- Lee, Shin, Micromachines 13, 1139 (2022), doi 10.3390/mi13071139. 3D-NAND TCAD: gate length 29.0 nm, space
  22.0 nm; taper angle 0.3 degrees, channel radius 56.5 nm top to 51.5 nm bottom over 1 um; bottom channel radius
  40-50 nm in mass production.
- Lee, Yoon, Sim, Song, Electronics 10, 2632 (2021), doi 10.3390/electronics10212632. Channel hole 120 nm at 150 nm
  string pitch; inter-word-line oxide 30, 25 and 19 nm.
- Miyaji et al., Japanese Journal of Applied Physics 53, 024201 (2014), doi 10.7567/JJAP.53.024201. 20 nm gate and
  20 nm space (40 nm pitch), 90 nm holes.
- Rachidi et al. (imec), IEEE International Memory Workshop 2022, doi 10.1109/IMW52921.2022.9779303. 25 nm z-pitch
  with 10 nm word lines.
- Kihara, Tomura, Sakamoto, Honda, Kojima (TEL), Symposium on VLSI Technology and Circuits 2023,
  doi 10.23919/VLSITechnologyandCir57934.2023.10185160. 10 um deep channel-hole etch for 400+ layers.
- Shen et al. (Lam), Japanese Journal of Applied Physics 62, SI0801 (2023), doi 10.35848/1347-4065/accbc7. About
  100 nm holes, aspect ratio above 100:1; bow-to-bottom CD difference as the profile metric.
- Oehrlein et al., Journal of Vacuum Science and Technology B 42, 041501 (2024), doi 10.1116/6.0003579. Aspect
  ratios approaching 100:1; 100 nm, 10 um holes through more than 400 layers.
- Huang, Huard, Shim, Nam, Song, Lu, Kushner, Journal of Vacuum Science and Technology A 37, 031304 (2019),
  doi 10.1116/1.5090606. Simulated profiles in 4.8 um SiO2: top / bow / bottom CD about 156 / 180 / 84 nm.
- Gao, Arpiainen, Puurunen, Journal of Vacuum Science and Technology A 33, 010601 (2015), doi 10.1116/1.4903941.
  PillarHall lateral channels, 200-1000 nm high, 20 um to 5 mm long (the measured analogue of a recess).
- Cremers, Puurunen, Dendooven, Applied Physics Reviews 6, 021302 (2019), doi 10.1063/1.5060967. Conformality review.
- Lam Research press release, 31 July 2024, "Lam Research's Cryo 3.0 ..." (no DOI): less than 0.1 % profile deviation
  ((CD_max - CD_min)/depth) at depths up to 10 um. EE Times Japan, 2 September 2024 (no DOI, secondary): top CD
  108 nm, bottom CD 99 nm with the cryogenic process; about 30 nm CD loss with the conventional process.
- Keller et al. (Onto Innovation), "Measuring the critical W-recess in 3D NAND", Semiconductor Digest, July 2023,
  p. 31 (no DOI): 2XX layer pairs, slit pitch above 2 um, slit width about 250 nm, W recess about 50 nm, channel
  holes about 100 nm.
- TechInsights / 3DInCites, Samsung 32-layer V-NAND teardown, December 2014 (no DOI): channel hole 120 nm before
  ONO; vertical tungsten slit-to-slit space 730 nm; hole pitch 245 nm between slits.
