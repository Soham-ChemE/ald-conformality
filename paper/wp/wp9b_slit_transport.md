# Work package 9b: the slit Knudsen diffusivity, corrected and refitted

## What was found
WP9 (f) showed the one-dimensional slit model 24 % slow against two ballistic solvers on a plain slit of aspect ratio 50
(t90 1.38 s against voxel 1.08-1.11 s and Machball Trench 1.12 s at 0.3 Pa). The cause is the slit's Knudsen
diffusivity: `transport.knudsen_diffusivity` uses Knudsen's long-channel form D_K = d_h vbar/3 with the hydraulic
diameter d_h = 4A/P, which is exact for a cylinder (transmission 4d/3L) and not for parallel plates. Between plates a
molecule that leaves one plate at a grazing angle travels w tan(theta) before the next collision, and the cosine law
gives <tan^2 theta> a logarithmic divergence: the mean-square axial step per wall collision, and with it the
diffusivity, grows as ln(L/w) until the channel length (or whatever else interrupts the flight) cuts it off. In the
random-walk form D = (1/2) nu <dx^2> with the wall-collision rate per molecule nu = (vbar/4)(P/A) = vbar/(2 w),
<dx^2> = w^2 [ln(L/w) + c] gives D_slit = (w vbar/4) [ln(L/w) + c]; the same result follows from Fick's law across
the channel, D = W L vbar/4, with the Clausing slit transmission W = (w/L) [ln(L/w) + c] (Berman, J. Appl. Phys. 36,
3356 (1965), tabulates W for slits; Steckelmacher, Rep. Prog. Phys. 49, 1083 (1986), reviews the Clausing and
Smoluchowski long-channel results; Cremers, Puurunen and Knoops, Appl. Phys. Rev. 6, 021302 (2019), discuss Knudsen
transport in lateral channels in the ALD context). Against the hydraulic-diameter value 2 w vbar/3 the ratio is

    D_slit / D_K = (3/8) [ln(L/w) + c].

The hydraulic-diameter value is therefore right only at one aspect ratio, L/w = e^(8/3) = 14.4, too high below it
and too low above it, by a factor 1.47 at AR 50 and 2.33 at AR 500.

### Transmission Monte Carlo (the constant c)
A test-particle Monte Carlo of the infinitely wide slit written for this check (cosine-law entry over the mouth,
diffuse re-emission at the plates, no gas-phase collisions; 4e5 particles per point) gives the transmission
probabilities below; W = 0.6844 at L/w = 1 is Clausing's value.

| L/w | 1 | 2 | 5 | 10 | 20 | 50 | 100 | 200 | 500 |
|---|---|---|---|---|---|---|---|---|---|
| W | 0.6844 | 0.5425 | 0.3572 | 0.2399 | 0.1545 | 0.0802 | 0.0472 | 0.0270 | 0.0124 |
| W L/w - ln(L/w) | 0.684 | 0.392 | 0.176 | 0.096 | 0.094 | 0.097 | 0.111 | 0.102 | -0.015 (noise +/- 0.08) |
| (3/8) W L/w | 0.257 | 0.407 | 0.670 | 0.899 | 1.159 | 1.503 | 1.768 | 2.025 | 2.325 |

So c = +0.10 for L/w = 10 to 200, and the leading term alone, (3/8) ln(L/w), is within 4 % of the full Monte Carlo
factor from L/w = 20 upwards.

### Machball calibration (task 1a-c)
Machball `Trench(AR, Nz = 8 AR)` (converged: AR 20 at Nz 80 / 160 / 320 gives 0.6242 / 0.6231 / 0.6226 s, AR 50 at
Nz 200 / 400 / 800 gives 1.1194 / 1.1170 / 1.1158 s) against the 1D slit model (`Slit(100 nm, AR x 100 nm)`,
200 cells, 200 substeps, no narrowing; 400 / 400 moves t90 at AR 50 by 0.2 %), s = 5e-3, N0 = 4e18 /m2, 573 K, 0.3 Pa
TMA. Both times are the time at which the bottom coverage itself reaches 0.9 (Machball's end criterion; the
repository's `dose_for_coverage` uses the step coverage, which is the same thing once the mouth is saturated, i.e.
for AR >= 50 here, and is why the earlier WP9 value 1.381 s agrees with the 1.384 s below). The matched factor is the
`knudsen_factor` at which the 1D t90 equals Machball's.

| AR | Machball t90 / t95 (s) | 1D t90 / t95, D_K = 2 w vbar/3 (s) | 1D / Machball (t90) | matched factor (t90 / t95) | (3/8) ln(AR) | 1D t90 with (3/8) ln(AR) | 1D t90 with (3/8)(ln AR + 0.10) |
|---|---|---|---|---|---|---|---|
| 10 | 0.5335 / 0.6882 | 0.5087 / 0.6533 | 0.954 | 0.59 / 0.52 (ill-conditioned, see below) | 0.863 | 0.5143 (-3.6 %) | 0.5126 (-3.9 %) |
| 20 | 0.6231 / 0.7782 | 0.6149 / 0.7647 | 0.987 | 0.945 / 0.917 | 1.123 | 0.5993 (-3.8 %) | 0.5952 (-4.5 %) |
| 50 | 1.1170 / 1.2929 | 1.3842 / 1.5585 | 1.239 | 1.405 / 1.387 | 1.467 | 1.0891 (-2.5 %) | 1.0735 (-3.9 %) |
| 100 | 2.6305 / 2.8373 | 4.2305 / 4.4503 | 1.608 | 1.721 / 1.715 | 1.727 | 2.6234 (-0.3 %) | 2.5767 (-2.0 %) |

Reading. At AR 10 and 20 the slit is reaction-limited: the kinetic floor ln(10) N0/(s n vbar/4) = 0.47 s is most of
the 0.53-0.62 s, the uncorrected model is already within 1-5 % of Machball, and the matched factor is a poorly
conditioned number (t90 moves 1 % for a 30 % change of D). At AR 50 and 100 the slit is transport-limited and the
matched factors, 1.405 and 1.721, follow the logarithm: (3/8) ln(AR) gives 1.467 and 1.727, 4 % and 0.3 % above
them. The 4 % at AR 50 is the continuum-versus-ballistic offset that the holes show with their exact diffusivity (1D
7 % fast at AR 50, 3 % at AR 100, WP9 (a)): the diffusion model with the right D is a few per cent fast against a
ballistic solver at moderate AR and converges to it with AR. For that reason the implemented factor drops the
Monte Carlo constant: with c = 0 the corrected 1D model is 0.3-3.8 % fast against Machball at AR 20-100 and within
5 % everywhere; with c = 0.10 it is 2.0-4.5 % fast. The factor, not the constant, is the point: the correction turns
a 24 % (AR 50) and 61 % (AR 100) error into 2.5 % and 0.3 %.

### Front shape, not only t90
The PillarHall fit is constrained by the shape of the saturation front, so the correction was also checked on the
coverage profile. Machball Trench against the 1D model at the same times, front position z50 (depth of half
coverage, in gaps) and front width (depth between coverage 0.8 and 0.2):
- AR 200, s = 5e-3, Nz 800, t = 0.78 / 1.99 s (front well away from the closed end): Machball z50 39.7 / 79.2, width
  48.2 / 53.4; 1D corrected z50 50.6 / 92.0, width 50.7 / 50.8; 1D uncorrected z50 35.9 / 65.2, width 36.0 / 36.0.
- AR 500, s = 8.5e-3, Nz 1000, 28 frames from z50 60 to 187: Machball width grows slowly from 39.9 to 48.4 (the
  cut-off length of the flights is the saturated length behind the front, which grows); 1D corrected width 42.1 at
  every depth (within 6 % of Machball at z50 = 60-100, 13 % narrower at z50 = 187), front 14-30 % ahead of
  Machball's; 1D uncorrected width 27.6 (31-43 % too narrow), front 15-25 % behind. A variant with the local
  cut-off ln(z50) (width 34-39) is worse than the full-length factor.
So the corrected model reproduces the ballistic front width to about 10 % and places the front 15-25 % too deep; the
uncorrected model has the front 30-40 % too narrow and 15-25 % too shallow. Since the front width fixes s (width ~
sqrt(D/s)) and the front position fixes the exposure (depth ~ sqrt(D x exposure)), a fit with the uncorrected model
returns an s that is low by about the factor D_slit/D_K and an exposure high by the same factor.

### Implementation (task 1d)
- `transport.slit_knudsen_factor(gap, length)` = (3/8) ln(L/w) for L > w, 1 otherwise, with the derivation and the
  references in its docstring; exported from `aldconf`.
- `geometry.Slit` takes a keyword-only field `knudsen_correction` ("none" default, "plates"), validated, and a method
  `knudsen_factor()`; `Slit(gap, length)` and `AreaEquivalentSlit(gap, length, area_factor)` are unchanged, so every
  earlier number is reproduced by default.
- `solver.simulate` and `simulate_driven` multiply `Process.knudsen_factor` by the geometry's factor (helper
  `_geometry_factor`, 1 for holes); `multispecies._diffusivities` does the same.
- `recess.simulate_recessed(..., knudsen_correction="none" | "plates")` (and `simulate_recessed_constant` through
  `**kw`) applies the factor to the slit (w_s, L_s) and to every recess (g_r, L_r) separately.
- tests/test_slit_transport.py (6 tests): the formula; the default reproduces the WP9 reference t90 1.3813 s at AR 50
  to 1 %; "plates" equals the explicit knudsen_factor to 1e-12; the recess model without recesses equals the corrected
  plain slit to 1e-9; the corrected slit matches Machball Trench within 5 % at AR 20 and 50 (and the uncorrected one is
  more than 20 % slow at AR 50). The two Machball tests skip when machball is absent, like test_machball_crosscheck.py.
  `PYTHONPATH=src .venv/bin/python -m pytest -q tests`: 124 passed in 32 s (118 before).

## PillarHall refit (task 2; coupling/fit_pillarhall_corrected.py; figures/pillarhall_fit_corrected.png, .csv)
Same data (Zenodo 3974438, 8 long-channel sheets averaged on a 2 um grid; plateau 51.54 nm, gpc 103.1 pm/cycle,
half-thickness depth 107.4 um = AR 215), same model (Slit 500 nm x 250 um, 500 cycles with channel narrowing, N0 =
4e18 fixed, 400 cells, 60 substeps per pulse), same procedure (nine Nelder-Mead starts on the mean profile, then one
refit per sheet from the joint optimum for the sheet-to-sheet spread). The only change is
`knudsen_correction="plates"`: factor (3/8) ln(500) = 2.330 on D_K.

| | uncorrected (coupling/fit_pillarhall.py, figures/pillarhall_fit.png) | corrected slit transport |
|---|---|---|
| s_TMA | 8.5e-3 | 1.96e-2 |
| TMA partial pressure at the mouth | 75 Pa | 32.3 Pa |
| exposure per 0.1 s pulse | 7.5 Pa s (56 000 L) | 3.23 Pa s (24 200 L) |
| rms of the fit (fraction of the plateau) | 0.038 | 0.0381 |
| sheet-to-sheet | s x/÷ 1.5 | s = 1.99e-2 x/÷ 1.46 (geometric mean of the eight sheets; range 7.9e-3 to 2.7e-2), p = 32.4 Pa x/÷ 1.02 |

All nine starts converge to the same optimum (s within 1 %, pressure within 0.1 %), and the rms is the same as the
uncorrected fit's: the data cannot tell the two transport models apart, because the fit has two free parameters and the
profile has two features (front position and front width). What changes is the meaning of the numbers. s rises by
the factor 2.3 by which D rose (front width ~ sqrt(D/s) at fixed measured width), and the pressure falls by the same
factor (penetration ~ sqrt(D x exposure) at fixed measured depth). The published fit of Ylilammi et al. 2018 to
their own 500 nm PillarHall profile returned 147 Pa; the uncorrected 75 Pa was a factor 2 below it and the corrected
32 Pa a factor 4.6.

### Against the literature (the WP1 reconciliation, rerun)
The 300 C values are 5.72e-3 (Ylilammi 2018, fit) and 4 +/- 1e-3 (Vandalon and Kessels 2016, SFG); the Arts et al.
2019 range at a 275 C set-point is 0.5-2e-3. The uncorrected 8.5e-3 was a factor 1.5 and 2.1 above the 300 C values
and 4 above the top of the Arts range. The corrected 1.96e-2 is a factor 3.4 and 4.9 above the 300 C values and 9.8
above the top of the Arts range. The corrected fit therefore lands farther from every literature value, not closer:
the 24 % transport error did not explain the high s, it hid part of it. Three remarks on what that means.
1. The direction is forced by the data. The measured front is sharper than a ballistic slit model with the literature
   s would make it: Machball Trench at AR 500 and s = 8.5e-3 gives a front 0.8 -> 0.2 width of 48 gaps at z50 = 187,
   against 27.6 for the uncorrected 1D model at the same s (which is what the data matched) and 42 for the corrected
   one. A fully ballistic fit would push s higher still, to roughly 8.5e-3 x (48/27.6)^2 = 2.6e-2.
2. Which transport the literature extraction used decides which of our numbers is comparable to it. The Arts slope
   relation dtheta/dz~ = sqrt(s/13.9) comes from the Yanguas-Gil and Elam continuum model, i.e. a Knudsen-diffusion
   slit model of the same family as our uncorrected one; a value extracted that way is a hydraulic-diffusivity value,
   and the uncorrected 8.5e-3 is the like-for-like number against the 2-7e-3 that Arts' Table I derives from the
   Ylilammi data by the same relation. Whether Ylilammi's own fit (5.72e-3) carried a logarithmic slit diffusivity
   could not be checked here (the full text is not readable through the tools of this session); WP1 records that its
   kinetic form (lumped coefficient with an equilibrium constant) differs anyway, and reversible adsorption broadens a
   front the way a lower s does.
3. The ideal-plate factor is an upper bound for this structure. The LHAR3 cavities are 70-100 um wide (sheet names
   W70-W100 in the workbook; 140-200 gaps) and carry support pillars (workbook metadata: "support pillars used",
   denser near the cavity mouth; the layout is in the SI of Yim et al. 2020, not readable here). Side walls and
   pillars interrupt the grazing flights that make the logarithm. The duct Monte Carlo (same code with diffuse side
   walls) gives, at L/w = 250, factors 2.17 (infinite width), 1.94 (b/w 200), 1.89 (140), 1.70 (50), 1.49 (20)
   against the formula's 2.07; at L/w = 500, 2.30 (infinite), 2.21 (b/w 200), 1.90 (140) and 1.88 (50), with a
   statistical noise of about +/- 0.05 at 2e5 particles. The finite width alone takes 8-18 % off the factor; a pillar pitch of a few tens of gaps would take it most of the way back to the hydraulic value.
   Without the pillar layout the honest statement is that the TMA sticking coefficient from this profile lies
   between the uncorrected 8.5e-3 (cut-off at about 14 gaps, i.e. the hydraulic value) and the ideal-slit 1.96e-2
   (cut-off at the channel length), with the ballistic-consistent ideal-slit value near 2.6e-2.
What the chamber conclusions inherit: nothing new for the patterned wafer (WP1: 14.00-14.20 s for s = 5e-4 to
8.5e-3; the wafer is supply-limited and a higher s only makes it more so; 2e-2 lies outside the swept range, so this
is an extrapolation of a monotone trend, not a run). The blanket-wafer times, which scale as about s^-0.4, would be
about 0.7 x the s = 8.5e-3 values at 2e-2 (0.33 s -> about 0.24 s worst face), extending the published-range spread
from a factor 3 to about 4.

## Word-line recesses restated (task 3; coupling/wordline_recess_corrected.py; figures/wordline_recess_corrected.csv, .png)
Same geometries, TMA parameters (s 5e-3, N0 4e18), chamber drive (centre-face density of
reactor/showerhead_pulse_patterned_AR100_long, one-way) and constant 0.3 Pa benchmark as coupling/wordline_recess.py,
with `knudsen_correction="plates"` on the nested model, the plain slit and the area-equivalent slit. Factors: slit
1.123 (250 nm, AR 20) and 1.467 (100 nm, AR 50); recesses 0.950 (29 nm x 365 nm, AR 12.6; below 1 because AR < 14.4)
and 1.328 (29 nm x 1 um, AR 34.5). Wall-area ratios and the rule are geometry only and unchanged (150-189; 38.3-48.0 s).

| geometry | chamber-driven slit bottom t90, old -> new | last recess bottom, old -> new (lag) | plain slit, old -> new | area-equivalent slit, old -> new | constant 0.3 Pa slit bottom / last recess, old -> new | exposure, old -> new (ratio to effusion limit) |
|---|---|---|---|---|---|---|
| w250_Lr365 | 5.37 -> 5.09 s | 5.47 -> 5.20 s (0.09 -> 0.10 s) | 2.30 -> 2.26 s | 5.39 -> 5.12 s | 2.70 / 2.79 -> 2.45 / 2.55 s | 0.84 -> 0.76 Pa s (4.6 -> 4.2) |
| w250_Lr1000 | 8.33 -> 7.96 s | 9.06 -> 8.49 s (0.72 -> 0.53 s) | 2.30 -> 2.26 s | 8.70 -> 8.20 s | 6.02 / 7.03 -> 5.54 / 6.24 s | 2.11 -> 1.87 Pa s (4.3 -> 3.8) |
| w100_Lr365 | 13.46 -> 10.99 s | 13.53 -> 11.07 s (0.07 -> 0.08 s) | 3.76 -> 3.27 s | 13.49 -> 11.03 s | 14.53 / 14.68 -> 10.03 / 10.17 s | 4.40 -> 3.05 Pa s (9.7 -> 6.7) |
| w100_Lr1000 | 16.52 -> 15.36 s | 16.73 -> 15.49 s (0.20 -> 0.13 s) | 3.76 -> 3.27 s | 16.66 -> 15.44 s | 37.2 / 39.2 -> 25.6 / 26.9 s | 11.7 -> 8.06 Pa s (9.6 -> 6.6) |

Against the voxel solver (WP9 (d), w100_Lr365 at 0.3 Pa: slit floor 10.70 s, deepest recess end 10.88 s, lag 0.18 s):
the corrected nested model gives 10.03 / 10.17 s, lag 0.14 s, i.e. 6.2 / 6.5 % fast where the uncorrected model was
36 / 35 % slow (14.53 / 14.68 s). The 6 % is the continuum-versus-ballistic offset of an AR-50 slit (2.5 % on the plain
slit against Machball) plus the recesses' own sub-problems and the voxel rounding of g_r and p_r to 30 / 50 nm (the
nested model with rounded dimensions was 1.7 % slower than nominal in WP9). The conclusions of 3.11 stand with
smaller margins: the structure saturates within 0.5-3 % of the plain slit with the same wall area (was 4 %), the last
recess bottom lags the slit bottom by 0.08-0.53 s (0.8-6.6 %; was 0.07-0.72 s, 0.5-8.6 %), and the chamber-driven
times move by 5-7 % for the 250 nm slits and 7-18 % for the 100 nm slits, less than the 24-36 % of the constant-density
benchmark because the AR-100 wafer's drive is supply-limited: the slit's conductance is not the bottleneck.

## Where else the slit transport was used (task 4)
`grep -rn "Slit(" coupling analysis tests scripts docs` (RecessedSlit and AreaEquivalentSlit included):
- coupling/fit_pillarhall.py, Slit(500 nm, 250 um): the quoted s = 8.5e-3, 75 Pa, 7.5 Pa s, rms 0.038 (README,
  CLAIMS row "PillarHall fit", paper 2.3, Discussion and the sticking table, docs/index.html). Refit above.
- coupling/wordline_recess.py: RecessedSlit (nested model), Slit (plain slit), AreaEquivalentSlit, all four
  geometries; figures/wordline_recess.csv and .png; CLAIMS rows "Word-line recesses ..." and "Rule for a word-line-
  slit wafer" (the rule row is geometry only and unaffected); paper 3.11 table. Restated above.
- coupling/voxel_mc_crosscheck.py part `slit`, Slit(100 nm, 5 um): the 1D reference 1.381 s in
  figures/voxel_mc_crosscheck.csv and CLAIMS row "Slit AR 50 at 0.3 Pa"; part `recess`, RecessedSlit w100_Lr365: the
  nested-model 14.53 / 14.68 s in CLAIMS row "Voxel word-line slit w100_Lr365". These rows state the uncorrected value
  as the thing being measured and stay correct as written; the new rows below add the corrected values.
- docs/build_assets.py: Slit(500 nm, 250 um) at s = 8.5e-3, 75 Pa to draw the page's PillarHall panel (not touched;
  the page quotes 8.5e-3 in two places and needs the integrator's edit).
- tests/test_scaling.py: Slit(500 nm, 200 um), scaling laws (sqrt(dose), AR^2, synthetic-profile recovery) that hold
  for any constant D and are unaffected; tests/test_recess.py: Slit(200 nm, 5 um) as the no-recess limit, unaffected
  (limits are tested against the same transport).
- Not slit-based and unaffected: every hole result (Hole, TaperedHole; d vbar/3 is exact), the chamber boundary
  conditions (holes), coupling/multispecies.py (holes), analysis/ (no Slit).

## Proposed paper text
2.3, replacing the PillarHall bullet (the last two sentences of the current bullet go):
```
- Validation against the Zenodo 3974438 PillarHall Al2O3 profiles (TMA/H2O, 300 C, 500 cycles x 0.1 s, 500 nm
  lateral channel, eight profiles averaged): the slit model with channel narrowing fits the saturation front to an
  rms of 0.038 of the plateau. The front position (AR 215 at half thickness) fixes the exposure and its width fixes
  s, but both depend on the slit's Knudsen diffusivity. With the hydraulic-diameter value D_K = 2 w vbar/3 the fit
  gives s_TMA = 8.5e-3 (profile-to-profile factor 1.5) and 7.5 Pa s per pulse (75 Pa TMA for 0.1 s). That value of
  D_K is exact for a cylinder but not for parallel plates, where grazing flights make the free-molecular diffusivity
  grow as (w vbar/4) ln(L/w) (Clausing; Berman 1965): the 1D slit model is 24 % slow against Machball's ballistic
  trench and the voxel solver at AR 50 and 61 % at AR 100, and with the factor (3/8) ln(L/w) on D_K it is within
  5 % of Machball from AR 20 to 100 and reproduces the ballistic front width to about 10 % at AR 500. With that
  correction (a factor 2.33 at AR 500) the same data give s_TMA = 1.96e-2 (sheet-to-sheet x/÷ 1.46) and 3.23 Pa s per pulse
  (32.3 Pa), at the same rms: the profile cannot distinguish the two transport models, it only re-labels the
  parameters, s by the factor D rose and the exposure by its inverse. Published values are 0.5-2e-3 (Arts 2019,
  275 C set-point, extracted with a continuum slit model of the hydraulic-diffusivity family), 4e-3 (Vandalon and
  Kessels 2016, SFG, 300 C) and 5.72e-3 (Ylilammi 2018, 300 C). The corrected value lies farther from all of them
  than the uncorrected one; the measured front is sharper than a ballistic slit with the literature s would make it.
  The ideal-plate factor is an upper bound for the real cavity, which is 140-200 gaps wide and carries support
  pillars that interrupt the grazing flights (finite width alone removes 8-13 % of the factor), so s from this
  profile lies between 8.5e-3 and 1.96e-2 depending on the cut-off length of the cavity; the chamber runs use 5e-3
  and the patterned-wafer conclusions are supply-limited and insensitive to it (Section 4).
```
3.11, replacing the last sentence of the word-line paragraph ("The voxel Monte Carlo solver ... follows."):
```
The 1D slit transport behind these numbers was 24 % slow against the ballistic solvers (Section 2.6); with the
parallel-plate correction of Section 2.3 the chamber-driven times become 5.1 / 8.0 s (250 nm slits) and 11.0 /
15.4 s (100 nm slits), the recess lags 0.1-0.5 s (1-7 %), the area-equivalent slit holds within 3 %, and the
constant-density w100_Lr365 structure saturates in 10.0 / 10.2 s against the voxel solver's 10.7 / 10.9 s (6 % fast,
the continuum-versus-ballistic offset of an AR-50 slit) instead of 14.5 / 14.7 s (36 % slow).
```
Discussion, sticking paragraph, replacing "and 8.5e-3 (this work, fit to the Zenodo 3974438 profiles at 300 C, factor
1.5 sheet-to-sheet)":
```
and 8.5e-3 to 1.96e-2 (this work, fit to the Zenodo 3974438 profiles at 300 C, factor 1.5 sheet-to-sheet; the lower
value with the hydraulic-diameter slit diffusivity of the literature extractions, the upper with the parallel-plate
diffusivity that matches the ballistic solvers, Section 2.3)
```
and, after "It sits a factor 1.5 to 2 above the two other 300 C values ...":
```
Correcting the slit transport does not close that gap: the 1D slit diffusivity was 24 % low at AR 50 and a factor
2.3 low at the channel's AR 500, and with the exact ballistic front width the same measured profile needs a larger
s, not a smaller one.
```

## CLAIMS rows (CLAIMS.md format)
| PillarHall fit, corrected slit transport (2026-11-11) | with the parallel-plate factor (3/8) ln(L/w) = 2.33 on the slit diffusivity: s_TMA = 1.96e-2 (sheet-to-sheet x/÷ 1.46), exposure 3.23 Pa s per 0.1 s pulse (32.3 Pa), rms 0.0381, same data and procedure as the 8.5e-3 / 75 Pa fit, which stays as the hydraulic-diffusivity value; s moves away from the 300 C literature values (3.4 x Ylilammi, 4.9 x Vandalon-Kessels); the ideal-plate factor is an upper bound for the pillar-supported cavity | figures/pillarhall_fit_corrected.png, .csv; paper/wp/wp9b_slit_transport.md | PYTHONPATH=src .venv/bin/python coupling/fit_pillarhall_corrected.py |
| Slit Knudsen diffusivity correction | D_slit / (2 w vbar/3) = (3/8) ln(L/w) (Clausing slit transmission W = (w/L)[ln(L/w) + 0.10], transmission Monte Carlo: W 0.6844 / 0.2399 / 0.0802 / 0.0124 at L/w 1 / 10 / 50 / 500); Machball Trench t90 0.6231 / 1.1170 / 2.6305 s at AR 20 / 50 / 100 vs 1D 0.6149 / 1.3842 / 4.2305 s uncorrected (-1 / +24 / +61 %) and 0.5993 / 1.0891 / 2.6234 s corrected (-3.8 / -2.5 / -0.3 %); matched factors 0.95 / 1.405 / 1.721 vs formula 1.123 / 1.467 / 1.727; front width at AR 500 within 6-13 % of Machball | src/aldconf/transport.py (slit_knudsen_factor), geometry.py (Slit.knudsen_correction), tests/test_slit_transport.py, paper/wp/wp9b_slit_transport.md | pytest tests/test_slit_transport.py (Machball within 5 % at AR 20 and 50) |
| Word-line recesses, corrected slit transport | chamber-driven slit bottom t90 5.09 / 7.96 s (250 nm slit, L_r 365 / 1000 nm), 10.99 / 15.36 s (100 nm slit); last recess bottom 0.10 / 0.53 / 0.08 / 0.13 s later (1-7 %); area-equivalent slit within 3 %; constant 0.3 Pa w100_Lr365 10.03 / 10.17 s vs voxel 10.70 / 10.88 s (6 % fast; uncorrected 14.53 / 14.68 s, 36 % slow); rule and wall ratios unchanged | figures/wordline_recess_corrected.csv, .png | PYTHONPATH=src .venv/bin/python coupling/wordline_recess_corrected.py |

## Run log (for paper/assumptions.md)
## Night 4, part 4: slit transport correction (src/aldconf/transport.py slit_knudsen_factor, geometry.py Slit.knudsen_correction, recess.py; coupling/fit_pillarhall_corrected.py, coupling/wordline_recess_corrected.py; figures/pillarhall_fit_corrected.png, .csv, wordline_recess_corrected.csv, .png; tests/test_slit_transport.py)
The 24 % slit error of WP9 (f) is the parallel-plate Knudsen diffusivity: D_K = d_h vbar/3 with d_h = 2 w is exact for
a cylinder and not for plates, where grazing flights give D_slit = (w vbar/4)[ln(L/w) + c], i.e. (3/8)[ln(L/w) + c]
times D_K (Clausing slit transmission; Berman 1965; Steckelmacher 1986). A transmission Monte Carlo of the infinite
slit (4e5 particles per point) gives W = 0.6844 / 0.5425 / 0.3572 / 0.2399 / 0.1545 / 0.0802 / 0.0472 / 0.0270 / 0.0124
at L/w = 1 / 2 / 5 / 10 / 20 / 50 / 100 / 200 / 500, so W L/w - ln(L/w) = +0.10 from L/w 10 to 200. Machball
Trench(AR, 8 AR) at s 5e-3, N0 4e18, 573 K, 0.3 Pa, 100 nm gap: t90 0.5335 / 0.6231 / 1.1170 / 2.6305 s at AR 10 / 20 /
50 / 100 (Nz convergence 0.2-0.3 %) against the 1D slit model 0.5087 / 0.6149 / 1.3842 / 4.2305 s (bottom coverage
0.9, 200 cells, 200 substeps; 400/400 moves AR 50 by 0.2 %); the knudsen_factor that matches Machball is 0.59 / 0.95 /
1.405 / 1.721 (AR 10-20 reaction-limited, kinetic floor 0.47 s, factor ill-conditioned) against (3/8) ln(AR) = 0.863 /
1.123 / 1.467 / 1.727; the 4 % at AR 50 is the continuum-vs-ballistic offset of WP9 (a). Implemented as
Slit(knudsen_correction="plates") (keyword-only, default "none" reproduces every earlier number), applied in simulate,
simulate_driven, multispecies and simulate_recessed (slit and recesses separately), factor (3/8) ln(L/w) without
the constant: corrected 1D t90 0.5143 / 0.5993 / 1.0891 / 2.6234 s, -3.6 / -3.8 / -2.5 / -0.3 % against Machball (with
c = 0.10: -3.9 / -4.5 / -3.9 / -2.0 %). Front shape (Machball Trench vs 1D at equal times): AR 500, s 8.5e-3, 0.8 -> 0.2
width 40-48 gaps ballistic over z50 60-187, 42.1 corrected (6-13 % narrow), 27.6 uncorrected (31-43 % narrow); front
position corrected 14-30 % deep, uncorrected 15-25 % shallow. PillarHall refit (same data, model, nine starts, per-sheet
refits; factor 2.330 at AR 500): s = 1.96e-2, 32.3 Pa, 3.23 Pa s per pulse, rms 0.0381 (unchanged), sheet-to-sheet
s = 1.99e-2 x/÷ 1.46 (geometric mean of the eight sheets; range 7.9e-3 to 2.7e-2), p = 32.4 Pa x/÷ 1.02; the nine starts agree within 1 %. s rises by the factor D rose and the pressure falls by it; against the 300 C
literature (5.72e-3, 4e-3) the corrected value is 3.4-4.9 x high (was 1.5-2.1 x), i.e. the correction moves the fit
away from the literature; the Arts slope relation itself is a hydraulic-diffusivity continuum result, so 8.5e-3 is the
like-for-like number against it; a fully ballistic fit would give about 2.6e-2; the ideal-plate factor is an upper bound
for the 140-200-gap-wide, pillar-supported cavity (duct Monte Carlo with side walls: 1.94 / 1.89 / 1.70 / 1.49 at b/w
200 / 140 / 50 / 20 against 2.17 infinite at L/w 250; 2.21 / 1.90 / 1.88 at b/w 200 / 140 / 50 vs 2.30 at L/w 500). Word-line recesses
restated (factors 1.123 / 1.467 slit, 0.950 / 1.328 recess): chamber-driven slit bottom 5.09 / 7.96 / 10.99 / 15.36 s
(was 5.37 / 8.33 / 13.46 / 16.52), last recess +0.10 / +0.53 / +0.08 / +0.13 s, plain slit 2.26 / 3.27 s, area-equivalent
within 0.5-3 %; constant 0.3 Pa w100_Lr365 10.03 / 10.17 s vs voxel 10.70 / 10.88 s (6 % fast; was 36 % slow); rule
and wall ratios unchanged. Everything else that used Slit: voxel_mc_crosscheck part slit/recess references (stay as
the uncorrected values they measure), docs/build_assets.py PillarHall panel (not touched), test_scaling and
test_recess (transport-independent). pytest 124 passed in 32 s (118 + 6). Wall time: Machball runs 1-20 s each,
transmission Monte Carlo 2 min, recess script 10 s, refit about 25 min single core. Caveats: the correction is for an
ideal, infinitely wide slit closed at one end with the cut-off at the channel length; finite width, pillars or an
open far end change the cut-off; it is a few per cent fast against ballistic solvers at AR 20-50 (the diffusion model's
own offset) and 14-30 % deep in front position at AR 500; the factor falls below 1 for AR < 14.4, where the slit is
reaction-limited and the factor hardly matters; no narrowing in the Machball comparisons.
