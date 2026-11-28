# Four checks on public data after the second blind test (post hoc, 2026-11-26)

None of this is part of a registered test. Every choice below was made after the results of both blind tests were known.
Numbers only; each line names the file that holds them.

## 1. Front shape without the optical tail: SEM-EDS line scans (Yim 2020, Fig. 5b and Fig. S8)
`eds_cycle_series.py`, `eds_cycle_series.csv`, `eds_cycle_series.png`; digitised traces and their provenance in `data/lhar_eds/`.
The model is the one registered for the second test (hydraulic kernel, s = 0.04, spot of two gaps), with and without narrowing.
Background 18 a.u. (the level far from the film in the unobstructed 1000-cycle trace), plateau = median at 10 to 60 um depth.

| cycles | fill | channel | rms with narrowing | rms without | measured slope (per gap height) |
|---|---|---|---|---|---|
| 250 | 0.10 | left / right | 0.039 / 0.051 | 0.045 / 0.056 | 0.040 / 0.040 |
| 500 | 0.22 | left / right | 0.038 / 0.038 | 0.060 / 0.061 | 0.020 / 0.021 |
| 1000 | 0.40 | left / right | 0.029 / 0.030 | 0.108 / 0.105 | 0.011 / 0.012 |

Mean rms 0.038 with narrowing, 0.072 without; six of six channels within 0.06 with narrowing, two of six without.
Limits: three samples from one laboratory and one process; the same samples as the reflectometry series, so this checks the
measurement, not a new process. The 250-cycle trace is partly hidden behind the others in the source figure (92 of its points
are flagged) and separates the two models weakly (fill 0.10). A first run with a background estimated per trace failed on the
250- and 500-cycle traces; the single background level above replaced it.

## 2. What sits behind the front in reflectometry profiles
`residual_nm.py`, `residual_nm.csv`. Thickness reported beyond 1.3 half depths: per-sample median 0 to 12.6 nm (single readings up to 33.7 nm),
with no scaling with the number of cycles (1000 cycles: 3.0 nm; 1080 cycles: 11.1 nm; 500 cycles: 0 to 12.6 nm) and local bumps at 127 to 129 um and
172 to 173 um, three and four times the 42.4 um pillar-row spacing. The EDS traces fall to a background that does not scale with
cycle number; a constant offset is removed with that background, so they do not by themselves exclude one. The three EDS samples
have little optical residual of their own (no points, no median, 3.0 nm).

## 3. Transport along the slit: Monte Carlo and hydraulic formula against kinetic theory
`docs/paper/benchmark_slit_diffusion.py`, `figures/benchmark_slit_diffusion.csv`, `figures/slit_diffusion_tma_n2.csv`;
`pillarmc slitdiffusion`; test `tests/test_pillarmc.py::test_slit_diffusion_limits`.
Reference: Siewert and Valougeorgis, Eur. J. Mech. B/Fluids 23, 645 (2004), heavy-species flow in the diffusion problem at
10 % heavy species. Ratio to the reference (slit width in viscosity mean free paths):

| width | Ne-Ar: Monte Carlo | Ne-Ar: hydraulic | He-Xe: Monte Carlo | He-Xe: hydraulic |
|---|---|---|---|---|
| 0.01 | 0.91 | 0.52 | 0.77 | 0.48 |
| 0.1 | 0.97 | 0.79 | 0.90 | 0.74 |
| 0.5 | 0.97 | 0.97 | 1.07 | 1.01 |
| 1 | 1.03 | 1.00 | 1.10 | 1.12 |
| 10 | 1.03 | 1.01 | 1.30 | 1.35 |
| 100 | 1.02 | 1.03 | 1.48 | 1.45 |

The Monte Carlo follows a dilute tracer; the reference has 10 % heavy species, whose counterflow at diffuse walls (Graham's
law) lowers the heavy flux in the continuum limit by 1/1.04 (Ne-Ar) and 1/1.47 (He-Xe) (column
`tracer_over_mixture_continuum`). That accounts for the wide-slit end. At the narrow end the Monte Carlo is 9 % (Ne-Ar) and
23 % (He-Xe) below the reference, and a longer run does not change it. The same walk with relaxation collisions (constant rate
kT / (m D12), full Maxwellian redraw: the kind of collision term the reference uses) gives 1.02, 1.00, 1.03 and 1.00 times the
reference for Ne-Ar at widths 0.01, 0.1, 0.5 and 1 (column `monte_carlo_relaxation_collisions`). So for Ne-Ar the narrow-end
difference is the collision model (hard spheres against relaxation), not the walk. For He-Xe the relaxation walk is 1.20 to
1.21 times the reference at all four widths: a 10 % xenon mixture is not a tracer problem, and He-Xe is not a clean test.
The hydraulic formula is 21 and 48 % below the reference at widths of 0.1 and 0.01.
TMA in N2 at 300 C, 500 nm (Monte Carlo over hydraulic): 2.34 at 30 Pa, 2.05 at 100 Pa, 1.62 at 300 Pa, 1.38 at 1000 Pa,
1.12 at 3000 Pa, 1.02 at 30000 Pa.

## 4. One process in two laboratories: ZnO from DEZ and H2O at 200 C
`docs/blindtest/posthoc/zno_cross_lab.csv`. Same model, sticking coefficient fitted to each profile: 0.004 at 1400 Pa
(Haimi 2026, sample B) and 0.020 to 0.025 at 86 Pa (Gonsalves 2026, left and right). No run with the gas-collision kernel was
made for these two profiles.

# Three follow-ups (2026-11-27)

## 5. Gap exponent of the narrowing law against pressure
`docs/paper/narrowing_exponent.py`, `figures/narrowing_exponent.csv`; `aldconf.lhar.gap_exponent`, `narrowing_limit(fill, n)`.
n = 1 + (1/2) d ln D / d ln p from the slit diffusivity of the Monte Carlo (TMA in N2, 300 C, 500 nm):
0.926 at 30 Pa, 0.916 at 100 Pa, 0.903 at 300 Pa, 0.881 at 1000 Pa, 0.877 at 1400 Pa, 0.844 at 3000 Pa, 0.775 at 10000 Pa.
Measured directly from pairs of front simulations at 400 and 500 nm: 0.917 (300 Pa), 0.892 (1400 Pa), 0.859 (3000 Pa).
A step front smeared over the stretch f^n has slope x depth = (1 - fill/2) / (n fill): the limit is 1/n times the gas-free one.

## 6. Slope tables over pressure, gap and spot size
`docs/paper/run_pillar_mc.py` (cases gas30, gas3000, gas10000 added), `docs/paper/build_slope_table.py`, `figures/slope_table.csv`
(8928 rows: hydraulic, gas-free and N2 at 30, 100, 300, 1000, 3000, 10000 Pa; spot 0 to 55 gaps), `aldconf.lhar.slope_curve_at`
and `sticking_from_slope_at` (interpolation in ln(pressure x gap), 15 to 5000 Pa um).
Check of the pressure x gap scaling, `docs/paper/slope_table_scaling_check.py`, `figures/slope_table_scaling_check.csv`: direct
simulations of the 2000 nm and 100 nm cavities at 300 Pa against the 500 nm tables, 198 points: rms difference in ln(slope)
0.049, largest 0.154; mean +0.03 to +0.07 for 100 nm and -0.01 to -0.03 for 2000 nm (the pillar spacing in gap units differs:
490 and 24.5 against 98). With a filled fraction of 0.1 or more and the spot: rms 0.014 to 0.035.
Example, slope 0.019 per gap at depth 208 gaps and filled fraction 0.19 (Yim 2020 sample 4, nominal gap): s = 0.022 at 100 Pa,
0.015 at 300 Pa, 0.0096 at 1000 Pa, 0.0068 at 3000 Pa; the slope method gives 0.0050.

## 7. Added after an independent review: the EDS front slopes
`eds_slopes.py`, `eds_slopes.csv`. The rms of section 1 is over the whole profile and does not show the front slope.
Between 25 and 75 % of the plateau the EDS fronts are steeper than the model with narrowing (EDS over model 1.21, 1.22 at
250 cycles; 1.41, 1.50 for the trace of sample 2; 1.16, 1.21 at 1000 cycles), steeper than the reflectometry slopes of the
same samples in Yim 2020 Table 1 (1.91, 1.92; 1.79, 1.93; 1.12, 1.17) and at 1.02 to 1.43 times the step-front value
(1 - f/2)/(f z50). In the model the spot size (1 um against 5.5 um) changes the slope by 3 % at 250 cycles and not at all
at 500 and 1000. The slope fits use 7, 14 and 23 digitised points. So the two measurements of the same samples disagree
on the slope by up to a factor of 1.9, this is not explained, and the statement in the first blind test's diagnosis that
it is not shown whether a clean profile can exceed the step-front value no longer holds: these do.
The comparison "without narrowing" in sections 1 and 5 is the same model with narrowing switched off at the same
s = 0.04, not a refitted model; s = 0.04 is 4.5 times the hydraulic value calibrated on the single-exposure profiles.

## 8. The slope-definition check (after section 7)
`slope_definition.py`, `slope_definition.csv`, `slope_definition_heldout.csv`. Section 7 compared the EDS slopes with the
slopes tabulated in Yim 2020 Table 1 (the authors' own line fits). Here the slope is read from the digitised reflectometry
profiles themselves with the definition used for the EDS scans (least squares between 25 and 75 % of the plateau, on the
contiguous front).
- For 8 of the 10 reflectometry profiles the profile slope is 0.97 to 1.10 times the tabulated one (median 1.08).
- Sample 2: profile slope 0.0169 against a tabulated 0.011 (1.54 times). Sample 1: 0.0094 against 0.013 (0.72), but its
  trace lies under the markers of sample 2 in the source figure and several points are missing (custodian's note), so its
  profile slope is not reliable.
- EDS against reflectometry for the same sample, same definition: 1.21 for sample 2 (not 1.8 to 1.9) and 1.07 for sample 11
  (1000 cycles). For the 250-cycle sample the reflectometry trace is hidden below 0.4 of the plateau and gives no slope, so
  the factor 1.9 against its tabulated slope cannot be checked.
- So most of the "factor 1.9" of section 7 is the tabulated slope of sample 2, not a disagreement between the two
  measurements. What remains: EDS fronts 7 to 21 % steeper than the reflectometry fronts, and both steeper than the model
  with narrowing; the EDS fronts still lie at 1.0 to 1.4 times the step-front value.
- The held-out slope test rescored on the profile slopes (9 samples at 500 nm with a digitised profile): complete models
  0.18 to 0.19 (25-75 % slope) and 0.22 to 0.23 (40-60 % secant) against 0.28 to 0.33 without narrowing; without sample 1:
  0.12 to 0.16 and 0.11 to 0.18 against 0.22 to 0.27. A baseline with no kinetics, k (1 - f/2)/(f z50), does as well on
  these samples (0.18 and 0.20; without sample 1 0.14 and 0.09): at about 20 % filling the slope is set by the filled
  fraction and the depth, and carries little information on the sticking coefficient.

## 9. Every front fitted on its own, with and without narrowing (final version, after two independent reviews)
`docs/paper/profile_fits.py`, `docs/paper/profile_fits_stats.py`, `figures/profile_fits*.csv`, `figures/profile_fits_stats.csv`.
Sixteen digitised TMA/water fronts of eleven samples (ten reflectometry, six EDS channels of three scans) and four ZnO fronts,
sticking coefficient and exposure free. Simulated range 0.0016 to 0.2 (Monte Carlo, 11 to 15 values) and 0.001 to 0.2
(hydraulic, 17 values). Model tables aligned at common front depths and interpolated in ln s; block bootstrap of the centred
residuals with the block length from their autocorrelation. Two fit windows (to 0.1 of the plateau, run first; to 0.3, run
second because several reflectometry curves keep a residual above 0.1 behind the front) and two film thicknesses (measured in
the channel, which is the appropriate one, and outside it, 5 to 9 % smaller). Earlier versions of this section quoted the
favourable window and the pooled set only, with intervals that were too narrow; two reviews found this.
`figures/profile_fits_stats.csv` holds every statistic. Checks on the final files: every estimate lies inside its interval,
no interval has zero width, the interpolation check is below 0.01 of the plateau for every fit.
- Without narrowing (film in channel), d ln s / d f over all 16 fronts: -6.5 +/- 1.0 (hydraulic), -6.2 +/- 0.9 (gas),
  -6.9 +/- 0.9 (gas-free) in the 0.3 window; -6.1 +/- 1.6, -6.0 +/- 1.6, -6.8 +/- 1.3 in the 0.1 window. With one sample
  left out: -4.9 to -12.0. Nine reflectometry fronts at 500 nm alone: -5.2 to -6.0 (0.3 window), -2.9 to -5.2 (0.1 window,
  0.8 to 1.7 standard errors). Hydraulic, EDS, 0.3 window: 0.012 and 0.0096 (250 cycles), 0.0054 and 0.0048 (22 % filled),
  0.0019 and 0.0020 (1000 cycles).
- With narrowing, film in channel, all 16: 0.9 +/- 2.5, -2.2 +/- 2.0, 0.7 +/- 2.1 (0.3 window); 2.5 +/- 3.8, 1.1 +/- 3.6,
  2.5 +/- 3.3 (0.1 window). Censored regression (fits without an upper limit as lower bounds): -0.2 +/- 1.9, -2.2 +/- 1.6,
  -0.3 +/- 2.8; 1.8 +/- 3.6, 0.2 +/- 3.3, 5 +/- 9. Ten reflectometry fronts: -4.2 +/- 2.8, -4.8 +/- 2.7, -1.6 +/- 3.7 (0.3).
- With narrowing, film outside the channel, all 16, 0.3 window: -2.7 +/- 2.0, -4.0 +/- 1.8, 0.3 +/- 2.4; ten reflectometry
  fronts: -5.8 +/- 2.6, -6.0 +/- 2.4, -3.8 +/- 3.9. So the correction accounts for between roughly half and all of the
  decline, depending on the film thickness to within a few per cent and on the set of fronts.
- The pooled zero with narrowing rests on the 1000-cycle sample (three fronts, the only ones above f = 0.22): without it the
  slope is -2.4 +/- 4.8, -2.3 +/- 4.3, -0.7 +/- 4.4 (0.3 window) and -6.3 +/- 7.7, -5.8 +/- 7.3, -4.0 +/- 7.1 (0.1 window). The
  uncorrected decline survives the same omission (-10.1 +/- 1.8, hydraulic, 0.3 window). Share of the uncorrected decline
  that remains with narrowing on the ten reflectometry fronts (0.3 window): 53 % and 68 % (hydraulic, gas; film in channel),
  70 % and 82 % (film outside). So the correction accounts for between about a third and all of the decline.
- Fits with narrowing whose interval has no upper limit: 3, 2, 10 of 16 (0.3 window); 4, 4, 12 (0.1 window).
- Scatter (sd of ln s, 16 fronts): 0.79 to 0.93 with narrowing against 0.71 to 0.77 without (0.3); 1.25 to 1.45 against
  0.85 to 0.87 (0.1). Window-to-window change of single values with narrowing: median factor 1.6 to 2.1, up to about 5;
  sample 16 by 10 to 38 because its 0.1-window fit fails (rms 0.12, with and without narrowing).
- Fit quality: mean rms 0.022 to 0.023 without and 0.023 to 0.026 with narrowing (0.3 window); the fit with narrowing has
  the lower rms for 5, 9 and 7 of 16 fronts; for the EDS scans of the 22 % filled sample it is worse by 0.004 to 0.017 in
  every kernel and window.
- Least-filled sample (2000 nm gap, 5 %): 0.032 (0.027 to 0.035) with and 0.027 (0.025 to 0.029) without narrowing
  (hydraulic), 0.035 and 0.033 (gas), 0.065 and 0.056 (gas-free).
- ZnO (Gonsalves 2026, about 10 % filled, hydraulic kernel): narrowing raises s by 1.9 and 2.1 at 200 C (0.3 window; 2.4
  and 2.2 in the 0.1 window) and by 1.1 and 1.2 at 175 C with overlapping intervals (1.3 in the 0.1 window).

## 10. A pair of line scans from a second reactor (2026-11-27)

Al2O3 from TMA and water at 300 C in a Beneq TFS-200 (2.1 to 2.2 mbar, 0.15 s pulses, 1.0 s purges), 250 and 500 cycles in a
500 nm gap: two reflectometry line scans digitised from Fig. 20 of the Chipmetrics LHAR4 analysis guide (V3, 2022), conditions
from Utriainen et al., Proc. SPIE 12008, 120080D (2022). Script `docs/paper/beneq_pair.py`; results
`figures/profile_fits_beneq.csv` and `figures/profile_fits_beneq_slopes.csv`. Made with all other data known.

- Without narrowing the fitted s falls from 250 to 500 cycles by a factor of 2.3 to 4.9 in all 12 combinations of model,
  window and reading (hydraulic, as read: 7.4e-3 and 3.2e-3); d ln s/d f = -7 to -13.
- With narrowing none of the 24 fits has the lower rms (500 cycles, 0.3 window: 0.032 to 0.104 against 0.017 to 0.025), 12
  of 24 have no upper limit, and the best values still fall in 10 of 12 combinations.
- Front slope times depth over the steepest front (1 - f/2)/(n f): 1.3 to 1.5 at 500 cycles, 0.8 to 1.1 at 250 cycles.
- Weight: lower than the main series. Growth per cycle 0.12 to 0.13 nm and 1.0 s purges; 3 to 6 nm behind the front;
  reflectometer objective not stated (5.5 um spot assumed); raster source with 0.56 um pixels.
- The reaction chamber of the Picosun R-150 was not pressure-monitored; about 3 mbar is the intermediate-space pressure
  (Utriainen et al. 2022, same reactor model and recipe as the main series).

## 11. How strongly narrowing acts: the exponent as an empirical dial (2026-11-27)

Script `docs/paper/exponent_scan.py`; results `figures/exponent_scan_{fits,summary,pair,bounds}.csv`. Made after all other
analyses.

- Upper bound on n from fronts steeper than (1 - f/2)/(n f): 0.70 and 0.75 (EDS, 22 % filled sample), 0.66 and 0.68
  (second reactor, 500 cycles); the other line-scan fronts 0.87 or more.
- Second-reactor pair, 500-cycle scan predicted from the 250-cycle scan (no kinetics, no transport model): floor subtracted,
  best n 0.65 at a depth scale of 1.02 (rms 0.022), 0.60 to 0.85 within 1.5 times that rms; at n = 1 and equal exposure rms
  0.15. As read: best n 0.4 (rms 0.066), 0 to 0.85.
- Main series, 16 fronts, hydraulic and carrier-gas models: mean rms has a shallow minimum at n = 0.7 (6 to 8 % below n = 1,
  1 to 3 % below n = 0). d ln s/d f at n = 0.7: -5.0 to -6.3 (83 to 96 % of the uncorrected slope); it is within about one
  standard error of zero only for n of 0.92 or more.
- So no single strength of narrowing both fits the three steep fronts and removes the decline of the main series. Not
  explained. A film thinner than the reflectometry value would do it for the second reactor (growth per cycle 20 to 30 %
  above the usual value) but not for the main series (two thickness measurements within 5 to 9 %).
- After an independent read (same day): the guide names neither reactor nor conditions; the link to the Beneq runs of
  Utriainen et al. is the identical table of thicknesses and penetration depths, so the conditions are inferred. Three steep
  fronts come from two samples. Without narrowing there are nine distinct combinations, not twelve. With narrowing 18 of the
  24 fits are open on one side, so they do not show whether the fall is removed. Pair prediction, like for like: n = 1 gives
  rms 0.050 with the depth scale free and 0.15 at equal exposure; equal exposure, best n 0.30 (rms 0.061). The prediction is
  approximate and returns n too high by 0.03 to 0.07 on synthetic fronts (`figures/exponent_scan_pair_bias.csv`). In the 0.1
  window the pooled slope cannot exclude n of 0.7 to 0.8.

## 12. Robustness of the profile fits to the gap and to reading errors (2026-11-27)

Script `docs/paper/robustness_fits.py`; results `figures/robustness_{gap,reading}.csv` and `_fits.csv`.

- Measured gap 460 nm in place of 500 nm (hydraulic and carrier-gas models, both windows): uncorrected d ln s/d f -5.7 to
  -6.3 (nominal -6.0 to -6.5); with narrowing 2.4 to 6.1 over 16 fronts (4 to 7 fits without an upper limit) and -1.6 to 2.3
  on the ten reflectometry fronts. Narrowing at its theoretical strength then tends to over-correct.
- Plateau 3 % higher or lower (hydraulic, 0.3 window): uncorrected -6.3 and -6.9; with narrowing -3.8 +/- 1.6 and
  4.9 +/- 2.8. Entrance shifted by 5 gap heights: uncorrected -6.6; with narrowing -0.8 and 2.0.
- So the uncorrected fall is insensitive to both; the corrected slope is not.
- Reference audit the same day: all registered references match Crossref/DataCite; Zenodo 3974438 is dated 2021.
