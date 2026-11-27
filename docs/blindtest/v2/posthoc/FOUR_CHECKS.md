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
