# Four checks on public data after the second blind test (post hoc, 2026-11-26)

None of this is part of a registered test. Every choice below was made after the results of both blind tests were known.
Numbers only; each line names the file that holds them.

## 1. Front shape without the optical tail: SEM-EDS line scans (Yim 2020, Fig. 5b and Fig. S8)
`eds_cycle_series.py`, `eds_cycle_series.csv`, `eds_cycle_series.png`; digitised traces and their provenance in `data/lhar_eds/`.
The model is the one registered for the second test (hydraulic kernel, s = 0.04, spot of two gaps), with and without narrowing.
Background 18 a.u. (the level far from the film in the unobstructed 1000-cycle trace), plateau = median at 10 to 60 um depth.

| cycles | fill | channel | rms with narrowing | rms without | measured slope (1/um) |
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
`residual_nm.py`, `residual_nm.csv`. Thickness reported beyond 1.3 half depths: 0 to 12.6 nm (median over samples about 5 nm),
with no scaling with the number of cycles (1000 cycles: 3.0 nm; 500 cycles: 0 to 12.6 nm) and local bumps at 127 to 129 um and
172 to 173 um, three and four times the 42.4 um pillar-row spacing. The EDS traces of the same samples fall to the background.

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
23 % (He-Xe) below the reference; a longer run does not change it (He-Xe, run length 16 and 400: 2.28 and 2.43 in the
reference's units); the cause is not established (hard spheres against the model collision term is one candidate).
The hydraulic formula is 21 to 52 % below the reference for widths of 0.1 and 0.01.
TMA in N2 at 300 C, 500 nm (Monte Carlo over hydraulic): 2.34 at 30 Pa, 2.05 at 100 Pa, 1.62 at 300 Pa, 1.38 at 1000 Pa,
1.13 at 3000 Pa, 1.03 at 30000 Pa.

## 4. One process in two laboratories: ZnO from DEZ and H2O at 200 C
`docs/blindtest/posthoc/zno_cross_lab.csv`. Same model, sticking coefficient fitted to each profile: 0.004 at 1400 Pa
(Haimi 2026, sample B) and 0.020 to 0.025 at 86 Pa (Gonsalves 2026, left and right). Replacing the hydraulic kernel by the
gas-collision kernel raises the low-pressure value more than the high-pressure one, so it widens the difference.
