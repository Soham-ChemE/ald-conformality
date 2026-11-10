# Work package 1: sticking-coefficient sensitivity of the saturating pulse

## What was done
Nine new two-way chamber runs at the baseline operating point (N2 300 Pa, 573 K, 300 sccm, U = 0.0501 m/s, 1 % TMA,
c0 = 3.79e20 /m3, N0 = 4e18 /m2, vbar 410 m/s, 15 mm gap, frozen carrier flow) with the TMA sticking coefficient
varied from 5e-4 to 8.5e-3, i.e. from the bottom of the Arts et al. 2019 range to the PillarHall fit of
coupling/fit_pillarhall.py; 5e-3 is the value used in every other chamber run.
- Patterned 3D-NAND-like wafer (100 nm holes, AR 100, 1.78e13 /m2, wall-area ratio 55.9): s = 5e-4, 1e-3, 2e-3, 8.5e-3,
  24 s pulse, endTime 26 s, writes every 0.05 s (copies of reactor/showerhead_pulse_patterned_AR100_long with only the
  wafer BC line, the pulse table and endTime/writeInterval changed). The existing long case (20 s pulse, endTime 22 s)
  is the s = 5e-3 point.
- Blanket wafer (langmuirSticking BC): s = 5e-4, 1e-3, 2e-3, 5e-3, 8.5e-3, 8 s pulse, endTime 9 s, writes every
  0.01 s (copies of reactor/showerhead_pulse_langmuir_p0.5).
Script reactor/sweep_sticking_sensitivity.py builds, runs (one scalarTransportFoam call per case through the OpenFOAM
launcher), post-processes and plots; it skips any case whose final time directory exists. Cases live under
reactor/sticking/. Post-processing: t90 = first write time at which the coverage reaches 0.9 (resolution one write
interval), at the centre face and at the worst face over faces 0-148; hole bottoms and flat surface for the patterned
wafer; supply-limit time (1 + 55.9) N0/(U c0) = 11.99 s (patterned) and N0/(U c0) = 0.2107 s (blanket);
tau_sat = t90_worst / supply-limit time; utilization up to saturation = 1/tau_sat (site-equivalents reacted per
site-equivalent delivered, as in reactor/doe_patterned.py and coupling/validate_supply_rule.py).
Literature on the TMA sticking probability was checked from the full texts of Arts et al. 2019 (JVST A) and
Ylilammi et al. 2018 (JAP) and the DOIs of the other sources were resolved through Crossref (section "References found").

## Results
All runs saturated well inside the pulse (final coverage 1.000 at every face in every case). Times in seconds.

| s | blanket t90 centre / worst | patterned flat t90 centre (worst) | patterned bottom t90 centre / worst | bottom t90 / supply time 11.99 s (centre, worst) | utilization patterned / blanket |
|---|---|---|---|---|---|
| 5e-4 | 0.92 / 1.01 | 6.00 (6.05) | 14.15 / 14.20 | 1.180, 1.184 | 84.4 % / 20.9 % |
| 1e-3 | 0.64 / 0.71 | 4.25 (4.30) | 14.05 / 14.10 | 1.172, 1.176 | 85.0 % / 29.7 % |
| 2e-3 | 0.47 / 0.53 | 3.00 (3.05) | 14.00 / 14.05 | 1.168, 1.172 | 85.3 % / 39.7 % |
| 5e-3 (reference) | 0.35 / 0.39 | 1.95 (2.00) | 14.05 / 14.10 | 1.172, 1.176 | 85.0 % / 54.0 % |
| 8.5e-3 | 0.30 / 0.33 | 1.55 (1.60) | 14.05 / 14.10 | 1.172, 1.176 | 85.0 % / 63.8 % |

Blanket wafer relative to its own supply-limit time 0.2107 s (worst face): 4.79, 3.37, 2.52, 1.85, 1.57.
Supporting numbers from the same fields (figures/sticking_sensitivity.csv and the case directories):
- Over the 17-fold range of s the patterned hole-bottom t90 moves by 0.15 s centre and 0.15 s worst, i.e. three write
  intervals, and tau_sat stays at 1.17-1.18. The flat surface of the patterned wafer scales as about s^-0.5
  (6.00 to 1.55 s) and the blanket wafer as about s^-0.4 (worst 1.01 to 0.33 s), not as 1/s: as the sticking falls,
  the TMA density at the wafer rises (centre, 0.1 s into the pulse: 0.229 c0 at s = 5e-4 against 0.065 c0 at 8.5e-3;
  blanket cases), which partly compensates. At s = 5e-4 the kinetic uptake coefficient s vbar/4 = 0.051 m/s equals
  the inlet speed 0.050 m/s, so the surface sees about half the feed density and the first-order Langmuir time
  ln(10) N0/(s vbar/4 c0) = 0.47 s stretches to the observed 0.92 s.
- Inside the holes the ordering is reversed: at 10 s the centre hole bottoms are at coverage 0.384 (s = 5e-4),
  0.262 (1e-3), 0.148 (2e-3), 0.045 (5e-3), 0.017 (8.5e-3). Low sticking lets TMA reach the bottom earlier (the
  classic conformality result), but the last 10 % of the bottom coverage waits for the supply in every case, so
  the finish time is the same. The flat-to-bottom lag shrinks from 12.5 s (8.5e-3) to 8.15 s (5e-4).
- The surface density over the patterned wafer at 10 s is 0.154-0.162 c0 for all five values of s: once the
  wafer is the sink for essentially the whole feed, the gas phase does not know the sticking coefficient.
- Blanket 5e-3 reproduces the reference numbers (0.345 / 0.385 s with 0.005 s writes) within one write interval.
- The blanket time at s = 5e-4 (1.01 s) is far below the 3.5 s that a 1/s extrapolation from 5e-3 would give; the
  8 s pulse covered it with a factor 8 of margin, and 2 s would have sufficed.
Wall times: blanket 84, 84, 94, 90, 96 s (9 s of simulated time, 900 writes each); patterned 206, 202, 200, 197 s
(26 s simulated, 520 writes each); 21 min in total for the nine new cases, run sequentially. The reference long case
took 171 s when it was run.

## Figure
figures/sticking_sensitivity.png (figures/sticking_sensitivity.csv holds the numbers).
Caption: Sensitivity of the saturating TMA pulse to the sticking coefficient at the baseline operating point (1 % TMA,
300 sccm, 300 Pa, 573 K, 15 mm gap). Left: time to 90 % coverage against s for the blanket wafer (worst face), the
flat surface between the holes of the AR-100 patterned wafer, and the hole bottoms of the patterned wafer (centre and
worst face); dash-dotted horizontal lines are the supply-limit times (1 + 55.9) N0/(U c0) = 11.99 s and
N0/(U c0) = 0.211 s. Right: the same times divided by their supply-limit time; the blue band is the 1.05-1.27 range of
the twelve patterned cases of the pulse-time map at s = 5e-3. In both panels the green band is the Arts et al. 2019
range (0.5-2e-3, lateral-high-aspect-ratio front slope, set-point 275 C), the vertical lines are the Ylilammi et al.
2018 fit (5.72e-3, 300 C), the Vandalon and Kessels 2016 sum-frequency-generation value (4e-3, 300 C), the PillarHall
fit of this work (8.5e-3, 300 C) and the value used in the chamber runs (5e-3).

## Literature and reconciliation
What the sources report for TMA on Al2O3 (initial sticking probability s0 unless stated):
- Arts, Vandalon, Puurunen, Utriainen, Gao, Kessels, Knoops, J. Vac. Sci. Technol. A 37, 030908 (2019): s0 =
  (0.5-2) x 10^-3 for Al(CH3)3 and (0.8-2) x 10^-4 for H2O at a set-point temperature of 275 C (estimated substrate
  temperature about 220 C), from the slope of the thickness-profile front at the half-coverage point in lateral
  high-aspect-ratio (LHAR) channels through their relation dtheta/dz~ at theta = 1/2 = sqrt(s0/13.9), with z~ the
  depth over the cavity height; the relation was computed with the Yanguas-Gil and Elam 2012 continuum model, in
  which the front slope is independent of the dose while the penetration depth grows with it. The stated limits come
  from assuming under 10 % variation in reactant dosing and about 10 nm uncertainty in cavity height. TMA doses were
  190 and 1200 mTorr s. Their Table I also lists: Ylilammi et al. 2018 data at 300 C read with the same relation,
  (2-7) x 10^-3; Ylilammi's own fit, 5.72 x 10^-3 at 300 C; Vandalon and Kessels 2016 (sum-frequency generation,
  300 C at the substrate), (4 +/- 1) x 10^-3 for TMA and (4 +/- 1) x 10^-4 for H2O. They report the H2O value to be
  temperature dependent ((1.5-2.3) x 10^-5 at 150 C set-point) and give TMA at the 275 C set-point only. They also
  note that natural pulse-to-pulse variation in the dose broadens the front and lowers the apparent sticking
  probability, most for sharp fronts.
- Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018): Picosun R-150, TMA/H2O at 300 C, 0.1 s pulse and
  4.0 s purge, 500 cycles, 150 sccm N2 per line, chamber about 300 Pa, PillarHall third-generation prototypes with a
  500 nm gap. Their Table I (parameters extracted from the measured profile): lumped sticking coefficient
  c = 0.00572, adsorption equilibrium constant K = 219 1/Pa, input partial pressure pA0 = 147 Pa, gpc at saturation
  105.6 pm, rms thickness error 0.974 nm. The partial pressure is a fitted quantity in their procedure, not a
  measured one; c is explicitly "a lumped sticking coefficient which combines the effects of all sequential and
  parallel surface processes" in a Langmuir model with desorption. The saturation adsorption density is computed from
  the film density and the saturated growth per cycle (their Eq. 9), not assumed.
- Arts et al., J. Phys. Chem. C 123, 27030 (2019) and J. Phys. Chem. C 125, 8244 (2021) treat recombination
  probabilities of O atoms in plasma ALD of SiO2, TiO2, Al2O3 and HfO2 from the same LHAR profile method; they do not
  report a TMA sticking probability. The TMA number attributed to "Arts et al." in this repository is the JVST A value.
- Rose and Bartha, Appl. Surf. Sci. 255, 6620 (2009) and Rose, Bartha, Endler, Appl. Surf. Sci. 256, 3778 (2010):
  the method papers for extracting sticking coefficients from profiles in nanoscale holes under free molecular flow;
  the 2009 abstract reports the method on tetrakis(dimethylamino)titanium with ozone, not TMA. Cited here for the
  method only.
- Yanguas-Gil and Elam, Chem. Vap. Deposition 18, 46 (2012): the reaction-diffusion model and its analytic saturation
  times; the Arts relation is derived with it. No TMA value is claimed from it here. The Cremers, Puurunen, Dendooven
  review (Appl. Phys. Rev. 6, 021302 (2019)) was not readable through the tools available in this session (the
  publisher and repository pages refused the request), so no number is taken from it.
- Yim et al., Phys. Chem. Chem. Phys. 22, 23107 (2020) is the source of the Zenodo 3974438 profiles used by
  coupling/fit_pillarhall.py (TMA/H2O at 300 C, 0.1 s TMA pulses, about 300 Pa, 500 nm LHAR gap per the workbook
  metadata); the full text was not readable in this session, so whether it states a TMA partial pressure is not
  confirmed here.

Reconciliation. The PillarHall fit of this work (s = 8.5e-3, sheet-to-sheet factor 1.5; mouth TMA pressure 75 Pa,
exposure 7.5 Pa s per 0.1 s pulse; rms 0.038 of the plateau) fits two parameters to one profile with the site
density fixed at N0 = 4e18 /m2 and the TMA partial pressure unknown. The two parameters are constrained by different
parts of the profile: the slope of the front at half thickness is set by the sticking coefficient (the Arts relation
makes this explicit, and it is why the s = 1e-3 and 1e-2 curves at the same exposure bracket the data visibly in
figures/pillarhall_fit.png), while the depth at which the front sits is set by the exposure, which enters only as
the ratio of delivered molecules to sites, n0 t_pulse / N0. A different N0 therefore moves the fitted pressure, not
the fitted s, to first order; N0 enters the front shape only through the channel-narrowing term over 500 cycles,
which is weak. The fit therefore gives s with the sheet-to-sheet scatter quoted and the exposure only as an
exposure-per-site estimate (the 75 Pa is 4e18 x 7.5 Pa s / 4e18, not a pressure measurement). Three things separate
8.5e-3 from the Arts range of 0.5-2e-3 and none of them is the chamber model: (i) temperature and data set, since
Arts report TMA at a 275 C set-point (about 220 C at the substrate) whereas the PillarHall profiles were grown at
300 C, and the 300 C values in Arts' own Table I are 2-7e-3 (Ylilammi data by the slope relation), 5.72e-3
(Ylilammi fit) and 4 +/- 1e-3 (sum-frequency generation), so the fit sits a factor 1.5 above the Ylilammi fit and 2
above the SFG value, within or at the edge of its own factor-1.5 scatter; (ii) the kinetic form, since the fit
uses first-order irreversible Langmuir uptake s (1 - theta) on 4e18 sites, Arts use an initial sticking probability
toward the bare surface in the Yanguas-Gil and Elam model, and Ylilammi use a lumped coefficient with an equilibrium
constant, so the three numbers are not the same quantity; (iii) the dose history, since cycle-to-cycle dose variation
and averaging eight profiles broaden the front, which Arts show biases a slope-derived value low rather than high.
The sensitivity runs show what the uncertainty costs: over 5e-4 to 8.5e-3, which covers every value above, the
AR-100 patterned wafer saturates in 14.00-14.20 s, 1.17-1.18 times the supply-limit time, with 84-85 % utilization,
so the design rule, the recipe-transfer factor and the pulse-time map do not depend on which literature value is
right; the blanket wafer moves from 0.33 to 1.01 s (1.6 to 4.8 times its supply limit) and the flat surface between
the holes from 1.55 to 6.00 s, so every blanket-wafer number in the paper carries the sticking uncertainty, about a
factor 3 across the published range. What would pin s down is a saturation profile grown with a known TMA partial
pressure at the channel mouth (or a measured dose, as Arts state theirs in mTorr s) at 300 C, together with the site
density taken from the measured growth per cycle and film density as Ylilammi do, so that the exposure is not a free
parameter and the fitted s is the only unknown; a sum-frequency-generation or QMS measurement at the same temperature
would give the independent check.

## Proposed paper text
Discussion (two paragraphs, ready to paste):

The TMA sticking coefficient is the least certain input of the chamber model: the published values at or near this
temperature span 0.5-2e-3 (Arts et al. 2019, lateral-high-aspect-ratio front slope at a 275 C set-point), 4 +/- 1e-3
(Vandalon and Kessels 2016, sum-frequency generation at 300 C), 5.72e-3 (Ylilammi et al. 2018, fit to a lateral
channel profile at 300 C) and 8.5e-3 (this work, fit to the Zenodo 3974438 profiles at 300 C, factor 1.5
sheet-to-sheet). We reran the baseline operating point (1 % TMA, 300 sccm, 300 Pa, 573 K) at s = 5e-4, 1e-3, 2e-3,
5e-3 and 8.5e-3 for the blanket wafer and the AR-100 patterned wafer. On the patterned wafer the hole bottoms reach
90 % coverage at 14.00-14.15 s (centre) and 14.05-14.20 s (worst face) for every value, 1.17-1.18 times the
supply-limit time of 11.99 s, with 84-85 % precursor utilization: a seventeen-fold change in s moves the saturating
pulse by three write intervals. Low sticking changes the order of events inside the hole (at 10 s the bottom coverage
is 0.38 at s = 5e-4 and 0.02 at 8.5e-3) but not the finish, because the last tenth of the bottom coverage waits for
molecules that the showerhead has not yet delivered. The design rule, the blanket-to-product transfer factor and the
pulse-time map are therefore independent of which published value is right.

The blanket wafer is not. Its saturation time runs from 0.33 s at s = 8.5e-3 to 1.01 s at 5e-4, about s^-0.4 rather
than 1/s, because the TMA density at the wafer rises as the uptake slows (0.065 c0 against 0.23 c0 at the centre
0.1 s into the pulse) and partly compensates; the flat surface between the holes of the patterned wafer behaves the
same way (1.55 to 6.00 s). Every blanket-wafer time in this paper therefore carries a factor-3 uncertainty from the
sticking coefficient, whereas the patterned-wafer times carry only the site-density uncertainty of the supply limit.
The 8.5e-3 of our PillarHall fit is a slope-constrained value: the front slope fixes s and the penetration depth fixes
the exposure per site, so the unknown TMA partial pressure trades against the assumed site density, not against s.
It sits a factor 1.5 to 2 above the two other 300 C values and a factor 4 above the upper bound of the 275 C
set-point data, a spread that also reflects three different kinetic definitions (initial sticking probability,
lumped Langmuir coefficient with an equilibrium constant, irreversible first-order Langmuir on 4e18 sites). A
profile grown at a known TMA partial pressure at 300 C, with the site density taken from the measured growth per
cycle, would remove the trade-off and fix s to the sheet-to-sheet scatter.

Abstract (one sentence, if wanted after "Chamber pressure, showerhead gap and carrier flow do not move that number;
precursor delivery rate does."): Nor does the TMA sticking coefficient: over the published range, 5e-4 to 8.5e-3,
the product-wafer pulse stays at 1.17-1.18 times the supply limit while the blanket-wafer pulse changes threefold.

## CLAIMS rows
| Sticking sensitivity, patterned AR-100 wafer | bottom t90 14.00-14.15 s centre, 14.05-14.20 s worst, for s = 5e-4 to 8.5e-3; 1.17-1.18 x the 11.99 s supply limit; utilization 84-85 % | figures/sticking_sensitivity.csv, figures/sticking_sensitivity.png | PYTHONPATH=src .venv/bin/python reactor/sweep_sticking_sensitivity.py (cases under reactor/sticking/) |
| Sticking sensitivity, blanket wafer | t90 0.92/1.01 s (s = 5e-4), 0.64/0.71, 0.47/0.53, 0.35/0.39, 0.30/0.33 s (8.5e-3), centre/worst; 1.6-4.8 x the 0.211 s supply limit, about s^-0.4 | same | same |
| Sticking sensitivity, flat surface of the patterned wafer | t90 6.00, 4.25, 3.00, 1.95, 1.55 s for s = 5e-4 to 8.5e-3, about s^-0.5 | same | same |
| Published TMA sticking probabilities on Al2O3 | 0.5-2e-3 at 275 C set-point (Arts 2019); 4 +/- 1e-3 at 300 C (Vandalon and Kessels 2016, in Arts 2019 Table I); 5.72e-3 at 300 C (Ylilammi 2018 Table I); PillarHall fit 8.5e-3 at 300 C (this work) | paper/wp/wp1_sticking.md, figures/pillarhall_fit.png | literature: doi 10.1116/1.5093620, 10.1063/1.4939654, 10.1063/1.5028178; coupling/fit_pillarhall.py |

## Run log (for paper/assumptions.md)
## Night 3, part 2: sticking-coefficient sensitivity (reactor/sweep_sticking_sensitivity.py; figures/sticking_sensitivity.png, .csv)
Nine new runs at the baseline point (300 Pa, 573 K, 300 sccm, 1 % TMA, 15 mm gap, N0 4e18, c0 3.79e20, frozen flow)
with s = 5e-4, 1e-3, 2e-3, 8.5e-3 on the AR-100 patterned wafer (24 s pulse, endTime 26 s, writes 0.05 s; the 20 s
long case is the 5e-3 point) and s = 5e-4, 1e-3, 2e-3, 5e-3, 8.5e-3 on the blanket wafer (8 s pulse, endTime 9 s,
writes 0.01 s). Cases are copies of showerhead_pulse_patterned_AR100_long and showerhead_pulse_langmuir_p0.5 with the
wafer BC line, pulse table and endTime/writeInterval substituted (each regex asserted to hit once), under
reactor/sticking/. Wall time 84-96 s per blanket case, 197-206 s per patterned case, 21 min in all.
| s | blanket t90 centre / worst | patterned flat t90 | patterned bottom t90 centre / worst | bottom t90 / 11.99 s (worst) | utilization patterned / blanket |
| 5e-4 | 0.92 / 1.01 s | 6.00 s | 14.15 / 14.20 s | 1.184 | 84.4 % / 20.9 % |
| 1e-3 | 0.64 / 0.71 | 4.25 | 14.05 / 14.10 | 1.176 | 85.0 / 29.7 |
| 2e-3 | 0.47 / 0.53 | 3.00 | 14.00 / 14.05 | 1.172 | 85.3 / 39.7 |
| 5e-3 | 0.35 / 0.39 | 1.95 | 14.05 / 14.10 | 1.176 | 85.0 / 54.0 |
| 8.5e-3 | 0.30 / 0.33 | 1.55 | 14.05 / 14.10 | 1.176 | 85.0 / 63.8 |
Reading:
- The patterned-wafer pulse is independent of s across the whole published range: 0.15 s (three write intervals)
  of movement for a seventeen-fold change, tau_sat 1.17-1.18 throughout, utilization 84-85 %. The surface density over
  the wafer at 10 s is 0.154-0.162 c0 for all five values: once the wafer takes essentially the whole feed, the gas
  phase does not know the sticking coefficient. Inside the holes the bottoms fill earlier at low s (coverage at 10 s:
  0.384 at 5e-4, 0.017 at 8.5e-3) and the flat-to-bottom lag shrinks from 12.5 to 8.15 s, but the last tenth waits for
  the supply in every case.
- The blanket wafer scales as about s^-0.4 (worst 1.01 to 0.33 s) and the flat surface of the patterned wafer as
  about s^-0.5 (6.00 to 1.55 s), weaker than 1/s because the TMA density at the wafer rises as the uptake slows
  (centre, 0.1 s: 0.229 c0 at 5e-4, 0.065 c0 at 8.5e-3). At s = 5e-4 the uptake coefficient s vbar/4 = 0.051 m/s
  equals the inlet speed, so the surface sees about half the feed and the first-order time ln(10) N0/(s vbar/4 c0)
  = 0.47 s becomes 0.92 s. The 1/s extrapolation (3.5 s at 5e-4) was wrong by a factor 3.5; the 8 s pulse was
  longer than needed by a factor 8.
- Blanket 5e-3 with 0.01 s writes gives 0.35 / 0.39 s against 0.345 / 0.385 s with 0.005 s writes: one write interval.
- Literature check (full texts of Arts 2019 JVST A and Ylilammi 2018 JAP): Arts give TMA s0 = (0.5-2)e-3 at a 275 C
  set-point (about 220 C substrate) from the front slope, and tabulate the 300 C values 2-7e-3 (Ylilammi data, slope
  relation), 5.72e-3 (Ylilammi fit; their fitted partial pressure 147 Pa, K = 219 1/Pa) and 4 +/- 1e-3 (Vandalon and
  Kessels 2016, SFG). Our PillarHall fit (8.5e-3 x/÷ 1.5 at 300 C) is a factor 1.5-2 above the 300 C values, not
  4x above the field; the front slope fixes s and the penetration depth fixes exposure per site, so the unknown TMA
  pressure trades against N0, not against s. Arts' J. Phys. Chem. C papers (2019, 2021) are on O-atom recombination
  and give no TMA value.
Caveats: single operating point (the DOE showed tau_sat 1.03-1.25 at p >= 316 Pa and up to 1.96 at 100 Pa with
600 sccm, where carrier bypass, not sticking, is the cause); s held constant with coverage (first-order Langmuir);
t90 resolution 0.05 s patterned, 0.01 s blanket; the Cremers 2019 review and Yim 2020 full texts were not readable
through the tools available, so no numbers were taken from them.

## References found
- K. Arts, V. Vandalon, R. L. Puurunen, M. Utriainen, F. Gao, W. M. M. Kessels, H. C. M. Knoops, "Sticking
  probabilities of H2O and Al(CH3)3 during atomic layer deposition of Al2O3 extracted from their impact on film
  conformality", J. Vac. Sci. Technol. A 37, 030908 (2019). doi 10.1116/1.5093620. Full text read (TU/e repository
  copy): TMA (0.5-2) x 10^-3 and H2O (0.8-2) x 10^-4 at Tset = 275 C; Table I as quoted above.
- M. Ylilammi, O. M. E. Ylivaara, R. L. Puurunen, "Modeling growth kinetics of thin films made by atomic layer
  deposition in lateral high-aspect-ratio structures", J. Appl. Phys. 123, 205301 (2018). doi 10.1063/1.5028178.
  Full text read (Aalto repository copy): Table I, Al2O3: c = 0.00572, K = 219 1/Pa, pA0 = 147 Pa, 300 C.
- V. Vandalon, W. M. M. Kessels, "What is limiting low-temperature atomic layer deposition of Al2O3? A vibrational
  sum-frequency generation study", Appl. Phys. Lett. 108, 011607 (2016). doi 10.1063/1.4939654. TMA (4 +/- 1) x 10^-3
  at 300 C as tabulated by Arts et al. 2019 (Table I); the APL text itself was not read.
- K. Arts, M. Utriainen, R. L. Puurunen, W. M. M. Kessels, H. C. M. Knoops, "Film conformality and extracted
  recombination probabilities of O atoms during plasma-assisted atomic layer deposition of SiO2, TiO2, Al2O3, and
  HfO2", J. Phys. Chem. C 123, 27030 (2019). doi 10.1021/acs.jpcc.9b08176 (correction: J. Phys. Chem. C 124, 1250,
  doi 10.1021/acs.jpcc.9b11082). O-atom recombination, no TMA value.
- K. Arts et al., "Oxygen recombination probability data for plasma-assisted atomic layer deposition of SiO2 and
  TiO2", J. Phys. Chem. C 125, 8244 (2021). doi 10.1021/acs.jpcc.1c01505. O-atom recombination, no TMA value.
- A. Yanguas-Gil, J. W. Elam, "Self-limited reaction-diffusion in nanostructured substrates: surface coverage
  dynamics and analytic approximations to ALD saturation times", Chem. Vap. Deposition 18, 46 (2012).
  doi 10.1002/cvde.201106938. Model used by Arts for the slope relation; no TMA value taken from it.
- M. Rose, J. W. Bartha, "Method to determine the sticking coefficient of precursor molecules in atomic layer
  deposition", Appl. Surf. Sci. 255, 6620 (2009). doi 10.1016/j.apsusc.2009.02.055. Method paper (TDMAT/ozone per
  the abstract).
- M. Rose, J. W. Bartha, I. Endler, "Temperature dependence of the sticking coefficient in atomic layer deposition",
  Appl. Surf. Sci. 256, 3778 (2010). doi 10.1016/j.apsusc.2010.01.025. Not read; cited by Arts 2019 as ref. 13.
- V. Cremers, R. L. Puurunen, J. Dendooven, "Conformality in atomic layer deposition: Current status overview of
  analysis and modelling", Appl. Phys. Rev. 6, 021302 (2019). doi 10.1063/1.5060967. Not readable in this session.
- J. Yim, O. M. E. Ylivaara, M. Ylilammi, V. Korpelainen, E. Haimi, E. Verkama, M. Utriainen, R. L. Puurunen,
  "Saturation profile based conformality analysis for atomic layer deposition: aluminum oxide in lateral
  high-aspect-ratio channels", Phys. Chem. Chem. Phys. 22, 23107 (2020). doi 10.1039/d0cp03358h. Source of the
  Zenodo 3974438 profiles; not readable in this session.
- "Molecular dynamics of atomic layer deposition: sticking coefficient investigation", Appl. Sci. 12, 2188 (2022).
  doi 10.3390/app12042188. A search snippet attributes to it a TMA initial sticking coefficient of 4.257e-3 at 573 K
  from a thickness-profile slope; the text was not readable in this session, so the number is not used above.
