# Work package 2: fab units (step coverage, thickness, purge, cycle time, wafers per hour)

2026-11-10. Script: `coupling/fab_units.py` (re-runnable; `PYTHONPATH=src .venv/bin/python coupling/fab_units.py`).
No OpenFOAM case was run; every number comes from time directories already on disk, from the Python feature
solver's transport module, or from a cited source.

## What was done

1. Step coverage and thickness. From the 1 % TMA map runs at AR 20, 50 and 100 (`reactor/map/patterned_AR*_x0.01`,
   pulses 4.61, 10.26 and 19.68 s, writes every 0.023, 0.051 and 0.098 s) and the 20 s AR-100 run
   (`reactor/showerhead_pulse_patterned_AR100_long`, writes every 0.05 s), the bottom/top step coverage
   SC(t) = thetaBottom/thetaFlat at the wafer centre and at the worst face (faces 0 to 148, as in the map), and the
   Al2O3 thickness after N cycles of a pulse of length t, N x GPC x theta, with GPC = 0.103 nm per cycle from the
   PillarHall record in this repository (51.5 nm plateau after 500 cycles, `coupling/fit_pillarhall.py`,
   `figures/pillarhall_fit.png`). In the irreversible Langmuir model the coverage at time t of a long pulse equals
   the coverage a pulse of length t leaves behind, up to what the gas left in the chamber adds after the valve
   closes (surface density 0.019 c0 over the patterned wafer, so negligible); the long run's bottom t90 (14.05 s)
   is the same as the map's (14.01 s) and the DOE's, which is the check. Tabulated t90, t95 and t99 per aspect ratio.
2. Purge, three parts. (a) The hole: Knudsen clearing of a 100 nm x 10 um hole with the mouth density set to zero,
   analytically (L^2/D_K, D_K from `aldconf.transport.knudsen_diffusivity`, and the slowest eigenmode) and
   numerically (implicit 1D transient, 200 cells, 20 ns steps). (b) The chamber: the outlet patch signal after the
   pulse ends (the `outlet()` reader of `coupling/postprocess_desorption.py`, copied), time to 10 %, 1 % and 0.1 %
   of the pulse-end value for the 20 s patterned TMA run, the two blanket TMA runs and the water kdes = 0 run, plus
   the exponential decay rate of the tail (log-linear fit over the last 40 % of each purge window) used to
   extrapolate levels the written window does not reach; compared with the gap transit 0.015 m / 0.050 m/s = 0.3 s
   and V/Q = 0.53 s. (c) Desorption: purge to 1 % for kdes = 0, 0.1, 1, 10 1/s with the same metric as
   `coupling/postprocess_desorption.py`; where the 3 s purge window does not reach 1 % (0.1 and 1 1/s) the tail is
   extrapolated as the desorption-fed outlet level kdes theta N0 A_wafer/(Q c0), which reproduces the simulated
   outlet at 3 s to within 7 % (0.0159 vs 0.0167 c0 at 0.1 1/s; 0.0238 vs 0.0254 at 1 1/s), with theta decaying at
   the rate fitted to the simulated coverage over the last 2 s (0.088 and 0.654 1/s; below kdes because desorbed
   water re-adsorbs). `figures/purge_times.csv`.
3. Cycle time and throughput. Full Al2O3 cycle (TMA, purge, H2O, purge) for the blanket wafer (0.345 s, 1.74 s), the
   patterned AR-100 wafer (14.05 s, 13.4 s) and the same-supply optimum from `figures/doe_gp_report.txt` (TMA
   9.83 s; H2O scaled by the same factor 9.83/14.05 to 9.38 s, because the water half-cycle was not in the DOE and
   on the patterned wafer both half-cycles are supply-limited, Section 3.5). Purges: TMA, chamber outlet to 1 % from
   the 20 s patterned run (1.90 s); H2O, the kdes = 0.1 1/s value (8.26 s, extrapolated as above). kdes is a
   parameter, not a measurement, so every cycle is also given with gas clearing only (kdes = 0, 2.20 s). Cycles for
   5 and 10 nm at 0.103 nm per cycle (49 and 98), minutes per wafer, wafers per hour (WPH) with no overhead and with
   60 s of wafer handling per wafer. The rule as a throughput formula, WPH against wall-area ratio for 1 % and 3 %
   TMA with the constants of `coupling/temporal_vs_spatial.py` (N0 4e18, U 0.050 m/s, c0 3.79e20 at 1 %, both
   half-cycles at the same supply flux) and the two purges above. `figures/cycle_time_wph.csv`,
   `figures/wph_vs_wall_ratio.csv`.

## Results

### Table 1a. Pulse for 90, 95 and 99 % bottom coverage, 1 % TMA, 300 sccm, 300 Pa, 15 mm (step coverage = thetaBottom/thetaFlat)

| wafer (run) | wall-area ratio | rule 1.2 (1 + A_wall/A_wafer) N0/(U c0) | t90 centre / worst | t95 centre / worst | t99 centre / worst | flat top t90 / t95 / t99 | t99 / rule |
|---|---|---|---|---|---|---|---|
| AR 20 (`map/patterned_AR20_x0.01`, 4.61 s pulse) | 11.2 | 3.08 s | 2.67 / 2.71 s | 2.78 / 2.81 s | 2.92 / 2.94 s | 1.84 / 2.14 / 2.60 s | 0.96 |
| AR 50 (`map/patterned_AR50_x0.01`, 10.26 s) | 28.0 | 7.32 s | 6.73 / 6.73 s | 6.83 / 6.83 s | 6.94 / 6.99 s | 1.94 / 2.35 / 3.06 s | 0.95 |
| AR 100 (`map/patterned_AR100_x0.01`, 19.68 s) | 55.9 | 14.38 s | 14.01 / 14.11 s | 14.21 / 14.21 s | 14.31 / 14.31 s | 1.96 / 2.35 / 3.04 s | 1.00 |
| AR 100 (`showerhead_pulse_patterned_AR100_long`, 20 s, 0.05 s writes) | 55.9 | 14.38 s | 14.05 / 14.10 s | 14.15 / 14.20 s | 14.25 / 14.30 s | 1.95 / 2.30 / 3.05 s | 0.99 |

The step-coverage and bottom-coverage times coincide to one write interval once the flat top has saturated
(1.84 to 1.96 s), which is long before any bottom. The saturation front is sharp: t99 - t90 is 0.25 s at AR 20,
0.26 s at AR 50 and 0.25 to 0.30 s at AR 100 (2 % of t90), so the rule, which is 2 to 12 % above t90, also covers
99 % step coverage at every aspect ratio (t99/rule 0.95 to 1.00). The worst face is within 0.1 s of the centre.
Time resolution: one write interval (0.023 to 0.098 s).

### Table 1b. Al2O3 after 50 cycles of a given pulse (0.103 nm per cycle; saturated film 5.15 nm), wafer centre

| pulse | AR 20: step coverage, bottom / top | AR 50: step coverage, bottom / top | AR 100: step coverage, bottom / top |
|---|---|---|---|
| 0.345 s (blanket t90) | 20 %, 0.26 / 1.29 nm | 0.5 %, 0.01 / 1.32 nm | 0 %, 0.00 / 1.44 nm |
| 1 s | 29 %, 0.95 / 3.23 nm | 0.8 %, 0.03 / 3.25 nm | 0 %, 0.00 / 3.16 nm |
| 2 s | 58 %, 2.76 / 4.79 nm | 1.7 %, 0.08 / 4.69 nm | 0 %, 0.00 / 4.67 nm |
| 5 s | 100 %, 5.15 / 5.15 nm | 25 %, 1.27 / 5.15 nm | 0.1 %, 0.00 / 5.15 nm |
| 10 s | 100 %, 5.15 / 5.15 nm | 100 %, 5.15 / 5.15 nm | 4.5 %, 0.23 / 5.15 nm |
| rule (3.08 / 7.32 / 14.38 s) | 100 % | 100 % | 100 % |

The blanket recipe (0.345 s) on an AR-100 wafer leaves 1.44 nm of the 5.15 nm target on the array top (the flat
top is itself at 28 % coverage, because it competes with the holes) and nothing at the bottom; on an AR-20 wafer
it leaves a 20 % step coverage. Full curves: `figures/step_coverage_vs_pulse.csv`.

### Table 2. Purge times (`figures/purge_times.csv`)

(a) The hole, 100 nm x 10 um, TMA at 573 K (vbar 410 m/s, D_K = d vbar/3 = 1.37e-5 m2/s), walls saturated, mouth
density zero:

| quantity | value |
|---|---|
| L^2/D_K | 7.3 us |
| slowest eigenmode 4 L^2/(pi^2 D_K) (Dirichlet mouth, Neumann bottom) | 3.0 us |
| bottom density to 1 %, first mode (4/pi) exp(-t/tau1) | 14.4 us |
| numeric: bottom density to 1/e | 3.7 us |
| numeric: bottom density to 1 % | 14.4 us |
| numeric: hole-mean density to 1 % | 13.1 us |

Analytic and numeric agree (L^2/D_K against the numeric 1/e time: factor 2.0; the eigenmode estimate of the 1 % time
against the numeric: 0.4 %). The hole empties in 15 us, 1.3e5 times faster than the chamber reaches 1 %; the hole
does not limit the purge. (A partly saturated hole empties faster still, because the walls still consume.)

(b) The chamber, outlet signal after the pulse ends (fraction of the pulse-end outlet value):

| run | pulse | outlet at pulse end | to 10 % | to 1 % | to 0.1 % | tail decay rate (time constant) |
|---|---|---|---|---|---|---|
| TMA, patterned AR 100, 20 s pulse, purge to 22 s | 20 s | 1.000 c0 | 1.10 s | 1.90 s | 2.73 s (extrapolated) | 2.75 1/s (0.36 s) |
| TMA, blanket, 0.5 s pulse, purge to 1.5 s | 0.5 s | 0.297 c0 (peak 0.367 at +0.11 s) | 1.25 s (extr.) | 2.14 s (extr.) | 3.04 s (extr.) | 2.57 1/s (0.39 s) |
| TMA, blanket, 1.26 s pulse, purge to 2.26 s (`map/blanket_x0.01`) | 1.26 s | 0.888 c0 | 1.07 s (extr.) | 1.95 s (extr.) | 2.83 s (extr.) | 2.63 1/s (0.38 s) |
| H2O, blanket, kdes 0, 5 s pulse, purge to 8 s | 5 s | 0.999 c0 | 1.14 s | 2.20 s | 3.25 s (extr.) | 2.18 1/s (0.46 s) |
| reference: gap transit 0.015 m / 0.050 m/s | | | | 0.30 s | | |
| reference: chamber volume / flow, pi 0.2^2 x 0.015 m3 / 3.54e-3 m3/s | | | | 0.53 s | | |

"extr." = extrapolated with the exponential tail fitted over the last 40 % of the written purge window; the blanket
0.5 s run has only a 1 s purge window and its outlet is still rising at pulse end (the gas has not crossed the
chamber), so its 10 % and 1 % times are extrapolations too. The outlet tail time constant, 0.36 to 0.39 s for TMA
and 0.46 s for water, lies between the gap transit (0.30 s) and V/Q (0.53 s); a decade of signal costs 0.8 to 1.1 s
of purge, and 1 % takes 1.9 s for TMA and 2.2 s for water (about 4 V/Q). The wafer flux clears much sooner
(below 1 % 0.05 s after the pulse, part 1 of the run log); the outlet is downstream of the whole chamber volume.

(c) Desorption, water on the blanket wafer, 5 s pulse at 1 %, outlet to 1 % of the pulse-end value:

| kdes | coverage at pulse end / 3 s after | outlet 3 s after the pulse | purge to 1 % | basis |
|---|---|---|---|---|
| 0 | 1.000 / 1.000 | 0.0017 c0 | 2.20 s | run (gas clearing) |
| 0.1 1/s | 0.951 / 0.755 | 0.0167 c0 (tail model 0.0159) | 8.26 s | extrapolated: tail kdes theta N0 A/(Q c0), theta decaying at 0.088 1/s |
| 1 1/s | 0.660 / 0.113 | 0.0254 c0 (tail model 0.0238) | 4.33 s | extrapolated, theta decaying at 0.654 1/s |
| 10 1/s | 0.163 / 0.001 | 0.0030 c0 | 2.40 s | run (surface empties within 0.3 s; gas clearing again) |

### Table 3a. Al2O3 cycle time and throughput (`figures/cycle_time_wph.csv`; GPC 0.103 nm per cycle; 49 cycles for 5 nm, 98 for 10 nm)

Purges: TMA 1.90 s (chamber outlet to 1 %, 20 s patterned run); H2O 8.26 s (kdes = 0.1 1/s, a stated choice, kdes
is not measured). In brackets: the same cycle with gas clearing only for water (kdes = 0, 2.20 s).

| recipe | TMA + purge + H2O + purge | cycle | 5 nm: min per wafer, WPH, WPH with 60 s handling | 10 nm: min per wafer, WPH, WPH with 60 s handling |
|---|---|---|---|---|
| blanket wafer, baseline (1 %, 300 sccm, 300 Pa, 15 mm) | 0.345 + 1.90 + 1.74 + 8.26 s | 12.2 s (6.2 s) | 10.0 min, 6.0, 5.5 (5.1 min, 11.9, 9.9) | 20.0 min, 3.0, 2.9 (10.1 min, 5.9, 5.4) |
| patterned AR 100, baseline | 14.05 + 1.90 + 13.4 + 8.26 s | 37.6 s (31.6 s) | 30.7 min, 1.95, 1.9 (25.8 min, 2.3, 2.2) | 61.4 min, 0.98, 0.96 (51.5 min, 1.2, 1.1) |
| patterned AR 100, same-supply optimum (2.51 %, 150 sccm, 631 Pa, 10 mm) | 9.83 + 1.90 + 9.38 + 8.26 s | 29.4 s (23.3 s) | 24.0 min, 2.5, 2.4 (19.0 min, 3.2, 3.0) | 48.0 min, 1.25, 1.2 (38.1 min, 1.6, 1.5) |

The patterned wafer costs 3.1 times the blanket cycle (37.6 against 12.2 s) and, at 5 nm, 1.9 wafers per hour
against 5.5; the pulses are 73 % of its cycle and 17 % of the blanket's. The same-supply optimum takes 22 % off the
patterned cycle (8.2 s of pulses). Handling at 60 s per wafer costs 3 to 9 % of throughput because a wafer takes
10 to 61 min. The purges for the optimum are kept at the baseline values because the purge was not in the DOE. This flatters
the optimum: it trades carrier gas for pressure, so at 150 sccm, 631 Pa and 10 mm the volumetric flow is
8.4e-4 m3/s (against 3.54e-3), the showerhead speed 0.012 m/s, the gap transit 0.84 s and V/Q 1.5 s, both 2.8
times the baseline. If the TMA purge scales with those (the baseline tail time constant sits between them) it
would be about 5.3 s instead of 1.9 s, the optimum's cycle about 32.8 s instead of 29.4 s, and its gain over the
baseline 4.8 s rather than 8.2 s. The water purge at kdes 0.1 is desorption-fed and scales with the precursor
supply, which the optimum holds within 25 % of the baseline, so it is left unchanged. This is a scaling estimate,
not a run.

### Table 3b. The rule as throughput: WPH for 5 nm (49 cycles) against wall-area ratio, both half-cycles at the same supply flux, purges 1.90 + 8.26 s (`figures/wph_vs_wall_ratio.csv`)

| A_wall/A_wafer | 1 % TMA: cycle, WPH, WPH with 60 s handling (kdes 0 purge) | 3 % TMA: cycle, WPH, WPH with 60 s handling (kdes 0 purge) |
|---|---|---|
| 1 | 11.2 s, 6.6, 5.9 (11.6) | 10.5 s, 7.0, 6.3 (13.0) |
| 11.2 (AR 20) | 16.3 s, 4.5, 4.2 (6.4) | 12.2 s, 6.0, 5.5 (10.0) |
| 28.0 (AR 50) | 25.1 s, 2.9, 2.8 (3.6) | 15.1 s, 4.9, 4.5 (7.1) |
| 55.9 (AR 100) | 39.3 s, 1.9, 1.8 (2.1) | 19.9 s, 3.7, 3.5 (4.9) |
| 112 | 67.6 s, 1.1, 1.1 (1.2) | 29.3 s, 2.5, 2.4 (3.0) |

The chamber run at AR 100 (37.6 s, 1.9 WPH) sits on the 1 % curve (39.3 s, 1.8 WPH; the rule is 2 % conservative).
Below a wall-area ratio of about 10 the purges set the cycle and precursor fraction buys little; above it the cycle
is 2 x 1.2 (1 + A_wall/A_wafer) N0/(U c0) and throughput falls as 1/(1 + A_wall/A_wafer): tripling the TMA fraction
doubles the AR-100 throughput (1.8 to 3.5 WPH) and does nothing for a blanket wafer (5.9 to 6.3).

## Figures

- `figures/step_coverage_vs_pulse.png`. Left: step coverage (hole bottom over flat top, wafer centre) against pulse
  time for AR 20, 50 and 100 at 1 % TMA, log x, with the 90, 95 and 99 % levels and the rule's pulse for each
  aspect ratio as dotted verticals (3.1, 7.3, 14.4 s); the 20 s AR-100 run is overlaid as the thin line. Right: Al2O3
  thickness at the hole bottom after 50 cycles of that pulse (0.103 nm per cycle), with the flat top of the array as
  the dotted black line and the saturated 5.15 nm marked.
- `figures/cycle_time_wph.png`. Left: the four segments of the Al2O3 cycle (TMA pulse, chamber purge to 1 %, H2O
  pulse, purge at kdes 0.1 1/s) stacked for the blanket baseline, the patterned AR-100 baseline and the same-supply
  optimum, with the cycle time and the 5 nm throughput (60 s handling) over each bar. Right: wafers per hour for
  5 nm against wall-area ratio from the rule at 1 % and 3 % TMA, with and without 60 s handling, and the AR-100
  chamber run as the point.

## Proposed paper text

### 3.10 Step coverage, purge and throughput

In fab units the rule reads as follows. Step coverage, the hole-bottom coverage over the flat top, follows the
supply integral: at 1 % TMA the bottoms of AR 20, 50 and 100 holes reach 90 % at 2.7, 6.7 and 14.1 s (worst face),
95 % at 2.8, 6.8 and 14.2 s and 99 % at 2.9, 7.0 and 14.3 s. The front is sharp, t99 - t90 is 0.25 to 0.3 s at
every aspect ratio, so the rule's pulse (3.1, 7.3, 14.4 s) covers 99 % step coverage as well as 90 % (t99/rule 0.95
to 1.00). The flat top between the holes saturates at 1.8 to 2.0 s regardless of aspect ratio. Converted to film at
0.103 nm per cycle (the PillarHall record's 51.5 nm after 500 cycles; the review value is about 0.1 nm per cycle
[Puurunen 2005]), 50 cycles of the blanket recipe (0.345 s) leave 1.4 nm of a 5.15 nm target on the top of an
AR-100 array and nothing at the bottom, and a 20 % step coverage on an AR-20 array; 50 cycles of the rule's pulse
leave 5.15 nm top and bottom on all three.

The purge has three parts and only one of them is set by the chamber. The hole itself empties by Knudsen
diffusion in microseconds: for 100 nm x 10 um and TMA at 573 K, D_K = 1.4e-5 m2/s, L^2/D_K = 7 us, and a transient
solution with the mouth density set to zero puts the bottom density at 1 % after 14 us. The chamber outlet signal
falls with a time constant of 0.36 to 0.39 s for TMA and 0.46 s for water, between the gap transit (0.3 s) and the
chamber volume over the flow (0.53 s): a decade costs about 1 s, 1 % of the pulse-end level takes 1.9 s (TMA) and
2.2 s (water), 0.1 % about 2.7 and 3.2 s. Desorbing water is the long part: at kdes = 0.1 1/s the exhaust is held
at 1.6 % of the feed level by the surface 3 s after the pulse and reaches 1 % only after about 8 s; at 1 1/s,
4.3 s; at 10 1/s the surface empties within 0.3 s and the purge is gas clearing again (2.4 s). kdes is a parameter,
not a measurement, and the water purge in the cycle below is the kdes = 0.1 1/s value by choice.

The full Al2O3 cycle, TMA pulse, purge to 1 %, water pulse, purge, is 12.2 s on the blanket wafer (0.35 + 1.9 +
1.7 + 8.3 s), 37.6 s on the AR-100 wafer (14.05 + 1.9 + 13.4 + 8.3 s) and 29.4 s at the same-supply optimum of
Section 3.8 (9.8 + 1.9 + 9.4 + 8.3 s, the water pulse scaled by the TMA ratio since the water half-cycle was not in
the design, and the purges held at the baseline values; the optimum's lower carrier flow raises the chamber
residence time 2.8-fold, so its TMA purge would scale to about 5 s and its cycle to about 33 s, a scaling estimate
rather than a run). At 0.103 nm per cycle, 5 nm is 49 cycles: 10 min per wafer and 5.5 wafers per hour with 60 s of
handling on the blanket wafer, 31 min and 1.9 wafers per hour on the AR-100 wafer, 24 min and 2.4 at the
optimum; 10 nm doubles the minutes. With gas clearing only for water (kdes = 0, 2.2 s) the three cycles are 6.2,
31.6 and 23.3 s. The rule turns this into a throughput formula, 3600 / (N [2 x 1.2 (1 + A_wall/A_wafer) N0/(U c0)
+ purges] + handling): below a wall-area ratio of about 10 the purges set the cycle and precursor fraction buys
little; above it throughput falls as 1/(1 + A_wall/A_wafer), and tripling the TMA fraction doubles the AR-100
throughput (1.8 to 3.5 wafers per hour) while leaving the blanket wafer where it was (5.9 to 6.3). The chamber run
at AR 100 sits on the 1 % curve within 2 %.

Sentence for the "What this means for a tool" box: on a 3D-NAND wafer the pulses are three quarters of the cycle
and the purges a quarter, the reverse of the blanket wafer, so the throughput lever is precursor delivery and the
purge lever is water desorption, not chamber flow.

## Proposed live-page text

Lead paragraph: In fab units: the rule's pulse gives 99 % step coverage, not just 90 %, because the saturation
front is 0.3 s wide at every aspect ratio. Fifty cycles of a blanket recipe leave 1.4 nm of a 5 nm film on top of
an AR-100 array and nothing at the bottom. The purge is not the chamber's doing: a 10 um hole empties in 14 us and
the chamber reaches 1 % in 1.9 s, while desorbing water can hold the exhaust up for 8 s. Put together, an Al2O3
cycle is 12 s on a blanket wafer and 38 s on the product wafer, 5.5 against 1.9 wafers per hour for 5 nm, and the
rule says how that falls with wall area and how much precursor fraction buys back.

Figure caption (cycle_time_wph.png): Al2O3 cycle on the blanket wafer, the AR-100 wafer and the same-supply
optimum (TMA pulse, purge to 1 %, H2O pulse, purge at kdes 0.1 1/s), and wafers per hour for 5 nm (49 cycles at
0.103 nm per cycle) against wall-area ratio from the rule at 1 % and 3 % TMA. On the product wafer the pulses are
73 % of the cycle; above a wall-area ratio of about 10 throughput falls as 1/(1 + A_wall/A_wafer) and precursor
fraction is the only knob that moves it. The point is the chamber run. Hover.

## CLAIMS rows

| Step coverage 90 / 95 / 99 % at 1 % TMA (worst face) | AR 20: 2.71 / 2.81 / 2.94 s; AR 50: 6.73 / 6.83 / 6.99 s; AR 100: 14.11 / 14.21 / 14.31 s; t99/rule 0.95-1.00 | figures/step_coverage_table.csv, figures/step_coverage_vs_pulse.png | PYTHONPATH=src .venv/bin/python coupling/fab_units.py |
| Film after 50 cycles of the blanket pulse (0.345 s) on the AR-100 wafer | 1.44 nm on the array top, 0.00 nm at the bottom, of 5.15 nm at 0.103 nm per cycle | figures/step_coverage_vs_pulse.csv | same |
| Growth per cycle used for thickness and cycle counts | 0.103 nm per cycle (PillarHall plateau 51.5 nm / 500 cycles); literature ~0.1 nm per cycle | figures/pillarhall_fit.png; literature: Puurunen, J. Appl. Phys. 97, 121301 (2005) | coupling/fit_pillarhall.py |
| Hole clearing by Knudsen diffusion, 100 nm x 10 um, TMA 573 K | L^2/D_K 7.3 us; bottom density to 1 % in 14.4 us (analytic eigenmode 14.4 us) | figures/purge_times.csv | coupling/fab_units.py |
| Chamber outlet clearing after the pulse | TMA to 10 / 1 / 0.1 %: 1.10 / 1.90 / 2.73 s (0.1 % extrapolated), tail time constant 0.36 s; H2O (kdes 0): 1.14 / 2.20 / 3.25 s; gap transit 0.30 s, V/Q 0.53 s | figures/purge_times.csv | same |
| Water purge to 1 % against desorption constant | kdes 0: 2.20 s; 0.1: 8.26 s (extrapolated); 1: 4.33 s (extrapolated); 10: 2.40 s | figures/purge_times.csv | same |
| Al2O3 cycle time (TMA 1 %, purge to 1 %, H2O, purge at kdes 0.1) | blanket 12.2 s; AR-100 37.6 s; same-supply optimum 29.4 s (H2O scaled by 9.83/14.05); with kdes 0 purge 6.2 / 31.6 / 23.3 s | figures/cycle_time_wph.csv, figures/cycle_time_wph.png | same |
| Wafers per hour, 5 nm (49 cycles), 60 s handling | blanket 5.5; AR-100 1.9; optimum 2.4 (no handling: 6.0 / 1.95 / 2.5); 10 nm: 2.9 / 0.96 / 1.2 | figures/cycle_time_wph.csv | same |
| Rule as throughput, 5 nm, 60 s handling | 1 % TMA: 5.9 WPH at wall ratio 1, 1.8 at 55.9; 3 %: 6.3 and 3.5; chamber run at 55.9 on the 1 % curve within 2 % | figures/wph_vs_wall_ratio.csv, figures/cycle_time_wph.png | same |

## Run log

### Night 3, part 3: fab units, purge and throughput (coupling/fab_units.py; figures/step_coverage_vs_pulse.png, .csv, step_coverage_table.csv, purge_times.csv, cycle_time_wph.png, .csv, wph_vs_wall_ratio.csv)
No new chamber runs. Inputs: the 1 % map runs at AR 20 / 50 / 100 (pulses 4.61 / 10.26 / 19.68 s, writes 0.023 /
0.051 / 0.098 s), the 20 s AR-100 run (0.05 s writes, purge to 22 s), the blanket TMA runs (0.5 s pulse to 1.5 s;
1.26 s pulse to 2.26 s), the water kdes 0 / 0.1 / 1 / 10 runs (5 s pulse to 8 s), `figures/pulse_time_map.csv`
(t_supply column, already (1 + A_wall/A_wafer) N0/(U c0)) and `figures/doe_gp_report.txt` (same-supply optimum
9.83 s). GPC 0.103 nm per cycle from the PillarHall plateau (night 1 part 10).
Step coverage SC = thetaBottom/thetaFlat at a face during the pulse, read as the coverage a pulse of that length
leaves (irreversible Langmuir; the chamber holds 0.019 c0 over the patterned wafer at the valve close, so the
post-pulse addition is negligible; the long run's 14.05 s and the map's 14.01 s agree to one write). Worst face over
faces 0-148 as in the map.
| AR | rule | t90 centre / worst | t95 | t99 | flat t90 / t95 / t99 | t99 / rule |
| 20 | 3.08 s | 2.67 / 2.71 | 2.78 / 2.81 | 2.92 / 2.94 | 1.84 / 2.14 / 2.60 | 0.96 |
| 50 | 7.32 s | 6.73 / 6.73 | 6.83 / 6.83 | 6.94 / 6.99 | 1.94 / 2.35 / 3.06 | 0.95 |
| 100 (map) | 14.38 s | 14.01 / 14.11 | 14.21 / 14.21 | 14.31 / 14.31 | 1.96 / 2.35 / 3.04 | 1.00 |
| 100 (20 s run) | 14.38 s | 14.05 / 14.10 | 14.15 / 14.20 | 14.25 / 14.30 | 1.95 / 2.30 / 3.05 | 0.99 |
Thickness after 50 cycles at the blanket pulse 0.345 s, bottom / top: AR 20 0.26 / 1.29 nm (SC 20 %), AR 50
0.01 / 1.32 nm, AR 100 0.00 / 1.44 nm; saturated 5.15 nm.
Purge (a), hole 100 nm x 10 um, TMA 573 K, vbar 410 m/s, D_K = d vbar/3 = 1.368e-5 m2/s: L^2/D_K 7.31 us; slowest
eigenmode 4 L^2/(pi^2 D_K) 2.96 us; bottom to 1 % by the first mode (4/pi) exp(-t/tau1): 14.36 us. Implicit 1D
transient (200 cells, 20 ns steps, n = 1 at t = 0, n(0) = 0, dn/dx(L) = 0): bottom to 1/e 3.70 us, to 1 % 14.42 us,
hole mean to 1 % 13.08 us. Hole clearing is 1.3e5 times faster than the chamber's 1 % time.
Purge (b), outlet patch mean after the pulse, fraction of the pulse-end value; tails fitted log-linearly over the
last 40 % of the written purge, extrapolated where the window ends first:
| run | outlet at pulse end | 10 % | 1 % | 0.1 % | tail decay (time constant) |
| TMA patterned AR 100, 20 s | 1.000 c0 | 1.10 s | 1.90 s | 2.73 s extr. | 2.75 1/s (0.36 s) |
| TMA blanket 0.5 s (1 s window; outlet still rising at pulse end, peak 0.367 c0 at +0.11 s) | 0.297 c0 | 1.25 extr. | 2.14 extr. | 3.04 extr. | 2.57 1/s (0.39 s) |
| TMA blanket 1.26 s (map, 1 s window) | 0.888 c0 | 1.07 extr. | 1.95 extr. | 2.83 extr. | 2.63 1/s (0.38 s) |
| H2O blanket kdes 0, 5 s | 0.999 c0 | 1.14 s | 2.20 s | 3.25 extr. | 2.18 1/s (0.46 s) |
References: gap transit 0.015 / 0.050 = 0.30 s; V/Q = pi 0.2^2 x 0.015 / 3.54e-3 = 0.53 s. The outlet time constant
sits between the two; the wafer flux itself clears within 0.05 s (part 1).
Purge (c), water outlet to 1 % of the pulse-end level (metric of coupling/postprocess_desorption.py): kdes 0 2.20 s
(run); 0.1 1/s 8.26 s and 1 1/s 4.33 s, extrapolated beyond the 3 s window as the desorption-fed level kdes theta
N0 A_wafer/(Q c0) (checked against the simulated outlet at 3 s: 0.0159 vs 0.0167 c0 and 0.0238 vs 0.0254 c0) with
theta decaying at the simulated rate over the last 2 s (0.088 and 0.654 1/s, below kdes because of re-adsorption);
10 1/s 2.40 s (run).
Cycle (TMA, purge 1.90 s, H2O, purge 8.26 s at kdes 0.1; in brackets kdes 0, 2.20 s): blanket 0.345 + 1.74 s pulses,
12.2 s (6.2); patterned AR 100 14.05 + 13.4 s, 37.6 s (31.6); same-supply optimum 9.83 + 9.38 s (H2O scaled by
9.83/14.05, not in the DOE; purges kept at baseline), 29.4 s (23.3). Caveat on the optimum: at 150 sccm, 631 Pa,
10 mm the volumetric flow is 8.4e-4 m3/s, the showerhead speed 0.012 m/s, gap transit 0.84 s and V/Q 1.5 s, 2.8
times the baseline; a TMA purge scaling with those would be about 5.3 s and the cycle about 32.8 s (estimate, not
run); the desorption-fed water purge scales with the precursor supply, held within 25 % of the baseline.
5 nm = 49 cycles, 10 nm = 98 at 0.103 nm per cycle. Minutes per wafer, 5 / 10 nm: 10.0 / 20.0, 30.7 / 61.4,
24.0 / 48.0. Wafers per hour at 5 nm, no overhead / 60 s handling: 6.0 / 5.5, 1.95 / 1.9, 2.5 / 2.4; at 10 nm
3.0 / 2.9, 0.98 / 0.96, 1.25 / 1.2.
Rule as throughput, 5 nm, both half-cycles at the same supply flux (constants of coupling/temporal_vs_spatial.py),
purges 1.90 + 8.26 s, 60 s handling: 1 % TMA 5.9 WPH at wall ratio 1, 4.2 at 11.2, 2.8 at 28, 1.8 at 55.9, 1.1 at
112; 3 % TMA 6.3, 5.5, 4.5, 3.5, 2.4. The AR-100 chamber run (37.6 s, 1.9 WPH) lies on the 1 % curve (39.3 s,
1.8 WPH).
Caveats: the purge criterion is 1 % of the pulse-end outlet level, a choice; kdes is a parameter; the blanket TMA
purge numbers are extrapolations from 1 s windows; the DOE optimum's purges are taken from the baseline chamber;
the time resolution of the step-coverage table is one write interval; no gas-phase reaction, no dead volumes, no
downstream plumbing, which in a real tool lengthen the purge.

## Sources

- R. L. Puurunen, "Surface chemistry of atomic layer deposition: A case study for the trimethylaluminum/water
  process", J. Appl. Phys. 97, 121301 (2005), doi:10.1063/1.1940727. Growth per cycle of the TMA/H2O process (about
  0.1 nm per cycle); the context for the 0.103 nm per cycle measured value used here.
- J. Yim, O. M. E. Ylivaara, M. Ylilammi, V. Korpelainen, E. Haimi, E. Verkama, M. Utriainen, R. L. Puurunen,
  "Saturation profile based conformality analysis for atomic layer deposition: aluminum oxide in lateral
  high-aspect-ratio channels", Phys. Chem. Chem. Phys. 22, 23107 (2020), doi:10.1039/D0CP03358H; data: Zenodo
  3974438, doi:10.5281/zenodo.3974438 (CC BY 4.0). The 51.5 nm / 500 cycles plateau (coupling/fit_pillarhall.py).
- M. Ylilammi, O. M. E. Ylivaara, R. L. Puurunen, "Modeling growth kinetics of thin films made by atomic layer
  deposition in lateral high-aspect-ratio structures", J. Appl. Phys. 123, 205301 (2018), doi:10.1063/1.5028178.
  Knudsen diffusivity d vbar/3 and the quasi-steady hole gas model.
- K. Arts, V. Vandalon, R. L. Puurunen, M. Utriainen, F. Gao, W. M. M. Kessels, H. C. M. Knoops, "Sticking
  probabilities of H2O and Al(CH3)3 during atomic layer deposition of Al2O3 extracted from their impact on film
  conformality", J. Vac. Sci. Technol. A 37, 030908 (2019), doi:10.1116/1.5093620. Sticking coefficients behind
  the chamber runs.
- S. M. George, "Atomic layer deposition: an overview", Chem. Rev. 110, 111 (2010), doi:10.1021/cr900056b. The
  four-step ALD cycle (reactant A, purge, reactant B, purge) as used for the cycle-time accounting.
- This repository: figures/pulse_time_map.csv (t_supply), figures/doe_gp_report.txt (same-supply optimum),
  paper/assumptions.md parts 1, 6, 7 and night 2 part 3 (baseline times, chamber flow, desorption runs).
