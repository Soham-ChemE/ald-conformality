# Blind test: analysis after scoring (post hoc)

Everything here was done after the registered predictions (tag `blind-predictions-v1`) were scored. None of it changes the
result of the test: 27 curves scored, 7 pass, 20 fail (`../results.csv`). It records why, as far as the released curves show.

## What the curves look like (`diagnosis.csv`, `overlay_*.png`)

`droop` is 1 minus the mean thickness/plateau between 0.5 and 0.8 of the half depth: how far the profile has already fallen
before the front. The model assumes a flat plateau (droop about 0).

| profile | pass | fail | no slope |
|---|---|---|---|
| flat plateau (droop below 0.05) | 3 | 3 | 0 |
| drooping plateau | 4 | 17 | 1 |

- 13 of the 20 failures have a droop of 0.10 or more: all three TiO2 curves of Gao 2015 (0.12 to 0.23; the 1000 nm curve
  passed only because its prediction was a wide range), the ZnO curves at 125 and 150 C of Haimi 2026 and Gonsalves 2026
  (0.13 to 0.30), and four of the Elam 2003 pore profiles (0.11 to 0.28). An irreversible Langmuir model cannot produce such
  profiles at any s, with or without narrowing. Choosing them for prediction from the conditions alone was the predictor's error.
- The clean ZnO profiles of Gonsalves 2026 at 200 C (droop 0.005) are fitted by the registered model with s free at an rms of
  0.022 to 0.025 and s = 0.020 to 0.025; at 175 C, s = 0.008. They failed (200 C) because the prediction used s = 0.002
  transferred from another paper and reactor. A sticking coefficient did not transfer between laboratories or temperatures.
- The six inhibited-ZnO fronts are shallower than the smallest simulated s (0.0025) allows; the registered range was too narrow.
- Elam 2003: pores open at both faces, probe 3 to 5 um, tails of 0.1 to 0.3 of the plateau beyond the front; the measured
  slopes are about half the predicted ones.

## Slopes above the narrowing limit

Seven curves scored a slope above the free-molecular narrowing limit (1 - f/2)/(f z50). With the registered slope measure the
steep fronts have 2 to 5 points in the band. With all raw points between 0.25 and 0.75 of the plateau:

| curve | scored / limit | raw 0.25-0.75 / limit | remark |
|---|---|---|---|
| Gonsalves 2026 ZnO 200 C, left / right | 1.20 / 1.08 | 0.98 / 0.94 | at the limit |
| Elam 2003 ZnO d = 46 nm | 1.11 | 0.88 | below |
| Gao 2015 TiO2 1000 nm | 1.03 | 0.32 | missing readings in the front |
| Haimi 2026 ZnO sample C (250 C) | 1.17 | 1.13 | plateau with bumps of 15 %; 1400 Pa |
| Gao 2015 TiO2 500 nm | 1.25 | 1.16 | droop 0.12, missing readings in the front |
| Gao 2015 TiO2 200 nm | 1.29 | 1.76 | droop 0.23; the curve follows the predicted ramp to half thickness, then drops by a quarter of the plateau between two scan points |

Does the limit rise at higher carrier pressure? The penetration scales with the gap to the power 0.92 at 300 Pa
(figures/mechanism_gap_scaling.csv), 0.89 at 1400 Pa and 0.86 at 3000 Pa (`gap_scaling_pressure.csv`), which raises the limit
by about 10 to 15 %: enough for Haimi sample C at 1400 Pa, not for the TiO2 curves. For those the limit is exceeded in data
that have a drooping plateau and gaps at the front; whether a clean profile can exceed it is not shown either way.

## What follows, as hypotheses for a later test (not results)
1. The model applies to profiles with a flat plateau (droop below about 0.05); a profile that fails this is outside it.
2. A front cannot be predicted without a sticking coefficient measured for the same process at the same temperature.
3. The narrowing limit should carry the gap exponent n(pressure): f is replaced by 1 - (1 - f)^n.
4. A slope measure needs at least about six points in the band; with fewer, fit the whole profile.
