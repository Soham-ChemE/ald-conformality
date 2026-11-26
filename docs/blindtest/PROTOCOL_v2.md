# Second blind test: protocol (written after the first test failed; committed before any new curve is opened)

The first test (tag `blind-predictions-v1`) scored 7 of 27. Its post hoc analysis (`posthoc/POSTHOC.md`) suggested three rules.
They are hypotheses, and this test is their registration. Roles and the order of events are as in `PROTOCOL.md`; the custody
folder is `~/Documents/MATLAB/ald-blind-custody-v2/`.

## Entry rule, applied by the custodian before anything is released
A curve enters only if all hold (the predictor never sees curves that fail, only their count and the reason):
1. measured thickness (or coverage) against depth in a slit or lateral channel, not used or seen by the predictor before
   (not in `data/lhar_profiles/`, `data/blind_revealed/`, Zenodo 3974438, or any figure the predictor digitised);
2. flat plateau: 1 - mean(thickness/plateau) between 0.5 and 0.8 of the half depth is below 0.05;
3. at least six measured points with thickness/plateau between 0.25 and 0.75;
4. half depth reached, and between 60 and 1000 gap heights;
5. gap, number of cycles (or film thickness) and spot or step size known;
6. it belongs to a series of at least two entering curves of the same process, temperature and reactor that differ in cycle
   number, gap height, pulse or exposure.

## Design: calibrate on one curve of a series, predict the others
For each series the custodian designates the calibration curve by a rule fixed here: the entering curve with the smallest
filled fraction 2t/H (ties: the smallest half depth). It is released in full. The predictor fits one sticking coefficient to
it with the complete model (hydraulic kernel, channel narrowing with exponent 1, top-hat spot; half depth fixed, amplitude
free; s over 0.0025 to 0.06) and then predicts every other curve of the series from its inputs (as listed in `PROTOCOL.md`).
If the calibration fit runs to an end of the s range, the series is predicted with that end value and flagged.

## What is predicted and how it is scored (fixed now)
For each predicted curve: the whole normalised profile at the measured half depth.
- Primary score: rms difference between the predicted and the measured thickness/plateau over all points from the entrance to
  1.6 times the half depth (measured values smoothed with a 3-point running median). Pass if rms <= 0.06, or <= 1.5 times the
  rms of that series' calibration fit, whichever is larger.
- Secondary (reported, not pass/fail): slope from all raw points between 0.25 and 0.75 against the predicted slope; depths at
  0.8 and 0.2 of the plateau over the half depth.
- Reference prediction, to show what narrowing adds: the same model without narrowing (fill = 0) at the same s. The test of
  the narrowing law is whether the complete model has the lower rms on the curves whose filled fraction exceeds the
  calibration curve's by more than 0.05.
All entering curves are scored and reported. Nothing is changed after the tag `blind-predictions-v2`.

## Known weaknesses
The predictor has seen 40 curves in the first test and knows the tabulated slopes of Yim 2020 Table 1 and the Haimi 2021
numbers; any curve from those sources is flagged by the custodian as "slope known to the predictor" and reported separately.
Published figures only. If fewer than six curves in two series enter, the test is reported as not informative.
