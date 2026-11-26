# Blind test of the lateral-channel model: protocol

Written and committed before any held-out curve was opened by the person or process that makes the predictions.

## Roles
- **Custodian**: a separate process that reads the source papers, digitises the held-out curves and keeps them in a folder
  outside this repository (`~/Documents/MATLAB/ald-blind-custody/`). It releases to the predictor only the *inputs* listed
  below. It does not make or see predictions before they are tagged.
- **Predictor**: writes the predictions from the inputs and the model in this repository, without opening the custody folder,
  the figures concerned, or any text that quotes their slopes or shapes.

## What the predictor may receive per curve (inputs)
Source and figure label; material and precursor; temperature; reactor pressure and carrier gas if stated; gap height; number
of cycles; plateau film thickness (or growth per cycle); measurement method and spot size; pulse and purge times; and the
half-thickness depth x50 (needed to set the exposure, which the papers do not give reliably). Nothing else: no slope, no
front width, no profile points, no penetration at other levels.

## What is predicted
For each curve, with the complete forward model (transport kernel, channel narrowing, measurement spot):
1. the front slope at half thickness, d(thickness/plateau)/d(x/H), as a value and an interval;
2. the 80 % and 20 % depths relative to x50;
3. where a prediction depends on a sticking coefficient that is not known for that precursor, either (a) the prediction is
   given as the narrowing-limited value with the range over s = 0.003 to 0.06, or (b) one curve of a series is designated in
   advance as the calibration curve, released in full, and the others are predicted with the s it gives.
Every prediction states in advance its pass limit: the predicted interval widened by 13 % (the scatter of repeated runs,
figures/mechanism_summary.json) and by the custodian's stated digitisation error.

## Order of events
1. This protocol is committed.
2. The custodian writes `inputs.json` (inputs only) into `docs/blindtest/` and keeps the curves sealed.
3. The predictor commits `predictions.csv` and the script that produced it, and the commit is tagged `blind-predictions-v1`
   and pushed. The tag time is the registration time.
4. The custodian releases the curves into `data/blind_revealed/` with provenance; the scoring script compares and writes
   `results.csv`. All results are reported, including failures. No prediction is changed after the tag; any later analysis
   is labelled post hoc.

## Known weaknesses
The curves are published figures, so "blind" means not seen by the predictor, not unpublished; half-thickness depths are
released as inputs; digitisation adds error; the predictor has read the text of some of the source papers for other purposes
(listed per curve in `inputs.json` by the custodian, with what was already known).
