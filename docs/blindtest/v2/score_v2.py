"""Scoring of the second blind test, fixed together with the predictions (PROTOCOL_v2.md). Reads data/blind_revealed_v2/<id>.csv
(x_over_H, thickness_over_plateau) and predictions_v2.json. Primary: rms between predicted and measured thickness/plateau
(3-point running median of the measurement) from the entrance to 1.6 half depths; pass if rms <= the series' pass limit.
Reference: the same without narrowing. Secondary: slope from all raw points between 0.25 and 0.75. Writes results_v2.csv.
  .venv/bin/python docs/blindtest/v2/score_v2.py"""
import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[3]; V2 = ROOT / "docs/blindtest/v2"; P = json.loads((V2 / "predictions_v2.json").read_text()); T = pd.read_csv(V2 / "predictions_v2.csv").set_index("id"); rows = []
for sid, S in P["series"].items():
    for cid, c in S["curves"].items():
        f = ROOT / "data/blind_revealed_v2" / f"{cid}.csv"
        if not f.exists(): rows.append(dict(id=cid, result="not released")); continue
        d = pd.read_csv(f).dropna(subset=["x_over_H", "thickness_over_plateau"]).sort_values("x_over_H"); x = d.x_over_H.values; y = d.thickness_over_plateau.values
        ys = y.copy(); ys[1:-1] = np.median(np.c_[y[:-2], y[1:-1], y[2:]], axis=1); z = T.loc[cid, "z50"]; m = (x >= 0) & (x <= 1.6 * z); xm = np.array(c["x"])
        ra = float(np.sqrt(np.mean((np.interp(x[m], xm, c["complete"], right=0.0) - ys[m]) ** 2))); rb = float(np.sqrt(np.mean((np.interp(x[m], xm, c["no_narrowing"], right=0.0) - ys[m]) ** 2)))
        b = (y <= 0.75) & (y >= 0.25) & (x >= 0); sl = float(-np.polyfit(x[b], y[b], 1)[0]) if b.sum() >= 3 else np.nan
        rows.append(dict(id=cid, series=sid, cycles=T.loc[cid, "cycles"], gap_nm=T.loc[cid, "gap_nm"], fill=T.loc[cid, "fill"], fill_minus_calibration=T.loc[cid, "fill_minus_calibration"], n_points=int(m.sum()), rms_complete=round(ra, 4), rms_no_narrowing=round(rb, 4), pass_limit=S["pass_limit"],
                         result="pass" if ra <= S["pass_limit"] else "fail", narrowing_helps=bool(ra < rb), slope_measured=round(sl, 5), slope_complete=T.loc[cid, "slope_complete"], slope_no_narrowing=T.loc[cid, "slope_no_narrowing"], n_points_25_75=int(b.sum())))
R = pd.DataFrame(rows); R.to_csv(V2 / "results_v2.csv", index=False); print(R.to_string(index=False))
s = R[R.result.isin(["pass", "fail"])]; print("pass", int((s.result == "pass").sum()), "of", len(s)); h = s[s.fill_minus_calibration > 0.05]; print("curves filled more than the calibration by > 0.05:", len(h), "| complete model better on", int(h.narrowing_helps.sum()), "| mean rms complete", round(float(h.rms_complete.mean()), 4), "no narrowing", round(float(h.rms_no_narrowing.mean()), 4))
