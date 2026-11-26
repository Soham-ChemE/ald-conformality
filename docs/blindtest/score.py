"""Scoring of the blind predictions (written and tagged together with the predictions, before the curves were released).
Reads data/blind_revealed/<id>.csv (custodian's columns x_over_H and thickness_over_plateau; x in units of the nominal gap or
pore diameter) and docs/blindtest/predictions.csv. For each curve: points sorted by x, 5-point running mean; slope = minus
the least-squares slope through the smoothed points with value between 0.35 and 0.65 that lie between the last crossing of
0.8 and the first crossing of 0.2 (at least three points, else the raw points in that band); r80 and r20 = first crossings of
0.8 and 0.2 over the first crossing of 0.5. A prediction passes if the measured slope lies in [pass_lo, pass_hi]; r80 and r20
are reported against their ranges and not used for pass or fail. Writes docs/blindtest/results.csv.
  .venv/bin/python docs/blindtest/score.py
"""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; P = pd.read_csv(ROOT / "docs/blindtest/predictions.csv"); rows = []
def cross(x, y, l):
    for i in range(1, len(y)):
        if y[i - 1] >= l > y[i]: return x[i - 1] + (y[i - 1] - l) / (y[i - 1] - y[i]) * (x[i] - x[i - 1])
    return np.nan
for r in P.itertuples():
    f = ROOT / "data/blind_revealed" / f"{r.id}.csv"
    if r.kind in ("excluded in advance", "not predicted") or not f.exists(): rows.append(dict(id=r.id, kind=r.kind, result="not scored")); continue
    d = pd.read_csv(f).dropna(subset=["x_over_H", "thickness_over_plateau"]).sort_values("x_over_H"); d = d[d.x_over_H >= 0]; x = d.x_over_H.values; y = d.thickness_over_plateau.values
    ys = np.convolve(y, np.ones(5) / 5, mode="same") if len(y) >= 9 else y.copy(); ys[:2], ys[-2:] = y[:2], y[-2:]
    x50 = cross(x, ys, 0.5); x80 = cross(x, ys, 0.8); x20 = cross(x, ys, 0.2); lo = x80 if np.isfinite(x80) else x.min(); hi = x20 if np.isfinite(x20) else x.max()
    m = (ys <= 0.65) & (ys >= 0.35) & (x >= lo) & (x <= hi)
    if m.sum() < 3: m = (y <= 0.65) & (y >= 0.35) & (x >= lo) & (x <= hi); ys = y
    slope = -np.polyfit(x[m], ys[m], 1)[0] if m.sum() >= 2 else np.nan
    ok = bool(r.pass_lo <= slope <= r.pass_hi) if np.isfinite(slope) else None
    rows.append(dict(id=r.id, kind=r.kind, n_points=int(len(x)), n_in_band=int(m.sum()), z50_measured=round(float(x50), 1), z50_input=r.z50, slope_measured=round(float(slope), 5), slope_pred=r.slope_pred, pass_lo=r.pass_lo, pass_hi=r.pass_hi,
                     result="pass" if ok else "fail" if ok is False else "no slope", measured_over_pred=round(float(slope / r.slope_pred), 3) if r.slope_pred == r.slope_pred else None, measured_over_narrowing_limit=round(float(slope / r.narrowing_limit_slope), 3),
                     r80_measured=round(float(x80 / x50), 3), r80_pred=r.r80_pred, r20_measured=round(float(x20 / x50), 3), r20_pred=r.r20_pred))
R = pd.DataFrame(rows); R.to_csv(ROOT / "docs/blindtest/results.csv", index=False)
print(R.to_string(index=False)); s = R[R.result.isin(["pass", "fail"])]; print(s.groupby("kind").result.value_counts().unstack(fill_value=0))
