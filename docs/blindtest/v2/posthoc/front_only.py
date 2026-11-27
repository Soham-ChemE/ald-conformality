"""POST HOC to the second blind test (after scoring; not a registered result). The registered score ran to 1.6 half depths
and so included the region behind the front, where the reflectometry profiles carry a residual signal and bumps of 0.1 to 0.5
of the plateau. Here the same predictions are compared with the measurement from the entrance to the first point where the
measured profile (3-point median) falls below 0.2 of the plateau. Writes front_only.csv."""
import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[4]; V2 = ROOT / "docs/blindtest/v2"; S = list(json.loads((V2 / "predictions_v2.json").read_text())["series"].values())[0]; rows = []
for cid, c in S["curves"].items():
    d = pd.read_csv(ROOT / "data/blind_revealed_v2" / f"{cid}.csv").dropna(subset=["x_over_H", "thickness_over_plateau"]).sort_values("x_over_H"); d = d[d.x_over_H >= 0]; x = d.x_over_H.values; y = d.thickness_over_plateau.values
    ys = y.copy(); ys[1:-1] = np.median(np.c_[y[:-2], y[1:-1], y[2:]], axis=1); k = int(np.argmax(ys < 0.2)) if (ys < 0.2).any() else len(x); m = np.arange(len(x)) <= k; xm = np.array(c["x"])
    ra = float(np.sqrt(np.mean((np.interp(x[m], xm, c["complete"]) - ys[m]) ** 2))); rb = float(np.sqrt(np.mean((np.interp(x[m], xm, c["no_narrowing"]) - ys[m]) ** 2)))
    tail = float(np.mean(ys[~m])) if (~m).sum() else np.nan
    rows.append(dict(id=cid, n_points=int(m.sum()), rms_complete=round(ra, 4), rms_no_narrowing=round(rb, 4), within_0p06=bool(ra <= 0.06), mean_signal_behind_front=round(tail, 3), max_signal_behind_front=round(float(np.max(ys[~m])) if (~m).sum() else np.nan, 3)))
D = pd.DataFrame(rows); D.to_csv(V2 / "posthoc/front_only.csv", index=False); print(D.to_string(index=False)); print("within 0.06:", int(D.within_0p06.sum()), "of", len(D), "| mean rms complete", round(D.rms_complete.mean(), 4), "no narrowing", round(D.rms_no_narrowing.mean(), 4))
