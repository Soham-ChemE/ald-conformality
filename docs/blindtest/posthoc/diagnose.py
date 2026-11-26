"""POST HOC diagnosis of the blind-test misses (run after scoring; nothing here is a prediction).
For every scored curve:
  droop      1 - mean thickness/plateau between 0.5 and 0.8 of the half depth: how far the 'plateau' has already fallen before
             the front (an ideal saturating ALD profile gives about 0; the model assumes it)
  fit        the registered model (hydraulic kernel, narrowing with exponent 1, spot) fitted to the whole curve with s free over
             the simulated range, half depth fixed, amplitude free: best s, rms, and whether s runs to the end of the range
  slopes     the scored slope and two alternatives (raw points 0.25-0.75; steepest 5-point secant) over the narrowing limit
Writes docs/blindtest/posthoc/diagnosis.csv.   PYTHONPATH=src .venv/bin/python docs/blindtest/posthoc/diagnose.py"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
P = pd.read_csv(ROOT / "docs/blindtest/predictions.csv"); P = P[P.kind.isin(["ceiling", "range", "transfer", "bound"])]; R = pd.read_csv(ROOT / "docs/blindtest/results.csv").set_index("id"); rows = []
def model(K, sk, z50, fill, spot):
    run = K[sk]; xx = np.arange(0.5, run[2][-1], 1.0)
    def prof(e): y = mt.cumulative(run, e, fill, xx, spot, 1.0); return y / y[:4].mean()
    def c50(y): i = int(np.argmax(y < 0.5)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i])
    g = lambda le: c50(prof(10 ** le)) - z50; lo, hi = np.log10(run[0][3]), np.log10(run[0][-2])
    if not (np.isfinite(g(lo)) and np.isfinite(g(hi))) or g(lo) * g(hi) > 0: return None
    return xx, prof(10 ** brentq(g, lo, hi, xtol=1e-4))
for r in P.itertuples():
    d = pd.read_csv(ROOT / "data/blind_revealed" / f"{r.id}.csv").dropna(subset=["x_over_H", "thickness_over_plateau"]).sort_values("x_over_H"); d = d[d.x_over_H >= 0]; x, y = d.x_over_H.values, d.thickness_over_plateau.values
    z50 = R.loc[r.id, "z50_measured"]; z50 = r.z50 if not np.isfinite(z50) else z50; m = (x > 0.5 * z50) & (x < 0.8 * z50); droop = float(1 - np.mean(y[m])) if m.sum() >= 2 else np.nan
    fill = r.fill if r.fill == r.fill else 0.184; spot = r.spot_units if r.spot_units == r.spot_units else 20.0; lim = (1 - fill / 2) / (fill * z50)
    b = (y <= 0.75) & (y >= 0.25); raw = -np.polyfit(x[b], y[b], 1)[0] if b.sum() >= 3 else np.nan
    sec = max(((y[i] - y[i + 4]) / (x[i + 4] - x[i]) for i in range(len(x) - 4) if x[i + 4] > x[i] and 0.2 < 0.5 * (y[i] + y[i + 4]) < 0.8), default=np.nan)
    best = (np.nan, np.nan, False); 
    if z50 / (1 - fill / 2) < 1150:
        K = mt.kernel("hyd" if z50 / (1 - fill / 2) < 300 else "hyd_deep"); res = []
        for sk in K:
            mm = model(K, sk, z50, fill, spot)
            if mm is None: continue
            ym = np.interp(x, mm[0], mm[1], right=0.0); a = float(ym @ y / (ym @ ym)); res.append((sk / (4 if r.geometry == "pore" else 1), float(np.sqrt(np.mean((a * ym - y) ** 2)))))
        if res: res = np.array(res); j = int(np.argmin(res[:, 1])); best = (res[j, 0], res[j, 1], j in (0, len(res) - 1))
    rows.append(dict(id=r.id, kind=r.kind, result=R.loc[r.id, "result"], fill=round(fill, 3), z50=round(float(z50), 1), droop=round(droop, 3), s_fit=best[0], rms_fit=round(best[1], 4) if best[1] == best[1] else None, s_fit_at_range_end=bool(best[2]),
                     scored_over_limit=round(float(R.loc[r.id, "slope_measured"] / lim), 2), raw_25_75_over_limit=round(float(raw / lim), 2), steepest_secant_over_limit=round(float(sec / lim), 2), n_points_25_75=int(b.sum())))
D = pd.DataFrame(rows); D.to_csv(ROOT / "docs/blindtest/posthoc/diagnosis.csv", index=False); print(D.to_string(index=False))
