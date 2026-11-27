"""POST HOC: the cycle series measured by SEM-EDS (Yim 2020 ESI Fig. S8: 250, 500, 1000 cycles; digitised in data/lhar_eds/),
which has no reflectometer offset behind the front (background 1 to 5 % of the plateau). Each scan crosses the open area and
both mirrored channels; the entrances are the steps at about 213 and 292 um. Per channel: background = 18 a.u. (the level far from the film in the
unobstructed trace), plateau = median between 10 and 60 um depth, thickness = (I - background) / (plateau - background).
The model is the one registered in the second blind test (hydraulic kernel, s = 0.04 from the 5 % filled calibration curve,
spot 1 um), with and without channel narrowing, at the half depth measured on each EDS front; rms from the entrance to 1.6
half depths (the registered window). The 1000-cycle and the V0069 samples had registered predictions (yim2020_s11, _s02); the
250-cycle sample had none (its reflectometry curve did not pass the entry rule). Writes eds_cycle_series.csv and .png."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[4]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
S_CAL, H, EL, ER, BG = 0.04, 0.5, 213.0, 292.0, 18.0     # background 18 a.u.: the level far from the film in the unobstructed 1000-cycle trace (data/lhar_eds/PROVENANCE.md)
TR = {"figS8_Al_V0057_250cyc": (250, 26.0, 10), "figS8_Al_V0069_500cyc_as_labelled": (500, 54.5, 2), "figS8_Al_V0031_1000cyc": (1000, 99.3, 11)}
def model(z50, fill, spot):
    K = mt.kernel("hyd" if z50 / (1 - fill / 2) < 300 else "hyd_deep"); run = K[min(K, key=lambda a: abs(np.log(a / S_CAL)))]; xx = np.arange(0.5, run[2][-1], 1.0)
    def prof(e): y = mt.cumulative(run, e, fill, xx, spot, 1.0); return y / y[:4].mean()
    def c50(y): i = int(np.argmax(y < 0.5)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i])
    g = lambda le: c50(prof(10 ** le)) - z50; return xx, prof(10 ** brentq(g, np.log10(run[0][3]), np.log10(run[0][-2]), xtol=1e-4))
rows = []; curves = []; fig, ax = plt.subplots(2, 3, figsize=(12, 6.4))
for j, (name, (cyc, film, sample)) in enumerate(TR.items()):
    d = pd.read_csv(ROOT / "data/lhar_eds" / f"{name}.csv")          # flagged points (partly hidden by another trace) are kept: dropping them removes whole stretches
    for i, (side, dep) in enumerate((("left", EL - d.x_um.values), ("right", d.x_um.values - ER))):
        m = dep >= 0; x = dep[m]; I = d.intensity.values[m]; o = np.argsort(x); x, I = x[o], I[o]
        pl = np.median(I[(x > 10) & (x < 60)]); bg = BG
        y = (I - bg) / (pl - bg); ys = y.copy(); ys[1:-1] = np.median(np.c_[y[:-2], y[1:-1], y[2:]], axis=1); k = int(np.argmax(ys < 0.5)); x50 = x[k - 1] + (ys[k - 1] - 0.5) / (ys[k - 1] - ys[k]) * (x[k] - x[k - 1])
        z = x / H; z50 = x50 / H; fill = 2 * film / 500.0; w = z <= 1.6 * z50; xa, ya = model(z50, fill, 2.0); xb, yb = model(z50, 0.0, 2.0)
        ra = float(np.sqrt(np.mean((np.interp(z[w], xa, ya, right=0) - ys[w]) ** 2))); rb = float(np.sqrt(np.mean((np.interp(z[w], xb, yb, right=0) - ys[w]) ** 2)))
        b = (y <= 0.75) & (y >= 0.25) & w; sl = float(-np.polyfit(z[b], y[b], 1)[0]) if b.sum() >= 3 else np.nan
        rows.append(dict(trace=name, cycles=cyc, yim_table1_sample=sample, channel=side, fill=round(fill, 3), x50_um=round(float(x50), 1), background_over_plateau=round(float(bg / pl), 3), n_points=int(w.sum()), rms_with_narrowing=round(ra, 4), rms_without=round(rb, 4), slope_measured=round(sl, 5), n_points_25_75=int(b.sum())))
        g = np.arange(0, 1.7 * z50, 4.0); curves.append(dict(cycles=cyc, channel=side, fill=round(fill, 3), z50=round(float(z50), 1), rms_with=round(ra, 4), rms_without=round(rb, 4), z=[round(float(v), 1) for v in z[z <= 1.7 * z50]], y=[round(float(v), 3) for v in ys[z <= 1.7 * z50]], grid=g.tolist(), with_narrowing=[round(float(v), 3) for v in np.interp(g, xa, ya, right=0)], without=[round(float(v), 3) for v in np.interp(g, xb, yb, right=0)]))
        a = ax[i, j]; a.plot(z, ys, "k.", ms=2.5); a.plot(xa, ya, "C3-", lw=1.3, label="with narrowing"); a.plot(xb, yb, "C0--", lw=1, label="without"); a.set_xlim(0, 1.7 * z50); a.set_ylim(-0.08, 1.25); a.set_title(f"{cyc} cycles, {side} channel: rms {ra:.3f} / {rb:.3f}", fontsize=8); a.legend(fontsize=6, frameon=False)
fig.tight_layout(); fig.savefig(ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series.png", dpi=100); import json; (ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").write_text(json.dumps(curves, separators=(",", ":"))); D = pd.DataFrame(rows); D.to_csv(ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series.csv", index=False); print(D.drop(columns="trace").to_string(index=False))
print("mean rms with narrowing", round(D.rms_with_narrowing.mean(), 4), "without", round(D.rms_without.mean(), 4), "| within 0.06:", int((D.rms_with_narrowing <= 0.06).sum()), "of", len(D))
