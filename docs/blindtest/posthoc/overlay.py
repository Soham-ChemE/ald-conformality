"""POST HOC (after the blind test was scored): overlay of the released curves and the registered model profile at the central
sticking coefficient, to see where the misses are. Writes docs/blindtest/posthoc/overlay_<group>.png. Not part of the test."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
P = pd.read_csv(ROOT / "docs/blindtest/predictions.csv"); P = P[P.kind.isin(["ceiling", "range", "transfer", "bound"])]; R = pd.read_csv(ROOT / "docs/blindtest/results.csv").set_index("id")
def model(r, s):
    K = mt.kernel("hyd" if r.z50 / (1 - r.fill / 2) < 300 else "hyd_deep"); sk = min(K, key=lambda a: abs(np.log(a / (4 * s if r.geometry == "pore" else s)))); run = K[sk]; xx = np.arange(0.5, run[2][-1], 1.0)
    def prof(e): y = mt.cumulative(run, e, r.fill, xx, r.spot_units, 1.0); return y / y[:4].mean()
    def c50(y): i = int(np.argmax(y < 0.5)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i])
    g = lambda le: c50(prof(10 ** le)) - r.z50
    try: return xx, prof(10 ** brentq(g, np.log10(run[0][3]), np.log10(run[0][-2]), xtol=1e-4))
    except Exception: return None, None
groups = {"gao": "gao2015", "elam": "elam2003", "arts": "arts2019", "gonsalves_plain": "gonsalves2026_fig2_ZnO_DEZ", "gonsalves_inhibited": "gonsalves2026_fig2_ZnO_MeOH", "haimi": "haimi2026"}
for gname, key in groups.items():
    g = P[P.id.str.startswith(key)]; n = len(g); fig, ax = plt.subplots(1, n, figsize=(3.6 * n, 3.2), squeeze=False)
    for a, r in zip(ax[0], g.itertuples()):
        d = pd.read_csv(ROOT / "data/blind_revealed" / f"{r.id}.csv"); a.plot(d.x_over_H, d.thickness_over_plateau, "k.", ms=3)
        if r.kind != "bound":
            for s, st in ((r.s_central, "C3-"), (r.s_low if r.s_low != r.s_central else r.s_high, "C0--")):
                x, y = model(r, s)
                if x is not None: a.plot(x, y, st, lw=1.2, label=f"s = {s:g}")
        res = R.loc[r.id]; a.set_title(f"{r.id.replace(key, '')[:26]}\n{res.result}: meas {res.slope_measured}, pred {r.slope_pred}", fontsize=7); a.set_xlim(0, max(r.z50 * 2.0, 50)); a.set_ylim(-0.05, 1.3); a.legend(fontsize=6, frameon=False); a.axvline(r.z50, color="grey", lw=.5)
    fig.tight_layout(); fig.savefig(ROOT / f"docs/blindtest/posthoc/overlay_{gname}.png", dpi=110); plt.close(fig)
print("ok")
