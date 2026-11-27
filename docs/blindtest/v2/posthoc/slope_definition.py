"""POST HOC: does the definition of the front slope explain why the EDS fronts are steeper than the tabulated reflectometry
slopes? For each digitised reflectometry profile of the TMA/water series (data/blind_revealed_v2/, the curves of the second
blind test) the slope is taken two ways on the front only (the contiguous run of readings around the first half-height crossing, so that
readings behind the front do not enter):
  ls_25_75   least squares over the points between 25 and 75 % of the plateau (the definition used for the EDS scans);
  secant     between the interpolated 60 % and 40 % crossings;
and compared with the slope the authors tabulate for that sample (Yim 2020 Table 1, their own linear fit) and, where an EDS
scan of the same sample exists, with the EDS slope of eds_cycle_series.csv. All slopes per gap height.
Writes slope_definition.csv.   .venv/bin/python docs/blindtest/v2/posthoc/slope_definition.py"""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[4]; HERE = Path(__file__).parent
T = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv").set_index("sample"); E = pd.read_csv(HERE / "eds_cycle_series.csv").groupby("yim_table1_sample").slope_measured.agg(["min", "max"]); rows = []
for f in sorted((ROOT / "data/blind_revealed_v2").glob("yim2020_s*.csv")):
    k = int(f.stem.split("_s")[1]); d = pd.read_csv(f); d = d[d.x_over_H >= 0].sort_values("x_over_H"); x, y = d.x_over_H.values, d.thickness_over_plateau.values
    ih = int(np.argmax(y < 0.5))                                    # first reading below half height
    def window(lo, hi):                                             # contiguous run of readings around the half-height crossing with lo <= y <= hi
        i = ih - 1
        while i > 0 and y[i] <= hi: i -= 1
        j = ih
        while j < len(y) - 1 and y[j] >= lo: j += 1
        return slice(i + (y[i] > hi), j + 1 - (y[j] < lo))
    def ls(lo, hi):
        w = window(lo, hi); return (float(-np.polyfit(x[w], y[w], 1)[0]) if (w.stop - w.start) >= 3 else np.nan), int(w.stop - w.start)
    def cr(l):
        w = window(min(l, .25), max(l, .75)); xs, ys = x[w], y[w]; k = int(np.argmax(ys < l)); k = max(k, 1); return float(xs[k - 1] + (ys[k - 1] - l) / (ys[k - 1] - ys[k]) * (xs[k] - xs[k - 1]))
    xf = x[window(.10, .90)]
    a, na = ls(.25, .75); sec = 0.2 / (cr(.4) - cr(.6)); tab = float(T.loc[k, "slope_dtheta_dxtilde"])
    rows.append(dict(sample=k, cycles=int(T.loc[k, "cycles"]), gap_nm=int(T.loc[k, "H_nm"]), objective=int(T.loc[k, "magnification"]), step_um=round(float(np.median(np.diff(d.x_um.values))), 2), points_front=len(xf), slope_table1=tab, slope_ls_25_75=round(a, 5), n_25_75=na, slope_secant_40_60=round(sec, 5),
                     ls_25_75_over_table=round(a / tab, 2), eds_slope_min=E["min"].get(k, np.nan), eds_slope_max=E["max"].get(k, np.nan), eds_over_optical_25_75=(round(float(E.loc[k].mean()) / a, 2) if k in E.index else np.nan)))
D = pd.DataFrame(rows); D.to_csv(HERE / "slope_definition.csv", index=False); print(D.to_string(index=False))
g = D[D.objective == 50]; print("50x samples: median ls_25_75 / table", round(g.ls_25_75_over_table.median(), 2), "range", g.ls_25_75_over_table.min(), g.ls_25_75_over_table.max(), "| secant / table median", round((g.slope_secant_40_60 / g.slope_table1).median(), 2))

# ---- the held-out slope test of the paper rescored on the slopes read from the profiles themselves ----
C = pd.read_csv(ROOT / "figures/mechanism_comparison.csv"); C = C[C["sample"].str.startswith("Yim")].copy(); C["k"] = C["sample"].str[3:].astype(int); C = C.merge(D[["sample", "slope_ls_25_75", "slope_secant_40_60"]].rename(columns={"sample": "k"}), on="k")
out = []; C["fill"] = 2 * C.film_nm / C.H_nm; C["step"] = (1 - C.fill / 2) / (C.fill * C.z50)
# sample 1 is listed separately: the custodian's note (written before scoring) says its trace lies under the markers of sample 2 and several points are missing
for sel, name in ((C.H_nm == 500, "500 nm gap, 9 samples with a digitised profile"), ((C.H_nm == 500) & (C.k != 1), "the same without sample 1 (trace partly hidden), 8 samples")):
    S = C[sel & C.predicted_slope.notna()]
    for m, g in S.groupby("model"):
        r = lambda col: round(float(np.sqrt(np.mean(np.log(g.predicted_slope / g[col]) ** 2))), 3)
        out.append(dict(samples=name, model=m, n=len(g), rms_ln_vs_table1=r("measured_slope"), rms_ln_vs_ls_25_75=r("slope_ls_25_75"), rms_ln_vs_secant_40_60=r("slope_secant_40_60")))
    g = S.drop_duplicates("k")
    def base(col): q = np.log(g.step / g[col]); return round(float(np.sqrt(((q - q.mean()) ** 2).mean())), 3)
    out.append(dict(samples=name, model="no kinetics: k (1 - f/2)/(f z50), k fitted", n=len(g), rms_ln_vs_table1=base("measured_slope"), rms_ln_vs_ls_25_75=base("slope_ls_25_75"), rms_ln_vs_secant_40_60=base("slope_secant_40_60")))
R = pd.DataFrame(out); R.to_csv(HERE / "slope_definition_heldout.csv", index=False); print(R.to_string(index=False))
