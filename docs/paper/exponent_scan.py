"""How strongly does channel narrowing act in the measured fronts? The closed form stretches the depth axis by xi^n; transport
theory gives n = 1 (no gas, hydraulic) or 0.92 (nitrogen, 300 Pa). Here n is varied as an empirical dial between 0 (no
narrowing) and 1, with everything else as in profile_fits.py, and three independent indications are collected:

 1. the 16 TMA/water fronts of the main series, fitted one by one for each n: the slope of ln s against the filled fraction
    (one process should give one s), the mean rms and the number of fits without an upper limit;
 2. the pair of scans from the second reactor: the 500-cycle scan predicted from the 250-cycle scan, with no kinetics and no
    transport model: theta_500(x) = [theta_250(x/a) + theta_250(x/(a (1 - f_250)^n))]/2, a being a depth scale between the
    two runs (1 if the exposure per cycle was the same);
 3. the upper bound on n from the fronts that are steeper than the steepest front (1 - f/2)/(n f).
Made after all other analyses, with all data known.
  .venv/bin/python docs/paper/exponent_scan.py   -> figures/exponent_scan_fits.csv, figures/exponent_scan_summary.csv,
                                                    figures/exponent_scan_pair.csv, figures/exponent_scan_bounds.csv
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(Path(__file__).parent))
import profile_fits as pf
import mechanism_test as mt
F = ROOT / "figures"; NS = (0.0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.92, 1.0)

# ---- 1. main series ----
if "pair-only" not in sys.argv:
    H = mt.heldout().set_index("sample"); profiles = []
    for f in sorted((ROOT / "data/blind_revealed_v2").glob("yim2020_s*.csv")):
        k = int(f.stem.split("_s")[1]); d = pd.read_csv(f); name = f"Yim{k}"
        if name in H.index: profiles.append(dict(id=f"opt_s{k:02d}", sample=name, method="reflectometry", z=d.x_over_H.values, y=d.thickness_over_plateau.values, spot_um=mt.SPOT_UM[int(H.loc[name, "magnification"])]))
    for c in json.loads((ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").read_text()):
        profiles.append(dict(id=f"eds_{c['cycles']}_{c['channel']}", sample={250: "Yim10", 500: "Yim2", 1000: "Yim11"}[c["cycles"]], method="EDS", z=np.array(c["z"]), y=np.array(c["y"]), spot_um=1.0))
    rows = []
    for win in (0.3, 0.1):
        pf.CUT = win
        for P in profiles:
            r = H.loc[P["sample"]]; fill = min(2 * r.film_nm / r.H_nm, 1.0); spot = P["spot_um"] * 1000 / r.H_nm
            for fam, kn in (("hydraulic", "hyd_deep" if r.H_nm == 500 else mt.kernel_for("hydraulic", r, fill)), ("gas", "gas300_deep" if (r.H_nm == 500 and r.z50 / (1 - fill / 2) > 250) else mt.kernel_for("gas", r, fill))):
                K = pf.kern(kn)
                for n in NS:
                    q = pf.fit(P["z"], P["y"], fill, spot, K, n, n > 0); rows.append(dict(window=win, id=P["id"], sample=P["sample"], method=P["method"], fill=round(fill, 3), transport=fam, n=n, s_fit=q["s_fit"], s_p16=q["s_p16"], s_p84=q["s_p84"], rms=q["rms"], no_upper_limit=q["no_upper_limit"]))
            print(win, P["id"], "done", flush=True)
    D = pd.DataFrame(rows); D.to_csv(F / "exponent_scan_fits.csv", index=False); out = []
    for (win, fam, n), g in D.groupby(["window", "transport", "n"]):
        for name, q in (("all 16 fronts", g), ("10 reflectometry fronts", g[g.method == "reflectometry"]), ("without the 1000-cycle sample", g[g["sample"] != "Yim11"])):
            c, cov = np.polyfit(q.fill, np.log(q.s_fit), 1, cov=True)
            out.append(dict(window=win, transport=fam, n=n, fronts=name, count=len(q), dlns_df=round(float(c[0]), 2), se=round(float(np.sqrt(cov[0, 0])), 2), sd_ln_s=round(float(np.log(q.s_fit).std(ddof=1)), 3), mean_rms=round(float(q.rms.mean()), 4), no_upper_limit=int(q.no_upper_limit.sum())))
    S = pd.DataFrame(out); S.to_csv(F / "exponent_scan_summary.csv", index=False); pd.set_option("display.width", 220); print(S[S.fronts == "all 16 fronts"].to_string(index=False))

# ---- 2. the pair from the second reactor ----
def load(c):
    d = pd.read_csv(ROOT / f"data/lhar_profiles/chipmetrics_guide_fig20_{c}cycles.csv"); d = d[d.x_um > 3].sort_values("x_um"); pl = d[(d.x_um > 40) & (d.x_um < 140)].t_nm.median(); fl = d[d.x_um > 220].t_nm.median(); return d.x_um.values, d.t_nm.values, float(pl), float(fl)
x1, t1, p1, f1 = load(250); x2, t2, p2, f2 = load(500); pr = []
for name, y1, y2, fa in (("as read", t1 / p1, t2 / p2, 2 * p1 / 500), ("floor subtracted", (t1 - f1) / (p1 - f1), (t2 - f2) / (p2 - f2), 2 * (p1 - f1) / 500)):
    g = lambda x: np.interp(x, x1, y1, left=y1[:5].mean(), right=y1[-20:].mean()); w = (x2 > 120) & (x2 < 215)
    for a in np.round(np.arange(0.94, 1.0601, 0.01), 2):
        for n in np.round(np.arange(0.0, 1.21, 0.05), 2):
            pr.append(dict(reading=name, filling_250=round(fa, 3), depth_scale=a, n=n, rms=round(float(np.sqrt(np.mean((0.5 * (g(x2 / a) + g(x2 / (a * (1 - fa) ** n)))[w] - y2[w]) ** 2))), 4)))
PR = pd.DataFrame(pr); PR.to_csv(F / "exponent_scan_pair.csv", index=False)
for name, g in PR.groupby("reading"):
    b = g.loc[g.rms.idxmin()]; ok = g[g.rms <= 1.5 * b.rms]; one = g[g.depth_scale == 1.0]; b1 = one.loc[one.rms.idxmin()]
    print(f"pair, {name}: best n {b.n} at depth scale {b.depth_scale} (rms {b.rms}); n from {ok.n.min()} to {ok.n.max()} within 1.5 times that rms (depth scale {ok.depth_scale.min()} to {ok.depth_scale.max()}); same exposure: n {b1.n} (rms {b1.rms}); rms at n = 1, same exposure: {float(one[one.n == 1.0].rms.iloc[0])}")

# ---- 3. bound from the fronts steeper than the steepest front ----
E = pd.read_csv(ROOT / "docs/blindtest/v2/posthoc/eds_slopes.csv"); B = pd.read_csv(F / "profile_fits_beneq_slopes.csv"); B = B[B.n == 1.0]
bd = [dict(front=f"EDS {int(r.cycles)} cycles, {r.channel}", fill=r.fill, slope_x_depth=round(r.slope_eds * r.z50_gaps, 2), n_upper_bound=round((1 - r.fill / 2) / (r.fill * r.slope_eds * r.z50_gaps), 2)) for r in E.itertuples()]
bd += [dict(front=f"second reactor {int(r.cycles)} cycles, {r.reading}", fill=r.fill, slope_x_depth=r.slope_x_z50, n_upper_bound=round(r.limit / r.slope_x_z50, 2)) for r in B.itertuples()]
BD = pd.DataFrame(bd); BD.to_csv(F / "exponent_scan_bounds.csv", index=False); print(BD.to_string(index=False))
