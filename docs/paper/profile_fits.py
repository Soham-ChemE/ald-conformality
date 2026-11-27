"""Full-profile test of channel narrowing (post hoc; made after the blind tests and after an independent review).
Every digitised front of the TMA/water PillarHall series is fitted on its own, with the sticking coefficient s and the
exposure free, with and without channel narrowing, for three transport kernels. If narrowing is the right physics the fitted
s should not depend on the number of cycles or on the filled fraction; without it the fitted s should fall as the channel fills.
Profiles: reflectometry curves of the second blind test (data/blind_revealed_v2/, Yim 2020 and Haimi 2021) and the EDS line
scans (docs/blindtest/v2/posthoc/eds_cycle_series_curves.json, two channels per scan). Fit window: from 5 gap heights to the
first reading below 0.1 of the plateau, so the residual behind the front of the reflectometry curves does not enter.
For each s of the kernel the exposure is fitted; s is then taken at the minimum of a parabola through the rms in ln s.
Writes figures/profile_fits.csv and figures/profile_fits_summary.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/profile_fits.py"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
H = mt.heldout().set_index("sample"); T = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv").set_index("sample")
profiles = []
for f in sorted((ROOT / "data/blind_revealed_v2").glob("yim2020_s*.csv")):
    k = int(f.stem.split("_s")[1]); d = pd.read_csv(f); name = f"Yim{k}"
    if name not in H.index: continue                                   # sample 8 / 9: calibration curves are not in the held-out table
    profiles.append(dict(id=f"opt_s{k:02d}", sample=name, method="reflectometry", z=d.x_over_H.values, y=d.thickness_over_plateau.values, spot_um=mt.SPOT_UM[int(H.loc[name, "magnification"])]))
for c in json.loads((ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").read_text()):
    name = {250: "Yim10", 500: "Yim2", 1000: "Yim11"}[c["cycles"]]; profiles.append(dict(id=f"eds_{c['cycles']}_{c['channel']}", sample=name, method="EDS", z=np.array(c["z"]), y=np.array(c["y"]), spot_um=1.0))
rows = []
for P in profiles:
    r = H.loc[P["sample"]]; fill = min(2 * r.film_nm / r.H_nm, 1.0); z, y = P["z"], P["y"]; o = np.argsort(z); z, y = z[o], y[o]
    iend = int(np.argmax((y < 0.1) & (z > 20))) if ((y < 0.1) & (z > 20)).any() else len(y) - 1; w = (z >= 5) & (np.arange(len(z)) <= iend); zw, yw = z[w], y[w]; spot = P["spot_um"] * 1000 / r.H_nm
    for fam in mt.FAMILY:
        K = mt.kernel(mt.kernel_for(fam, r, fill)); n = mt.GAP_EXPONENT[fam]
        for narrow in (True, False):
            S, RMS, EXP = [], [], []
            for s in sorted(K):
                run = K[s]; E = run[0]
                def cost(le):
                    m = mt.cumulative(run, 10 ** le, fill if narrow else 0.0, np.concatenate([[1.0, 2.0, 3.0, 4.0], zw]), spot, n); return float(np.sqrt(np.mean((m[4:] / m[:4].mean() - yw) ** 2)))
                g = np.linspace(np.log10(E[1]), np.log10(E[-1]), 25); c = [cost(v) for v in g]; j = int(np.argmin(c)); res = minimize_scalar(cost, bounds=(g[max(j - 1, 0)], g[min(j + 1, len(g) - 1)]), method="bounded", options=dict(xatol=1e-3))
                S.append(s); RMS.append(res.fun); EXP.append(10 ** res.x)
            S, RMS = np.array(S), np.array(RMS); j = int(np.argmin(RMS)); edge = j in (0, len(S) - 1); sb, rb = S[j], RMS[j]
            if not edge:
                a = np.polyfit(np.log(S[j - 1:j + 2]), RMS[j - 1:j + 2], 2)
                if a[0] > 0: sb = float(np.exp(-a[1] / (2 * a[0]))); rb = float(np.polyval(a, np.log(sb)))
            ok = S[RMS <= 1.10 * RMS[j]]
            rows.append(dict(id=P["id"], sample=P["sample"], method=P["method"], cycles=int(r.cycles), H_nm=int(r.H_nm), fill=round(fill, 3), n_points=int(w.sum()), transport=fam, narrowing=narrow, s_fit=float(f"{sb:.3g}"), rms=round(rb, 4), s_lo_10pct=float(ok.min()), s_hi_10pct=float(ok.max()), at_grid_edge=bool(edge), exposure_Pa_s=round(EXP[j], 3), s_grid_min=float(S.min()), s_grid_max=float(S.max())))
    print(P["id"], "done", flush=True)
D = pd.DataFrame(rows); D.to_csv(ROOT / "figures/profile_fits.csv", index=False)
out = []
for (fam, narrow), g in D.groupby(["transport", "narrowing"]):
    for name, q in (("all profiles", g), ("500 nm gap", g[g.H_nm == 500]), ("500 nm gap, EDS only", g[(g.H_nm == 500) & (g.method == "EDS")]), ("500 nm gap, reflectometry only", g[(g.H_nm == 500) & (g.method == "reflectometry")])):
        if len(q) < 3: continue
        ls = np.log(q.s_fit); sl = np.polyfit(q.fill, ls, 1)[0] if q.fill.nunique() > 2 else np.nan
        out.append(dict(transport=fam, narrowing=narrow, profiles=name, n=len(q), n_at_grid_edge=int(q.at_grid_edge.sum()), s_geomean=float(f"{np.exp(ls.mean()):.3g}"), s_min=q.s_fit.min(), s_max=q.s_fit.max(), sd_ln_s=round(float(ls.std(ddof=1)), 3), dlns_dfill=round(float(sl), 2), mean_rms=round(float(q.rms.mean()), 4)))
Sm = pd.DataFrame(out); Sm.to_csv(ROOT / "figures/profile_fits_summary.csv", index=False); pd.set_option("display.width", 250)
print(D[D.transport == "gas"].pivot_table(index=["id", "cycles", "fill"], columns="narrowing", values=["s_fit", "rms", "at_grid_edge"], aggfunc="first").to_string()); print(Sm.to_string(index=False))
