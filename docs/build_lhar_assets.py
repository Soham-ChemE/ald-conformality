#!/usr/bin/env python3
"""Data for the page's re-analysis scenes (section 03b), from committed result files only:
figures/ballistic_slope_calibration.csv, figures/lhar_slope_table.csv, figures/lhar_single_exposure_fit.csv,
figures/lhar_single_curve_*.csv, figures/lhar_refit.csv, figures/pillarhall_fit_cutoff.csv, data/lhar_profiles/, the Zenodo
workbook, and figures/pillar_mc_*.csv when present. Writes docs/assets/data-lhar.js (window.ALD.lhar).
  PYTHONPATH=src .venv/bin/python docs/build_lhar_assets.py
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
DOCS = Path(__file__).resolve().parent; ROOT = DOCS.parent; FIG = ROOT / "figures"; sys.path.insert(0, str(ROOT / "src"))
from aldconf.data import load_pillarhall, default_path
def r(x, d=4):
    if isinstance(x, (list, tuple, np.ndarray)): return [r(v, d) for v in x]
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else float(f"{float(x):.{d}g}")
out = {}
# ---- slope relation ----
cal = pd.read_csv(FIG / "ballistic_slope_calibration.csv"); b = cal[cal.model == "ballistic"]; h = cal[cal.model == "1d_hydraulic"]
curves = []
for s0 in sorted(b.s.unique()):
    d = b[b.s == s0].sort_values("z50"); zg = np.geomspace(max(70, d.z50.min()), d.z50.max(), 36)
    curves.append(dict(s=float(s0), z50=r(zg, 4), A=r(np.interp(zg, d.z50, d.A), 4)))
bb = b[(b.z50 >= 80) & (b.z50 <= 500) & (b.s >= 1e-3)]; X = np.c_[np.ones(len(bb)), np.log(bb.z50), np.log(bb.s)]; k, *_ = np.linalg.lstsq(X, np.log(bb.A), rcond=None)
out["slope"] = dict(curves=curves, hydraulic=dict(z50=r(h.z50.values, 4), A=r(h.A.values, 4), median=r(h.A.median(), 4)), arts=13.9, gonsalves=11.1,
                    relation=dict(pre=r(np.exp(k[0]), 3), ez=r(k[1], 3), es=r(k[2], 3)), A_range=r([bb.A.min(), bb.A.max()], 3))
# ---- sticking coefficient by method ----
st = pd.read_csv(FIG / "lhar_slope_table.csv"); sf = pd.read_csv(FIG / "lhar_single_exposure_fit.csv"); rf = pd.read_csv(FIG / "lhar_refit.csv"); cut = pd.read_csv(FIG / "pillarhall_fit_cutoff.csv")
ref = st[(st.source == "Yim2020") & (st.H_nm == 500) & (st.magnification == 50) & (st.cycles == 500) & (st.pulse.astype(str) == "0.1")]
gm = lambda v: float(np.exp(np.log(v).mean()))
def g1(df, **kw):
    m = df
    for a, v in kw.items(): m = m[m[a] == v]
    return m.iloc[0]
ar = st[st["sample"] == 101].iloc[0]
methods = ["hydraulic (Arts relation)", "corrected, cut-off at the front", "corrected, cut-off at pillars or 1 mm", "ideal-slit ballistic"]
rows = [dict(name="Yim 2020, eight samples, tabulated slopes, 300 °C", key="yim", s=r([gm(ref.s_arts), gm(ref.s_corr_front), gm(ref.s_corr_pillar), gm(ref.s_ballistic)], 3),
             samples=[dict(sample=int(q["sample"]), s=r([q.s_arts, q.s_corr_front, q.s_corr_pillar, q.s_ballistic], 3), slope=float(q.slope), z50=float(q.x50_tilde)) for _, q in ref.iterrows()]),
        dict(name="Zenodo 3974438, single-exposure fit, 300 °C", key="zen", s=r([g1(sf, dataset="zenodo", model="1d_hydraulic").s, None, g1(sf, dataset="zenodo", model="1d_corr_1mm").s, g1(sf, dataset="zenodo", model="ballistic").s], 3)),
        dict(name="Ylilammi 2018 Fig. 6, single-exposure fit, 300 °C", key="yli", s=r([g1(sf, dataset="ylilammi_al2o3", model="1d_hydraulic").s, None, g1(sf, dataset="ylilammi_al2o3", model="1d_corr_1mm").s, g1(sf, dataset="ylilammi_al2o3", model="ballistic").s], 3)),
        dict(name="Zenodo 3974438, all cycles with narrowing, 300 °C", key="zenm", s=r([8.5e-3, None, g1(cut, L_phys_m=0.001).s, None], 3)),
        dict(name="Ylilammi 2018 Fig. 6, all cycles with narrowing, 300 °C", key="ylim", s=r([g1(rf, dataset="ylilammi_al2o3", transport="hyd").s, None, g1(rf, dataset="ylilammi_al2o3", transport="plates:1e-3").s, None], 3)),
        dict(name="Arts 2019 Fig. 2, TMA-limited, about 220 °C, digitised", key="arts", s=r([ar.s_arts, ar.s_corr_front, None, ar.s_ballistic], 3))]
out["methods"] = dict(names=methods, rows=rows, sfg=[3e-3, 4e-3, 5e-3], gakis=3.44e-3)
# ---- full-profile fits ----
def measured(ds):
    if ds == "zenodo":
        df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
        sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
        for sh in sheets:
            g = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
        mean = np.nanmean(np.array(prof), 0); y = mean / np.nanmean(mean[(grid > 5) & (grid < 60)]); ok = ~np.isnan(y); return grid[ok] * 2, y[ok]
    d = pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv"); return d.x_um.values * 2, d.thickness_nm.values / d.thickness_nm[(d.x_um > 5) & (d.x_um < 60)].mean()
fits = {}
for ds, title, xmax_fit in (("zenodo", "Zenodo 3974438 (Yim 2020 sample 8), mean of 8 scans", 300), ("ylilammi_al2o3", "Ylilammi et al. 2018, Fig. 6", None)):
    xm, ym = measured(ds); g = np.arange(0, 402, 2.0); e = dict(title=title, AR=r(g, 4), measured=dict(AR=r(xm, 4), y=r(ym, 4)), fit_to=xmax_fit, models={})
    for mdl in ("1d_hydraulic", "1d_corr_1mm", "ballistic"):
        c = pd.read_csv(FIG / f"lhar_single_curve_{ds}_{mdl}.csv"); q = g1(sf, dataset=ds, model=mdl)
        e["models"][mdl] = dict(y=r(np.interp(g, c.x_um * 2, c.model_norm), 4), s=r(q.s, 3), exposure=r(q.exposure_Pa_s, 3), rms=r(q.rms_all, 3), rms_front=r(q.rms_front, 3), rms_tail=r(q.rms_tail, 3))
    fits[ds] = e
out["fits"] = fits
# ---- pillar-layout Monte Carlo (added when the result file exists) ----
pm = FIG / "pillar_mc_summary.json"
if pm.exists(): out["pillar_mc"] = json.loads(pm.read_text())
(DOCS / "assets" / "data-lhar.js").write_text("window.ALD=window.ALD||{};window.ALD.lhar=" + json.dumps(out, separators=(",", ":")) + ";\n")
print("wrote docs/assets/data-lhar.js", (DOCS / "assets" / "data-lhar.js").stat().st_size, "bytes;", len(curves), "ballistic curves;", [x["key"] for x in rows])
