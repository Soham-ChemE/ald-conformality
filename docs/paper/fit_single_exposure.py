"""Like-for-like transport comparison on measured LHAR profiles: fit the single-exposure coverage profile theta(x) (no
narrowing, no multi-cycle averaging) to the normalised thickness with
  (a) ballistic: Machball trench AR 800 frames from docs/paper/ballistic_frames.py (figures/ballistic_frames/s_*_Nz1200.npz);
      exposure fitted per sticking value by interpolating the frames in time, s by a parabola through the three best grid
      values in log s; the curve saved is interpolated between the two neighbouring grid profiles;
  (b) 1D slit, hydraulic-diameter D_K;  (c) 1D slit, D_K x (3/8) ln(L_c/w) with L_c = 1 mm.
Fit window: the Zenodo mean profile is fitted for x <= 150 um only. Beyond that, six of the eight scans rise again (mean
normalised thickness 0.02 to 0.14 between 160 and 178 um), a feature no transport model here describes; fitting it would
credit whichever model has the highest tail. Residuals are reported for the plateau (y > 0.9), the front (0.1 < y < 0.9)
and the tail (y < 0.1) separately. No interval is quoted from the residual surface (the residuals are systematic); the
spread is taken from fitting the eight Zenodo scans one by one.
Writes figures/lhar_single_exposure_fit.csv, figures/lhar_single_exposure_scans.csv, figures/lhar_single_curve_<dataset>_<model>.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/fit_single_exposure.py
"""
import sys, csv, glob
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize, minimize_scalar
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Slit, Precursor, Process, simulate, slit_knudsen_factor
from aldconf.data import load_pillarhall, default_path
KB, GAP, T, MASS, N0, TPULSE, XMAX_ZENODO = 1.380649e-23, 500e-9, 573.15, 72.09, 4e18, 0.1, 150.0
def zenodo_scans():
    df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
    sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
    for sh in sheets:
        g = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
    return sheets, grid, np.array(prof)
def load(dataset):
    if dataset == "zenodo":
        _, grid, prof = zenodo_scans(); mean = np.nanmean(prof, 0); y = mean / np.nanmean(mean[(grid > 5) & (grid < 60)]); ok = ~np.isnan(y) & (grid <= XMAX_ZENODO)
    else:
        d = pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv"); grid = d.x_um.values.astype(float); y = d.thickness_nm.values / d.thickness_nm[(d.x_um > 5) & (d.x_um < 60)].mean(); ok = ~np.isnan(y)
    return grid[ok] * 1e-6, y[ok]
runs = []
for f in sorted(glob.glob(str(ROOT / "figures/ballistic_frames/s_*_Nz1200.npz"))):
    d = np.load(f); runs.append((float(d["s"]), d["exposure_Pa_s"], d["z"] * GAP, d["cov"]))
runs.sort(key=lambda r: r[0]); print("ballistic runs:", [r[0] for r in runs])
def frame(run, t):
    s, expo, z, cov = run; i = min(max(np.searchsorted(expo, t), 1), len(expo) - 1); w = (t - expo[i - 1]) / (expo[i] - expo[i - 1]); return z, (1 - w) * cov[i - 1] + w * cov[i]
def fit_ballistic(xg, yy):
    tab = []
    for run in runs:
        expo = run[1]
        def cost(le): z, th = frame(run, 10 ** le); return float(np.mean((np.interp(xg, z, th) - yy) ** 2))
        r = minimize_scalar(cost, bounds=(np.log10(expo[1]), np.log10(expo[-2])), method="bounded", options={"xatol": 1e-4}); tab.append((run[0], r.fun, 10 ** r.x))
    tab = np.array(tab); j = min(max(int(np.argmin(tab[:, 1])), 1), len(tab) - 2)
    p = np.polyfit(np.log(tab[j - 1:j + 2, 0]), tab[j - 1:j + 2, 1], 2); ls = -p[1] / (2 * p[0]) if p[0] > 0 else np.log(tab[j, 0]); ls = min(max(ls, np.log(tab[j - 1, 0])), np.log(tab[j + 1, 0]))
    k = j if ls >= np.log(tab[j, 0]) else j - 1; w = (ls - np.log(tab[k, 0])) / (np.log(tab[k + 1, 0]) - np.log(tab[k, 0]))
    z, a = frame(runs[k], tab[k, 2]); _, b = frame(runs[k + 1], tab[k + 1, 2]); th = (1 - w) * a + w * b
    return float(np.exp(ls)), float(np.exp((1 - w) * np.log(tab[k, 2]) + w * np.log(tab[k + 1, 2]))), z, th, tab
def fit_1d(xg, yy, factor):
    def model(log_s, log_n0):
        pr = Precursor(sticking=10 ** log_s, site_density=N0, mass_amu=MASS, gpc=1e-10)
        res = simulate(Process(Slit(GAP, 400 * GAP, knudsen_correction="none"), pr, T=T, p_partial=10 ** log_n0 * KB * T, t_dose=TPULSE, cycles=1, n_cells=400, n_substeps=200, knudsen_factor=factor, narrowing=False)); return res
    def cost(pp): res = model(*pp); return float(np.mean((np.interp(xg, res.x, res.theta) - yy) ** 2))
    best = None
    for s0 in (-2.6, -2.0):
        for n0 in (19.7, 20.3):
            r = minimize(cost, [s0, n0], method="Nelder-Mead", options={"xatol": 0.003, "fatol": 1e-9, "maxiter": 300})
            if best is None or r.fun < best.fun: best = r
    res = model(*best.x); return 10 ** best.x[0], 10 ** best.x[1] * KB * T * TPULSE, res.x, res.theta
def regions(xg, yy, xm, ym):
    m = np.interp(xg, xm, ym); out = {}
    for name, sel in (("all", np.ones(len(yy), bool)), ("plateau", yy > 0.9), ("front", (yy > 0.1) & (yy <= 0.9)), ("tail", yy <= 0.1)):
        out[f"rms_{name}"] = float(np.sqrt(np.mean((m[sel] - yy[sel]) ** 2))) if sel.any() else float("nan"); out[f"n_{name}"] = int(sel.sum())
    return out
rows = []
for dataset in ("zenodo", "ylilammi_al2o3"):
    xg, yy = load(dataset); half = float(np.interp(0.5, yy[::-1], xg[::-1])); print(dataset, len(xg), "points, half depth", round(half * 1e6, 1), "um")
    def curve(tag, xm, ym): np.savetxt(ROOT / "figures" / f"lhar_single_curve_{dataset}_{tag}.csv", np.c_[xm * 1e6, ym], delimiter=",", header="x_um,model_norm", comments="")
    s, e, z, th, tab = fit_ballistic(xg, yy); curve("ballistic", z, th); rows.append(dict(dataset=dataset, model="ballistic", factor=float("nan"), s=s, exposure_Pa_s=e, n_points=len(xg), **regions(xg, yy, z, th))); print("  ", rows[-1], flush=True)
    for tag, factor in (("1d_hydraulic", 1.0), ("1d_corr_1mm", float(slit_knudsen_factor(GAP, 1e-3)))):
        s, e, xm, ym = fit_1d(xg, yy, factor); curve(tag, xm, ym); rows.append(dict(dataset=dataset, model=tag, factor=factor, s=s, exposure_Pa_s=e, n_points=len(xg), **regions(xg, yy, xm, ym))); print("  ", rows[-1], flush=True)
with open(ROOT / "figures/lhar_single_exposure_fit.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
# scan-to-scan spread, Zenodo
sheets, grid, prof = zenodo_scans(); srows = []
for sh, pr_ in zip(sheets, prof):
    y = pr_ / np.nanmean(pr_[(grid > 5) & (grid < 60)]); ok = ~np.isnan(y) & (grid <= XMAX_ZENODO); xg, yy = grid[ok] * 1e-6, y[ok]
    if yy.min() > 0.3: continue
    sb, eb, z, th, _ = fit_ballistic(xg, yy); sh_, eh, xm, ym = fit_1d(xg, yy, 1.0)
    srows.append(dict(scan=sh, n_points=len(xg), s_hydraulic=sh_, s_ballistic=sb, ratio=sb / sh_, rms_hydraulic=regions(xg, yy, xm, ym)["rms_all"], rms_ballistic=regions(xg, yy, z, th)["rms_all"])); print("  ", srows[-1], flush=True)
with open(ROOT / "figures/lhar_single_exposure_scans.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(srows[0])); w.writeheader(); w.writerows(srows)
print("done")
