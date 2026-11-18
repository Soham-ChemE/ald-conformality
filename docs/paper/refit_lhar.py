"""Refit a published lateral-channel (LHAR) saturation profile with the 1D slit model under a chosen slit transport:
  dataset: zenodo | ylilammi_al2o3 | ylilammi_tio2 | gao_al2o3
  transport: hyd            hydraulic-diameter Knudsen diffusivity D_K = 2 w vbar / 3 (as in the published fits)
             plates:<L_m>   D_K multiplied by (3/8) ln(L/w) with the cut-off length L in metres (e.g. plates:1e-3)
Same procedure as coupling/fit_pillarhall_corrected.py: modelled domain 250 um closed, 400 cells, 60 substeps, narrowing on,
N0 = 4e18 fixed, nine Nelder-Mead starts over (log s, log n0); the mouth density and s are fitted, the pulse time is the
published one. Appends one row to figures/lhar_refit.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/refit_lhar.py ylilammi_al2o3 plates:1e-3
"""
import sys, csv, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Slit, Precursor, Process, simulate, slit_knudsen_factor
from aldconf.data import load_pillarhall, default_path
KB = 1.380649e-23
dataset, transport = sys.argv[1], sys.argv[2]
DOMAIN = 250e-6
if dataset == "zenodo":
    df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
    sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
    for sh in sheets:
        g = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
    mean = np.nanmean(np.array(prof), 0); plateau = np.nanmean(mean[(grid > 5) & (grid < 60)]); y = mean / plateau
    GAP, CYCLES, TPULSE, T, MASS = 500e-9, 500, 0.1, 573.15, 72.09; L_PHYS_NOTE = "0.5-5 mm"
elif dataset == "ylilammi_al2o3":
    d = pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv"); grid = d.x_um.values; th = d.thickness_nm.values
    plateau = th[(grid > 5) & (grid < 60)].mean(); y = th / plateau
    GAP, CYCLES, TPULSE, T, MASS = 500e-9, 500, 0.1, 573.15, 72.09; L_PHYS_NOTE = "1 mm"
elif dataset == "ylilammi_tio2":
    d = pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig7_tio2.csv"); grid = d.x_um.values; th = d.thickness_nm.values
    plateau = th[(grid >= 0) & (grid < 20)].mean(); y = th / plateau
    GAP, CYCLES, TPULSE, T, MASS = 500e-9, 1000, 0.1, 383.15, 189.68; L_PHYS_NOTE = "1 mm"
elif dataset == "gao_al2o3":
    d = pd.read_csv(ROOT / "data/lhar_profiles/gao2015_fig5a_al2o3.csv"); grid = d.x_um.values; th = d.thickness_nm.values
    plateau = th[(grid >= 10) & (grid <= 50)].mean(); y = th / plateau
    # pulse time not stated in Gao 2015: 0.1 s is nominal, so only the exposure (p x t) is meaningful, not p
    GAP, CYCLES, TPULSE, T, MASS = 500e-9, 1109, 0.1, 573.15, 72.09; L_PHYS_NOTE = "first-generation LHAR"
else:
    raise SystemExit("unknown dataset")
GPC = plateau * 1e-9 / CYCLES
if transport == "hyd": factor = 1.0; L_cut = float("nan")
else: L_cut = float(transport.split(":")[1]); factor = float(slit_knudsen_factor(GAP, L_cut))
print(f"{dataset}: plateau {plateau:.2f} nm, gpc {GPC*1e12:.1f} pm, {len(grid)} points; transport {transport}, factor {factor:.3f}", flush=True)
ok = ~np.isnan(y)
def model(log_s, log_n0):
    pr = Precursor(sticking=10 ** log_s, site_density=4e18, mass_amu=MASS, gpc=GPC)
    res = simulate(Process(Slit(GAP, DOMAIN, knudsen_correction="none"), pr, T=T, p_partial=10 ** log_n0 * KB * T, t_dose=TPULSE, cycles=CYCLES, n_cells=400, n_substeps=60, knudsen_factor=factor, narrowing=True))
    th = res.thickness / res.thickness[:10].mean(); return np.interp(grid * 1e-6, res.x, th)
def cost(p): m = model(*p); return float(np.mean((m[ok] - y[ok]) ** 2))
t0 = time.time(); best = None
for s0 in (-3.0, -2.3, -1.7):
    for n0 in (19.5, 20.0, 20.5):
        r = minimize(cost, [s0, n0], method="Nelder-Mead", options={"xatol": 0.01, "fatol": 1e-7, "maxiter": 200})
        print(f"   start ({s0}, {n0}): s = {10**r.x[0]:.3e}, p = {10**r.x[1]*KB*T:.2f} Pa, rms {np.sqrt(r.fun):.4f}", flush=True)
        if best is None or r.fun < best.fun: best = r
s_fit, p_fit = 10 ** best.x[0], 10 ** best.x[1] * KB * T
half = float(np.interp(0.5, y[ok][::-1], grid[ok][::-1])) if y[ok].min() < 0.5 else float("nan")
row = dict(dataset=dataset, transport=transport, L_cut_m=L_cut, factor=factor, plateau_nm=plateau, gpc_pm=GPC * 1e12, half_depth_um=half, s=s_fit, p_Pa=p_fit, exposure_Pa_s=p_fit * TPULSE, rms=float(np.sqrt(best.fun)), n_points=int(ok.sum()), wall_s=time.time() - t0)
out = ROOT / "figures" / "lhar_refit.csv"; new = not out.exists()
with open(out, "a", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(row))
    if new: w.writeheader()
    w.writerow(row)
np.savetxt(ROOT / "figures" / f"lhar_refit_curve_{dataset}_{transport.replace(':','_')}.csv", np.c_[grid, y, model(*best.x)], delimiter=",", header="x_um,measured_norm,model_norm", comments="")
print(row, flush=True)
