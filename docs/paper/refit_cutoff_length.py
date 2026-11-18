"""Sensitivity of the corrected PillarHall fit to the cut-off length in the slit factor (3/8) ln(L/w).
Same data, grid, model, domain (250 um, 400 cells) and nine-start procedure as coupling/fit_pillarhall_corrected.py,
but the slit diffusivity is multiplied by the factor evaluated at a chosen physical channel length (0.5, 1 or 5 mm)
through Process(knudsen_factor=...), instead of at the 250 um modelled domain. Writes figures/pillarhall_fit_cutoff.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/refit_cutoff_length.py 1e-3
"""
import sys, csv, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Slit, Precursor, Process, simulate, slit_knudsen_factor
from aldconf.data import load_pillarhall, default_path
L_PHYS = float(sys.argv[1])
df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
for sh in sheets:
    g = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
prof = np.array(prof); mean = np.nanmean(prof, 0); plateau = np.nanmean(mean[(grid > 5) & (grid < 60)]); y = mean / plateau
GAP, L, CYCLES, TPULSE, T, KB = 500e-9, 250e-6, 500, 0.1, 573.15, 1.380649e-23; GPC = plateau * 1e-9 / CYCLES
factor = float(slit_knudsen_factor(GAP, L_PHYS)); print(f"L_phys {L_PHYS*1e3:.1f} mm: factor {factor:.3f}", flush=True)
def model(log_s, log_n0):
    pr = Precursor(sticking=10 ** log_s, site_density=4e18, mass_amu=72.09, gpc=GPC)
    res = simulate(Process(Slit(GAP, L, knudsen_correction="none"), pr, T=T, p_partial=10 ** log_n0 * KB * T, t_dose=TPULSE, cycles=CYCLES, n_cells=400, n_substeps=60, knudsen_factor=factor, narrowing=True))
    th = res.thickness / res.thickness[:10].mean(); return np.interp(grid * 1e-6, res.x, th)
ok = ~np.isnan(y)
def cost(p): m = model(*p); return float(np.mean((m[ok] - y[ok]) ** 2))
t0 = time.time(); best = None
for s0 in (-2.3, -1.7, -1.4):
    for n0 in (19.5, 20.0, 20.5):
        r = minimize(cost, [s0, n0], method="Nelder-Mead", options={"xatol": 0.01, "fatol": 1e-7, "maxiter": 200})
        print(f"   start ({s0}, {n0}): s = {10**r.x[0]:.3e}, p = {10**r.x[1]*KB*T:.2f} Pa, rms {np.sqrt(r.fun):.4f}", flush=True)
        if best is None or r.fun < best.fun: best = r
s_fit, p_fit = 10 ** best.x[0], 10 ** best.x[1] * KB * T
row = dict(L_phys_m=L_PHYS, factor=factor, s=s_fit, p_TMA_Pa=p_fit, exposure_Pa_s=p_fit * TPULSE, rms=float(np.sqrt(best.fun)), wall_s=time.time() - t0)
out = ROOT / "figures" / "pillarhall_fit_cutoff.csv"; new = not out.exists()
with open(out, "a", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(row)); 
    if new: w.writeheader()
    w.writerow(row)
print(row, flush=True)
