"""Fit the TMA sticking coefficient and the exposure (mouth density) to the Zenodo 3974438 PillarHall profiles.

Experiment (workbook METADATA): Al2O3, TMA/H2O, 300 C, 500 cycles, 0.1 s TMA pulse, 4 s purge, ~300 Pa N2, LHAR3
lateral channel of 500 nm gap. The TMA partial pressure is not given, so the mouth density n0 is fitted together
with s; the front slope constrains s, the penetration depth constrains the exposure n0 * t_pulse (Arts 2019).
Model: aldconf.simulate, Slit(500 nm), 500 cycles with channel narrowing (gpc from the measured plateau thickness),
N0 = 4e18 fixed (n0 and N0 enter only as a ratio in the dose)."""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Slit, Precursor, Process, simulate
from aldconf.data import load_pillarhall, default_path

df = load_pillarhall(default_path())
good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
# average the long-channel sheets (0.5, 1, 5 mm) on a common depth grid
sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]
grid = np.arange(0, 180, 2.0)
prof = []
for sh in sheets:
    g = good[good.sheet == sh].sort_values("depth_um")
    prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
prof = np.array(prof); mean = np.nanmean(prof, 0); std = np.nanstd(prof, 0)
plateau = np.nanmean(mean[(grid > 5) & (grid < 60)])
y = mean / plateau
print(f"{len(sheets)} sheets averaged; plateau {plateau:.2f} nm -> gpc {plateau/500*1000:.1f} pm/cycle; half-thickness depth {np.interp(0.5, y[::-1], grid[::-1]):.1f} um (AR {np.interp(0.5, y[::-1], grid[::-1])*2:.0f})")

GAP, L, CYCLES, TPULSE, T = 500e-9, 250e-6, 500, 0.1, 573.15
GPC = plateau * 1e-9 / CYCLES
def model(log_s, log_n0, n_cells=400, substeps=60):
    pr = Precursor(sticking=10 ** log_s, site_density=4e18, mass_amu=72.09, gpc=GPC)
    n0 = 10 ** log_n0
    p_partial = n0 * 1.380649e-23 * T
    proc = Process(Slit(GAP, L), pr, T=T, p_partial=p_partial, t_dose=TPULSE, cycles=CYCLES, n_cells=n_cells, n_substeps=substeps, narrowing=True)
    res = simulate(proc)
    th = res.thickness / res.thickness[: 10].mean()
    return np.interp(grid * 1e-6, res.x, th)
def cost(p):
    m = model(*p); ok = ~np.isnan(y)
    return float(np.mean((m[ok] - y[ok]) ** 2))
best = None
for s0 in (-3.0, -2.3, -1.7):
    for n0 in (19.5, 20.0, 20.5):
        r = minimize(cost, [s0, n0], method="Nelder-Mead", options={"xatol": 0.01, "fatol": 1e-7, "maxiter": 200})
        if best is None or r.fun < best.fun: best = r
log_s, log_n0 = best.x
s_fit, n0_fit = 10 ** log_s, 10 ** log_n0
p_fit = n0_fit * 1.380649e-23 * T
print(f"fit: s = {s_fit:.2e}, n0 = {n0_fit:.2e} m^-3 (TMA partial pressure {p_fit:.2f} Pa), exposure {p_fit*TPULSE:.3f} Pa s = {p_fit*TPULSE/1.333e-4:.0f} L per pulse, rms {np.sqrt(best.fun):.4f}")
# crude uncertainty: profile-to-profile scatter -> refit on each sheet
fits = []
for i, sh in enumerate(sheets):
    yi = prof[i] / np.nanmean(prof[i][(grid > 5) & (grid < 60)])
    def ci(p):
        m = model(*p); ok = ~np.isnan(yi); return float(np.mean((m[ok] - yi[ok]) ** 2))
    ri = minimize(ci, best.x, method="Nelder-Mead", options={"xatol": 0.02, "fatol": 1e-6, "maxiter": 120})
    fits.append(ri.x); print(f"   {sh}: s = {10**ri.x[0]:.2e}, p = {10**ri.x[1]*1.380649e-23*T:.2f} Pa")
fits = np.array(fits)
print(f"sheet-to-sheet: s = {10**fits[:,0].mean():.2e} (x/÷ {10**fits[:,0].std():.2f}), p = {10**fits[:,1].mean()*1.380649e-23*T:.2f} Pa (x/÷ {10**fits[:,1].std():.2f})")
fig, ax = plt.subplots(figsize=(7, 4))
ax.errorbar(grid * 2 / 1.0, y, yerr=std / plateau, fmt=".", ms=3, color="k", alpha=0.6, label=f"Zenodo 3974438, mean of {len(sheets)} profiles")
ax.plot(grid * 2, model(*best.x), "r-", label=f"fit: s = {s_fit:.1e}, p_TMA = {p_fit:.2f} Pa")
for ls, lab in ((-3, "s = 1e-3"), (-2, "s = 1e-2")):
    ax.plot(grid * 2, model(ls, log_n0), "--", lw=0.8, label=f"{lab}, same exposure")
ax.set_xlabel("depth / gap (aspect ratio)"); ax.set_ylabel("thickness / plateau"); ax.set_ylim(-0.05, 1.1)
ax.set_title("Al2O3 TMA/H2O 300 C, 500 cycles x 0.1 s: saturation profile in a 500 nm lateral channel", fontsize=9); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(ROOT / "figures/pillarhall_fit.png", dpi=150); print("saved figures/pillarhall_fit.png")
