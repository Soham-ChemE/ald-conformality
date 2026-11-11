"""Refit of the TMA sticking coefficient and mouth exposure to the Zenodo 3974438 PillarHall profiles with the
corrected slit transport (Slit(knudsen_correction="plates"): D_slit = (3/8) ln(L/w) x 2 w vbar/3, paper/wp/wp9b_slit_transport.md).

Same data, grid, model and fitting procedure as coupling/fit_pillarhall.py (Slit 500 nm x 250 um, 500 cycles with
channel narrowing, N0 = 4e18 fixed, nine Nelder-Mead starts, then one refit per sheet for the sheet-to-sheet spread);
only the slit diffusivity differs. The uncorrected fit (s = 8.5e-3, 75 Pa) is kept in figures/pillarhall_fit.png.
Output: figures/pillarhall_fit_corrected.png, figures/pillarhall_fit_corrected.csv."""
import sys, csv, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Slit, Precursor, Process, simulate, slit_knudsen_factor
from aldconf.data import load_pillarhall, default_path

CORRECTION = "plates"
df = load_pillarhall(default_path())
good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]
grid = np.arange(0, 180, 2.0)
prof = []
for sh in sheets:
    g = good[good.sheet == sh].sort_values("depth_um")
    prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
prof = np.array(prof); mean = np.nanmean(prof, 0); std = np.nanstd(prof, 0)
plateau = np.nanmean(mean[(grid > 5) & (grid < 60)])
y = mean / plateau
half = np.interp(0.5, y[::-1], grid[::-1])
print(f"{len(sheets)} sheets averaged; plateau {plateau:.2f} nm -> gpc {plateau/500*1000:.1f} pm/cycle; half-thickness depth {half:.1f} um (AR {half*2:.0f})")

GAP, L, CYCLES, TPULSE, T = 500e-9, 250e-6, 500, 0.1, 573.15
GPC = plateau * 1e-9 / CYCLES
KB = 1.380649e-23
print(f"slit AR {L/GAP:.0f}: parallel-plate factor {float(slit_knudsen_factor(GAP, L)):.3f} on D_K = 2 w vbar/3")
def model(log_s, log_n0, n_cells=400, substeps=60, correction=CORRECTION):
    pr = Precursor(sticking=10 ** log_s, site_density=4e18, mass_amu=72.09, gpc=GPC)
    n0 = 10 ** log_n0
    proc = Process(Slit(GAP, L, knudsen_correction=correction), pr, T=T, p_partial=n0 * KB * T, t_dose=TPULSE, cycles=CYCLES,
                   n_cells=n_cells, n_substeps=substeps, narrowing=True)
    res = simulate(proc)
    th = res.thickness / res.thickness[: 10].mean()
    return np.interp(grid * 1e-6, res.x, th)
def cost(p):
    m = model(*p); ok = ~np.isnan(y)
    return float(np.mean((m[ok] - y[ok]) ** 2))
t0 = time.time(); best = None; n_eval = 0
for s0 in (-3.0, -2.3, -1.7):
    for n0 in (19.5, 20.0, 20.5):
        r = minimize(cost, [s0, n0], method="Nelder-Mead", options={"xatol": 0.01, "fatol": 1e-7, "maxiter": 200})
        n_eval += r.nfev
        print(f"   start ({s0}, {n0}): s = {10**r.x[0]:.2e}, p = {10**r.x[1]*KB*T:.2f} Pa, rms {np.sqrt(r.fun):.4f} ({r.nfev} evaluations)", flush=True)
        if best is None or r.fun < best.fun: best = r
log_s, log_n0 = best.x
s_fit, n0_fit = 10 ** log_s, 10 ** log_n0
p_fit = n0_fit * KB * T
print(f"fit (corrected transport): s = {s_fit:.2e}, n0 = {n0_fit:.2e} m^-3 (TMA partial pressure {p_fit:.2f} Pa), exposure {p_fit*TPULSE:.3f} Pa s = {p_fit*TPULSE/1.333e-4:.0f} L per pulse, rms {np.sqrt(best.fun):.4f}  [{n_eval} evaluations, {time.time()-t0:.0f} s]", flush=True)
rows = [dict(fit="joint_corrected", correction=CORRECTION, s=s_fit, n0_m3=n0_fit, p_TMA_Pa=p_fit, exposure_Pa_s=p_fit * TPULSE, rms=np.sqrt(best.fun))]
fits = []
for i, sh in enumerate(sheets):
    yi = prof[i] / np.nanmean(prof[i][(grid > 5) & (grid < 60)])
    def ci(p):
        m = model(*p); ok = ~np.isnan(yi); return float(np.mean((m[ok] - yi[ok]) ** 2))
    ri = minimize(ci, best.x, method="Nelder-Mead", options={"xatol": 0.02, "fatol": 1e-6, "maxiter": 120})
    fits.append(ri.x)
    rows.append(dict(fit=sh, correction=CORRECTION, s=10 ** ri.x[0], n0_m3=10 ** ri.x[1], p_TMA_Pa=10 ** ri.x[1] * KB * T, exposure_Pa_s=10 ** ri.x[1] * KB * T * TPULSE, rms=np.sqrt(ri.fun)))
    print(f"   {sh}: s = {10**ri.x[0]:.2e}, p = {10**ri.x[1]*KB*T:.2f} Pa, rms {np.sqrt(ri.fun):.4f}", flush=True)
fits = np.array(fits)
print(f"sheet-to-sheet: s = {10**fits[:,0].mean():.2e} (x/÷ {10**fits[:,0].std():.2f}), p = {10**fits[:,1].mean()*KB*T:.2f} Pa (x/÷ {10**fits[:,1].std():.2f})")
rows.append(dict(fit="sheet_to_sheet_geometric_mean", correction=CORRECTION, s=10 ** fits[:, 0].mean(), n0_m3=10 ** fits[:, 1].mean(), p_TMA_Pa=10 ** fits[:, 1].mean() * KB * T,
                 exposure_Pa_s=10 ** fits[:, 1].mean() * KB * T * TPULSE, rms=float("nan")))
rows.append(dict(fit="sheet_to_sheet_spread_factor", correction=CORRECTION, s=10 ** fits[:, 0].std(), n0_m3=10 ** fits[:, 1].std(), p_TMA_Pa=10 ** fits[:, 1].std(), exposure_Pa_s=10 ** fits[:, 1].std(), rms=float("nan")))
# the uncorrected model at the corrected fit's parameters and at its own published fit, for the figure
rows.append(dict(fit="uncorrected_published", correction="none", s=8.5e-3, n0_m3=75.0 / (KB * T), p_TMA_Pa=75.0, exposure_Pa_s=7.5, rms=0.038))
with open(ROOT / "figures/pillarhall_fit_corrected.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"wall time {time.time()-t0:.0f} s")
fig, ax = plt.subplots(figsize=(7, 4))
ax.errorbar(grid * 2, y, yerr=std / plateau, fmt=".", ms=3, color="k", alpha=0.6, label=f"Zenodo 3974438, mean of {len(sheets)} profiles")
ax.plot(grid * 2, model(*best.x), "r-", label=f"corrected slit transport: s = {s_fit:.1e}, p_TMA = {p_fit:.1f} Pa")
ax.plot(grid * 2, model(np.log10(8.5e-3), np.log10(75.0 / (KB * T)), correction="none"), "b--", lw=1.0, label="uncorrected fit (D_K = 2 w vbar/3): s = 8.5e-3, 75 Pa")
for ls, lab in ((-3, "s = 1e-3"), (-2, "s = 1e-2")):
    ax.plot(grid * 2, model(ls, log_n0), ":", lw=0.8, label=f"{lab}, corrected, same exposure")
ax.set_xlabel("depth / gap (aspect ratio)"); ax.set_ylabel("thickness / plateau"); ax.set_ylim(-0.05, 1.1)
ax.set_title("Al2O3 TMA/H2O 300 C, 500 cycles x 0.1 s, 500 nm channel: refit with the parallel-plate slit diffusivity", fontsize=8); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(ROOT / "figures/pillarhall_fit_corrected.png", dpi=150); print("saved figures/pillarhall_fit_corrected.png, .csv")
