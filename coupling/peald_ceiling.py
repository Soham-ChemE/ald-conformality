"""Plasma ALD: radical recombination sets a ceiling on reachable aspect ratio.
Feature solver at constant mouth density (0.3 Pa-equivalent radical density, 573 K, O-atom mass 16), sticking s on open
sites, loss probability r on every wall collision. Metrics: time to 90 % bottom coverage vs r for AR 20/50/100, and the
largest AR whose bottom reaches 90 % within a fixed 5 s exposure, vs r. r values from Arts et al. 2019/2021: O atoms
recombine with r ~ 6e-5 on SiO2 and TiO2, higher on Al2O3 and HfO2 (material-dependent)."""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Hole, Precursor, Process, simulate, step_coverage, dose_for_coverage

D, T, P, S = 100e-9, 573.15, 0.3, 0.05        # radicals stick readily on the methyl-terminated surface; s = 0.05 assumed
RS = np.array([0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
ARS = (20, 50, 100)
def proc(AR, r, t):
    pr = Precursor(name="O", sticking=S, recombination=r, mass_amu=16.0, site_density=4e18)
    return Process(Hole(D, AR * D), pr, T=T, p_partial=P, t_dose=t, n_cells=200, n_substeps=150, narrowing=False)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
rows = []
for AR in ARS:
    ts = []
    for r in RS:
        t90 = dose_for_coverage(proc(AR, r, 1.0), target=0.90, t_lo=1e-3, t_hi=300.0)
        ts.append(t90); rows.append((AR, r, t90))
    ts = np.array(ts)
    axes[0].plot(np.where(RS > 0, RS, 3e-6), ts, "o-", label=f"AR {AR}")
    print(f"AR {AR}: t90 bottom (s) vs r:", " ".join(f"{r:g}:{t:.2f}" for r, t in zip(RS, ts)))
axes[0].set_xscale("log"); axes[0].set_yscale("log"); axes[0].set_xlabel("radical recombination probability r (leftmost point: r = 0)")
axes[0].set_ylabel("time to 90 % bottom coverage (s)"); axes[0].set_title("Plasma ALD: pulse time vs wall recombination, 100 nm hole", fontsize=9); axes[0].legend(fontsize=8); axes[0].grid(True, which="both", alpha=0.3)
# largest AR reachable within a fixed exposure
TEXP = 5.0
def reach(r):
    def f(AR): return step_coverage(simulate(proc(AR, r, TEXP))) - 0.9
    lo, hi = 5.0, 2000.0
    if f(hi) > 0: return hi
    if f(lo) < 0: return np.nan
    return brentq(f, lo, hi, xtol=1.0)
rr = np.logspace(-5, -1.5, 8)
armax = [reach(r) for r in rr]
axes[1].plot(rr, armax, "s-", color="tab:red"); axes[1].axvspan(5e-5, 8e-5, color="0.9", label="O on SiO2 / TiO2 (Arts 2019)")
axes[1].set_xscale("log"); axes[1].set_yscale("log"); axes[1].set_xlabel("radical recombination probability r"); axes[1].set_ylabel(f"largest AR with 90 % bottom coverage in {TEXP:g} s")
axes[1].set_title("The recombination ceiling", fontsize=9); axes[1].legend(fontsize=8); axes[1].grid(True, which="both", alpha=0.3)
print("largest AR reachable in 5 s vs r:", " ".join(f"{r:.1e}:{a:.0f}" for r, a in zip(rr, armax)))
fig.tight_layout(); fig.savefig(ROOT / "figures/peald_ceiling.png", dpi=150); print("saved figures/peald_ceiling.png")
