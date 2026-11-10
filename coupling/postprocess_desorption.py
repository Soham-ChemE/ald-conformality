"""Water half-cycle with reversible adsorption: coverage and outlet water vs time for kdes = 0, 0.1, 1, 10 1/s
(blanket wafer, 5 s H2O pulse at 1 %, then purge). Outlet signal = what a downstream sensor sees during the purge."""
import re, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path[:0] = [str(ROOT / "coupling"), str(ROOT / "reactor")]
from postprocess_langmuir import load as load_flat
from patterned_lib import _block, _field
def outlet(case):
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times]); o = np.zeros(len(t))
    for i, (_, d) in enumerate(times): o[i] = _field(_block((d / "T").read_text(), "outlet"), "value", 30).mean()
    return t, o
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2)); rows = []
for k in (0, 0.1, 1, 10):
    case = ROOT / f"reactor/showerhead_h2o_blanket_kdes{k:g}"
    t, th, face, s, vbar, c0, N0, pulse = load_flat(case); to, o = outlet(case)
    i_end = np.argmin(abs(t - pulse))
    th_end = th[i_end, 0]; th_final = th[-1, 0]
    t90 = t[np.argmax(th[:, 0] >= 0.9)] if th[:, 0].max() >= 0.9 else np.nan
    # purge metric: time after the pulse for the outlet water to fall below 1 % of its pulse-end level
    o_end = o[np.argmin(abs(to - pulse))]; after = to > pulse + 0.02
    below = np.where(after & (o < 0.01 * max(o_end, 1e-9)))[0]; t_purge = (to[below[0]] - pulse) if len(below) else np.nan
    rows.append((k, t90, th_end, th_final, t_purge, o_end))
    ax[0].plot(t, th[:, 0], label=f"kdes = {k:g} 1/s"); ax[1].plot(to, o, label=f"kdes = {k:g} 1/s"); ax[2].semilogy(to - pulse, np.maximum(o, 1e-6), label=f"kdes = {k:g} 1/s")
    print(f"kdes={k:5g} 1/s: t90 {t90:.2f} s, coverage at pulse end {th_end:.3f}, 3 s after {th_final:.3f}, outlet at pulse end {o_end:.3f} c0, purge to 1 % of that: {t_purge:.2f} s")
for a in ax[:2]: a.axvspan(0, 5, color="0.92")
ax[0].set_xlabel("time (s)"); ax[0].set_ylabel("surface coverage, centre"); ax[0].set_title("Water coverage: adsorption vs desorption", fontsize=9); ax[0].legend(fontsize=8)
ax[1].set_xlabel("time (s)"); ax[1].set_ylabel("outlet H2O / feed"); ax[1].set_title("Outlet water during pulse and purge", fontsize=9); ax[1].legend(fontsize=8)
ax[2].set_xlabel("time after pulse end (s)"); ax[2].set_ylabel("outlet H2O / feed (log)"); ax[2].set_xlim(0, 3); ax[2].set_ylim(1e-5, 1.5); ax[2].set_title("Purge tail: desorbing water keeps the exhaust wet", fontsize=9); ax[2].legend(fontsize=8); ax[2].grid(True, which="both", alpha=0.3)
fig.tight_layout(); fig.savefig(ROOT / "figures/h2o_desorption.png", dpi=150); print("saved figures/h2o_desorption.png")
