"""Patterned-wafer (two-way) results versus the one-way coupling and the flat wafer."""
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "coupling"), str(ROOT / "reactor"), str(ROOT / "src")]
from patterned_lib import load_patterned
from postprocess_langmuir import load as load_flat
OUT = ROOT / "figures"
r_mm = np.arange(150) + 0.5
U, A = 0.050, np.pi * 0.15 ** 2

cases = {"AR 100": ROOT / "reactor/showerhead_pulse_patterned_AR100", "AR 50": ROOT / "reactor/showerhead_pulse_patterned_AR50"}
tf, thf, facef, sf, vbarf, c0f, N0f, pulsef = load_flat(ROOT / "reactor/sweep/base/pulse")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
axes[0].plot(tf, facef[:, 0], "k:", label="unpatterned wafer")
for lab, cd in cases.items():
    t, f, prm, pulse = load_patterned(cd)
    phi = prm["holeDensity"] * np.pi * prm["d"] ** 2 / 4
    wall_ratio = prm["holeDensity"] * np.pi * prm["d"] * prm["AR"] * prm["d"]
    def t90(arr, j): 
        k = np.argmax(arr[:, j] >= 0.9); return t[k] if arr[:, j].max() >= 0.9 else np.nan
    # TMA captured by the wafer up to each time (per unit wafer area), from Keff * c_face
    uptake = f["Keff"] * f["value"] * prm["c0"]
    absorbed = 2 * np.pi * np.sum((uptake * np.diff(t, prepend=0)[:, None]).sum(0) * r_mm * 1e-3) * 1e-3
    injected = U * prm["c0"] * A * pulse
    print(f"== {lab}: open area {phi*100:.1f} %, wall area / wafer area = {wall_ratio:.0f}, pulse {pulse:g} s, run to {t[-1]:g} s")
    print(f"   flat-area t90: centre {t90(f['thetaFlat'],0):.2f} s, r=140 {t90(f['thetaFlat'],140):.2f} s   (unpatterned wafer: 0.35 / 0.38 s)")
    print(f"   hole bottom t90: centre {t90(f['thetaBottom'],0):.2f} s, r=100 {t90(f['thetaBottom'],100):.2f}, r=140 {t90(f['thetaBottom'],140):.2f}, r=149 {t90(f['thetaBottom'],149):.2f} s")
    print(f"   hole bottom coverage at end of pulse ({pulse:g} s): centre {f['thetaBottom'][np.argmin(abs(t-pulse)),0]:.3f}, r=140 {f['thetaBottom'][np.argmin(abs(t-pulse)),140]:.3f}; at end of run {f['thetaBottom'][-1,0]:.3f}")
    print(f"   surface TMA density / c0 at 0.3 s: {f['value'][np.argmin(abs(t-0.3)),0]:.3f} (unpatterned {facef[np.argmin(abs(tf-0.3)),0]:.3f}); captured {absorbed/injected*100:.0f} % of the injected TMA")
    axes[0].plot(t, f["value"][:, 0], label=f"patterned, {lab}")
    axes[1].plot(t, f["thetaBottom"][:, 0], label=f"{lab} hole bottom, centre"); axes[1].plot(t, f["thetaBottom"][:, 140], "--", label=f"{lab} hole bottom, r = 140 mm")
    axes[2].plot(r_mm, [t90(f["thetaBottom"], j) for j in range(150)], label=f"{lab} hole bottom")
axes[1].plot(tf, thf[:, 0], "k:", label="unpatterned flat wafer")
axes[0].axvspan(0, 1.0, color="0.92"); axes[0].set_xlim(0, 3); axes[0].set_xlabel("time (s)"); axes[0].set_ylabel("TMA density at wafer surface / c0")
axes[0].set_title("Loading: a patterned wafer depresses the surface precursor density", fontsize=9); axes[0].legend(fontsize=7)
axes[1].axvspan(0, 1.0, color="0.92"); axes[1].set_xlabel("time (s)"); axes[1].set_ylabel("coverage"); axes[1].set_ylim(0, 1.05)
axes[1].set_title("Hole-bottom coverage, two-way coupled", fontsize=9); axes[1].legend(fontsize=7)
axes[2].set_xlabel("radius (mm)"); axes[2].set_ylabel("time to 90 % coverage at the hole bottom (s)"); axes[2].set_title("Saturation time of the hole bottom across the wafer", fontsize=9); axes[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "coupling_patterned_wafer.png", dpi=150); print("saved", OUT / "coupling_patterned_wafer.png")
