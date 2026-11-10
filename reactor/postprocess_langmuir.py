"""Coverage and uptake on the self-limiting wafer (langmuirSticking BC) versus radius and time."""
import re, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / "figures"
NR = 150
r_mm = np.arange(NR) + 0.5

def patch_block(txt, name="wafer"):
    m = re.search(r"\n\s+%s\s*\{(.*?)\n\s+\}" % name, txt, re.S)
    if m: return m.group(1)
    m = re.search(r"\n\s+%s\s*\{([^}]*)\}" % name, txt, re.S)   # one-line block in the initial 0/T
    return m.group(1) if m else ""

def list_or_uniform(blk, key, n):
    m = re.search(r"\b%s\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)" % key, blk, re.S)
    if m: return np.fromstring(m.group(2), sep="\n")
    m = re.search(r"\b%s\s+uniform\s+([0-9.eE+-]+)" % key, blk)
    return np.full(n, float(m.group(1))) if m else None

def scalar_entry(blk, key):
    return float(re.search(r"\b%s\s+([0-9.eE+-]+);" % key, blk).group(1))

def load(case):
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times]); theta = np.zeros((len(t), NR)); face = np.zeros((len(t), NR))
    s = vbar = c0 = N0 = None
    for i, (_, d) in enumerate(times):
        blk = patch_block((d / "T").read_text())
        if s is None and "langmuirSticking" in blk:
            s, vbar, c0, N0 = (scalar_entry(blk, k) for k in ("s", "vbar", "c0", "N0"))
        th = list_or_uniform(blk, "theta", NR); va = list_or_uniform(blk, "value", NR)
        theta[i] = 0 if th is None else th; face[i] = 0 if va is None else va
    pulse = float(re.search(r"\(\(0 1\) \(([0-9.]+) 1\)", (case / "0/T").read_text()).group(1))
    return t, theta, face, s, vbar, c0, N0, pulse

cases = sorted(ROOT.glob("showerhead_pulse_langmuir_p*"))
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
for case in cases:
    t, theta, face, s, vbar, c0, N0, pulse = load(case)
    uptake = (vbar / 4) * s * (1 - theta) * face * c0          # molecules / m2 / s
    lab = f"pulse {pulse:g} s"
    axes[0].plot(t, theta[:, 0], label=f"{lab}, centre"); axes[0].plot(t, theta[:, 140], "--", label=f"{lab}, r = 140 mm")
    axes[1].plot(r_mm, theta[-1], label=lab)
    # time to reach 90 % coverage per radius (nan where never reached)
    t90 = np.array([t[np.argmax(theta[:, j] >= 0.9)] if theta[:, j].max() >= 0.9 else np.nan for j in range(NR)])
    axes[2].plot(r_mm, t90, label=lab)
    supply = 0.050 * c0                                            # molecules/m2/s delivered per unit area at the showerhead
    print(f"{case.name}: s={s:g} c0={c0:.3g} N0={N0:.3g}; supply flux U c0 = {supply:.3g} /m2/s -> monolayer in {N0/supply:.3f} s")
    print(f"   final coverage: centre {theta[-1,0]:.3f}, r=100 {theta[-1,100]:.3f}, r=140 {theta[-1,140]:.3f}, r=149 {theta[-1,149]:.3f}")
    print(f"   t90 centre {t90[0]:.3f} s, r=140 {t90[140]:.3f} s, worst on wafer {np.nanmax(t90[:149]):.3f} s; peak uptake at centre {uptake[:,0].max():.3g} /m2/s")
axes[0].set_xlabel("time (s)"); axes[0].set_ylabel("coverage theta"); axes[0].set_title("Surface coverage vs time", fontsize=9); axes[0].legend(fontsize=7)
axes[1].set_xlabel("radius (mm)"); axes[1].set_ylabel("final coverage"); axes[1].set_ylim(0, 1.05); axes[1].set_title("Coverage across the wafer after the pulse", fontsize=9); axes[1].legend(fontsize=8)
axes[2].set_xlabel("radius (mm)"); axes[2].set_ylabel("time to 90 % coverage (s)"); axes[2].set_title("Saturation time across the wafer", fontsize=9); axes[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "reactor_langmuir_coverage.png", dpi=150); print("saved", OUT / "reactor_langmuir_coverage.png")
