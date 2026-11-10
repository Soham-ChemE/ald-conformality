"""Compare wafer dose profiles and capture fractions across the sticking-coefficient sweep."""
import re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent / "figures"
NX, NR, DELTA, PULSE = 30, 150, 0.25e-3, 0.1
VBAR = 410.0          # m/s, TMA at 573 K
U_IN, A_SH = 0.050, np.pi * 0.15 ** 2

def read_internal(path):
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)", path.read_text(), re.S)
    return None if not m else np.fromstring(m.group(2), sep="\n")

def read_wafer_face(path):
    """Wafer patch face values (the mixed BC's c_face), needed for the kinetic flux j = k c_face."""
    txt = path.read_text()
    m = re.search(r"\n\s+wafer\s*\{(.*?)\n\s+\}", txt, re.S)
    blk = m.group(1)
    mv = re.search(r"value\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)", blk, re.S)
    if mv:
        return np.fromstring(mv.group(2), sep="\n")
    mu = re.search(r"value\s+uniform\s+([0-9.eE+-]+)", blk)
    return np.full(NR, float(mu.group(1))) if mu else None

def load_case(case, s):
    DT = float(re.search(r"^DT\s+([0-9.eE+-]+);", (case / "constant/transportProperties").read_text(), re.M).group(1))
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times]); flux = np.zeros((len(t), NR))
    k = (s or 0.0) * VBAR / 4
    for i, (_, d) in enumerate(times):
        T = read_internal(d / "T")
        if T is None: continue
        cell = T[np.arange(NR) * NX]
        if s is None:                      # absorbing wall: c_face = 0
            flux[i] = DT * cell / DELTA
        else:
            face = read_wafer_face(d / "T")
            flux[i] = k * face if face is not None else DT * cell / DELTA
    dose = np.trapezoid(flux, t, axis=0)
    r = (np.arange(NR) + 0.5) * 1e-3
    captured = 2 * np.pi * np.sum(dose * r) * 1e-3 / (U_IN * A_SH * PULSE)
    return t, flux, dose, captured

cases = [("absorbing (s -> 1)", ROOT / "showerhead_pulse", None)] + \
        [(f"s = {s:g}", ROOT / f"showerhead_pulse_s{s}", s) for s in (1, 0.1, 0.01, 0.001)]
r_mm = np.arange(NR) + 0.5
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
rows = []
for lab, case, s in cases:
    t, flux, dose, cap = load_case(case, s)
    c = dose[:50].mean()
    tail = t[np.where(flux[:, 0] > 0.01 * flux[:, 0].max())[0][-1]]
    spread140 = (dose[:140].max() - dose[:140].min()) / dose[:140].mean() * 100
    rows.append((lab, cap * 100, c, spread140, tail))
    axes[0].plot(t, flux[:, 0], label=lab)
    axes[1].plot(r_mm, dose / c, label=lab)
    axes[2].plot(r_mm, dose, label=lab)
axes[0].axvspan(0, PULSE, color="0.92"); axes[0].set_xlim(0, 0.6); axes[0].set_xlabel("time (s)"); axes[0].set_ylabel("centre wafer flux / c0 (m/s)")
axes[0].set_title("Centre flux during pulse and purge", fontsize=9); axes[0].legend(fontsize=7)
axes[1].set_ylim(0.85, 1.15); axes[1].axhline(1, color="0.7", lw=0.8); axes[1].set_xlabel("radius (mm)"); axes[1].set_ylabel("dose / centre dose")
axes[1].set_title("Within-wafer dose shape", fontsize=9)
axes[2].set_yscale("log"); axes[2].set_xlabel("radius (mm)"); axes[2].set_ylabel("dose per pulse / c0 (m)"); axes[2].set_title("Absolute dose", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "reactor_sticking_sweep.png", dpi=150)
print(f"{'case':20s} {'captured %':>10s} {'centre dose/c0 (m)':>20s} {'spread r<140 %':>15s} {'flux<1% at (s)':>15s}")
for lab, cap, c, sp, tail in rows:
    print(f"{lab:20s} {cap:10.1f} {c:20.3e} {sp:15.1f} {tail:15.3f}")
print("saved", OUT / "reactor_sticking_sweep.png")
