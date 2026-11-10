"""Dose delivered to the wafer, per radius, from the showerhead_pulse OpenFOAM run.

Mesh facts (blockMeshDict): block 0 is the wafer region, 30 axial x 150 radial cells, cells ordered axial-fastest,
so cell i*30 is the wafer-adjacent cell at radial index i (r = (i + 0.5) mm). Axial cell height 0.5 mm, so the
wall-normal distance from the wafer face to the first cell centre is delta = 0.25 mm. With T = 0 on the wafer,
the finite-volume wall flux is  j = DT * T_cell / delta  (per unit showerhead concentration).
"""
import re, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CASE = Path(__file__).resolve().parent / "showerhead_pulse"
OUT = Path(__file__).resolve().parents[1] / "figures"
NX, NR = 30, 150
DELTA = 0.25e-3
DT = float(re.search(r"^DT\s+([0-9.eE+-]+);", (CASE / "constant/transportProperties").read_text(), re.M).group(1))
PULSE = 0.1

def read_internal(path):
    txt = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)", txt, re.S)
    if not m:
        return None
    return np.fromstring(m.group(2), sep="\n")

times = sorted((float(d.name), d) for d in CASE.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
t = np.array([tt for tt, _ in times])
Tw = np.empty((len(times), NR))
for k, (_, d) in enumerate(times):
    T = read_internal(d / "T")
    Tw[k] = 0.0 if T is None else T[np.arange(NR) * NX]
r_mm = np.arange(NR) + 0.5
flux = DT * Tw / DELTA                       # m/s * (dimensionless concentration): molar flux per unit c0
dose = np.trapezoid(flux, t, axis=0)             # m  (flux-time per unit c0); multiply by c0 [mol/m3] for mol/m2

# mass balance: TMA entering during the pulse through the showerhead (per unit c0), full 360 deg
U_in = 0.050; A_sh = np.pi * 0.15 ** 2
n_in = U_in * A_sh * PULSE                   # m3 per unit c0
n_wafer = 2 * np.pi * np.sum(dose * (r_mm * 1e-3)) * 1e-3   # integrate dose over the wafer area (dr = 1 mm)
print(f"DT = {DT:.2e} m2/s, {len(t)} time samples, t_end = {t[-1]:.3f} s")
print(f"peak flux at centre {flux[:, 0].max():.3e}, at r=149 mm {flux[:, -1].max():.3e}  (m/s per unit c0)")
print(f"dose centre {dose[0]:.3e}, mid {dose[NR//2]:.3e}, edge {dose[-1]:.3e}  (m per unit c0)")
print(f"within-wafer dose non-uniformity (max-min)/mean = {(dose.max()-dose.min())/dose.mean()*100:.1f} %")
print(f"TMA injected during pulse {n_in:.3e}, absorbed by wafer {n_wafer:.3e}  -> wafer capture fraction {n_wafer/n_in*100:.1f} %")
i10 = np.argmax(flux[:, 0] > 0.1 * flux[:, 0].max())
tail = t[np.where(flux[:, 0] > 0.01 * flux[:, 0].max())[0][-1]]
print(f"centre flux rises past 10% of peak at t = {t[i10]:.3f} s; falls below 1% of peak at t = {tail:.3f} s (purge needed)")

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ri, lab in ((0, "centre"), (NR // 2, "r = 75 mm"), (NR - 1, "r = 149 mm")):
    axes[0].plot(t, flux[:, ri], label=lab)
axes[0].axvspan(0, PULSE, color="0.9", label="TMA pulse")
axes[0].set_xlabel("time (s)"); axes[0].set_ylabel("wafer flux / c0  (m/s)"); axes[0].set_xlim(0, 1.0)
axes[0].set_title("TMA flux to the wafer: 0.1 s pulse, 300 Pa N2, 15 mm gap, absorbing wafer", fontsize=9)
axes[0].legend(fontsize=8)
axes[1].plot(r_mm, dose / dose.mean(), "k")
axes[1].set_xlabel("radius (mm)"); axes[1].set_ylabel("dose / mean dose")
axes[1].set_title("Within-wafer dose after pulse + purge", fontsize=9)
axes[1].set_ylim(0.9, 1.1); axes[1].axhline(1, color="0.7", lw=0.8)
fig.tight_layout(); fig.savefig(OUT / "reactor_dose_vs_r.png", dpi=150)
np.savetxt(OUT / "reactor_dose_vs_r.csv", np.c_[r_mm, dose], delimiter=",", header="r_mm,dose_m_per_unit_c0", comments="")
print("saved", OUT / "reactor_dose_vs_r.png")
