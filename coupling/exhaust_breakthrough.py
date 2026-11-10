"""Exhaust precursor signal as a saturation end-point, and pattern-aware prediction of it.

From the two-way runs: the precursor mole fraction at the chamber outlet versus time (what a downstream QMS/RGA or
absorption sensor would see), the wafer-average coverage, and the hole-bottom coverage. Questions:
(1) Does the outlet signal show a breakthrough when the wafer saturates, and how far ahead of the hole bottoms?
(2) Can the early-pulse outlet level be used to estimate the wafer's effective area ratio A (loading), so that the
    design rule t = 1.2 (1 + A) N0/(U c0) predicts the pulse end before it happens?
Outlet value: area-weighted mean of T on the outlet patch (read from the written T files).
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "coupling"), str(ROOT / "reactor"), str(ROOT / "src")]
from patterned_lib import load_patterned, _block, _field
N0, U = 4e18, 0.0501
KB, T = 1.380649e-23, 573.15

def outlet_series(case):
    case = Path(case)
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times]); out = np.zeros(len(t))
    for i, (_, d) in enumerate(times):
        txt = (d / "T").read_text()
        blk = _block(txt, "outlet")
        # inletOutlet: 'value' holds the face values (outflow = cell values)
        v = _field(blk, "value", 30)
        out[i] = v.mean()
    return t, out

rows = []
fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
for case in sorted((ROOT / "reactor/map").glob("patterned_AR*_x*")) + sorted((ROOT / "reactor/map").glob("blanket_x*")):
    name = case.name
    m = re.search(r"AR(\d+)_x([0-9.]+)", name); mb = re.search(r"blanket_x([0-9.]+)", name)
    if m: AR, x = int(m.group(1)), float(m.group(2)); A = 1.78e13 * np.pi * 100e-9 * AR * 100e-9
    else: AR, x = 0, float(mb.group(1)); A = 0.0
    c0 = x * 300.0 / (KB * T)
    t, out = outlet_series(case)
    pulse = float(re.search(r"\(\(0 1\) \(([0-9.]+) 1\)", (case / "0/T").read_text()).group(1))
    if m:
        tt, f, prm, _ = load_patterned(case)
        theta_b = f["thetaBottom"]; theta_flat = f["thetaFlat"]
        t90b = tt[np.argmax(theta_b[:, 0] >= 0.9)] if theta_b[:, 0].max() >= 0.9 else np.nan
        wafer_mean_cov = (0.14 * theta_b.mean(1) * 0 + theta_flat.mean(1))  # flat only, display
    else:
        t90b = np.nan
    # breakthrough: time when the outlet signal passes 50 % of its pulse-end plateau
    i_end = np.argmin(abs(t - pulse)); plateau = out[max(0, i_end - 3):i_end + 1].mean()
    t50 = t[np.argmax(out >= 0.5 * plateau)] if plateau > 0 else np.nan
    # early-pulse loading estimate: at t_e the outlet level relative to the blanket outlet at the same time ratio ~ supply split
    t_e = min(0.5, pulse / 3); i_e = np.argmin(abs(t - t_e)); early = out[i_e]
    rows.append(dict(name=name, AR=AR, x=x, A_true=A, pulse=pulse, t90_bottom=t90b, t_break50=t50, early_outlet=early, plateau=plateau,
                     t_rule=1.2 * (1 + A) * N0 / (U * c0)))
    if m and x == 0.01:
        axes[0].plot(t, out, label=f"AR {AR} (A = {A:.0f})"); axes[0].axvline(t90b, color="0.6", lw=0.6)
    if x == 0.01 and not m:
        axes[0].plot(t, out, "k:", label="blanket")
df = pd.DataFrame(rows).sort_values(["x", "AR"])
axes[0].set_xlabel("time (s)"); axes[0].set_ylabel("precursor mole fraction at the outlet / feed"); axes[0].set_xlim(0, 20)
axes[0].set_title("Exhaust breakthrough at 1 % TMA (vertical lines: hole bottoms reach 90 %)", fontsize=9); axes[0].legend(fontsize=7)
p = df[df.AR > 0]
axes[1].plot(p.t90_bottom, p.t_break50, "o", label="outlet 50 % breakthrough")
lim = [0, 1.1 * p.t90_bottom.max()]; axes[1].plot(lim, lim, "k--", lw=0.8); axes[1].set_xlabel("time hole bottoms reach 90 % (s)"); axes[1].set_ylabel("exhaust breakthrough time (s)")
axes[1].set_title("Does the exhaust see the hole bottoms saturate?", fontsize=9); axes[1].legend(fontsize=8)
# early outlet level vs true loading
for x, mk in zip(sorted(df.x.unique()), "osd^"):
    d = df[df.x == x]; axes[2].plot(1 + d.A_true, d.early_outlet, mk, label=f"{x*100:g} % TMA")
axes[2].set_xscale("log"); axes[2].set_xlabel("1 + A_wall/A_wafer (true)"); axes[2].set_ylabel("outlet level early in the pulse / feed")
axes[2].set_title("Early exhaust level encodes the wafer loading", fontsize=9); axes[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(ROOT / "figures/exhaust_breakthrough.png", dpi=150)
df.to_csv(ROOT / "figures/exhaust_breakthrough.csv", index=False)
print(df[["name", "A_true", "t_rule", "t90_bottom", "t_break50", "early_outlet", "plateau"]].round(3).to_string())
print("saved figures/exhaust_breakthrough.png")
