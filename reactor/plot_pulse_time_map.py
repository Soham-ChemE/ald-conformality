"""Headline figure: saturating pulse time versus TMA supply and aspect ratio, blanket vs patterned wafer, with the
pure supply-limit line (1 + wall area ratio) x N0 / (U c0)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = Path(__file__).resolve().parents[1] / "figures"
df = pd.read_csv(OUT / "pulse_time_map.csv")
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
ax = axes[0]
b = df[df.kind == "blanket"].sort_values("x")
ax.plot(b.x * 100, b.t90_worst, "ko-", label="blanket wafer (flat surface)")
ax.plot(b.x * 100, b.t_supply, "k:", lw=1, label="supply limit, blanket")
for AR, col in zip((20, 50, 100), ("tab:green", "tab:orange", "tab:red")):
    p = df[(df.kind == "patterned") & (df.AR == AR)].sort_values("x")
    ax.plot(p.x * 100, p.t90_worst, "s-", color=col, label=f"AR {AR} hole bottoms, worst on wafer ({p.wall_ratio.iloc[0]:.0f}x area)")
    ax.plot(p.x * 100, p.t_supply, ":", color=col, lw=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("TMA fraction in the feed (%)"); ax.set_ylabel("pulse time to 90 % coverage (s)")
ax.set_title("Saturating TMA pulse: blanket vs 3D-NAND-like wafer (dotted: pure supply limit)", fontsize=9); ax.legend(fontsize=7); ax.grid(True, which="both", alpha=0.3)
ax = axes[1]
for AR, col in zip((20, 50, 100), ("tab:green", "tab:orange", "tab:red")):
    p = df[(df.kind == "patterned") & (df.AR == AR)].sort_values("x")
    ax.plot(p.x * 100, p.t90_worst / p.t_supply, "s-", color=col, label=f"AR {AR}")
ax.plot(b.x * 100, b.t90_worst / b.t_supply, "ko-", label="blanket")
ax.set_xscale("log"); ax.set_xlabel("TMA fraction in the feed (%)"); ax.set_ylabel("t90 / supply-limit time"); ax.set_ylim(bottom=0)
ax.set_title("How far above the supply limit each wafer sits", fontsize=9); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(OUT / "pulse_time_map.png", dpi=150)
print(df.round(3).to_string()); print("saved", OUT / "pulse_time_map.png")
