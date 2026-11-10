"""Plot the one-factor-at-a-time operating sweep: saturation time (centre, worst on wafer), capture and the
supply-limited monolayer time, versus pressure, flow, gap and TMA fraction."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "figures"
df = pd.read_csv(OUT / "reactor_sweep_operating.csv").set_index("name")
base = df.loc["base"]
factors = [("p", "chamber pressure (Pa)", ["p100", "base", "p1000"], True),
           ("Q", "N2 flow (sccm)", ["Q150", "base", "Q600"], True),
           ("g", "gap (mm)", ["g10", "base", "g25"], False),
           ("x", "TMA fraction", ["x0.005", "base", "x0.03"], True)]
fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=False)
for ax, (col, label, names, logx) in zip(axes, factors):
    names = [n for n in names if n in df.index]
    d = df.loc[names]
    ax.plot(d[col], d["t90_centre"], "o-", label="t90 centre")
    ax.plot(d[col], d["t90_worst"], "s--", label="t90 worst on wafer")
    ax.plot(d[col], d["t_mono"], "^:", color="0.5", label="monolayer supply time N0/(U c0)")
    if logx: ax.set_xscale("log")
    ax.set_xlabel(label); ax.set_ylabel("time (s)"); ax.set_ylim(bottom=0)
    ax2 = ax.twinx(); ax2.plot(d[col], d["captured_pct"], "d-", color="tab:red", alpha=0.6, label="captured %")
    ax2.set_ylim(0, 100); ax2.set_ylabel("TMA captured (%)", color="tab:red")
    ax.set_title(f"vary {col}", fontsize=10)
axes[0].legend(fontsize=7, loc="upper left")
fig.suptitle("Self-limiting wafer: time to 90 % coverage and TMA capture vs operating point (baseline 300 Pa, 300 sccm, 15 mm, 1 % TMA; s = 5e-3)", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "reactor_sweep_operating.png", dpi=150)
cols = ["p", "Q", "g", "x", "Pe", "t_mono", "t90_centre", "t90_worst", "t90_spread_pct", "captured_pct", "theta_at_0p1s"]
print(df[cols].round(4).to_string())
print("saved", OUT / "reactor_sweep_operating.png")
