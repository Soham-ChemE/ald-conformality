"""DOE figure: (a) saturating pulse vs precursor supply rate with the rule; (b) utilization vs pressure and flow;
(c) time-utilization trade-off coloured by pressure. From figures/doe_results.csv (25 runs, two-way multiscale)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; d = pd.read_csv(ROOT / "figures/doe_results.csv")
KB, T, N0, WALL = 1.380649e-23, 573.15, 4e18, 55.9
d["supply"] = d.U * d.x * d.p / (KB * T)          # molecules / m2 / s at the wafer
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
sc = ax[0].scatter(d.supply, d.t90_worst, c=np.log10(d.p), cmap="viridis", s=60, edgecolor="k")
s = np.logspace(np.log10(d.supply.min() * 0.8), np.log10(d.supply.max() * 1.2), 50)
ax[0].plot(s, 1.0 * (1 + WALL) * N0 / s, "k:", label="supply limit (1 + A) N0 / (U c0)"); ax[0].plot(s, 1.2 * (1 + WALL) * N0 / s, "k--", lw=1, label="design rule, 1.2 x")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("precursor supply at the wafer, U c0 (molecules / m2 / s)"); ax[0].set_ylabel("saturating pulse, worst point (s)")
ax[0].set_title("25 operating points: pulse time vs supply (colour: log10 pressure)", fontsize=9); ax[0].legend(fontsize=8); plt.colorbar(sc, ax=ax[0], label="log10 p (Pa)")
for p, mk in ((100, "o"), (316.228, "s"), (1000, "^")):
    dd = d[np.isclose(d.p, p, rtol=1e-3)].sort_values("Q"); ax[1].plot(dd.Q, dd.utilization * 100, mk, ms=8, label=f"{p:.0f} Pa")
ax[1].set_xscale("log"); ax[1].set_xlabel("N2 carrier flow (sccm)"); ax[1].set_ylabel("precursor utilization up to saturation (%)"); ax[1].set_ylim(40, 100)
ax[1].set_title("Utilization rises with pressure and falls with carrier flow", fontsize=9); ax[1].legend(fontsize=8); ax[1].grid(True, alpha=0.3)
sc2 = ax[2].scatter(d.t90_worst / d.t_supply, d.utilization * 100, c=np.log10(d.p), cmap="viridis", s=60, edgecolor="k")
ax[2].set_xlabel("saturating pulse / supply-limit time"); ax[2].set_ylabel("utilization (%)"); ax[2].set_title("The two are the same number: utilization = 1 / (t90 / t_supply)", fontsize=9); ax[2].grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(ROOT / "figures/doe_results.png", dpi=150); print("saved figures/doe_results.png")
print(d.groupby("p")[["utilization", "tau_sat"]].mean().round(3)); print(d.groupby("Q")[["utilization", "tau_sat"]].mean().round(3)); print(d.groupby("g")[["utilization", "tau_sat", "edge_lag_pct"]].mean().round(3))
