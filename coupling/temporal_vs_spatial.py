"""Temporal vs spatial ALD from the supply-limit design rule.

Design rule (this work): molecules that must reach the wafer per unit area per half-cycle,
    Phi_req = 1.2 (1 + A_wall/A_wafer) N0.
Temporal tool: the whole wafer sees the showerhead supply flux J_T = U c0 for a pulse t_T = Phi_req / J_T, then purges.
    Cycle time  tau_T = t_TMA + t_H2O + 2 t_purge  (two half-cycles with their own J; purge 0.5 s each, a fixed cost
    independent of pattern, from the chamber clearing times in part 1).
Spatial tool (separated-chemistry, Olympia-type architecture): the wafer moves through a precursor zone with flux
J_S for a dwell t_d per pass; exposure per pass J_S t_d. Passes needed per half-cycle  n = ceil(Phi_req / (J_S t_d)).
    Cycle time  tau_S = n * tau_rot, with the rotation period tau_rot = t_d / f_zone (f_zone = zone fraction of the
    orbit). No purge steps; the zone separation does the purging.
Assumptions, stated, not fitted: J_S = 30 x J_T (the precursor zone runs near-pure precursor at higher local pressure;
a factor to sweep), f_zone = 0.25, tau_rot = 1 s. Output: cycle time vs A_wall/A_wafer for both, and the crossover.
"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
N0, U, C0_1PCT = 4e18, 0.050, 3.79e20
J_T = U * C0_1PCT                      # 1 % precursor: 1.9e19 /m2/s
wall = np.logspace(0, 2.2, 40)         # A_wall/A_wafer from 1 to 160
phi_req = 1.2 * (1 + wall) * N0
t_purge = 0.5
tau_T = 2 * phi_req / J_T + 2 * t_purge          # same J for both half-cycles (1 % each)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(wall, tau_T, "k-", lw=2, label="temporal showerhead, 1 % precursor, 300 sccm (this work's rule)")
ax.plot(wall, 2 * phi_req / (3 * J_T) + 2 * t_purge, "k--", lw=1, label="temporal, 3 % precursor")
for ratio, c in ((10, "tab:blue"), (30, "tab:orange"), (100, "tab:green")):
    J_S = ratio * J_T; tau_rot = 1.0; f_zone = 0.25; t_d = f_zone * tau_rot
    n = np.ceil(phi_req / (J_S * t_d))
    tau_S = 2 * n * tau_rot
    ax.step(wall, tau_S, where="post", color=c, label=f"spatial, zone flux {ratio}x showerhead, 1 s rotation, 25 % zone")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("patterned wall area / wafer area"); ax.set_ylabel("cycle time, TMA + H2O (s)")
ax.axvline(11, color="0.7", lw=0.8); ax.axvline(56, color="0.7", lw=0.8); ax.text(11, 1.1, "AR 20", fontsize=7, ha="center"); ax.text(56, 1.1, "AR 100", fontsize=7, ha="center")
ax.set_title("Temporal vs spatial ALD cycle time from the supply-limit rule (assumed zone fluxes)", fontsize=9); ax.legend(fontsize=7); ax.grid(True, which="both", alpha=0.3)
fig.tight_layout(); fig.savefig(ROOT / "figures/temporal_vs_spatial.png", dpi=150)
for w in (1, 11, 28, 56, 112):
    i = np.argmin(abs(wall - w)); print(f"wall {w:3d}x: temporal 1 % {tau_T[i]:6.1f} s, 3 % {2*phi_req[i]/(3*J_T)+1:6.1f} s | spatial 30x: {2*np.ceil(phi_req[i]/(30*J_T*0.25))*1.0:5.0f} s")
print("saved figures/temporal_vs_spatial.png")
