"""Validation of the supply-limit rule against measured high-surface-area ALD and against the Argonne plug-flow theory.

Normalised saturation time  tau_sat = (measured or simulated time to saturation) / t0,
t0 = (reactive sites) / (precursor molecules supplied per second)   [Lu, Yanguas-Gil, ..., Elam, JVST A 40 (2022);
Yanguas-Gil and Elam, arXiv 2408.13116]. Supply-limited means tau_sat -> 1 from above as the Damkoehler number grows.

Measured point, Argonne fluidized-bed Al2O3 ALD (OSTI 1902316 / JVST A 2022), all numbers theirs:
  10 g Silicycle silica gel, 940 m2 total; TMA supplied 0.0115 +/- 0.0020 mol during the 24 min TMA dose; TMA needed for
  saturation (from the deposited Al2O3 mass) 0.0103 mol; QMS breakthrough (saturation) at 24 min; utilization 91.8 +/- 14.6 %.
  -> t0 = 0.0103 / (0.0115/24 min) = 21.5 min, tau_sat = 24 / 21.5 = 1.12 (range 1.00-1.27 from the +/-0.002 mol supply).
  Their Damkoehler number: ~1e4.
Our chamber model: tau_sat = t90 / [(1 + A) N0/(U c0)] from figures/pulse_time_map.csv; the 1.2 in the design rule is this
ratio. Damkoehler number for the chamber, their definition Da = (S/V) beta (vbar/4) t_res with S the total reactive area,
V the chamber volume (1.88e-3 m3), beta = 5e-3, vbar/4 = 102.5 m/s, t_res = 0.53 s.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = Path(__file__).resolve().parents[1] / "figures"
df = pd.read_csv(OUT / "pulse_time_map.csv")
A_WAFER, V, BETA, VQ, TRES = np.pi * 0.15 ** 2, 1.88e-3, 5e-3, 102.5, 0.53
df["tau_sat"] = df.t90_worst / df.t_supply
df["Da"] = ((1 + df.wall_ratio) * A_WAFER / V) * BETA * VQ * TRES
# Argonne point
t0_min = 0.0103 / (0.0115 / 24.0); tau_arg = 24.0 / t0_min
tau_lo = 24.0 / (0.0103 / ((0.0115 + 0.0020) / 24.0)); tau_hi = 24.0 / (0.0103 / ((0.0115 - 0.0020) / 24.0))
print(f"Argonne fluidized bed: t0 = {t0_min:.1f} min, measured 24 min -> tau_sat = {tau_arg:.2f} (range {min(tau_lo,tau_hi):.2f}-{max(tau_lo,tau_hi):.2f}), Da ~ 1e4")
print(df[["kind", "x", "AR", "wall_ratio", "Da", "t_supply", "t90_worst", "tau_sat"]].round(3).to_string())
fig, ax = plt.subplots(figsize=(7.5, 4.5))
for kind, mk, lab in (("blanket", "o", "this work: blanket wafer (chamber model)"), ("patterned", "s", "this work: patterned wafers AR 20-100 (two-way model)")):
    d = df[df.kind == kind]; ax.plot(d.Da, d.tau_sat, mk, ms=7, label=lab)
lo, hi = min(tau_lo, tau_hi), max(tau_lo, tau_hi)
ax.errorbar([1e4], [tau_arg], yerr=[[tau_arg - lo], [hi - tau_arg]], fmt="D", color="tab:red", ms=8, capsize=4,
            label="measured: Argonne fluidized-bed Al2O3 ALD, 940 m2 silica gel (Lu et al. 2022)")
Da = np.logspace(0.5, 4.5, 100)
ax.plot(Da, 1 / (1 - np.log(2) / Da).clip(0.05, None), "k:", lw=1, label="plug-flow theory: tau for 90 % coverage*")
ax.axhline(1, color="0.6", lw=0.8)
ax.set_xscale("log"); ax.set_xlabel("Damkoehler number  (S/V) beta (vbar/4) t_res"); ax.set_ylabel("saturation time / supply-limit time t0")
ax.set_ylim(0.8, 4.5); ax.set_title("Supply-limited ALD: simulations, theory and a measurement on one axis", fontsize=9)
ax.legend(fontsize=7); ax.grid(True, which="both", alpha=0.3)
fig.text(0.01, 0.01, "*rough: inverted from theta(t0) = 1 - ln2/Da (Lu et al. Eq. 9), i.e. the coverage reached when exactly one site-equivalent has been supplied", fontsize=6)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(OUT / "validation_supply_rule.png", dpi=150); print("saved", OUT / "validation_supply_rule.png")
