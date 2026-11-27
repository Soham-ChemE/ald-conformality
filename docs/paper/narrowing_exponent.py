"""Gap exponent n of the channel-narrowing law as a function of carrier pressure. The penetration of a saturating front goes as
(D H)^(1/2) at fixed exposure, and in a slit the diffusivity is D = H vbar g(p H), so the depth scales with the gap as H^n with
    n = 1 + (1/2) d ln D / d ln p     (at fixed gap).
D comes from the Monte Carlo of tracer diffusion along an open slit (`pillarmc slitdiffusion`, TMA in N2, 300 C, 500 nm), the
derivative from a quadratic fit to ln D against ln p over a window of +-1 decade. Compared with the exponents measured directly
from pairs of full front simulations at 400 and 500 nm (figures/mechanism_gap_scaling.csv, docs/blindtest/posthoc/
gap_scaling_pressure.csv). Writes figures/narrowing_exponent.csv; pressures scale to other gaps as p x gap / 500 nm.
  PYTHONPATH=src .venv/bin/python docs/paper/narrowing_exponent.py
"""
import subprocess
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; EXE = ROOT / "lib/pillarmc/pillarmc"
P = [10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000]; D, ERR = [], []
for p in P:
    o = subprocess.check_output([str(EXE), "slitdiffusion", str(p), "0.5", "10000"], text=True).split(); v = {o[i]: float(o[i + 1]) for i in range(0, len(o), 2)}
    D.append(v["D_slit"]); ERR.append(v["rel_err"]); print(p, v["D_slit"], v["rel_err"], flush=True)
lp, lD = np.log(P), np.log(D); rows = []
g1 = pd.read_csv(ROOT / "figures/mechanism_gap_scaling.csv"); g2 = pd.read_csv(ROOT / "docs/blindtest/posthoc/gap_scaling_pressure.csv")
direct = {300: g1[g1.transport == "gas 300 Pa"].exponent.mean(), 1400: g2[g2.pressure_Pa == 1400].exponent.mean(), 3000: g2[g2.pressure_Pa == 3000].exponent.mean()}
for p in (30, 100, 300, 1000, 1400, 3000, 10000, 30000, 100000):
    k = np.abs(lp - np.log(p)) <= np.log(10) + 1e-9; c = np.polyfit(lp[k] - np.log(p), lD[k], 2)
    rows.append(dict(pressure_Pa_at_500nm=p, pressure_x_gap_Pa_um=p * 0.5, D_m2_s=float(np.exp(c[2])), exponent_n=round(1 + c[1] / 2, 3), exponent_direct=round(direct[p], 3) if p in direct else None))
d = pd.DataFrame(rows); d.to_csv(ROOT / "figures/narrowing_exponent.csv", index=False); print(d.to_string(index=False))
