"""Observables that would separate the slit transport models without knowing the sticking coefficient or the dose
(Phase 2). Complete model (narrowing + spot) at the same per-cycle exposure:
  (1) half-thickness depth at carrier pressures 0, 100, 300, 1000 Pa in a 500 nm gap after 500 cycles;
  (2) reduced depth (x50 / gap) in a 2000 nm and in a 100 nm gap over that in a 500 nm gap, with the measured gaps and film
      thicknesses of Yim 2020 samples 9, 7 and 8, for the Monte Carlo with N2 at 300 Pa and without gas;
  and the hydraulic model with Bosanquet for comparison. Writes figures/identifiability_experiments.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/identifiability_experiments.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
def depth(kn, s, E, fill, n, spot):
    K = mt.kernel(kn); ss = np.array(list(K)); lo, hi = ss[ss <= s].max(), ss[ss >= s].min(); z = lambda sv: mt.front(K[sv], E, fill, spot, n)[0]
    if lo == hi: return z(lo)
    w = (np.log(s) - np.log(lo)) / (np.log(hi) - np.log(lo)); return (1 - w) * z(lo) + w * z(hi)
kB, T, amu = 1.380649e-23, 573.15, 1.66053906660e-27; v = np.sqrt(8 * kB * T / (np.pi * 72.09 * amu)); mu = 72.09 * 28.0134 / (72.09 + 28.0134) * amu
Dg = lambda p: 3 / (8 * p / (kB * T) * (4.825e-10) ** 2) * np.sqrt(kB * T / (2 * np.pi * mu)); De = lambda H, p: 1 / (1 / (2 * H * v / 3) + 1 / Dg(p))
G500, G2000, G100, F500, F2000, F100 = 519.1, 2034.7, 108.6, 51.1, 48.8, 47.8; rows = []
for s in (0.008, 0.0158):
    for E in (3.0, 4.0):
        d = {p: depth(k, s, E, 0.2, 0.92 if p else 1.0, 11.0) for p, k in ((0, "walls90_pillars49"), (100, "gas100"), (300, "gas300"), (1000, "gas1000"))}
        rows.append(dict(test="depth against carrier pressure, 500 nm gap, 500 cycles", s=s, exposure=E, **{f"z50_{p}Pa": round(v_, 1) for p, v_ in d.items()}, monte_carlo=round(d[1000] / d[100], 3), hydraulic_bosanquet=round(float(np.sqrt(De(5e-7, 1000) / De(5e-7, 100))), 3), quantity="depth at 1000 Pa / depth at 100 Pa"))
for s in (0.0125, 0.02):
    for E in (3.5, 4.5):
        b_g = depth("gas300", s, E, 2 * F500 / G500, 0.92, 5500 / G500); b_f = depth("walls90_pillars49", s, E, 2 * F500 / G500, 1.0, 5500 / G500)
        rows.append(dict(test="reduced depth, 2000 nm gap over 500 nm gap", s=s, exposure=E, monte_carlo=round(depth("h2000_gas300", s, E, 2 * F2000 / G2000, 0.92, 5500 / G2000) / b_g, 3), monte_carlo_no_gas=round(depth("h2000_free", s, E, 2 * F2000 / G2000, 1.0, 5500 / G2000) / b_f, 3),
                         hydraulic_bosanquet=round(float(np.sqrt((De(G2000 * 1e-9, 300) / G2000) / (De(G500 * 1e-9, 300) / G500))), 3), measured=round((375 / G2000) / (106 / G500), 3), quantity="Yim 2020 samples 9 and 8; not stated to share a run"))
        rows.append(dict(test="reduced depth, 100 nm gap over 500 nm gap", s=s, exposure=E, monte_carlo=round(depth("h100_gas300", s, E, 2 * F100 / G100, 0.92, 5500 / G100) / b_g, 3), monte_carlo_no_gas=round(depth("h100_free", s, E, 2 * F100 / G100, 1.0, 5500 / G100) / b_f, 3),
                         measured=round((10 / G100) / (106 / G500), 3), quantity="Yim 2020 samples 7 and 8; depth of sample 7 tabulated as 10 um; membrane sag not modelled"))
d = pd.DataFrame(rows); d.to_csv(ROOT / "figures/identifiability_experiments.csv", index=False); print(d.drop(columns=[c for c in d.columns if c.startswith("z50_")]).to_string(index=False))
