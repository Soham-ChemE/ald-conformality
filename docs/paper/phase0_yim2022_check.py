"""Phase 0, check 1: the narrowing law against an independent multi-cycle model.
Yim et al., Phys. Chem. Chem. Phys. 24, 8645 (2022), ESI Table S1, cycle series: Ylilammi diffusion-reaction model with
narrowing, H = 500 nm, c = 0.01, Kn = 7.6, N = 5 to 500 cycles. Their slope is of the scaled profile (nm per cycle per gap
height), i.e. gpc x d(theta)/d(x/H) with gpc = q M / (rho N_A) = 0.0949 nm for their baseline (which reproduces their
'fraction filled' 0.0949 at 250 cycles), so d(theta)/d(x/H) = slope / 0.0949.
Test: take a single-exposure front with their N = 5 depth and slope (hydraulic kernel of mechanism_test.py, s chosen to give
that slope at that depth), apply the law (average over the gap stretch f^n, f uniform on [1 - fill, 1]) with their fill
fractions, and compare the half-thickness depth and the local slope at half thickness with their table. Nothing is fitted
to the N > 5 rows. Writes figures/phase0_yim2022_check.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/phase0_yim2022_check.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
GPC = 4e18 * 0.050 / (3500 * 6.02214076e23) * 1e9                         # nm per cycle: q = 4 nm^-2, M = 0.050 kg/mol, rho = 3500 kg/m3
TABLE = [(5, 0.0019, 261.5, 0.0029), (10, 0.0038, 261.3, 0.0029), (20, 0.0076, 260.8, 0.0030), (50, 0.0190, 259.4, 0.0029), (100, 0.0380, 257.0, 0.0029), (250, 0.0950, 250.1, 0.0026), (500, 0.1898, 238.5, 0.0018)]
K = mt.kernel("hyd"); xx = np.arange(0.5, 399, 0.25)
def prof(run, e, fill, n):
    y = mt.cumulative(run, e, fill, xx, 0.0, n); return y / y[:4].mean()
def x50_slope(y):
    i = int(np.argmax(y < 0.5)); x50 = xx[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i]) * 0.25
    j = slice(i - 8, i + 8); return x50, -np.polyfit(xx[j], y[j], 1)[0]                          # local slope over +-2 gaps
def single(s):                                                                                  # exposure putting the N = 5 front at 261.5
    run = K[s]; g = lambda le: x50_slope(prof(run, 10 ** le, TABLE[0][1], 1.0))[0] - TABLE[0][2]
    le = brentq(g, np.log10(run[0][3]), np.log10(run[0][-2]), xtol=1e-5); return run, 10 ** le
target = TABLE[0][3] / GPC; ss = np.array(list(K)); sl = np.array([x50_slope(prof(*single(s), TABLE[0][1], 1.0))[1] for s in ss])
s0 = float(np.exp(np.interp(np.log(target), np.log(sl), np.log(ss)))); lo, hi = ss[ss <= s0].max(), ss[ss >= s0].min(); w = (np.log(s0) - np.log(lo)) / ((np.log(hi) - np.log(lo)) or 1)
rows = []
for n in (1.0, 0.95):
    for N, fill, x50_t, slope_t in TABLE:
        y = (1 - w) * prof(*single(lo), fill, n) + w * prof(*single(hi), fill, n); x50, slope = x50_slope(y)
        rows.append(dict(gap_exponent=n, cycles=N, fill=fill, x50_yim=x50_t, x50_law=round(x50, 1), x50_error_pct=round(100 * (x50 / x50_t - 1), 2), slope_yim_scaled=slope_t,
                         slope_law_scaled=round(slope * GPC, 5), slope_yim_theta=round(slope_t / GPC, 4), slope_law_theta=round(slope, 4), slope_error_pct=round(100 * (slope * GPC / slope_t - 1), 1)))
d = pd.DataFrame(rows); d.to_csv(ROOT / "figures/phase0_yim2022_check.csv", index=False)
print(f"gpc {GPC:.4f} nm; single-exposure front: hydraulic kernel, s = {s0:.4f} (their model: c = 0.01 with gas-phase diffusion, Kn = 7.6)"); print(d.to_string(index=False))
