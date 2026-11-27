"""Check that the 500 nm slope tables serve other gaps through pressure x gap. Direct simulations of the 2000 nm and 100 nm
PillarHall cavities with N2 at 300 Pa (figures/pillar_mc/h2000_gas300_*, h100_gas300_*; pressure x gap = 600 and 30 Pa um,
pillar spacing 24.5 and 490 gaps instead of 98) against `aldconf.lhar.slope_curve_at`, which interpolates the 500 nm tables.
Spot 5.5 um (2.75 and 55 gaps) and no spot; filled fraction 0, 0.1, 0.2. Writes figures/slope_table_scaling_check.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/slope_table_scaling_check.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
from aldconf.lhar import slope_curve_at, gap_exponent
rows = []
for name, gap in (("h2000_gas300", 2000.0), ("h100_gas300", 100.0)):
    K = mt.kernel(name); n = gap_exponent(300, gap)
    for z50 in (150, 200, 250):
        for fill in (0.0, 0.1, 0.2):
            for spot_um in (0.0, 5.5):
                try: st, slt = slope_curve_at(z50, fill, 300, gap, spot_um)
                except ValueError: continue
                for s, run in K.items():
                    if not st[0] <= s <= st[-1]: continue
                    sl, _ = mt.slope_at_depth(run, z50, fill, spot_um * 1000 / gap, n)
                    if not np.isfinite(sl): continue
                    tab = float(np.exp(np.interp(np.log(s), np.log(st), np.log(slt))))
                    rows.append(dict(gap_nm=gap, pressure_x_gap_Pa_um=300 * gap / 1000, z50=z50, fill=fill, spot_um=spot_um, s=s, slope_direct=round(float(sl), 6), slope_table=round(tab, 6), ln_ratio=round(float(np.log(tab / sl)), 4)))
d = pd.DataFrame(rows); d.to_csv(ROOT / "figures/slope_table_scaling_check.csv", index=False)
print(d.groupby(["gap_nm", "fill", "spot_um"]).ln_ratio.agg(["count", "mean", lambda v: np.sqrt(np.mean(v ** 2)), lambda v: np.abs(v).max()]).round(3).to_string())
print("all: rms", round(float(np.sqrt(np.mean(d.ln_ratio ** 2))), 3), "max", round(float(d.ln_ratio.abs().max()), 3), "n", len(d))
