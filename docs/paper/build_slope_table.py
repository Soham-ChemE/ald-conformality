"""Tabulate the front slope of the complete forward model (transport kernel + channel narrowing + measurement spot) so that a
measured slope can be turned into a sticking coefficient with its sensitivity, or flagged as not identifiable.
Grid: transport (hydraulic 1D; Monte Carlo of the real 500 nm cavity without gas and with N2 at 100, 300, 1000 Pa),
half-thickness depth 150 to 300 gaps, filled fraction 0 to 0.4, spot 0 or 11 gaps (5.5 um on a 500 nm gap), every simulated s.
Writes figures/slope_table.csv (used by aldconf.lhar). About ten minutes.
  PYTHONPATH=src .venv/bin/python docs/paper/build_slope_table.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
KERNELS = {"hydraulic": ("hyd_deep", 1.0), "free": ("walls90_pillars49", 1.0), "gas100": ("gas100", 0.92), "gas300": ("gas300", 0.92), "gas1000": ("gas1000", 0.92)}
rows = []
for name, (kn, n) in KERNELS.items():
    K = mt.kernel(kn)
    for z50 in (150, 200, 250, 300):
        for fill in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4):
            for spot in (0.0, 11.0):
                for s, run in K.items():
                    sl, e = mt.slope_at_depth(run, z50, fill, spot, n)
                    rows.append(dict(transport=name, z50=z50, fill=fill, spot_gaps=spot, s=s, slope=None if not np.isfinite(sl) else round(float(sl), 6), exposure_Pa_s=None if not np.isfinite(sl) else round(float(e), 3)))
    print(name, "done", flush=True)
pd.DataFrame(rows).to_csv(ROOT / "figures/slope_table.csv", index=False); pd.DataFrame(rows).to_csv(ROOT / "src/aldconf/lhar_slope_table.csv", index=False); print(len(rows), "rows")
