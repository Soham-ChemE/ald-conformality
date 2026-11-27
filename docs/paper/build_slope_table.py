"""Tabulate the front slope of the complete forward model (transport kernel + channel narrowing + measurement spot) so that a
measured slope can be turned into a sticking coefficient with its sensitivity, or flagged as not identifiable.
Grid: transport (hydraulic 1D; Monte Carlo of the real 500 nm cavity without gas and with N2 at 30 to 10000 Pa),
half-thickness depth 150 to 300 gaps, filled fraction 0 to 0.4, spot 0 to 55 gaps (5.5 um is 11 gaps at 500 nm, 2.75 at
2000 nm, 55 at 100 nm), every simulated s. The gap exponent of the narrowing average is 1 without gas and comes from
figures/narrowing_exponent.csv with gas (0.92 at 300 Pa, the value measured directly and used throughout).
Writes figures/slope_table.csv (used by aldconf.lhar). About an hour.
  PYTHONPATH=src .venv/bin/python docs/paper/build_slope_table.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
_n = pd.read_csv(ROOT / "figures/narrowing_exponent.csv").set_index("pressure_Pa_at_500nm").exponent_n
KERNELS = {"hydraulic": ("hyd_deep", 1.0), "free": ("walls90_pillars49", 1.0), "gas300": ("gas300", 0.92)}
for _p in (30, 100, 1000, 3000, 10000): KERNELS[f"gas{_p}"] = (f"gas{_p}", float(_n[_p]))
if len(sys.argv) > 1: KERNELS = {k: v for k, v in KERNELS.items() if k in sys.argv[1:]}
rows = []
for name, (kn, n) in KERNELS.items():
    K = mt.kernel(kn)
    for z50 in (150, 200, 250, 300):
        for fill in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4):
            for spot in (0.0, 2.75, 5.5, 11.0, 22.0, 55.0):
                for s, run in K.items():
                    sl, e = mt.slope_at_depth(run, z50, fill, spot, n)
                    rows.append(dict(transport=name, gap_exponent=n, z50=z50, fill=fill, spot_gaps=spot, s=s, slope=None if not np.isfinite(sl) else round(float(sl), 6), exposure_Pa_s=None if not np.isfinite(sl) else round(float(e), 3)))
    print(name, "done", flush=True)
d = pd.DataFrame(rows)
if len(sys.argv) > 1: d = pd.concat([pd.read_csv(ROOT / "figures/slope_table.csv").query("transport not in @sys.argv[1:]"), d])
d.to_csv(ROOT / "figures/slope_table.csv", index=False); d.to_csv(ROOT / "src/aldconf/lhar_slope_table.csv", index=False); print(len(d), "rows")
pd.read_csv(ROOT / "figures/narrowing_exponent.csv")[["pressure_x_gap_Pa_um", "exponent_n"]].to_csv(ROOT / "src/aldconf/lhar_gap_exponent.csv", index=False)
