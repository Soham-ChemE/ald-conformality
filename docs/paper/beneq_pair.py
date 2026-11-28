"""Two reflectometry line scans of Al2O3 from TMA and water in a second reactor (Beneq TFS-200, 300 C, 2.1 to 2.2 mbar,
0.15 s pulses and 1.0 s purges; 250 and 500 cycles in a PillarHall LHAR4 chip, gap 500 nm), fitted one by one with the
sticking coefficient and the exposure free, with and without channel narrowing, by the functions of profile_fits.py.

Source of the curves: Chipmetrics, PillarHall LHAR4 Test Chip Analysis Guide V3 (13.04.2022), Fig. 20, a raster image of 663 x 414
pixels (1 pixel = 0.56 um = 1.1 gap heights; 0.25 nm). Process conditions: Utriainen et al., Proc. SPIE 12008, 120080D (2022),
which reports the same thickness series. The figure does not state the reflectometer objective; the 5.5 um spot of the 50x
objective of the same instrument family (Yim et al. 2020) is assumed. The pillar layout of the LHAR4 chip is taken as the standard one.

The scans keep a level of 3 to 6 nm behind the front and the film thickness at the opening is given as 33.5 and 64.2 nm, so three
readings are fitted:
  as read            thickness over the plateau of the scan; filling from that plateau
  floor subtracted   the level behind the front taken off both; filling from the difference
  opening thickness  as read, with the filling from the tabulated thickness at the opening
  .venv/bin/python docs/paper/beneq_pair.py      -> figures/profile_fits_beneq.csv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(Path(__file__).parent))
import profile_fits as pf
H_NM = 500.0; SPOT = 5.5 * 1000 / H_NM; OPENING = {250: 33.5, 500: 64.2}
KERNEL = {"hydraulic": "hyd_deep", "gas": "gas300_deep", "ballistic": "free_deep"}; EXPONENT = {"hydraulic": 1.0, "gas": 0.92, "ballistic": 1.0}
rows = []
for cyc in (250, 500):
    d = pd.read_csv(ROOT / f"data/lhar_profiles/chipmetrics_guide_fig20_{cyc}cycles.csv"); d = d[d.x_um > 3]          # the first 3 um are the edge of the opening
    plateau = float(d[(d.x_um > 40) & (d.x_um < 140)].t_nm.median()); floor = float(d[d.x_um > 220].t_nm.median()); z = d.x_um.values * 1000 / H_NM
    readings = {"as read": (d.t_nm.values / plateau, plateau, (0.3,)), "floor subtracted": ((d.t_nm.values - floor) / (plateau - floor), plateau - floor, (0.3, 0.1)), "opening thickness": (d.t_nm.values / plateau, OPENING[cyc], (0.3,))}
    for name, (y, film, windows) in readings.items():
        fill = 2 * film / H_NM
        for win in windows:
            pf.CUT = win
            for fam in ("hydraulic", "gas", "ballistic"):
                K = pf.kern(KERNEL[fam])
                for narrow in (False, True):
                    rows.append(dict(cycles=cyc, reading=name, window=win, plateau_nm=round(plateau, 1), floor_nm=round(floor, 1), film_nm=round(film, 1), fill=round(fill, 3), transport=fam, narrowing=narrow, **pf.fit(z, y, fill, SPOT, K, EXPONENT[fam], narrow)))
            print(cyc, name, win, "done", flush=True)
D = pd.DataFrame(rows); D.to_csv(ROOT / "figures/profile_fits_beneq.csv", index=False)
P = D.pivot_table(index=["reading", "window", "transport", "narrowing"], columns="cycles", values="s_fit"); P["s500_over_s250"] = (P[500] / P[250]).round(2)
R = D.pivot_table(index=["reading", "window", "transport", "narrowing"], columns="cycles", values="rms"); P["rms_250"] = R[250]; P["rms_500"] = R[500]
pd.set_option("display.width", 200); print(P.to_string())
