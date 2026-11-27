"""POST HOC, after an independent review of the EDS comparison: the front slopes of the EDS scans (eds_cycle_series.csv,
least squares between 25 and 75 % of the plateau, per gap height) against (a) the model of eds_cycle_series.py with channel
narrowing at a 1 um and at a 5.5 um spot, (b) the tabulated reflectometry slope of the same sample (Yim 2020 Table 1) and
(c) the step-front value (1 - f/2)/(f z50). The profile rms of eds_cycle_series.py does not show this comparison.
Writes eds_slopes.csv.   PYTHONPATH=src .venv/bin/python docs/blindtest/v2/posthoc/eds_slopes.py"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[4]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
HERE = Path(__file__).parent; S_CAL, H = 0.04, 0.5
def model_slope(z50, fill, spot):
    K = mt.kernel("hyd" if z50 / (1 - fill / 2) < 300 else "hyd_deep"); run = K[min(K, key=lambda a: abs(np.log(a / S_CAL)))]; xx = np.arange(0.5, run[2][-1], 1.0)
    def prof(e): y = mt.cumulative(run, e, fill, xx, spot, 1.0); return y / y[:4].mean()
    def c50(y): i = int(np.argmax(y < 0.5)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i])
    y = prof(10 ** brentq(lambda le: c50(prof(10 ** le)) - z50, np.log10(run[0][3]), np.log10(run[0][-2]), xtol=1e-4)); cr = lambda l: np.interp(l, y[::-1], xx[::-1]); return 0.5 / (cr(.25) - cr(.75))
E = pd.read_csv(HERE / "eds_cycle_series.csv"); T = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv").set_index("sample"); rows = []
for _, r in E.iterrows():
    z50 = r.x50_um / H; f = r.fill; refl = float(T.loc[int(r.yim_table1_sample), "slope_dtheta_dxtilde"]); m1, m5 = model_slope(z50, f, 2.0), model_slope(z50, f, 11.0); step = (1 - f / 2) / (f * z50)
    rows.append(dict(cycles=int(r.cycles), yim_table1_sample=int(r.yim_table1_sample), channel=r.channel, fill=f, z50_gaps=round(z50, 1), n_points_25_75=int(r.n_points_25_75), slope_eds=r.slope_measured, slope_model_narrowing_spot_1um=round(m1, 5), slope_model_narrowing_spot_5p5um=round(m5, 5), slope_reflectometry_table1=refl, step_front_value=round(step, 5),
                     eds_over_model=round(r.slope_measured / m1, 3), eds_over_reflectometry=round(r.slope_measured / refl, 3), eds_over_step_value=round(r.slope_measured / step, 3), model_spot_1um_over_5p5um=round(m1 / m5, 3)))
D = pd.DataFrame(rows); D.to_csv(HERE / "eds_slopes.csv", index=False); print(D.to_string(index=False))
