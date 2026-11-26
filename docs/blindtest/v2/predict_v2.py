"""Second blind test, predictions (PROTOCOL_v2.md). Uses docs/blindtest/v2/inputs.json and the released calibration curve only.
One sticking coefficient is fitted to the calibration curve with the complete model (hydraulic kernel, channel narrowing with
exponent 1, top-hat spot; half depth fixed at the released value, amplitude free; s over the simulated 0.0025 to 0.06), then
the whole profile of every other curve of the series is predicted at its released half depth, with and without narrowing.
Writes predictions_v2.json (profiles on a grid of 1 gap) and predictions_v2.csv.
  PYTHONPATH=src .venv/bin/python docs/blindtest/v2/predict_v2.py
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[3]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
V2 = ROOT / "docs/blindtest/v2"; INP = json.loads((V2 / "inputs.json").read_text()); out = {"series": {}}; rows = []; ALLC = next(v for v in INP.values() if isinstance(v, list) and v and isinstance(v[0], dict) and "id" in v[0])
def model(sk, z50, fill, spot):
    K = mt.kernel("hyd" if z50 / (1 - fill / 2) < 300 else "hyd_deep"); run = K[min(K, key=lambda a: abs(np.log(a / sk)))]; xx = np.arange(0.5, run[2][-1], 1.0)
    def prof(e): y = mt.cumulative(run, e, fill, xx, spot, 1.0); return y / y[:4].mean()
    def c50(y): i = int(np.argmax(y < 0.5)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i])
    g = lambda le: c50(prof(10 ** le)) - z50; lo, hi = np.log10(run[0][3]), np.log10(run[0][-2])
    if not (np.isfinite(g(lo)) and np.isfinite(g(hi))) or g(lo) * g(hi) > 0: return None, None
    return xx, prof(10 ** brentq(g, lo, hi, xtol=1e-4))
def slope(x, y):
    b = (y <= 0.75) & (y >= 0.25); return float(-np.polyfit(x[b], y[b], 1)[0])
for S in INP["series"]:
    curves = [c for c in ALLC if c["series_id"] == S["series_id"]]; cal = next(c for c in curves if c["role"] == "calibration"); H = cal["gap_height_nm"]; d = pd.read_csv(V2 / f"calibration_{cal['id']}.csv"); x = d.x_over_H.values; y = d.thickness_over_plateau.values
    fill = cal["filled_fraction_2t_over_H"]; z50 = cal["half_thickness_depth_x50_over_H"]; spot = cal["spot_size_um"] * 1000 / H; m = (x >= 0) & (x <= 1.6 * z50); tab = []
    for sk in mt.S_HYD:
        xm, ym = model(sk, z50, fill, spot)
        if xm is None: continue
        yi = np.interp(x[m], xm, ym, right=0.0); a = float(yi @ y[m] / (yi @ yi)); tab.append((sk, float(np.sqrt(np.mean((a * yi - y[m]) ** 2)))))
    tab = np.array(tab); j = int(np.argmin(tab[:, 1])); s_cal, rms_cal, edge = float(tab[j, 0]), float(tab[j, 1]), bool(j in (0, len(tab) - 1))
    out["series"][S.get("series_id", cal["series_id"])] = dict(calibration=cal["id"], s=s_cal, rms_calibration=round(rms_cal, 4), s_at_range_end=edge, s_grid=[float(v) for v in tab[:, 0]], rms_grid=[round(float(v), 4) for v in tab[:, 1]], pass_limit=round(max(0.06, 1.5 * rms_cal), 4), curves={})
    for c in curves:
        if c["role"] != "predict": continue
        Hc = c["gap_height_nm"]; f = min(c["filled_fraction_2t_over_H"], 0.98); z = c["half_thickness_depth_x50_over_H"]; sp = (c.get("spot_size_um") or 5.5) * 1000 / Hc
        xa, ya = model(s_cal, z, f, sp); xb, yb = model(s_cal, z, 0.0, sp)
        if xa is None: rows.append(dict(id=c["id"], note="half depth outside the simulated exposures")); continue
        k = xa <= 1.7 * z; out["series"][cal["series_id"]]["curves"][c["id"]] = dict(x=[float(v) for v in xa[k]], complete=[round(float(v), 4) for v in ya[k]], no_narrowing=[round(float(v), 4) for v in np.interp(xa[k], xb, yb)])
        cr = lambda yy, l: float(xa[int(np.argmax(yy < l)) - 1])
        rows.append(dict(id=c["id"], series=cal["series_id"], s=s_cal, gap_nm=Hc, cycles=c["cycles"], fill=round(f, 3), z50=z, spot_gaps=round(sp, 1), slope_complete=round(slope(xa, ya), 5), slope_no_narrowing=round(slope(xb, yb), 5),
                         r80_complete=round(cr(ya, 0.8) / z, 3), r20_complete=round(cr(ya, 0.2) / z, 3), fill_minus_calibration=round(f - fill, 3), slope_known_to_predictor=c.get("slope_known_to_predictor")))
    print("series", cal["series_id"], "calibration", cal["id"], "s =", s_cal, "rms", round(rms_cal, 4), "edge", edge, "| grid", [(float(a), round(float(b), 4)) for a, b in tab])
(V2 / "predictions_v2.json").write_text(json.dumps(out)); D = pd.DataFrame(rows); D.to_csv(V2 / "predictions_v2.csv", index=False); print(D.to_string(index=False))
