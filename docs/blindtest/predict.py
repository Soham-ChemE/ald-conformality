"""Blind predictions (see PROTOCOL.md). Written from docs/blindtest/inputs.json only; the custody folder was not opened.

Model: single-exposure coverage of the 1D hydraulic-diameter kernel in reduced units (docs/paper/mechanism_test.py), which
holds for any precursor because the exposure is set by the measured half-thickness depth; channel narrowing by the averaging
law with gap exponent 1; top-hat spot. A cylindrical pore of diameter d is the same problem with distance in units of d and
sticking 4 s (wall area per volume 4/d against 2/H, hydraulic diameter d against 2H), and the film fills f = 2 t / d.
Slope definition used for prediction and for scoring: least-squares line through the points with thickness/plateau between
0.35 and 0.65, in units of gap (or diameter). r80 and r20 are the depths at 0.8 and 0.2 of the plateau over the half depth.

Kinds of prediction
  ceiling      the filled fraction is so large that the slope hardly depends on s: value at s = 0.06, range over s = 0.0025-0.06
  range        slope depends on s and s is not known for this precursor: range over s = 0.0025-0.06 (a weak test)
  transfer     s taken from another published sample or laboratory, stated per curve
  bound        only an upper bound: the narrowing limit (1 - f/2) / (f z50)
Pass limits: the predicted interval widened by 13 % and by the custodian's digitisation error (10 % on the slope where the
curve was read from a raster image, 5 % for vector data).
  PYTHONPATH=src .venv/bin/python docs/blindtest/predict.py
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
INP = json.loads((ROOT / "docs/blindtest/inputs.json").read_text()); INP = INP if isinstance(INP, list) else INP["curves"]; BY = {c["id"]: c for c in INP}
val = lambda v: v["value"] if isinstance(v, dict) and "value" in v else v

def shape(run, z50, fill, spot):
    """Slope (regression 0.35-0.65), r80, r20 of the model profile whose half depth is z50."""
    zmax = run[2][-1]; xx = np.arange(0.5, zmax, 0.5 if zmax < 500 else 1.0)
    def prof(e): y = mt.cumulative(run, e, fill, xx, spot, 1.0); return y / y[:4].mean()
    def cross(y, l): i = int(np.argmax(y < l)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - l) / (y[i - 1] - y[i]) * (xx[i] - xx[i - 1])
    g = lambda le: cross(prof(10 ** le), 0.5) - z50; lo, hi = np.log10(run[0][3]), np.log10(run[0][-2])
    if not (np.isfinite(g(lo)) and np.isfinite(g(hi))) or g(lo) * g(hi) > 0: return None
    y = prof(10 ** brentq(g, lo, hi, xtol=1e-4)); m = (y <= 0.65) & (y >= 0.35) & (xx > 2)
    return -np.polyfit(xx[m], y[m], 1)[0], cross(y, 0.8) / z50, cross(y, 0.2) / z50

def predict(cid, unit_nm, kind, spot_um, s_values, s_central, raster, note, tube=False, film_nm=None):
    c = BY[cid]; x50 = val(c["x50_um"]); t = film_nm if film_nm is not None else val(c["plateau_thickness_nm"]); fill = min(2 * t / unit_nm, 0.98); z50 = x50 * 1000 / unit_nm; spot = spot_um * 1000 / unit_nm
    K = mt.kernel("hyd" if z50 / (1 - fill / 2) < 300 else "hyd_deep"); out = {}
    for s in sorted(set(list(s_values) + [s_central])):
        sk = min(K, key=lambda a: abs(np.log(a / ((4 * s) if tube else s)))); r = shape(K[sk], z50, fill, spot)
        if r is not None: out[s] = r
    if s_central not in out: return dict(id=cid, kind="not predicted", note="half depth outside the simulated exposures")
    sl = np.array([v[0] for v in out.values()]); r8 = np.array([v[1] for v in out.values()]); r2 = np.array([v[2] for v in out.values()]); w = 1.13 * (1.10 if raster else 1.05)
    return dict(id=cid, kind=kind, unit_nm=unit_nm, geometry="pore" if tube else "slit", fill=round(fill, 3), z50=round(z50, 1), spot_units=round(spot, 1), s_central=s_central, s_low=min(out), s_high=max(out),
                slope_pred=round(float(out[s_central][0]), 5), slope_lo=round(float(sl.min()), 5), slope_hi=round(float(sl.max()), 5), pass_lo=round(float(sl.min() / w), 5), pass_hi=round(float(sl.max() * w), 5),
                r80_pred=round(float(out[s_central][1]), 3), r80_lo=round(float(r8.min()) - 0.03, 3), r80_hi=round(float(r8.max()) + 0.03, 3), r20_pred=round(float(out[s_central][2]), 3), r20_lo=round(float(r2.min()) - 0.03, 3), r20_hi=round(float(r2.max()) + 0.03, 3),
                narrowing_limit_slope=round((1 - fill / 2) / (fill * z50), 5), note=note)

ALL = (0.0025, 0.004, 0.0063, 0.01, 0.0158, 0.025, 0.04, 0.06); rows = []
# 1. Gao 2015 Fig. 5(b): TiO2 from TiCl4/H2O, three gaps, about 100 nm of film. s unknown here (Ylilammi 2018 fitted c = 0.10 for this chemistry on another chip).
rows.append(predict("gao2015_fig5b_TiO2_110C_gap200nm", 200, "ceiling", 5, ALL, 0.06, False, "film closes the gap: slope set by narrowing"))
rows.append(predict("gao2015_fig5b_TiO2_110C_gap500nm", 500, "ceiling", 5, ALL, 0.06, False, "41 % filled"))
rows.append(predict("gao2015_fig5b_TiO2_110C_gap1000nm", 1000, "range", 5, ALL, 0.06, False, "20 % filled: depends on s; central value at the top of the simulated range because the published TiCl4 value is 0.10"))
# 2. Elam 2003: ZnO from DEZ/H2O in alumina pores, cycles chosen to fill 40 % of the diameter; EPMA resolution 3-5 um (4 um used); s unknown
for cid, d in (("elam2003_fig6a_ZnO_d19nm_5s", 19), ("elam2003_fig6a_ZnO_d46nm_5s", 46), ("elam2003_fig6a_ZnO_d65nm_5s", 65), ("elam2003_fig5a_ZnO_d65nm_exposure1s", 65), ("elam2003_fig5a_ZnO_d65nm_exposure3s", 65), ("elam2003_fig5a_ZnO_d65nm_exposure10s", 65)):
    rows.append(predict(cid, d, "ceiling", 4, [s / 4 for s in ALL], 0.06 / 4, True, "pore, 40 % filled; spot 4 um; the pore is open at both faces 50 um apart, the far front is ignored", tube=True))
# 3. Arts 2019 Fig. 3(a): H2O-limited, s of water not transferable: upper bound only
for cid in ("arts2019_fig3a_Al2O3_H2Olim_Tset150C", "arts2019_fig3a_Al2O3_H2Olim_Tset275C", "arts2019_fig3a_Al2O3_H2Olim_Tset400C"):
    c = BY[cid]; f = 2 * 46 / 500; z50 = val(c["x50_um"]) * 2; lim = (1 - f / 2) / (f * z50)
    rows.append(dict(id=cid, kind="bound", unit_nm=500, geometry="slit", fill=round(f, 3), z50=round(z50, 1), spot_units=20.0, slope_pred=None, slope_lo=0.0, slope_hi=round(lim, 5), pass_lo=0.0, pass_hi=round(lim * 1.13 * 1.10, 5), narrowing_limit_slope=round(lim, 5), note="upper bound only: slope cannot exceed the narrowing limit"))
# 4. ZnO in PillarHall: s for DEZ taken from one published fit (Haimi 2026 sample B, 200 C: c = 0.002), transferred to other samples and to another laboratory
ZS = (0.0025,)   # 0.002 is below the simulated range; the nearest simulated value is used and the range covers a factor of 2.5 either side
for cid, raster in (("gonsalves2026_fig2_ZnO_DEZ_H2O_200C_left", False), ("gonsalves2026_fig2_ZnO_DEZ_H2O_200C_right", False), ("gonsalves2026_fig2_ZnO_DEZ_H2O_175C_left", False), ("gonsalves2026_fig2_ZnO_DEZ_H2O_175C_right", False), ("gonsalves2026_fig2_ZnO_DEZ_H2O_150C_left", False), ("gonsalves2026_fig2_ZnO_DEZ_H2O_150C_right", False),
                    ("haimi2026_fig5_ZnO_sampleA_125C_125cyc", True), ("haimi2026_fig5_ZnO_sampleC_250C_125cyc", True), ("haimi2026_fig6_ZnO_sampleD_125C_56cyc_DEZ1s", True), ("haimi2026_fig6_ZnO_sampleE_125C_56cyc_DEZ0.5s", True)):
    rows.append(predict(cid, 500, "transfer", 1, (0.0025, 0.004, 0.0063), 0.0025, raster, "s = 0.002 from Haimi 2026 sample B assumed to transfer (nearest simulated 0.0025; range 0.0025-0.0063); hydraulic transport"))
for cid in ("gonsalves2026_fig2_ZnO_MeOH_DEZ_H2O_150C_left", "gonsalves2026_fig2_ZnO_MeOH_DEZ_H2O_150C_right", "gonsalves2026_fig2_ZnO_MeOH_DEZ_H2O_175C_left", "gonsalves2026_fig2_ZnO_MeOH_DEZ_H2O_175C_right", "gonsalves2026_fig2_ZnO_MeOH_DEZ_H2O_200C_left", "gonsalves2026_fig2_ZnO_MeOH_DEZ_H2O_200C_right"):
    rows.append(predict(cid, 500, "range", 1, ALL, 0.0025, False, "inhibited process: s unknown, range over the simulated s (weak test)"))
D = pd.DataFrame(rows); done = set(D.id)
EXC = {"dendooven": "macroscopic slot 5 mm wide with gap 0.1-0.5 mm inferred from the aspect ratio, free-molecular at 0.3 Pa: outside the geometries simulated; s unknown", "rose": "holes with four TEM points per curve and, for three of five, no half depth in range",
       "haimi2026_fig5_ZnO_sampleB": "the predictor knows this sample's fitted c and depth (it supplies s for the transfer predictions)", "elam2003_fig5a_ZnO_d65nm_exposure30s": "half depth not reached"}
ex = [dict(id=c["id"], kind="excluded in advance", note=next(v for k, v in EXC.items() if k in c["id"])) for c in INP if c["id"] not in done and any(k in c["id"] for k in EXC)]
left = [c["id"] for c in INP if c["id"] not in done and not any(k in c["id"] for k in EXC)]; assert not left, left
D = pd.concat([D, pd.DataFrame(ex)], ignore_index=True); D.to_csv(ROOT / "docs/blindtest/predictions.csv", index=False)
print(D[["id", "kind", "fill", "z50", "slope_pred", "slope_lo", "slope_hi", "pass_lo", "pass_hi", "r80_pred", "r20_pred", "narrowing_limit_slope"]].to_string(index=False)); print(D.kind.value_counts().to_dict())
