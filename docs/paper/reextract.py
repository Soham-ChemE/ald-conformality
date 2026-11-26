"""Phase 3: every lateral-channel TMA sticking coefficient re-extracted with one physics.
For each published front (Yim 2020 Table 1, Haimi 2021, Arts 2019) the measured slope, with a 13 % one-sigma error (9 % for the
digitised Arts front), is inverted through the complete forward model (transport + channel narrowing + spot) of
mechanism_test.py with the measured gaps (figures/mechanism_slope_curves_measgap.json); the two full profiles (Zenodo 3974438,
Ylilammi 2018 Fig. 6) come from the profile calibration (s range within 5 % of the best rms). Three transport families.
Then: does the sample-to-sample spread shrink, and is there one s inside every sample's interval?
Writes figures/reextraction.csv and figures/reextraction_summary.json.
  PYTHONPATH=src .venv/bin/python docs/paper/reextract.py
"""
import sys, json, os
from pathlib import Path
import numpy as np, pandas as pd
os.environ["MECH_GAPS"] = "measured"
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
FIG = ROOT / "figures"; cur = json.loads((FIG / "mechanism_slope_curves_measgap.json").read_text()); cal = pd.read_csv(FIG / "mechanism_calibration_measgap.csv").set_index("model"); H = mt.heldout().set_index("sample")
FAM = {"hydraulic": "hydraulic+narrowing+spot", "gas300": "gas+narrowing+spot", "free": "ballistic+narrowing+spot"}
PUB = {"HaimiA": ("1.6e-3 (slope method, Haimi 2021)", 1.6e-3), "Arts": ("(0.5-2)e-3 (slope method, Arts 2019)", 1.0e-3)}
def invert(c, slope, err):
    ls, lsl = np.log(c["s"]), np.maximum.accumulate(np.log(c["slope"]))
    def inv(x): return None if (x >= lsl[-1] or x <= lsl[0]) else float(np.exp(np.interp(x, lsl, ls)))
    x = np.log(slope); return inv(x), inv(x - err), inv(x + err), float(slope / c["slope"][-1]), (c["s"][0], c["s"][-1])
rows = []
for smp, r in H.iterrows():
    err = 0.09 if smp == "Arts" else 0.13; row = dict(dataset=("Yim 2020 sample " + smp[3:]) if smp.startswith("Yim") else {"HaimiA": "Haimi 2021 specimen A (EDS)", "Arts": "Arts 2019 TMA-limited (about 220 C)"}[smp], sample=smp, kind="front slope", gap_nm=round(r.H_nm * r.gs, 1), cycles=int(r.cycles),
                                                      fill=round(float(min(2 * r.film_nm / r.H_nm, 1.0)), 3), z50=round(float(r.z50), 1), slope=round(float(r.measured_slope), 5), objective=int(r.magnification), published=PUB.get(smp, ("4e-3 (slope method, mean of samples 2, 4, 8, 12; Yim 2020)", 4e-3))[0], slope_method=float(f"{13.9 * r.measured_slope ** 2:.3g}"))
    for fam, model in FAM.items():
        c = cur.get(model, {}).get(smp)
        if c is None or len(c["s"]) < 3: row.update({f"s_{fam}": None, f"lo_{fam}": None, f"hi_{fam}": None, f"ceiling_{fam}": None, f"note_{fam}": "no kernel for this structure"}); continue
        s, lo, hi, frac, rng = invert(c, r.measured_slope, err)
        note = "above ceiling" if (s is None and frac >= 1) else "below simulated range" if s is None else "; ".join(n for n, t in (("no upper bound", hi is None), ("no lower bound", lo is None)) if t)
        row.update({f"s_{fam}": None if s is None else float(f"{s:.3g}"), f"lo_{fam}": None if lo is None else float(f"{lo:.3g}"), f"hi_{fam}": None if hi is None else float(f"{hi:.3g}"), f"ceiling_{fam}": round(frac, 2), f"note_{fam}": note})
    rows.append(row)
for ds, col, name, pub in (("zenodo", "zenodo", "Zenodo 3974438 (Yim 2020 sample 8), full profile", "5.0e-3 (slope) / 5.4e-3 (Klabbers 2026)"), ("ylilammi", "ylilammi", "Ylilammi 2018 Fig. 6, full profile", "5.72e-3 (fit, Ylilammi 2018)")):
    row = dict(dataset=name, sample=ds, kind="full profile", gap_nm=519.1 if ds == "zenodo" else 500.0, cycles=500, fill=0.197 if ds == "zenodo" else 0.208, published=pub, slope_method=float(f"{cal.loc['hydraulic', 's_' + col]:.3g}"))
    for fam, model in FAM.items():
        g = cal.loc[model]; s = np.array([float(v) for v in g.s_grid.split()]); rms = np.array([float(v) for v in g["rms_" + col + "_grid"].split()]); ok = rms <= 1.05 * rms.min()
        row.update({f"s_{fam}": float(f"{g['s_' + col]:.3g}"), f"lo_{fam}": None if ok[0] else float(s[ok][0]), f"hi_{fam}": None if ok[-1] else float(s[ok][-1]), f"ceiling_{fam}": None, f"note_{fam}": "; ".join(n for n, t in (("no upper bound", ok[-1]), ("no lower bound", ok[0])) if t)})
    rows.append(row)
D = pd.DataFrame(rows); D.to_csv(FIG / "reextraction.csv", index=False)
# ---- does the spread collapse? 300 C TMA samples only (Arts is at about 220 C) ----
T = D[D["sample"] != "Arts"]; summ = {"n_300C": int(len(T)), "slope_method": {}, "complete": {}}
sl = T[T.kind == "front slope"]; v = np.log(sl.slope_method.astype(float)); summ["slope_method"] = dict(n=int(len(v)), min=float(f"{np.exp(v.min()):.3g}"), max=float(f"{np.exp(v.max()):.3g}"), max_over_min=round(float(np.exp(v.max() - v.min())), 1), sd_ln=round(float(v.std(ddof=1)), 3))
for fam in FAM:
    b = T[T[f"s_{fam}"].notna()]; v = np.log(b[f"s_{fam}"].astype(float)); lo = T[f"lo_{fam}"].dropna().astype(float); hi = T[f"hi_{fam}"].dropna().astype(float)
    both = sl[sl[f"s_{fam}"].notna()]; v0 = np.log(both.slope_method.astype(float)); v1 = np.log(both[f"s_{fam}"].astype(float))
    summ["complete"][fam] = dict(n_with_value=int(len(b)), n_above_ceiling=int((T[f"note_{fam}"] == "above ceiling").sum()), n_no_kernel=int((T[f"note_{fam}"] == "no kernel for this structure").sum()), min=float(f"{np.exp(v.min()):.3g}"), max=float(f"{np.exp(v.max()):.3g}"), max_over_min=round(float(np.exp(v.max() - v.min())), 1), sd_ln=round(float(v.std(ddof=1)), 3),
                                 same_samples_slope_method_sd_ln=round(float(v0.std(ddof=1)), 3), same_samples_complete_sd_ln=round(float(v1.std(ddof=1)), 3), same_samples_n=int(len(both)),
                                 largest_lower_bound=float(f"{lo.max():.3g}") if len(lo) else None, smallest_upper_bound=float(f"{hi.min():.3g}") if len(hi) else None,
                                 one_s_fits_all=bool(len(lo) == 0 or len(hi) == 0 or lo.max() <= hi.min()),
                                 samples_setting_the_bounds=dict(lower=T.loc[T[f"lo_{fam}"].astype(float).idxmax(), "sample"] if len(lo) else None, upper=T.loc[T[f"hi_{fam}"].astype(float).idxmin(), "sample"] if len(hi) else None))
(FIG / "reextraction_summary.json").write_text(json.dumps(summ, indent=1))
show = ["dataset", "fill", "slope_method"] + [f"{k}_{f}" for f in FAM for k in ("s", "lo", "hi")]
print(D[show].to_string(index=False)); print(json.dumps(summ, indent=1))
