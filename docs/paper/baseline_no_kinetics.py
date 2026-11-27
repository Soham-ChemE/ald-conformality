"""A baseline without kinetics for the held-out slope test of Section 3.7: slope = k (1 - f/2) / (f z50), the step-front form
with one constant k fitted in ln(slope) (k = 1 is the step-front value itself). No transport model and no sticking coefficient.
Samples and slopes as in figures/mechanism_comparison.csv (complete gas model rows); f = 2 film / H with the gap of that file.
Writes figures/baseline_no_kinetics.csv.   .venv/bin/python docs/paper/baseline_no_kinetics.py"""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; C = pd.read_csv(ROOT / "figures/mechanism_comparison.csv")
M = C[(C.model == C.model.unique()[-1]) & (C["sample"] != "Arts")].drop_duplicates("sample").copy()   # the 15 scored held-out samples
M["fill"] = 2 * M.film_nm / M.H_nm; M["step"] = (1 - M.fill / 2) / (M.fill * M.z50); M["ln_ratio_k1"] = np.log(M.step / M.measured_slope)
rows = []
for name, S in (("15 scored", M), ("13 that every model reaches (500 nm gap)", M[M.H_nm == 500]), ("500 nm gap, 50x objective or EDS", M[(M.H_nm == 500) & ((M.magnification == 50) | M["sample"].str.startswith("Haimi"))])):
    lnk = -S.ln_ratio_k1.mean(); r = S.ln_ratio_k1 + lnk; rows.append(dict(set=name, n=len(S), k_fitted=round(float(np.exp(lnk)), 3), rms_ln_slope_k_fitted=round(float(np.sqrt((r ** 2).mean())), 3), rms_ln_slope_k_1=round(float(np.sqrt((S.ln_ratio_k1 ** 2).mean())), 3)))
D = pd.DataFrame(rows); D.to_csv(ROOT / "figures/baseline_no_kinetics.csv", index=False); print(C.model.unique()); print(M[["sample", "H_nm", "magnification", "fill", "z50", "measured_slope", "step"]].round(4).to_string(index=False)); print(D.to_string(index=False))
