"""POST HOC: is the signal behind the front a film or an offset of the reflectometer? For every reflectometry curve released in
the second test (entering and not entered), the level behind the front in nanometres: median and maximum of the thickness
between 1.25 and 1.6 half depths, with the plateau thickness of inputs.json or Yim 2020 Table 1; and the positions of the
local maxima behind the front (to compare with the pillar rows: inter-distance 49 um, rows 42.4 um apart for layout v1b).
A film left by the process would scale with the number of cycles; an instrument offset would be the same in nanometres.
Writes residual_nm.csv."""
import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[4]; INP = json.loads((ROOT / "docs/blindtest/v2/inputs.json").read_text()); C = {c["id"]: c for c in INP["curves"]}
T = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv").set_index("sample"); rows = []
for f in sorted(list((ROOT / "data/blind_revealed_v2").glob("*.csv")) + list((ROOT / "data/blind_revealed_v2/not_entered").glob("yim*.csv")) + list((ROOT / "data/blind_revealed_v2/not_entered").glob("haimi*opt*.csv"))):
    cid = f.stem; d = pd.read_csv(f).dropna(subset=["x_um", "thickness_over_plateau"]).sort_values("x_um"); x = d.x_um.values; y = d.thickness_over_plateau.values
    if cid in C: t = C[cid].get("plateau_thickness_nm_published") or 108.0; cyc = C[cid]["cycles"]; gap = C[cid]["gap_height_nm"]; x50 = C[cid]["half_thickness_depth_x50_um"]
    else: n = int(cid.split("_s")[1]); t = T.loc[n, "S_IIb_nm"]; cyc = int(T.loc[n, "cycles"]); gap = int(T.loc[n, "H_nm"]); x50 = float(T.loc[n, "x50_um"])
    m = (x >= 1.25 * x50) & (x <= 1.6 * x50)
    if m.sum() < 3: rows.append(dict(id=cid, cycles=cyc, gap_nm=gap, plateau_nm=t, n_points_behind=int(m.sum()))); continue
    ym = y[m]; xm = x[m]; pk = [float(xm[i]) for i in range(1, len(ym) - 1) if ym[i] > ym[i - 1] and ym[i] >= ym[i + 1] and ym[i] > np.median(ym) + 0.05]
    rows.append(dict(id=cid, cycles=cyc, gap_nm=gap, plateau_nm=t, n_points_behind=int(m.sum()), median_fraction=round(float(np.median(ym)), 3), median_nm=round(float(np.median(ym) * t), 1), max_nm=round(float(ym.max() * t), 1), peaks_um=" ".join(f"{p:.0f}" for p in pk[:4])))
D = pd.DataFrame(rows); D.to_csv(ROOT / "docs/blindtest/v2/posthoc/residual_nm.csv", index=False); print(D.to_string(index=False))
g = D.dropna(subset=["median_nm"]); print("median level behind the front, nm, by cycles:"); print(g.groupby("cycles").median_nm.agg(["count", "median", "min", "max"]))
print("correlation of median_nm with cycles:", round(float(np.corrcoef(g.cycles, g.median_nm)[0, 1]), 2), "| of median_fraction with cycles:", round(float(np.corrcoef(g.cycles, g.median_fraction)[0, 1]), 2))
