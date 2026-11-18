"""Ballistic counterpart of the Arts et al. (2019) slope relation, from saved Machball trench frames
(docs/paper/ballistic_frames.py). For each sticking coefficient and output frame: the half-coverage position z50 (gap
units) and the secant slope |d theta / d(x/H)| between theta = 0.6 and 0.4, both with linearly interpolated crossings.
End-wall control: the trenches are long enough that the end wall is at least 500 heights beyond every front used
  s = 2e-4 ... 1e-3 : AR 1700, 1700 sections, fronts to z50 = 1150 (and at least 1.5 front widths from the end wall;
                      the s = 1e-4 run is not used: its front is wider than the trench)
  s = 2e-3, 4e-3    : AR 1100, 1100 sections, fronts to 560
  s = 8e-3 ... 6.4e-2: AR 1000, 1000 sections, fronts to 500
(an earlier AR-600 calibration was contaminated by the end wall for fronts beyond about 350 heights; section size 1 H
changes A by about 1 % against 0.25 H). The continuum constant is checked with this repository's 1D slit model.
Writes figures/ballistic_slope_calibration.csv.   PYTHONPATH=src .venv/bin/python docs/paper/ballistic_slope_calibration.py
"""
import sys, csv, glob
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
from aldconf import Slit, Precursor, Process, simulate
N0, T, MASS, GAP = 4e18, 573.15, 72.09, 500e-9
def cross(z, th, level):
    i = int(np.argmax(th < level))
    if i == 0: return None
    return z[i - 1] + (th[i - 1] - level) / (th[i - 1] - th[i]) * (z[i] - z[i - 1])
def front(z, th):
    c = [cross(z, th, l) for l in (0.6, 0.5, 0.4)]
    if any(v is None for v in c) or c[2] <= c[0]: return None
    return c[1], 0.2 / (c[2] - c[0])
rows = []
files = [f for nz in (1700, 1100, 1000) for f in glob.glob(str(ROOT / f"figures/ballistic_frames/s_*_Nz{nz}.npz"))]
for f in sorted(files, key=lambda f: float(np.load(f)["s"])):
    d = np.load(f); s = float(d["s"]); AR = float(d["AR"]); z = d["z"]
    for e, th in zip(d["exposure_Pa_s"], d["cov"]):
        fr = front(z, th)
        if fr is None or fr[0] < 60 or AR - fr[0] < 500 or s < 2e-4 or AR - fr[0] < 1.5 / fr[1]: continue   # end wall at least 500 H and 1.5 front widths (1/slope) away
        rows.append(dict(model="ballistic", s=s, AR_trench=AR, exposure_Pa_s=float(e), z50=fr[0], slope=fr[1], A=s / fr[1] ** 2))
    print(f"s {s:g} (AR {AR:g}): {sum(1 for r in rows if r['s'] == s and r['model'] == 'ballistic')} frames", flush=True)
for s in (1e-3, 4e-3, 1.6e-2, 6.4e-2):
    for expo in (0.3, 1.0, 3.0, 8.0, 20.0):
        pr = Precursor(sticking=s, site_density=N0, mass_amu=MASS, gpc=1e-10)
        res = simulate(Process(Slit(GAP, 1000 * GAP, knudsen_correction="none"), pr, T=T, p_partial=expo / 0.1, t_dose=0.1, cycles=1, n_cells=2000, n_substeps=300, narrowing=False))
        fr = front(res.x / GAP, res.theta)
        if fr is not None and 60 <= fr[0] <= 500: rows.append(dict(model="1d_hydraulic", s=s, AR_trench=1000, exposure_Pa_s=expo, z50=fr[0], slope=fr[1], A=s / fr[1] ** 2))
with open(ROOT / "figures" / "ballistic_slope_calibration.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("done", len(rows))
