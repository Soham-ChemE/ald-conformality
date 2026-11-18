"""Monte Carlo of one ALD exposure in the real PillarHall cavity geometry (lib/pillarmc: exact plates, cylindrical pillars
on a hexagonal lattice, optional diffuse side walls), to test whether the pillars or the finite width remove the
ideal-slit ballistic tail and the logarithmic transport enhancement.
  PYTHONPATH=src .venv/bin/python docs/paper/run_pillar_mc.py run     (builds lib/pillarmc if needed; minutes per case)
  PYTHONPATH=src .venv/bin/python docs/paper/run_pillar_mc.py fit     (fits every geometry to the two measured profiles)
Cases (gap 0.5 um, cavity 400 um deep, pillar diameter 4 um, Yim 2020 ESI Table S2):
  ideal               no pillars, unbounded width (the Machball trench; validation)
  pillars49           standard layout v1b, 49 um inter-distance, unbounded (periodic) width
  walls90             no pillars, diffuse side walls 90 um apart (the W90 L1mm channel opening)
  walls90_pillars49   both: the nominal real cavity
  pillars28           densest layout v1a
  pillars10           10 um inter-distance (first-generation structures of Gao 2015; pillar diameter assumed 4 um)
Raw frames go to figures/pillar_mc/ (not tracked); the summary is figures/pillar_mc_summary.json and pillar_mc_fits.csv.
"""
import sys, subprocess, json, csv
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
EXE = ROOT / "lib/pillarmc/pillarmc"; OUT = ROOT / "figures/pillar_mc"; OUT.mkdir(exist_ok=True)
CASES = {"ideal": (0, 2, 0), "pillars49": (49, 2, 0), "walls90": (0, 2, 90), "walls90_pillars49": (49, 2, 90), "pillars28": (28, 2, 0), "pillars10": (10, 2, 0)}
S_GRID = (0.004, 0.0063, 0.01, 0.0158, 0.025, 0.04)
E_END, NSTEP, M, H, L = 6.0, 300, 400000, 0.5, 400
def run():
    if not EXE.exists(): subprocess.check_call(["clang++", "-O3", "-std=c++17", "-o", str(EXE), str(ROOT / "lib/pillarmc/pillarmc.cpp")])
    for name, (a, R, W) in CASES.items():
        for s in S_GRID:
            f = OUT / f"{name}_s{s:g}.csv"
            if f.exists(): continue
            # the hydraulic-like cases need more exposure for the same depth: give every run the same generous end exposure scaled by 1/s^0
            subprocess.check_call([str(EXE), f"{s:g}", f"{E_END * (2.5 if a and a <= 28 else 1.0):g}", str(NSTEP), str(M), str(H), str(L), str(a), str(R), str(W), "1", "8", str(f)], stderr=subprocess.DEVNULL)
            print("done", f.name, flush=True)
def load(f):
    d = np.loadtxt(f, delimiter=",", comments="#"); return d[:, 0], d[:, 1:]
def front(z, th):
    def cross(l):
        i = int(np.argmax(th < l)); return None if i == 0 else z[i - 1] + (th[i - 1] - l) / (th[i - 1] - th[i]) * (z[i] - z[i - 1])
    c = [cross(l) for l in (0.6, 0.5, 0.4)]
    return None if any(v is None for v in c) else (c[1], 0.2 / (c[2] - c[0]))
def at(E, TH, e):
    i = min(max(np.searchsorted(E, e), 1), len(E) - 1); w = (e - E[i - 1]) / (E[i] - E[i - 1]); return (1 - w) * TH[i - 1] + w * TH[i]
def measured(ds):
    from aldconf.data import load_pillarhall, default_path
    if ds == "zenodo":
        df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
        sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
        for sh in sheets:
            g = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
        mean = np.nanmean(np.array(prof), 0); y = mean / np.nanmean(mean[(grid > 5) & (grid < 60)]); ok = ~np.isnan(y) & (grid <= 150); return grid[ok] * 2, y[ok]
    d = pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv"); return d.x_um.values * 2, d.thickness_nm.values / d.thickness_nm[(d.x_um > 5) & (d.x_um < 60)].mean()
def fit():
    from scipy.optimize import minimize_scalar
    rows = []; summary = {"cases": {}, "s_grid": list(S_GRID)}
    for name in CASES:
        runs = []
        for s in S_GRID:
            f = OUT / f"{name}_s{s:g}.csv"
            if f.exists(): E, TH = load(f); runs.append((s, E, TH))
        if not runs: continue
        z = np.arange(runs[0][2].shape[1]) + 0.5
        # transport diagnostics at s = 0.0158: exposure needed to put the front at 215 heights, slope constant A and tail there
        diag = {}
        for s, E, TH in runs:
            z50 = np.array([(front(z, th) or (np.nan, np.nan))[0] for th in TH]); ok = np.isfinite(z50)
            if ok.sum() < 5 or np.nanmax(z50) < 215: continue
            e215 = float(np.interp(215, z50[ok], E[ok])); fr = front(z, at(E, TH, e215)); th = at(E, TH, e215)
            diag[f"{s:g}"] = dict(exposure_for_z50_215=round(e215, 3), A=round(s / fr[1] ** 2, 1), tail_at_z50_plus_60=round(float(np.interp(275, z, th)), 4), tail_at_z50_plus_100=round(float(np.interp(315, z, th)), 4))
        summary["cases"][name] = dict(geometry=dict(a_um=CASES[name][0], R_um=CASES[name][1], W_um=CASES[name][2]), diagnostics=diag, fits={})
        for ds in ("zenodo", "ylilammi_al2o3"):
            xg, yy = measured(ds); tab = []
            for s, E, TH in runs:
                def cost(le): return float(np.mean((np.interp(xg, z, at(E, TH, 10 ** le)) - yy) ** 2))
                r = minimize_scalar(cost, bounds=(np.log10(E[2]), np.log10(E[-2])), method="bounded", options={"xatol": 1e-4}); tab.append((s, r.fun, 10 ** r.x))
            tab = np.array(tab); j = min(max(int(np.argmin(tab[:, 1])), 1), len(tab) - 2); p = np.polyfit(np.log(tab[j - 1:j + 2, 0]), tab[j - 1:j + 2, 1], 2)
            ls = -p[1] / (2 * p[0]) if p[0] > 0 else np.log(tab[j, 0]); ls = min(max(ls, np.log(tab[j - 1, 0])), np.log(tab[j + 1, 0])); fmin = float(np.polyval(p, ls)) if p[0] > 0 else tab[j, 1]
            edge = bool(int(np.argmin(tab[:, 1])) in (0, len(tab) - 1))
            sn = tab[int(np.argmin(tab[:, 1]))]; srun = [r_ for r_ in runs if r_[0] == sn[0]][0]; th = at(srun[1], srun[2], sn[2]); m = np.interp(xg, z, th)
            reg = {k: float(np.sqrt(np.mean((m[sel] - yy[sel]) ** 2))) for k, sel in (("plateau", yy > 0.9), ("front", (yy > 0.1) & (yy <= 0.9)), ("tail", yy <= 0.1))}
            rows.append(dict(geometry=name, dataset=ds, s=float(np.exp(ls)), rms=float(np.sqrt(max(fmin, 0))), s_best_grid=float(sn[0]), rms_best_grid=float(np.sqrt(sn[1])), exposure_Pa_s=float(sn[2]), rms_plateau=reg["plateau"], rms_front=reg["front"], rms_tail=reg["tail"], at_grid_edge=edge))
            summary["cases"][name]["fits"][ds] = dict(s=float(f"{np.exp(ls):.3g}"), rms=round(float(np.sqrt(max(fmin, 0))), 4), exposure_Pa_s=round(float(sn[2]), 3), rms_front=round(reg["front"], 4), rms_tail=round(reg["tail"], 4), at_grid_edge=edge,
                                                       curve=dict(AR=[float(v) for v in np.arange(0, 402, 4.0)], y=[round(float(v), 4) for v in np.interp(np.arange(0, 402, 4.0), z, th)]))
            print(rows[-1], flush=True)
    with open(ROOT / "figures/pillar_mc_fits.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    (ROOT / "figures/pillar_mc_summary.json").write_text(json.dumps(summary, indent=1)); print(json.dumps({k: v["diagnostics"] for k, v in summary["cases"].items()}, indent=1))
if __name__ == "__main__":
    {"run": run, "fit": fit}[sys.argv[1]]()
