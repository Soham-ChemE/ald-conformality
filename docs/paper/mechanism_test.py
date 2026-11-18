"""Which transport physics explains the published lateral-channel (PillarHall) Al2O3 profiles? Calibrate on two full
profiles, freeze, and predict the held-out front slopes.

Transport kernels (single-exposure coverage theta(x/H; exposure, s) of an irreversible Langmuir film):
  hydraulic     1D diffusion with the hydraulic-diameter Knudsen diffusivity (the model behind the published slope method)
  ballistic     Monte Carlo, real cavity (plates, pillars, side walls), free molecular            (figures/pillar_mc)
  gas           the same Monte Carlo with hard-sphere TMA-N2 collisions at the stated reactor pressure
Each kernel is used with and without channel narrowing. Narrowing is not fitted: the film of every earlier cycle has
narrowed the gap, the per-cycle exposure is fixed, and in a diffusion-limited slit the penetration scales with the gap, so
the profile after N cycles is the average of the single-exposure profile stretched by f = H_k/H_0, f uniform between
1 - 2 t_film/H_0 and 1 (checked against the full multi-cycle continuum model: max difference 0.011, 0.024 at pinch-off).

  calibrate   one sticking coefficient per model from the Zenodo 3974438 mean profile and Ylilammi 2018 Fig. 6 (exposures free)
  predict     front slope at half coverage of each held-out sample at its measured penetration, with the calibrated s frozen;
              writes figures/mechanism_predictions.csv (no measured slopes in it)
  compare     joins the predictions with the measured slopes (Yim 2020 Table 1, Haimi 2021) and scores each model;
              writes figures/mechanism_comparison.csv and figures/mechanism_summary.json
  PYTHONPATH=src .venv/bin/python docs/paper/mechanism_test.py calibrate|predict|compare
"""
import sys, json, csv
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar, brentq
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
MC = ROOT / "figures/pillar_mc"; FIG = ROOT / "figures"
KB = 1.380649e-23
S_HYD = (0.0025, 0.004, 0.005, 0.0063, 0.008, 0.01, 0.0125, 0.0158, 0.02, 0.025, 0.03, 0.04, 0.06)


def at(E, TH, e):
    i = min(max(np.searchsorted(E, e), 1), len(E) - 1); w = (e - E[i - 1]) / (E[i] - E[i - 1]); return (1 - w) * TH[i - 1] + w * TH[i]


def mc_kernel(case):
    out = {}
    for f in sorted(MC.glob(f"{case}_s*.csv")):
        lines = [l for l in f.read_text().split("\n") if l and not l.startswith("#")]
        if len(lines) < 100: continue
        n = lines[0].count(","); d = np.array([np.fromstring(l, sep=",") for l in lines if l.count(",") == n and l[-1].isdigit()])
        if len(d) < 100 or len(d) not in (150, 200, 300): continue                       # an unfinished run is not used
        out[float(f.stem.split("_s")[-1])] = (d[:, 0], d[:, 1:], np.arange(d.shape[1] - 1) + 0.5)
    return dict(sorted(out.items()))


def hydraulic_kernel(deep=False):
    """Single-exposure coverage of the 1D hydraulic-diameter model on the same exposure ladder as the Monte Carlo
    (deep: to 1200 gaps and 90 Pa s, for the deep front of Arts 2019)."""
    cache = FIG / ("mechanism_hydraulic_kernel_deep.npz" if deep else "mechanism_hydraulic_kernel.npz")
    if cache.exists():
        z = np.load(cache); return {float(s): (z["E"], z[f"th_{s}"], z["z"]) for s in z["s"]}
    from aldconf import Slit, Precursor, Process, simulate
    E = (90.0 if deep else 16.0) * (np.arange(1, 121) / 120.0) ** 2; zz = np.arange(1200 if deep else 400) + 0.5; out = {}
    for s in S_HYD:
        TH = []
        for e in E:
            pr = Precursor(sticking=s, site_density=4e18, mass_amu=72.09, gpc=1e-10)
            r = simulate(Process(Slit(500e-9, 600e-6 if deep else 200e-6, knudsen_correction="none"), pr, T=573.15, p_partial=e / 0.1, t_dose=0.1, cycles=1, n_cells=1200 if deep else 400, n_substeps=60, knudsen_factor=1.0, narrowing=False))
            TH.append(np.interp(zz, r.x / 500e-9, r.thickness / 1e-10))
        out[s] = (E, np.array(TH), zz)
    np.savez_compressed(cache, E=E, z=zz, s=np.array(list(out)), **{f"th_{s}": v[1] for s, v in out.items()})
    return out


def cumulative(run, e, fill, x, spot=0.0):
    """Normalised thickness after many cycles at reduced distance x (in units of the initial gap), averaged over a
    measurement spot of full width `spot` (same units; top-hat)."""
    E, TH, z = run; th = at(E, TH, e)
    if spot > 0: return np.mean([cumulative(run, e, fill, np.abs(x + d)) for d in np.linspace(-0.5, 0.5, 9) * spot], 0)
    if fill <= 0: return np.interp(x, z, th, right=0.0)
    fs = np.linspace(max(1 - fill, 0.02), 1.0, 21)
    return np.mean([np.interp(x / f, z, th, right=0.0) for f in fs], 0)


def front(run, e, fill, spot=0.0):
    xx = np.arange(0.5, run[2][-1], 0.5 if run[2][-1] < 500 else 1.0); y = cumulative(run, e, fill, xx, spot); y = y / y[:4].mean()
    def cross(l):
        i = int(np.argmax(y < l)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - l) / (y[i - 1] - y[i]) * (xx[i] - xx[i - 1])
    c = [cross(l) for l in (0.6, 0.5, 0.4)]; return c[1], 0.2 / (c[2] - c[0])


def slope_at_depth(run, z50, fill, spot=0.0):
    """Front slope when the exposure is chosen so that the half-thickness depth equals the measured one."""
    E = run[0]; g = lambda le: front(run, 10 ** le, fill, spot)[0] - z50
    lo, hi = np.log10(E[3]), np.log10(E[-2])
    if not (np.isfinite(g(lo)) and np.isfinite(g(hi))) or g(lo) * g(hi) > 0: return np.nan, np.nan
    le = brentq(g, lo, hi, xtol=1e-4); return front(run, 10 ** le, fill, spot)[1], 10 ** le


def measured_profile(ds):
    import run_pillar_mc as r
    x, y = r.measured(ds)
    fill = 2 * (51.1 if ds == "zenodo" else float(pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv").query("x_um > 5 and x_um < 60").thickness_nm.mean())) / 500.0
    return x, y, fill


# model name: (kernel for the calibration structure, narrowing)
MODELS = {"hydraulic": ("hyd", False), "hydraulic+narrowing": ("hyd", True), "ballistic": ("walls90_pillars49", False), "ballistic+narrowing": ("walls90_pillars49", True),
          "gas": ("gas300", False), "gas+narrowing": ("gas300", True)}
SPOT_UM = {50: 5.5, 10: 25.0, 5: 50.0, 0: 1.0, 100: 10.0}   # reflectometer objective -> spot (Yim 2020: 5-6, 25, 50 um); 0 = EDS line scan, 1 um step; 100 = Arts 2019 (10 um)
FREE_KERNEL = {2000: "h2000_free", 100: "h100_free"}  # gas-free Monte Carlo of the other gap heights


def kernel(name):
    return hydraulic_kernel() if name == "hyd" else hydraulic_kernel(True) if name == "hyd_deep" else mc_kernel(name)


def calibrate():
    rows = []; data = {ds: measured_profile(ds) for ds in ("zenodo", "ylilammi_al2o3")}
    for model, (kname, narrow) in MODELS.items():
        K = kernel(kname); tab = []
        for s, run in K.items():
            mse = {}
            for ds, (x, y, fill) in data.items():
                c = lambda le: float(np.mean((cumulative(run, 10 ** le, fill if narrow else 0, x, 11.0 if narrow else 0.0) - y) ** 2))
                q = minimize_scalar(c, bounds=(np.log10(run[0][3]), np.log10(run[0][-2])), method="bounded", options={"xatol": 1e-4}); mse[ds] = (q.fun, 10 ** q.x)
            tab.append((s, mse["zenodo"][0], mse["ylilammi_al2o3"][0], mse["zenodo"][1], mse["ylilammi_al2o3"][1]))
        tab = np.array(tab)
        def best(col):
            j = int(np.argmin(col)); edge = j in (0, len(col) - 1); j = min(max(j, 1), len(col) - 2); p = np.polyfit(np.log(tab[j - 1:j + 2, 0]), col[j - 1:j + 2], 2)
            ls = np.clip(-p[1] / (2 * p[0]), np.log(tab[j - 1, 0]), np.log(tab[j + 1, 0])) if p[0] > 0 else np.log(tab[j, 0])
            return float(np.exp(ls)), float(np.sqrt(max(np.polyval(p, ls), 0))) if p[0] > 0 else float(np.sqrt(col[j])), bool(edge)
        joint = best(tab[:, 1] + tab[:, 2]); zen = best(tab[:, 1]); yli = best(tab[:, 2])
        rows.append(dict(model=model, kernel=kname, narrowing=narrow, s_joint=joint[0], rms_joint=float(np.sqrt(joint[1] ** 2 / 2)), s_zenodo=zen[0], rms_zenodo=zen[1], s_ylilammi=yli[0], rms_ylilammi=yli[1],
                         at_grid_edge=joint[2] or zen[2] or yli[2], n_s=len(tab), s_grid=" ".join(f"{v:g}" for v in tab[:, 0]),
                         rms_zenodo_grid=" ".join(f"{np.sqrt(v):.4f}" for v in tab[:, 1]), rms_ylilammi_grid=" ".join(f"{np.sqrt(v):.4f}" for v in tab[:, 2])))
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rows[-1].items() if "grid" not in k or k == "at_grid_edge"}, flush=True)
    pd.DataFrame(rows).to_csv(FIG / "mechanism_calibration.csv", index=False)


# held-out samples: Yim 2020 Table 1 (sample 8 is the Zenodo calibration profile and is excluded) and Haimi 2021 specimen A
def heldout():
    t = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv"); rows = []
    for r in t.itertuples():
        if r.sample in (8, 14): continue      # 8 = calibration profile; 14 is the water-limited run (Yim 2020 use it for H2O) and lies beyond the simulated depth
        rows.append(dict(sample=f"Yim{r.sample}", series=r.series, H_nm=r.H_nm, layout=r.pillar_layout, cycles=r.cycles, magnification=r.magnification, film_nm=r.S_IIb_nm, z50=r.x50_tilde, measured_slope=r.slope_dtheta_dxtilde))
    rows.append(dict(sample="HaimiA", series="EDS", H_nm=500, layout="v1b", cycles=1080, magnification=0, film_nm=108.0, z50=195.2, measured_slope=0.0108))
    # Arts 2019 Fig. 2, TMA-limited profile (digitised: data/lhar_profiles/arts2019_fig2_tma_limited.csv; figures/lhar_slope_table.csv).
    # Another reactor (Oxford OpAL, pressure not stated) and 220 C instead of 300 C, so s need not be the same: scored separately.
    rows.append(dict(sample="Arts", series="Arts", H_nm=500, layout="arts", cycles=400, magnification=100, film_nm=46.0, z50=486.1, measured_slope=0.00768))
    return pd.DataFrame(rows)


def gas_kernel_for(row):
    if row.layout == "arts": return "arts_gas27"
    if row.H_nm == 2000: return "h2000_gas300"
    if row.H_nm == 100: return "h100_gas300"
    return {"v1a": "v1a_gas300", "v2a": "v2a_gas300"}.get(row.layout, "gas300")


def predict():
    cal = pd.read_csv(FIG / "mechanism_calibration.csv").set_index("model"); H = heldout(); rows = []; cache = {}; curves = {}
    for model, (kname, narrow) in MODELS.items():
        s_star = float(cal.loc[model, "s_joint"])
        for r in H.itertuples():
            kn = gas_kernel_for(r) if kname == "gas300" else (FREE_KERNEL.get(r.H_nm, kname) if kname != "hyd" else kname)
            if r.layout == "arts": kn = {"hyd": "hyd_deep", "walls90_pillars49": "arts_free"}.get(kname, kn)
            spot = SPOT_UM[r.magnification] * 1000.0 / r.H_nm if narrow else 0.0          # the complete models include the spot; the bare kernels do not
            if kn not in cache: cache[kn] = kernel(kn)
            K = cache[kn]; fill = 2 * r.film_nm / r.H_nm if narrow else 0.0
            if len(K) < 2: rows.append(dict(model=model, sample=r.sample, kernel=kn, s=s_star, predicted_slope=np.nan, exposure_Pa_s=np.nan)); continue
            ss = np.array(list(K)); sl = np.array([slope_at_depth(K[s], r.z50, fill, spot) for s in ss]); ok = np.isfinite(sl[:, 0])
            curves.setdefault(model, {})[r.sample] = dict(kernel=kn, s=[float(v) for v in ss[ok]], slope=[round(float(v), 5) for v in sl[ok, 0]])
            pred = float(np.exp(np.interp(np.log(s_star), np.log(ss[ok]), np.log(sl[ok, 0])))) if ok.sum() >= 2 else np.nan
            rows.append(dict(model=model, sample=r.sample, kernel=kn, s=s_star, predicted_slope=pred, exposure_Pa_s=float(np.interp(np.log(s_star), np.log(ss[ok]), sl[ok, 1])) if ok.sum() >= 2 else np.nan,
                             extrapolated=bool(s_star < ss[ok].min() or s_star > ss[ok].max()) if ok.any() else True))
        print(model, "done", flush=True)
    pd.DataFrame(rows).to_csv(FIG / "mechanism_predictions.csv", index=False)
    (FIG / "mechanism_slope_curves.json").write_text(json.dumps(curves))


def compare():
    P = pd.read_csv(FIG / "mechanism_predictions.csv"); H = heldout().set_index("sample"); P = P.join(H, on="sample")
    P["log_ratio"] = np.log(P.predicted_slope / P.measured_slope)
    # replicate scatter of the measured slope: the nominally identical 50x runs (samples 2, 4, 12, 15, 16; purge has no effect)
    rep = np.log(H.loc[["Yim2", "Yim4", "Yim12", "Yim15", "Yim16"], "measured_slope"].values); sig = float(np.std(rep, ddof=1))
    P.to_csv(FIG / "mechanism_comparison.csv", index=False); summ = {"replicate_sigma_ln_slope": round(sig, 3), "models": {}}
    for m, g in P.groupby("model", sort=False):
        fair = g[g.log_ratio.notna() & (g["sample"] != "Arts")]
        summ["models"][m] = dict(s=float(g.s.iloc[0]), n=int(len(fair)), rms_ln=round(float(np.sqrt(np.mean(fair.log_ratio ** 2))), 3), mean_ln=round(float(fair.log_ratio.mean()), 3),
                                 chi2_per_point=round(float(np.mean((fair.log_ratio / sig) ** 2)), 2), worst=fair.loc[fair.log_ratio.abs().idxmax(), "sample"],
                                 per_sample={r.sample: round(float(r.predicted_slope), 4) for r in g.itertuples()})
        print(m, {k: v for k, v in summ["models"][m].items() if k != "per_sample"})
    # one s for every held-out slope: chi2(s) per model; and how much each slope still depends on s (0.5 = slope ~ sqrt(s), 0 = narrowing-limited)
    curves = json.loads((FIG / "mechanism_slope_curves.json").read_text()); sgrid = np.geomspace(0.003, 0.06, 40); summ["global"] = {}; summ["identifiability"] = {}
    for mname, per in curves.items():
        chi = np.zeros_like(sgrid); n = 0
        for smp, c in per.items():
            if smp == "Arts" or len(c["s"]) < 3: continue
            pred = np.exp(np.interp(np.log(sgrid), np.log(c["s"]), np.log(c["slope"]))); chi += (np.log(pred / H.loc[smp, "measured_slope"]) / sig) ** 2; n += 1
        if n == 0: continue
        j = int(np.argmin(chi)); inside = sgrid[chi <= chi[j] + 1.0]
        summ["global"][mname] = dict(n=n, s_best=float(f"{sgrid[j]:.3g}"), chi2_min_per_point=round(float(chi[j] / n), 2), s_lo=float(f"{inside.min():.3g}"), s_hi=float(f"{inside.max():.3g}"), at_scan_edge=bool(j in (0, len(sgrid) - 1) or inside.min() == sgrid[0] or inside.max() == sgrid[-1]))
        if "narrowing" in mname:
            s0 = float(P[P.model == mname].s.iloc[0]); summ["identifiability"][mname] = {}
            for smp, c in per.items():
                if len(c["s"]) < 3: continue
                ls = np.log(c["s"]); f = lambda x: np.interp(x, ls, np.log(c["slope"])); h = 0.35
                r = H.loc[smp]; fill = 2 * r.film_nm / r.H_nm; z1 = r.z50 / (1 - min(fill, 1) / 2)
                summ["identifiability"][mname][smp] = dict(dlnslope_dlns=round(float((f(np.log(s0) + h) - f(np.log(s0) - h)) / (2 * h)), 3), fill=round(float(fill), 3), narrowing_number=round(float(r.measured_slope * min(fill, 1) * z1), 2))
    for k, v in summ["global"].items(): print("global", k, v)
    print(pd.DataFrame(summ["identifiability"]).to_string() if summ["identifiability"] else "")
    (FIG / "mechanism_summary.json").write_text(json.dumps(summ, indent=1))
    print(P.pivot(index="sample", columns="model", values="predicted_slope").join(H[["measured_slope", "z50", "H_nm", "cycles"]]).round(4).to_string())


if __name__ == "__main__":
    {"calibrate": calibrate, "predict": predict, "compare": compare}[sys.argv[1]]()
