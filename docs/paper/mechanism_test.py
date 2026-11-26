"""What sets the front slope of the published lateral-channel (PillarHall) Al2O3 profiles? Calibrate one sticking
coefficient on two full profiles, freeze it, and predict the front slopes of the other samples.

Transport kernels (single-exposure coverage theta(x/H; exposure, s) of an irreversible Langmuir film):
  hydraulic     1D diffusion with the hydraulic-diameter Knudsen diffusivity (the model behind the published slope method)
  ballistic     Monte Carlo of the real cavity (plates, pillars, side walls), free molecular            (figures/pillar_mc)
  gas           the same Monte Carlo with hard-sphere TMA-N2 collisions at the stated reactor pressure (300 Pa)
Each kernel is run bare, with the measurement spot, with channel narrowing, and with both (2 x 2).
  narrowing   not fitted. The film of every earlier cycle has narrowed the gap; the per-cycle exposure is fixed and the
              penetration scales with the gap (exponent n: 1 for diffusion, 1.01 measured for the gas-free Monte Carlo, 0.92
              with gas, figures/mechanism_gap_scaling.csv), so the N-cycle profile is the single-exposure profile averaged
              over a stretch f^n, f uniform on [1 - 2 t_film/H, 1]. Against the full multi-cycle hydraulic model the law is
              within 0.011 (0.021 at closure), figures/mechanism_narrowing_check.csv. The mechanism itself is in the
              Ylilammi 2018 model; what is added here is the closed form and the test.
  spot        top-hat average over the stated spot width (Yim 2020: 5-6, 25, 50 um; Haimi 2021 EDS: 1 um step; Arts 2019: 10 um)
The slope is the secant between 0.6 and 0.4 of the plateau; Yim 2020 and Arts 2019 fit a line through the front, which on
the Zenodo profile gives 0.019 against 0.0206 for the secant (8 %).

  calibrate   s per model from the Zenodo 3974438 mean profile and Ylilammi 2018 Fig. 6 (exposures free), with the s range
              whose rms is within 5 % of the best
  predict     slope of every other sample at its measured half-thickness depth, for all s on the kernel grid and at the
              calibrated s (figures/mechanism_predictions.csv, mechanism_slope_curves.json)
  compare     scores against the measured slopes (figures/mechanism_comparison.csv, mechanism_summary.json)
The test is out of sample but not blind and not pre-registered: the held-out table was known when the models were built.
  PYTHONPATH=src .venv/bin/python docs/paper/mechanism_test.py calibrate|predict|compare
"""
import sys, json, os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar, brentq
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
MC = ROOT / "figures/pillar_mc"; FIG = ROOT / "figures"
S_HYD = (0.0025, 0.004, 0.005, 0.0063, 0.008, 0.01, 0.0125, 0.0158, 0.02, 0.025, 0.03, 0.04, 0.06)
FAMILY = {"hydraulic": "hyd", "ballistic": "walls90_pillars49", "gas": "gas300"}
GAP_EXPONENT = {"hydraulic": 1.0, "ballistic": 1.0, "gas": 0.92}
MODELS = {f + ("+narrowing" if n else "") + ("+spot" if sp else ""): (f, n, sp) for f in FAMILY for n in (False, True) for sp in (False, True)}
COMPLETE = [f + "+narrowing+spot" for f in FAMILY]
# MECH_GAPS=measured uses the sacrificial-oxide thicknesses of Yim 2020 ESI Table S5 (wafers M5, M8/M9, M15: 108.6, 519.1, 2034.7 nm)
# instead of the nominal 100 / 500 / 2000 nm for the Yim 2020 chips and the Zenodo profile; outputs get the suffix _measgap
MEAS = os.environ.get("MECH_GAPS") == "measured"; SUF = "_measgap" if MEAS else ""; GAP_MEAS = {100: 108.6, 500: 519.1, 2000: 2034.7}
SPOT_UM = {50: 5.5, 10: 25.0, 5: 50.0, 0: 1.0, 100: 10.0}   # objective -> spot width; 0 = EDS line scan; 100 = Arts 2019 reflectometer


def at(E, TH, e):
    i = min(max(np.searchsorted(E, e), 1), len(E) - 1); w = (e - E[i - 1]) / (E[i] - E[i - 1]); return (1 - w) * TH[i - 1] + w * TH[i]


def mc_kernel(case):
    out = {}
    for f in sorted(MC.glob(f"{case}_s*.csv")):
        lines = [l for l in f.read_text().split("\n") if l and not l.startswith("#")]
        if len(lines) < 100: continue
        n = lines[0].count(","); d = np.array([np.fromstring(l, sep=",") for l in lines if l.count(",") == n and l[-1].isdigit()])
        if len(d) not in (150, 200, 300): continue                                       # an unfinished run is not used
        out[float(f.stem.split("_s")[-1])] = (d[:, 0], d[:, 1:], np.arange(d.shape[1] - 1) + 0.5)
    return dict(sorted(out.items()))


def hydraulic_kernel(deep=False):
    """Single-exposure coverage of the 1D hydraulic-diameter model (deep: to 1200 gaps and 90 Pa s, for deep fronts)."""
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


_KC = {}
def kernel(name):
    if name not in _KC: _KC[name] = hydraulic_kernel() if name == "hyd" else hydraulic_kernel(True) if name == "hyd_deep" else mc_kernel(name)
    return _KC[name]


def cumulative(run, e, fill, x, spot=0.0, n=1.0):
    """Normalised thickness after many cycles at reduced distance x (units of the initial gap): average over the gap
    stretch f^n (channel narrowing) and over a top-hat measurement spot of full width `spot`."""
    E, TH, z = run; th = at(E, TH, e)
    if spot > 0: return np.mean([cumulative(run, e, fill, np.abs(x + d), 0.0, n) for d in ((np.arange(9) + 0.5) / 9 - 0.5) * spot], 0)
    if fill <= 0: return np.interp(x, z, th, right=0.0)
    fs = np.linspace(max(1 - fill, 0.02), 1.0, 21) ** n
    return np.mean([np.interp(x / f, z, th, right=0.0) for f in fs], 0)


def front(run, e, fill, spot=0.0, n=1.0):
    xx = np.arange(0.5, run[2][-1], 0.5 if run[2][-1] < 500 else 1.0); y = cumulative(run, e, fill, xx, spot, n); y = y / y[:4].mean()
    def cross(l):
        i = int(np.argmax(y < l)); return np.nan if i == 0 else xx[i - 1] + (y[i - 1] - l) / (y[i - 1] - y[i]) * (xx[i] - xx[i - 1])
    c = [cross(l) for l in (0.6, 0.5, 0.4)]; return c[1], 0.2 / (c[2] - c[0])


def slope_at_depth(run, z50, fill, spot=0.0, n=1.0):
    """Front slope and exposure when the exposure is chosen so that the half-thickness depth equals the measured one."""
    E = run[0]; g = lambda le: front(run, 10 ** le, fill, spot, n)[0] - z50
    lo, hi = np.log10(E[3]), np.log10(E[-2])
    if not (np.isfinite(g(lo)) and np.isfinite(g(hi))) or g(lo) * g(hi) > 0: return np.nan, np.nan
    le = brentq(g, lo, hi, xtol=1e-4); return front(run, 10 ** le, fill, spot, n)[1], 10 ** le


def measured_profile(ds):
    import run_pillar_mc as r
    x, y = r.measured(ds); gap = GAP_MEAS[500] if (MEAS and ds == "zenodo") else 500.0        # the Zenodo chip is wafer M08; the Ylilammi 2018 chip is an earlier generation
    fill = 2 * (51.1 if ds == "zenodo" else float(pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv").query("x_um > 5 and x_um < 60").thickness_nm.mean())) / gap
    return x * 500.0 / gap, y, fill


def calibrate():
    rows = []; data = {ds: measured_profile(ds) for ds in ("zenodo", "ylilammi_al2o3")}
    for model, (fam, narrow, spot) in MODELS.items():
        K = kernel(FAMILY[fam]); tab = []
        for s, run in K.items():
            mse = {}
            for ds, (x, y, fill) in data.items():
                c = lambda le: float(np.mean((cumulative(run, 10 ** le, fill if narrow else 0, x, 11.0 if spot else 0.0, GAP_EXPONENT[fam]) - y) ** 2))
                q = minimize_scalar(c, bounds=(np.log10(run[0][3]), np.log10(run[0][-2])), method="bounded", options={"xatol": 1e-4}); mse[ds] = q.fun
            tab.append((s, mse["zenodo"], mse["ylilammi_al2o3"]))
        tab = np.array(tab); joint = tab[:, 1] + tab[:, 2]
        def best(col):
            j = int(np.argmin(col)); edge = j in (0, len(col) - 1); j = min(max(j, 1), len(col) - 2); p = np.polyfit(np.log(tab[j - 1:j + 2, 0]), col[j - 1:j + 2], 2)
            ls = np.clip(-p[1] / (2 * p[0]), np.log(tab[j - 1, 0]), np.log(tab[j + 1, 0])) if p[0] > 0 else np.log(tab[j, 0])
            return float(np.exp(ls)), float(np.sqrt(max(np.polyval(p, ls), 0))) if p[0] > 0 else float(np.sqrt(col[j])), bool(edge)
        jb, zb, yb = best(joint), best(tab[:, 1]), best(tab[:, 2])
        # s range whose joint rms is within 5 % of the best (interpolated in ln s); open if it reaches the end of the grid
        lsf = np.linspace(np.log(tab[0, 0]), np.log(tab[-1, 0]), 400); rf = np.sqrt(np.interp(lsf, np.log(tab[:, 0]), joint) / 2); ok = rf <= 1.05 * rf.min()
        rows.append(dict(model=model, kernel=FAMILY[fam], narrowing=narrow, spot=spot, s_joint=jb[0], rms_joint=float(np.sqrt(jb[1] ** 2 / 2)), s_lo=float(np.exp(lsf[ok][0])), s_hi=float(np.exp(lsf[ok][-1])),
                         open_below=bool(ok[0]), open_above=bool(ok[-1]), s_zenodo=zb[0], rms_zenodo=zb[1], s_ylilammi=yb[0], rms_ylilammi=yb[1], n_s=len(tab), s_grid=" ".join(f"{v:g}" for v in tab[:, 0]),
                         rms_zenodo_grid=" ".join(f"{np.sqrt(v):.4f}" for v in tab[:, 1]), rms_ylilammi_grid=" ".join(f"{np.sqrt(v):.4f}" for v in tab[:, 2])))
        print({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in rows[-1].items() if "grid" not in k and k not in ("kernel", "narrowing", "spot", "n_s")}, flush=True)
    pd.DataFrame(rows).to_csv(FIG / f"mechanism_calibration{SUF}.csv", index=False)


def heldout():
    """Samples not used in the calibration. Yim 2020 Table 1 without sample 8 (the Zenodo calibration profile) and sample
    14 (0.4 s pulse on layout v2a: the run Yim 2020 use for the water step, front at 386 gaps); Haimi 2021 specimen A
    (EDS: 1080 cycles, 500 nm gap, about 300 Pa; half-thickness depth 195.2, slope 0.0108; film 108 nm from the 0.1 nm
    growth per cycle of the series); Arts 2019 Fig. 2 TMA-limited profile (digitised, figures/lhar_slope_table.csv; 400
    cycles, about 46 nm; another reactor and 220 C, so its s need not be the same: reported separately, not scored)."""
    t = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv"); rows = []
    for r in t.itertuples():
        if r.sample in (8, 14): continue
        g = GAP_MEAS[r.H_nm] / r.H_nm if MEAS else 1.0      # the table's reduced depth and slope use the nominal gap
        rows.append(dict(sample=f"Yim{r.sample}", series=r.series, H_nm=r.H_nm, layout=r.pillar_layout, cycles=r.cycles, magnification=r.magnification, film_nm=r.S_IIb_nm / g, z50=r.x50_tilde / g, measured_slope=r.slope_dtheta_dxtilde * g, gs=g))
    rows.append(dict(sample="HaimiA", series="EDS", H_nm=500, layout="v1b", cycles=1080, magnification=0, film_nm=108.0, z50=195.2, measured_slope=0.0108, gs=1.0))
    rows.append(dict(sample="Arts", series="Arts", H_nm=500, layout="arts", cycles=400, magnification=100, film_nm=46.0, z50=486.1, measured_slope=0.00768, gs=1.0))
    return pd.DataFrame(rows)


def kernel_for(fam, r, fill):
    deep = r.z50 / (1 - min(fill, 1.0) / 2) > 300          # the uncoated-gap front the exposure ladder must reach
    if fam == "hydraulic": return "hyd_deep" if (deep or r.layout == "arts") else "hyd"
    if r.layout == "arts": return {"ballistic": "arts_free", "gas": "arts_gas27"}[fam]
    if r.H_nm == 2000: return {"ballistic": "h2000_free", "gas": "h2000_gas300"}[fam]
    if r.H_nm == 100: return {"ballistic": "h100_free", "gas": "h100_gas300"}[fam]
    if fam == "gas": return {"v1a": "v1a_gas300", "v2a": "v2a_gas300"}.get(r.layout, "gas300_deep" if deep else "gas300")
    return {"v1a": "pillars28", "v2a": "walls90"}.get(r.layout, "free_deep" if deep else "walls90_pillars49")   # gas-free: v1a without side walls, v2a without its 98 um pillars


def predict():
    cal = pd.read_csv(FIG / f"mechanism_calibration{SUF}.csv").set_index("model"); H = heldout(); rows = []; curves = {}
    for model, (fam, narrow, spot) in MODELS.items():
        s_star = float(cal.loc[model, "s_joint"])
        for r in H.itertuples():
            fill_true = 2 * r.film_nm / r.H_nm; kn = kernel_for(fam, r, fill_true if narrow else 0.0); K = kernel(kn)
            fill = fill_true if narrow else 0.0; sp = SPOT_UM[r.magnification] * 1000.0 / (r.H_nm * r.gs) if spot else 0.0
            ss = np.array(list(K)); sl = np.array([slope_at_depth(K[s], r.z50, fill, sp, GAP_EXPONENT[fam]) for s in ss]) if len(ss) else np.zeros((0, 2)); ok = np.isfinite(sl[:, 0]) if len(ss) else np.array([], bool)
            if ok.sum() < 2: rows.append(dict(model=model, sample=r.sample, kernel=kn, s=s_star, predicted_slope=np.nan, exposure_Pa_s=np.nan, extrapolated=True)); continue
            curves.setdefault(model, {})[r.sample] = dict(kernel=kn, s=[float(v) for v in ss[ok]], slope=[round(float(v), 5) for v in sl[ok, 0]], exposure=[round(float(v), 3) for v in sl[ok, 1]])
            ex = bool(s_star < ss[ok].min() or s_star > ss[ok].max())
            rows.append(dict(model=model, sample=r.sample, kernel=kn, s=s_star, predicted_slope=float(np.exp(np.interp(np.log(s_star), np.log(ss[ok]), np.log(sl[ok, 0])))),
                             exposure_Pa_s=float(np.interp(np.log(s_star), np.log(ss[ok]), sl[ok, 1])), extrapolated=ex))
        print(model, "done", flush=True)
    pd.DataFrame(rows).to_csv(FIG / f"mechanism_predictions{SUF}.csv", index=False); (FIG / f"mechanism_slope_curves{SUF}.json").write_text(json.dumps(curves))


def compare():
    P = pd.read_csv(FIG / f"mechanism_predictions{SUF}.csv"); H = heldout().set_index("sample"); P = P.join(H, on="sample")
    P.loc[P.extrapolated.astype(bool), "predicted_slope"] = np.nan                         # a prediction outside the simulated s range is not scored
    P["log_ratio"] = np.log(P.predicted_slope / P.measured_slope); P.to_csv(FIG / f"mechanism_comparison{SUF}.csv", index=False)
    scored = P[P["sample"] != "Arts"]; piv = scored.pivot(index="sample", columns="model", values="log_ratio")
    common = list(piv.dropna().index); summ = {"n_heldout": int(len(H) - 1), "common_samples": common, "models": {}}
    # scatter of nominally identical runs (50x objective, layout v1b, 500 nm, 500 cycles, 0.1 s pulse): the slope falls with depth, so the scatter is taken on slope x depth
    rep = H.loc[["Yim2", "Yim4", "Yim12", "Yim15", "Yim16"]]; summ["replicate_sigma_ln_slope_times_depth"] = round(float(np.std(np.log(rep.measured_slope * rep.z50), ddof=1)), 3)
    summ["replicate_sigma_ln_slope"] = round(float(np.std(np.log(rep.measured_slope), ddof=1)), 3)
    rms = lambda v: round(float(np.sqrt(np.mean(np.square(v)))), 3)
    for m in MODELS:
        g = scored[(scored.model == m) & scored.log_ratio.notna()].set_index("sample"); c = g.loc[common]; c9 = c.drop(index="Yim9", errors="ignore")
        summ["models"][m] = dict(s=float(f"{g.s.iloc[0]:.3g}"), n=int(len(g)), rms_ln=rms(g.log_ratio), mean_ln=round(float(g.log_ratio.mean()), 3), n_common=len(c), rms_ln_common=rms(c.log_ratio),
                                 rms_ln_common_without_Yim9=rms(c9.log_ratio), worst=str(g.log_ratio.abs().idxmax()), missing=sorted(set(H.index) - set(g.index) - {"Arts"}))
        print(m, {k: v for k, v in summ["models"][m].items()})
    # ceiling: the slope of the complete model for s -> large (largest simulated s), how close the measurement is to it, and the sensitivity d ln slope / d ln s at three s
    curves = json.loads((FIG / f"mechanism_slope_curves{SUF}.json").read_text()); summ["ceiling"] = {}; summ["single_s_scan"] = {}; sig = summ["replicate_sigma_ln_slope_times_depth"]
    for m in COMPLETE:
        per = curves.get(m, {}); summ["ceiling"][m] = {}
        for smp, c in per.items():
            ls, lsl = np.log(c["s"]), np.log(c["slope"]); r = H.loc[smp]; fill = min(2 * r.film_nm / r.H_nm, 1.0)
            def sens(s0):
                if not (c["s"][0] * 1.3 <= s0 <= c["s"][-1] / 1.3): return None
                return round(float((np.interp(np.log(s0 * 1.3), ls, lsl) - np.interp(np.log(s0 / 1.3), ls, lsl)) / (2 * np.log(1.3))), 3)
            inv = float(np.exp(np.interp(np.log(r.measured_slope), lsl, ls))) if lsl[0] < np.log(r.measured_slope) < lsl.max() and np.all(np.diff(lsl) > 0) else None
            summ["ceiling"][m][smp] = dict(fill=round(float(fill), 3), s_max=c["s"][-1], slope_at_s_max=c["slope"][-1], measured_over_ceiling=round(float(r.measured_slope / c["slope"][-1]), 2),
                                           sensitivity_at_0p004=sens(0.004), sensitivity_at_0p01=sens(0.01), sensitivity_at_0p02=sens(0.02), s_matching_measured_slope=None if inv is None else float(f"{inv:.3g}"))
        # one s for all scored slopes, on the s range every sample's curve covers (no extrapolation), with and without the 2000 nm sample
        for tag, drop in (("all", ()), ("without_Yim9", ("Yim9",))):
            use = [k for k in per if k != "Arts" and k not in drop and len(per[k]["s"]) >= 3]
            if not use: continue
            lo = max(per[k]["s"][0] for k in use); hi = min(per[k]["s"][-1] for k in use); sg = np.geomspace(lo, hi, 60); chi = np.zeros_like(sg)
            for k in use: chi += (np.interp(np.log(sg), np.log(per[k]["s"]), np.log(per[k]["slope"])) - np.log(H.loc[k, "measured_slope"])) ** 2
            j = int(np.argmin(chi)); inside = sg[chi <= chi[j] + sig ** 2]
            summ["single_s_scan"].setdefault(m, {})[tag] = dict(n=len(use), s_range_scanned=[float(f"{lo:.3g}"), float(f"{hi:.3g}")], s_best=float(f"{sg[j]:.3g}"), rms_ln=round(float(np.sqrt(chi[j] / len(use))), 3),
                                                              s_lo=float(f"{inside.min():.3g}"), s_hi=float(f"{inside.max():.3g}"), open_below=bool(inside.min() == sg[0]), open_above=bool(inside.max() == sg[-1]))
    # the cycle series: with a fixed recipe the per-cycle exposure implied by each front should not depend on the cycle number
    summ["implied_exposure_cycle_series_Pa_s"] = {}
    for m in ("gas", "gas+narrowing+spot", "hydraulic", "hydraulic+narrowing+spot"):
        g = P[(P.model == m) & (P.series == "D")].sort_values("cycles"); g8 = P[(P.model == m) & P["sample"].isin(["Yim4", "Yim12"])]
        summ["implied_exposure_cycle_series_Pa_s"][m] = {"250": round(float(g[g.cycles == 250].exposure_Pa_s.iloc[0]), 2), "500 (mean of samples 4, 12)": round(float(g8.exposure_Pa_s.mean()), 2), "1000": round(float(g[g.cycles == 1000].exposure_Pa_s.iloc[0]), 2)}
    for k in ("single_s_scan", "implied_exposure_cycle_series_Pa_s"): print(k, json.dumps(summ[k]))
    print(pd.DataFrame(summ["ceiling"]["gas+narrowing+spot"]).T.to_string())
    (FIG / f"mechanism_summary{SUF}.json").write_text(json.dumps(summ, indent=1))
    print(P.pivot(index="sample", columns="model", values="predicted_slope")[["hydraulic", "ballistic", "gas"] + COMPLETE].join(H[["measured_slope", "z50", "H_nm", "cycles"]]).round(4).to_string())


if __name__ == "__main__":
    {"calibrate": calibrate, "predict": predict, "compare": compare}[sys.argv[1]]()
