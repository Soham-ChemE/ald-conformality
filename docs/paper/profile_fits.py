"""Full-profile test of channel narrowing (post hoc; made after the blind tests and after an independent review).
Every digitised front is fitted on its own, with the sticking coefficient s and the exposure free, with and without channel
narrowing. If narrowing is the right physics the fitted s of one process should not depend on the number of cycles or on the
filled fraction; without it the fitted s should fall as the channel fills.
  TMA/water (three transport kernels): reflectometry curves of the second blind test (data/blind_revealed_v2/, Yim 2020) and
      the EDS line scans (docs/blindtest/v2/posthoc/eds_cycle_series_curves.json, two channels per scan).
  ZnO from diethylzinc and water (hydraulic kernel only; no Monte Carlo was run for this gas pair): the profiles of the
      first blind test with a flat plateau (fall of less than 0.05 before the front), 7 to 11 % filled.
The filled fraction uses the thicker of the two film thicknesses of Yim 2020 Table 1 (S_IIb); FIT_FILM=S_I repeats the
TMA/water fits with the thinner one (files *_filmSI.csv). Fronts deeper than 250 gap heights (uncoated-gap depth) use the
1200-gap kernels of the standard layout. The plateau is fixed at the published normalisation, not fitted.
Fit window: from 5 gap heights to the first reading below 0.1 of the plateau that lies beyond 20 gap heights (the whole
curve if there is none); the exposure search covers front depths of 0.8 to 1.25 times the measured half depth. Several reflectometry curves never fall below
0.1 (a residual of up to 0.13 of the plateau stays behind the front), so the fits are repeated with the window closed at
0.3 of the plateau (FIT_CUT=0.3; files *_cut0.3.csv), which keeps the upper part of the front only.
For every simulated s the model is tabulated on a fine exposure ladder; the tables are interpolated in ln s to six steps
between simulated values (`interp_rms_change` is the largest change of the best rms when a simulated s is replaced by the
interpolation between its neighbours), and s is the value of smallest rms on that ladder. Uncertainty: moving-block bootstrap of
the centred residuals of the best fit (block length = the lag at which their autocorrelation falls below 1/e, between 5
and 40 readings; 300 replicates, refitted on the same ladder), which keeps the correlation of the residuals; the 16th to 84th percentile is reported, with flags when that interval reaches an end of the simulated range.
The kernels include runs added for this test and not used by the earlier analyses: Monte Carlo at s = 0.0016, 0.0025, 0.1 and 0.2
(files *hi_s*.csv) and the hydraulic model at s = 0.001, 0.0016, 0.1 and 0.2.
Writes figures/profile_fits.csv, figures/profile_fits_summary.csv and figures/profile_fits_zno.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/profile_fits.py"""
import sys, json, os
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
CUT = float(os.environ.get("FIT_CUT", 0.1)); FILM = os.environ.get("FIT_FILM", "S_IIb"); SUF = ("" if CUT == 0.1 else f"_cut{CUT:g}") + ("" if FILM == "S_IIb" else "_filmSI")     # FIT_CUT=0.3 repeats the fits on the upper part of the front only
NB, BLOCK, NE = 300, 5, 160; rng = np.random.default_rng(20261127)
HI = {"gas300": "gas300hi", "gas300_deep": "gas300deephi", "walls90_pillars49": "freehi", "free_deep": "freedeephi"}
S_HYD_EXT = (0.001, 0.0016, 0.1, 0.2)
def hyd_ext(deep):
    """The hydraulic kernel of mechanism_test at four more sticking coefficients (same solver, grid and exposure ladder)."""
    cache = mt.FIG / ("profile_fits_hydraulic_ext_deep.npz" if deep else "profile_fits_hydraulic_ext.npz")
    if cache.exists():
        z = np.load(cache); return {float(s): (z["E"], z[f"th_{s}"], z["z"]) for s in z["s"]}
    from aldconf import Slit, Precursor, Process, simulate
    E = (90.0 if deep else 16.0) * (np.arange(1, 121) / 120.0) ** 2; zz = np.arange(1200 if deep else 400) + 0.5; out = {}
    for s in S_HYD_EXT:
        TH = []
        for e in E:
            pr = Precursor(sticking=s, site_density=4e18, mass_amu=72.09, gpc=1e-10)
            r = simulate(Process(Slit(500e-9, 600e-6 if deep else 200e-6, knudsen_correction="none"), pr, T=573.15, p_partial=e / 0.1, t_dose=0.1, cycles=1, n_cells=1200 if deep else 400, n_substeps=60, knudsen_factor=1.0, narrowing=False))
            TH.append(np.interp(zz, r.x / 500e-9, r.thickness / 1e-10))
        out[s] = (E, np.array(TH), zz); print("hydraulic kernel", "deep" if deep else "", s, flush=True)
    np.savez_compressed(cache, E=E, z=zz, s=np.array(list(out)), **{f"th_{s}": v[1] for s, v in out.items()})
    return out
def kern(name):
    K = dict(mt.kernel(name))
    if name in HI: K.update(mt.mc_kernel(HI[name]))
    elif name not in ("hyd", "hyd_deep"): K.update(mt.mc_kernel(name + "hi"))      # the other layouts: files <case>hi_s*.csv
    if name in ("hyd", "hyd_deep"): K.update(hyd_ext(name == "hyd_deep"))
    return dict(sorted(K.items()))
ND = 81
def tables(z, y, fill, spot, K, n, narrow):
    """Model on the data grid for every simulated s, tabulated at ND common half-depths around the measured one (0.8 to 1.25
    times it). For each s the model is computed on a fine exposure ladder and interpolated between neighbouring exposures to
    each target depth, so that entry k of every s has the same front position and only the shape differs."""
    k50 = int(np.argmax(y < 0.5)); z50 = z[k50 - 1] + (y[k50 - 1] - 0.5) / (y[k50 - 1] - y[k50]) * (z[k50] - z[k50 - 1]) if k50 > 0 else z[0]; DT = np.linspace(0.8, 1.25, ND) * z50; S = sorted(K); out = []; zz = np.concatenate([[1.0, 2.0, 3.0, 4.0], z]); xd = np.arange(2.0, 1.6 * z50, 2.0)
    for s in S:
        run = K[s]; E = run[0]; f = fill if narrow else 0.0
        d = np.array([mt.front(run, e, f, spot, n)[0] for e in E[1::4]]); Ec = E[1::4]; ok = np.isfinite(d)
        if ok.sum() < 2 or np.nanmax(d) < 0.7 * z50: out.append(np.full((ND, len(z)), np.nan)); continue
        lo = Ec[ok][np.searchsorted(d[ok], 0.7 * z50).clip(0, ok.sum() - 1)]; hi = Ec[ok][np.searchsorted(d[ok], 1.4 * z50).clip(0, ok.sum() - 1)]
        lo, hi = min(lo, hi) * 0.8, min(max(lo, hi) * 1.25, E[-1]); G = np.geomspace(max(lo, E[1]), hi, NE)
        M = np.array([mt.cumulative(run, e, f, np.concatenate([zz, xd]), spot, n) for e in G]); M = M / M[:, :4].mean(1)[:, None]; P = M[:, 4:4 + len(z)]; Q = M[:, 4 + len(z):]
        d50 = np.array([xd[max(int(np.argmax(q < 0.5)), 1) - 1] + (q[max(int(np.argmax(q < 0.5)), 1) - 1] - 0.5) / max(q[max(int(np.argmax(q < 0.5)), 1) - 1] - q[max(int(np.argmax(q < 0.5)), 1)], 1e-9) * 2.0 if (q < 0.5).any() else np.nan for q in Q])
        T = np.full((ND, len(z)), np.nan); g = np.isfinite(d50) & (np.concatenate([[True], np.diff(d50) > 0]) if len(d50) > 1 else True); dd, PP = d50[g], P[g]
        for k, t in enumerate(DT):
            if len(dd) < 2 or t < dd[0] or t > dd[-1]: continue
            a = int(np.searchsorted(dd, t)); a = min(max(a, 1), len(dd) - 1); w = (t - dd[a - 1]) / (dd[a] - dd[a - 1]); T[k] = (1 - w) * PP[a - 1] + w * PP[a]
        out.append(T)
    return np.array(S), np.array(out)                                   # (nS,), (nS, ND, nz); NaN where a depth is not reached
NSUB = 6
def refine(S, T):
    """Tables on a finer ladder of s: linear interpolation in ln s between neighbouring simulated values, at the same index of
    the exposure ladder (each ladder spans the same range of front depths around the measured one)."""
    Sf, Tf = [], []
    for i in range(len(S) - 1):
        for w in np.arange(NSUB) / NSUB: Sf.append(float(np.exp((1 - w) * np.log(S[i]) + w * np.log(S[i + 1])))); Tf.append((1 - w) * T[i] + w * T[i + 1])
    Sf.append(float(S[-1])); Tf.append(T[-1]); return np.array(Sf), np.array(Tf)
def pick(S, T, y):
    r = np.sqrt(((T - y) ** 2).mean(2)); r = np.where(np.isfinite(r), r, np.inf); per = r.min(1); j = int(np.argmin(per)); return float(S[j]), float(per[j]), j in (0, len(S) - 1), T[j, int(np.argmin(r[j]))]
def interp_check(S, T, y):
    """Largest change of the best rms when a simulated s is replaced by the interpolation between its two neighbours."""
    worst = 0.0
    for i in range(1, len(S) - 1):
        w = (np.log(S[i]) - np.log(S[i - 1])) / (np.log(S[i + 1]) - np.log(S[i - 1])); Ti = (1 - w) * T[i - 1] + w * T[i + 1]
        ra = np.sqrt(((Ti - y) ** 2).mean(1)); rb = np.sqrt(((T[i] - y) ** 2).mean(1)); ra, rb = ra[np.isfinite(ra)], rb[np.isfinite(rb)]
        if len(ra) and len(rb): worst = max(worst, abs(ra.min() - rb.min()))
    return worst
def fit(z, y, fill, spot, K, n, narrow):
    o = np.argsort(z); z, y = z[o], y[o]; m = (y < CUT) & (z > 20); iend = int(np.argmax(m)) if m.any() else len(y) - 1; w = (z >= 5) & (np.arange(len(z)) <= iend); zw, yw = z[w], y[w]
    S0, T0 = tables(zw, yw, fill, spot, K, n, narrow); reach = np.isfinite(T0).all(2).mean(1); s_eff = float(S0[reach >= 0.9].max()) if (reach >= 0.9).any() else np.nan
    chk = interp_check(S0, T0, yw); S, T = refine(S0, T0); s, rms, edge, best = pick(S, T, yw); res = yw - best; res = res - res.mean(); nz = len(yw)
    ac = np.array([np.corrcoef(res[:-k], res[k:])[0, 1] if nz - k > 3 else 0.0 for k in range(1, min(41, nz // 3))]); below = np.where(ac < np.exp(-1))[0]; blk = int(np.clip((below[0] + 1) if len(below) else 40, BLOCK, max(BLOCK, nz // 3))); nb = int(np.ceil(nz / blk)); B = []
    for _ in range(NB):
        st = rng.integers(0, nz - blk + 1, nb); rs = np.concatenate([res[i:i + blk] for i in st])[:nz]; B.append(pick(S, T, best + rs)[0])
    B = np.array(B); lo, hi = np.percentile(B, [16, 84]); top, bot = S.max(), S.min()
    return dict(n_points=int(nz), s_fit=float(f"{s:.3g}"), s_p16=float(f"{lo:.3g}"), s_p84=float(f"{hi:.3g}"), s_p2_5=float(f"{np.percentile(B, 2.5):.3g}"), s_p97_5=float(f"{np.percentile(B, 97.5):.3g}"), rms=round(rms, 4), at_grid_edge=bool(edge),
                no_upper_limit=bool(hi >= 0.999 * top), no_lower_limit=bool(lo <= 1.001 * bot), boot_frac_at_edge=round(float(np.mean((B >= 0.999 * top) | (B <= 1.001 * bot))), 2), s_grid_min=float(bot), s_grid_max=float(top), s_max_reaching_all_depths=s_eff, n_s_simulated=len(S0), block_length=blk, interp_rms_change=round(chk, 4))
if __name__ == "__main__":
    H = mt.heldout().set_index("sample"); TB = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv").set_index("sample"); profiles = []
    def kname(fam, r, fill):                 # deep fronts use the 1200-gap kernels of the standard layout, whatever the pillar layout
        if r.H_nm == 500 and fam == "hydraulic": return "hyd_deep"          # the 1200-gap hydraulic kernel reaches every depth the narrowed model needs
        if r.H_nm == 500 and r.z50 / (1 - fill / 2) > 250: return {"hydraulic": "hyd_deep", "gas": "gas300_deep", "ballistic": "free_deep"}[fam]
        return mt.kernel_for(fam, r, fill)
    for f in sorted((ROOT / "data/blind_revealed_v2").glob("yim2020_s*.csv")):
        k = int(f.stem.split("_s")[1]); d = pd.read_csv(f); name = f"Yim{k}"
        if name in H.index: profiles.append(dict(id=f"opt_s{k:02d}", sample=name, method="reflectometry", z=d.x_over_H.values, y=d.thickness_over_plateau.values, spot_um=mt.SPOT_UM[int(H.loc[name, "magnification"])]))
    for c in json.loads((ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").read_text()):
        profiles.append(dict(id=f"eds_{c['cycles']}_{c['channel']}", sample={250: "Yim10", 500: "Yim2", 1000: "Yim11"}[c["cycles"]], method="EDS", z=np.array(c["z"]), y=np.array(c["y"]), spot_um=1.0))
    rows = []
    for P in profiles:
        r = H.loc[P["sample"]]; film = float(TB.loc[int(P["sample"][3:]), "S_I_nm"]) if FILM == "S_I" else r.film_nm; fill = min(2 * film / r.H_nm, 1.0); spot = P["spot_um"] * 1000 / r.H_nm
        for fam in mt.FAMILY:
            K = kern(kname(fam, r, fill))
            for narrow in (True, False): rows.append(dict(id=P["id"], sample=P["sample"], method=P["method"], cycles=int(r.cycles), H_nm=int(r.H_nm), fill=round(fill, 3), kernel=kname(fam, r, fill), transport=fam, narrowing=narrow, **fit(P["z"], P["y"], fill, spot, K, mt.GAP_EXPONENT[fam], narrow)))
        print(P["id"], "done", flush=True)
    D = pd.DataFrame(rows); D.to_csv(ROOT / f"figures/profile_fits{SUF}.csv", index=False); out = []
    for (fam, narrow), g in D.groupby(["transport", "narrowing"]):
        for name, q in (("all profiles", g), ("EDS", g[g.method == "EDS"]), ("reflectometry, 500 nm gap", g[(g.H_nm == 500) & (g.method == "reflectometry")])):
            ls = np.log(q.s_fit); c = np.polyfit(q.fill, ls, 1, cov=True) if q.fill.nunique() > 2 else None
            out.append(dict(transport=fam, narrowing=narrow, profiles=name, n=len(q), n_at_grid_edge=int(q.at_grid_edge.sum()), s_geomean=float(f"{np.exp(ls.mean()):.3g}"), s_min=q.s_fit.min(), s_max=q.s_fit.max(), sd_ln_s=round(float(ls.std(ddof=1)), 3),
                            dlns_dfill=round(float(c[0][0]), 2) if c else np.nan, dlns_dfill_se=round(float(np.sqrt(c[1][0, 0])), 2) if c else np.nan, mean_rms=round(float(q.rms.mean()), 4)))
    pd.DataFrame(out).to_csv(ROOT / f"figures/profile_fits_summary{SUF}.csv", index=False)
    if FILM != "S_IIb": sys.exit(0)                                   # the thin-film repeat covers the TMA/water fronts only
    # ---- ZnO: flat-plateau profiles of the first blind test ----
    Dg = pd.read_csv(ROOT / "docs/blindtest/posthoc/diagnosis.csv"); Pr = pd.read_csv(ROOT / "docs/blindtest/predictions.csv").set_index("id"); zr = []
    for q in Dg[Dg.id.str.contains("ZnO_DEZ_H2O|haimi2026_fig5") & (Dg.droop.abs() < 0.05)].itertuples():
        d = pd.read_csv(ROOT / "data/blind_revealed" / f"{q.id}.csv").dropna(subset=["x_over_H", "thickness_over_plateau"]); d = d[d.x_over_H >= 0]; spot = Pr.loc[q.id, "spot_units"] if "spot_units" in Pr.columns and np.isfinite(Pr.loc[q.id, "spot_units"]) else 20.0
        K = kern("hyd" if q.z50 / (1 - q.fill / 2) < 250 else "hyd_deep")
        for narrow in (True, False): zr.append(dict(id=q.id, fill=q.fill, z50=q.z50, transport="hydraulic", narrowing=narrow, **fit(d.x_over_H.values, d.thickness_over_plateau.values, q.fill, spot, K, 1.0, narrow)))
        print(q.id, "done", flush=True)
    Z = pd.DataFrame(zr); Z.to_csv(ROOT / f"figures/profile_fits_zno{SUF}.csv", index=False); pd.set_option("display.width", 250)
    print(pd.DataFrame(out).to_string(index=False)); print("largest interpolation check", D.interp_rms_change.max(), Z.interp_rms_change.max())
