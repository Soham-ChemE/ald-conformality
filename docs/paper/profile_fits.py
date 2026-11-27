"""Full-profile test of channel narrowing (post hoc; made after the blind tests and after an independent review).
Every digitised front is fitted on its own, with the sticking coefficient s and the exposure free, with and without channel
narrowing. If narrowing is the right physics the fitted s of one process should not depend on the number of cycles or on the
filled fraction; without it the fitted s should fall as the channel fills.
  TMA/water (three transport kernels): reflectometry curves of the second blind test (data/blind_revealed_v2/, Yim 2020) and
      the EDS line scans (docs/blindtest/v2/posthoc/eds_cycle_series_curves.json, two channels per scan).
  ZnO from diethylzinc and water (hydraulic kernel only; no Monte Carlo was run for this gas pair): the profiles of the
      first blind test with a flat plateau (fall of less than 0.05 before the front), 7 to 11 % filled.
Fit window: from 5 gap heights to the first reading below 0.1 of the plateau. Several reflectometry curves never fall below
0.1 (a residual of up to 0.13 of the plateau stays behind the front), so the fits are repeated with the window closed at
0.3 of the plateau (FIT_CUT=0.3; files *_cut0.3.csv), which keeps the upper part of the front only.
For every s of the kernel the model is tabulated on a fine exposure grid; s is taken at the minimum of a parabola through the
per-s minimum rms in ln s. Uncertainty: moving-block bootstrap of the residuals of the best fit (blocks of 5 readings, 300
replicates, refitted on the same tables), which keeps the correlation of the residuals; the 16th to 84th percentile of the
refitted s is reported, with the fraction of replicates that end at either end of the simulated range.
The kernels include runs added for this test and not used by the earlier analyses: Monte Carlo at s = 0.0016, 0.0025, 0.1 and 0.2
(files *hi_s*.csv) and the hydraulic model at s = 0.001, 0.0016, 0.1 and 0.2.
Writes figures/profile_fits.csv, figures/profile_fits_summary.csv and figures/profile_fits_zno.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/profile_fits.py"""
import sys, json, os
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
CUT = float(os.environ.get("FIT_CUT", 0.1)); SUF = "" if CUT == 0.1 else f"_cut{CUT:g}"     # FIT_CUT=0.3 repeats the fits on the upper part of the front only
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
    if name in ("hyd", "hyd_deep"): K.update(hyd_ext(name == "hyd_deep"))
    return dict(sorted(K.items()))
def tables(z, y, fill, spot, K, n, narrow):
    """Model on the data grid for every s and a fine exposure ladder around the measured half depth."""
    z50 = z[int(np.argmax(y < 0.5))]; S = sorted(K); out = []
    for s in S:
        run = K[s]; E = run[0]; f = fill if narrow else 0.0
        d = np.array([mt.front(run, e, f, spot, n)[0] for e in E[1::4]]); Ec = E[1::4]; ok = np.isfinite(d)
        lo = Ec[ok][np.searchsorted(d[ok], 0.7 * z50).clip(0, ok.sum() - 1)] if ok.any() else E[1]; hi = Ec[ok][np.searchsorted(d[ok], 1.35 * z50).clip(0, ok.sum() - 1)] if ok.any() else E[-1]
        lo, hi = min(lo, hi) * 0.8, min(max(lo, hi) * 1.25, E[-1]); G = np.geomspace(max(lo, E[1]), hi, NE)
        M = np.array([mt.cumulative(run, e, f, np.concatenate([[1.0, 2.0, 3.0, 4.0], z]), spot, n) for e in G]); out.append(M[:, 4:] / M[:, :4].mean(1)[:, None])
    return np.array(S), np.array(out)                                   # (nS,), (nS, NE, nz)
def pick(S, T, y):
    r = np.sqrt(((T - y) ** 2).mean(2)); per = r.min(1); j = int(np.argmin(per)); s = S[j]; rb = per[j]; edge = j in (0, len(S) - 1)
    if not edge:
        a = np.polyfit(np.log(S[j - 1:j + 2]), per[j - 1:j + 2], 2)
        if a[0] > 0: s = float(np.exp(-a[1] / (2 * a[0])))
    return s, float(rb), edge, T[j, int(np.argmin(r[j]))]
def fit(z, y, fill, spot, K, n, narrow):
    o = np.argsort(z); z, y = z[o], y[o]; m = (y < CUT) & (z > 20); iend = int(np.argmax(m)) if m.any() else len(y) - 1; w = (z >= 5) & (np.arange(len(z)) <= iend); zw, yw = z[w], y[w]
    S, T = tables(zw, yw, fill, spot, K, n, narrow); s, rms, edge, best = pick(S, T, yw); res = yw - best; nz = len(yw); nb = int(np.ceil(nz / BLOCK)); B, E = [], 0
    for _ in range(NB):
        st = rng.integers(0, nz - BLOCK + 1, nb); rs = np.concatenate([res[i:i + BLOCK] for i in st])[:nz]; sb, _, eb, _ = pick(S, T, best + rs); B.append(sb); E += eb
    B = np.array(B); return dict(n_points=int(nz), s_fit=float(f"{s:.3g}"), s_p16=float(f"{np.percentile(B, 16):.3g}"), s_p84=float(f"{np.percentile(B, 84):.3g}"), s_p2_5=float(f"{np.percentile(B, 2.5):.3g}"), s_p97_5=float(f"{np.percentile(B, 97.5):.3g}"),
                         rms=round(rms, 4), at_grid_edge=bool(edge), boot_frac_at_edge=round(E / NB, 2), s_grid_min=float(S.min()), s_grid_max=float(S.max()))
if __name__ == "__main__":
    H = mt.heldout().set_index("sample"); profiles = []
    for f in sorted((ROOT / "data/blind_revealed_v2").glob("yim2020_s*.csv")):
        k = int(f.stem.split("_s")[1]); d = pd.read_csv(f); name = f"Yim{k}"
        if name in H.index: profiles.append(dict(id=f"opt_s{k:02d}", sample=name, method="reflectometry", z=d.x_over_H.values, y=d.thickness_over_plateau.values, spot_um=mt.SPOT_UM[int(H.loc[name, "magnification"])]))
    for c in json.loads((ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").read_text()):
        profiles.append(dict(id=f"eds_{c['cycles']}_{c['channel']}", sample={250: "Yim10", 500: "Yim2", 1000: "Yim11"}[c["cycles"]], method="EDS", z=np.array(c["z"]), y=np.array(c["y"]), spot_um=1.0))
    rows = []
    for P in profiles:
        r = H.loc[P["sample"]]; fill = min(2 * r.film_nm / r.H_nm, 1.0); spot = P["spot_um"] * 1000 / r.H_nm
        for fam in mt.FAMILY:
            K = kern(mt.kernel_for(fam, r, fill))
            for narrow in (True, False): rows.append(dict(id=P["id"], sample=P["sample"], method=P["method"], cycles=int(r.cycles), H_nm=int(r.H_nm), fill=round(fill, 3), transport=fam, narrowing=narrow, **fit(P["z"], P["y"], fill, spot, K, mt.GAP_EXPONENT[fam], narrow)))
        print(P["id"], "done", flush=True)
    D = pd.DataFrame(rows); D.to_csv(ROOT / f"figures/profile_fits{SUF}.csv", index=False); out = []
    for (fam, narrow), g in D.groupby(["transport", "narrowing"]):
        for name, q in (("all profiles", g), ("EDS", g[g.method == "EDS"]), ("reflectometry, 500 nm gap", g[(g.H_nm == 500) & (g.method == "reflectometry")])):
            ls = np.log(q.s_fit); c = np.polyfit(q.fill, ls, 1, cov=True) if q.fill.nunique() > 2 else None
            out.append(dict(transport=fam, narrowing=narrow, profiles=name, n=len(q), n_at_grid_edge=int(q.at_grid_edge.sum()), s_geomean=float(f"{np.exp(ls.mean()):.3g}"), s_min=q.s_fit.min(), s_max=q.s_fit.max(), sd_ln_s=round(float(ls.std(ddof=1)), 3),
                            dlns_dfill=round(float(c[0][0]), 2) if c else np.nan, dlns_dfill_se=round(float(np.sqrt(c[1][0, 0])), 2) if c else np.nan, mean_rms=round(float(q.rms.mean()), 4)))
    pd.DataFrame(out).to_csv(ROOT / f"figures/profile_fits_summary{SUF}.csv", index=False)
    # ---- ZnO: flat-plateau profiles of the first blind test ----
    Dg = pd.read_csv(ROOT / "docs/blindtest/posthoc/diagnosis.csv"); Pr = pd.read_csv(ROOT / "docs/blindtest/predictions.csv").set_index("id"); zr = []
    for q in Dg[Dg.id.str.contains("ZnO_DEZ_H2O|haimi2026_fig5") & (Dg.droop.abs() < 0.05)].itertuples():
        d = pd.read_csv(ROOT / "data/blind_revealed" / f"{q.id}.csv").dropna(subset=["x_over_H", "thickness_over_plateau"]); d = d[d.x_over_H >= 0]; spot = Pr.loc[q.id, "spot_units"] if "spot_units" in Pr.columns and np.isfinite(Pr.loc[q.id, "spot_units"]) else 20.0
        K = kern("hyd" if q.z50 / (1 - q.fill / 2) < 300 else "hyd_deep")
        for narrow in (True, False): zr.append(dict(id=q.id, fill=q.fill, z50=q.z50, transport="hydraulic", narrowing=narrow, **fit(d.x_over_H.values, d.thickness_over_plateau.values, q.fill, spot, K, 1.0, narrow)))
        print(q.id, "done", flush=True)
    Z = pd.DataFrame(zr); Z.to_csv(ROOT / f"figures/profile_fits_zno{SUF}.csv", index=False); pd.set_option("display.width", 250)
    print(D[["id", "fill", "transport", "narrowing", "s_fit", "s_p16", "s_p84", "rms", "at_grid_edge", "boot_frac_at_edge"]].to_string(index=False)); print(pd.DataFrame(out).to_string(index=False)); print(Z[["id", "fill", "z50", "narrowing", "s_fit", "s_p16", "s_p84", "rms", "at_grid_edge", "boot_frac_at_edge"]].to_string(index=False))
