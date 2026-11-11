#!/usr/bin/env python3
"""Build the live page's data file and static assets from this repository's own results.

  python docs/build_assets.py --data       docs/assets/data.js (window.ALD): chamber fields and wafer coverage from the OpenFOAM
                                           cases, the pulse-time map, the DOE table, validation points, water, plasma and
                                           robustness series, the PillarHall fit, all recomputed or read from figures/*.csv
  python docs/build_assets.py --preview    docs/assets/social-preview.png
  python docs/build_assets.py --snap       static PNG fallbacks for every [data-viz] scene in both themes (Playwright)
  python docs/build_assets.py --hero-png   docs/assets/hero.png, the static hero fallback
  python docs/build_assets.py --qa         both-theme screenshots and console check into docs/qa/
  python docs/build_assets.py --lighthouse desktop Lighthouse run (needs the lighthouse CLI on PATH)
  python docs/build_assets.py              everything

The OpenFOAM time directories are not committed (see .gitignore); --data needs the cases to have been run with reactor/Allrun,
reactor/sweep_pattern_map.py, reactor/make_h2o_desorption_cases.py and the edge-ring / hole-size cases. Everything else reads
figures/*.csv, which are committed, so the page can be rebuilt from the repository alone except for the hero field frames.
"""
from __future__ import annotations
import json, re, sys, subprocess
from pathlib import Path
import numpy as np

DOCS = Path(__file__).resolve().parent; ROOT = DOCS.parent; OUT = DOCS / "assets"; FIG = ROOT / "figures"
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "coupling"), str(ROOT / "reactor")]
KB = 1.380649e-23
# the operating point of every chamber run on the page unless stated otherwise (reactor/sweep_pattern_map.py)
P_TOT, T_K, U, N0, HOLES, D_HOLE, AR = 300.0, 573.15, 0.0501, 4e18, 1.78e13, 100e-9, 100
C0_1PCT = 0.01 * P_TOT / (KB * T_K)                      # 3.79e20 /m3 at 1 % TMA
A_WAFER = np.pi * 0.15 ** 2
WALL = HOLES * np.pi * D_HOLE * AR * D_HOLE               # 55.9: patterned wall area / wafer area


def r(x, d=4):
    if hasattr(x, 'tolist') and not np.isscalar(x): x = x.tolist()
    if isinstance(x, (list, tuple)): return [r(v, d) for v in x]
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), d)


def t90_of(t, series, level=0.9):
    i = np.argmax(series >= level); return float(t[i]) if series.max() >= level else float("nan")


# ---------------- chamber fields and wafer coverage ----------------
def internal_field(path):
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)", path.read_text(), re.S)
    return np.fromstring(m.group(2), sep="\n")


def pool(a, fr, fz):
    nr, nz = a.shape; return a[: nr // fr * fr, : nz // fz * fz].reshape(nr // fr, fr, nz // fz, fz).mean(axis=(1, 3))


def hero_data():
    from patterned_lib import load_patterned
    from postprocess_langmuir import load as load_flat
    case = ROOT / "reactor/showerhead_pulse_patterned_AR100_long"
    t, f, params, pulse = load_patterned(case)
    frames_t = [0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 20.3, 20.6, 21, 22]
    frames = []
    for tt in frames_t:
        d = case / (f"{tt:g}"); c = internal_field(d / "T")
        b0 = pool(c[:4500].reshape(150, 30), 5, 3); b1 = pool(c[4500:].reshape(50, 30), 5, 3)       # [r, z] -> 30x10 and 10x10
        frames.append(r(np.vstack([b0, b1]).ravel(), 3))                                            # 40 radii (5 mm bins) x 10 heights (1.5 mm)
    keep = slice(0, None, 2); rad = slice(0, 150, 5)
    cov = {k: r(f[k][keep][:, rad], 3) for k in ("thetaFlat", "thetaMid", "thetaBottom")}
    tb, tf = f["thetaBottom"], f["thetaFlat"]
    hero = {"pulse_s": pulse, "times": r(t[keep], 3), "radii_mm": r(np.arange(150)[rad] + 0.5, 1), "cov": cov,
            "field_times": frames_t, "field": frames, "field_shape": [40, 10], "field_r_mm": 5.0, "field_z_mm": 1.5,
            "t90_flat_centre": t90_of(t, tf[:, 0]), "t90_bottom_centre": t90_of(t, tb[:, 0]), "t90_bottom_worst": float(np.nanmax([t90_of(t, tb[:, j]) for j in range(149)])),
            "t_supply": (1 + WALL) * N0 / (U * C0_1PCT), "rule_s": 1.2 * (1 + WALL) * N0 / (U * C0_1PCT), "params": {k: params[k] for k in params}}
    # blanket wafer at the same operating point (langmuirSticking, 0.5 s pulse)
    tb_, th, face, s, vbar, c0, N0_, pulse_b = load_flat(ROOT / "reactor/showerhead_pulse_langmuir_p0.5")
    hero["blanket"] = {"times": r(tb_[::2], 3), "centre": r(th[::2, 0], 3), "edge": r(th[::2, 140], 3), "t90_centre": t90_of(tb_, th[:, 0]), "t90_worst": float(np.nanmax([t90_of(tb_, th[:, j]) for j in range(149)])), "pulse_s": pulse_b, "s": s}
    # one hole at the wafer centre, driven by the chamber's face concentration: coverage along the depth vs time (aldconf feature solver)
    from aldconf import Hole, Precursor, Process, simulate_driven
    pr = Precursor(sticking=params["s"], site_density=params["N0"], mass_amu=72.09)
    proc = Process(Hole(D_HOLE, AR * D_HOLE), pr, T=T_K, p_partial=1.0, n_cells=200)
    depths = np.linspace(0, 1, 21)
    x, theta, ht, hist = simulate_driven(proc, t, f["value"][:, 0] * params["c0"], substeps_per_sample=5, record_depths=depths)
    sel = [i for i, tt in enumerate(ht) if any(abs(tt - q) < 1e-6 for q in (0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))]
    hero["hole"] = {"depths": r(depths, 2), "times": r(ht[sel], 3), "theta": r(hist[sel], 3), "t90_bottom": t90_of(ht, hist[:, -1]), "t90_mouth": t90_of(ht, hist[:, 0])}
    print(f"hero: flat t90 {hero['t90_flat_centre']:.2f} s, bottom t90 {hero['t90_bottom_centre']:.2f} s (worst {hero['t90_bottom_worst']:.2f}), standalone hole {hero['hole']['t90_bottom']:.2f} s; rule {hero['rule_s']:.2f} s; blanket {hero['blanket']['t90_centre']:.2f}/{hero['blanket']['t90_worst']:.2f} s")
    return hero


# ---------------- tables from figures/*.csv ----------------
def csv_rows(path):
    import pandas as pd
    return pd.read_csv(path)


def map_data():
    df = csv_rows(FIG / "pulse_time_map.csv")
    rows = [{"kind": a.kind, "x": r(a.x, 4), "AR": int(a.AR), "wall": r(a.wall_ratio, 2), "t_supply": r(a.t_supply, 4), "t90_centre": r(a.t90_centre, 3), "t90_worst": r(a.t90_worst, 3), "t90_flat": r(a.t90_flat, 3), "ratio": r(a.t90_worst / a.t_supply, 3)} for a in df.itertuples()]
    pat = [q["ratio"] for q in rows if q["kind"] == "patterned"]; bl = [q["ratio"] for q in rows if q["kind"] == "blanket"]
    print(f"map: patterned ratio {min(pat):.2f}-{max(pat):.2f}, blanket {min(bl):.2f}-{max(bl):.2f}")
    return {"rows": rows, "patterned_ratio": [r(min(pat), 2), r(max(pat), 2)], "blanket_ratio": [r(min(bl), 2), r(max(bl), 2)]}


def validation_data():
    df = csv_rows(FIG / "pulse_time_map.csv"); V, BETA, VQ, TRES = 1.88e-3, 5e-3, 102.5, 0.53
    tau = df.t90_worst / df.t_supply; Da = ((1 + df.wall_ratio) * A_WAFER / V) * BETA * VQ * TRES; S = (1 + df.wall_ratio) * A_WAFER
    pts = [{"kind": k, "AR": int(a), "x": r(x, 4), "Da": r(d, 2), "tau": r(t, 3), "S_m2": r(s, 3), "util": r(1 / t, 3)} for k, a, x, d, t, s in zip(df.kind, df.AR, df.x, Da, tau, S)]
    t0 = 0.0103 / (0.0115 / 24.0); tau_arg = 24.0 / t0; lo = 24.0 / (0.0103 / ((0.0115 + 0.0020) / 24.0)); hi = 24.0 / (0.0103 / ((0.0115 - 0.0020) / 24.0))
    drum_mass = [2, 5, 15, 30, 75]; drum_util = [0.247, 0.479, 0.467, 0.553, 0.707]
    DaT = np.logspace(0.5, 4.5, 60); theory = 1 / np.clip(1 - np.log(2) / DaT, 0.05, None)
    return {"points": pts, "fluidized": {"Da": 1e4, "tau": r(tau_arg, 3), "tau_lo": r(min(lo, hi), 3), "tau_hi": r(max(lo, hi), 3), "S_m2": 940, "util": 0.918, "util_lo": 0.146, "util_hi": 0.082, "source": "Lu et al., JVST A 40 (2022), OSTI 1902316"},
            "drum": {"S_m2": [m * 20.0 for m in drum_mass], "util": drum_util, "mass_g": drum_mass, "source": "Coile et al., JVST A 38, 052403 (2020), Table 2"},
            "theory": {"Da": r(DaT, 2), "tau": r(theory, 3)}}


def doe_data():
    df = csv_rows(FIG / "doe_results.csv")
    runs = [{"run": int(a.run), "x": r(a.x, 4), "Q": r(a.Q, 1), "p": r(a.p, 1), "g": r(a.g, 1), "cx": int(a.c_x), "cQ": int(a.c_Q), "cp": int(a.c_p), "cg": int(a.c_g), "t90": r(a.t90_worst, 3), "edge_lag": r(a.edge_lag_pct, 2), "util": r(a.utilization, 3), "tau": r(a.tau_sat, 3)} for a in df.itertuples()]
    rsm = {}
    for blk in re.split(r"\n== ", (FIG / "doe_rsm_report.txt").read_text())[1:]:
        name = blk.split(":")[0]; coef = {m.group(1): float(m.group(2)) for m in re.finditer(r"\n\s+(\w+) = ([-0-9.eE]+)", blk)}; r2 = float(re.search(r"RSquare ([0-9.]+)", blk).group(1)); rsm[name] = {"coef": coef, "r2": r2}
    gp = (FIG / "doe_gp_report.txt").read_text()
    print("doe:", len(runs), "runs; rsm", {k: v["r2"] for k, v in rsm.items()})
    return {"runs": runs, "rsm": rsm, "gp_report": gp.strip(), "ranges": {"x": [0.005, 0.03], "Q": [150, 600], "p": [100, 1000], "g": [10, 25]}, "util_range": [r(df.utilization.min(), 3), r(df.utilization.max(), 3)]}


# ---------------- water half-cycle with reversible adsorption ----------------
def water_data():
    from postprocess_langmuir import load as load_flat
    from patterned_lib import _block, _field, load_patterned
    out = []
    for k in (0, 0.1, 1, 10):
        case = ROOT / f"reactor/showerhead_h2o_blanket_kdes{k:g}"
        t, th, face, s, vbar, c0, N0_, pulse = load_flat(case)
        times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
        o = np.array([_field(_block((d / "T").read_text(), "outlet"), "value", 30).mean() for _, d in times])
        i_end = np.argmin(abs(t - pulse)); out.append({"kdes": k, "times": r(t[::4], 3), "centre": r(th[::4, 0], 3), "outlet": r(o[::4], 4), "t90": r(t90_of(t, th[:, 0]), 3), "theta_end": r(th[i_end, 0], 3), "theta_final": r(th[-1, 0], 3), "pulse_s": pulse})
    tp, fp, pp, pulse_p = load_patterned(ROOT / "reactor/showerhead_h2o_patterned_AR100")
    tb, fb, _, _ = load_patterned(ROOT / "reactor/showerhead_pulse_patterned_AR100_long")
    res = {"kdes": out, "s_h2o": out[0] and 1e-4, "patterned_h2o_t90_bottom": r(t90_of(tp, fp["thetaBottom"][:, 0]), 3), "patterned_tma_t90_bottom": r(t90_of(tb, fb["thetaBottom"][:, 0]), 3),
           "patterned_h2o_profile": {"times": r(tp[::3], 2), "bottom": r(fp["thetaBottom"][::3, 0], 3), "flat": r(fp["thetaFlat"][::3, 0], 3)}}
    print("water:", [(q["kdes"], q["t90"]) for q in out], "patterned H2O bottom t90", res["patterned_h2o_t90_bottom"], "TMA", res["patterned_tma_t90_bottom"])
    return res


# ---------------- plasma ALD: recombination ceiling (coupling/peald_ceiling.py, recomputed) ----------------
def peald_data():
    from scipy.optimize import brentq
    from aldconf import Hole, Precursor, Process, simulate, step_coverage, dose_for_coverage
    D_, T_, P_, S_ = 100e-9, 573.15, 0.3, 0.05
    RS = [0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]; ARS = (20, 50, 100)
    def proc(ar, rr, t):
        return Process(Hole(D_, ar * D_), Precursor(name="O", sticking=S_, recombination=rr, mass_amu=16.0, site_density=4e18), T=T_, p_partial=P_, t_dose=t, n_cells=200, n_substeps=150, narrowing=False)
    t90 = {str(ar): [r(dose_for_coverage(proc(ar, rr, 1.0), target=0.90, t_lo=1e-3, t_hi=300.0), 3) for rr in RS] for ar in ARS}
    TEXP = 5.0
    def reach(rr):
        f = lambda ar: step_coverage(simulate(proc(ar, rr, TEXP))) - 0.9
        lo, hi = 5.0, 2000.0
        if f(hi) > 0: return hi
        if f(lo) < 0: return float("nan")
        return brentq(f, lo, hi, xtol=1.0)
    rr_ = np.logspace(-5, -1.5, 8); armax = [r(reach(x), 1) for x in rr_]
    print("peald: t90 AR100", t90["100"], "ceiling", armax)
    return {"r": RS, "t90": t90, "ceiling_r": r(rr_, 7), "ceiling_AR": armax, "t_exp": TEXP, "s": S_, "p_Pa": P_}


# ---------------- robustness: edge ring, temperature, hole size ----------------
def robust_data():
    from patterned_lib import load_patterned
    df = csv_rows(FIG / "edgering_t90_profiles.csv")
    ring = {c: {"centre": r(df[c].iloc[0], 3), "worst": r(df[c].iloc[:149].max(), 3)} for c in df.columns[1:]}
    tdf = csv_rows(FIG / "temperature_sweep.csv")
    temp = [{"T_C": int(a.T_C), "blanket_centre": r(a.blanket_t90_centre, 3), "blanket_worst": r(a.blanket_t90_worst, 3), "bottom_centre": r(a.patterned_bottom_t90_centre, 3), "bottom_worst": r(a.patterned_bottom_t90_worst, 3)} for a in tdf.itertuples()]
    holes = []
    for d_nm, case in ((50, "holesize_d50nm_AR100"), (200, "holesize_d200nm_AR100")):
        t, f, p, _ = load_patterned(ROOT / "reactor" / case); holes.append({"d_nm": d_nm, "centre": r(t90_of(t, f["thetaBottom"][:, 0]), 3), "worst": r(np.nanmax([t90_of(t, f["thetaBottom"][:, j]) for j in range(149)]), 3), "holeDensity": p["holeDensity"]})
    print("robust: ring", ring, "temp", [(q["T_C"], q["bottom_centre"]) for q in temp], "holes", [(q["d_nm"], q["centre"]) for q in holes])
    return {"edge_ring": ring, "temperature": temp, "hole_size": holes}


# ---------------- PillarHall fit (coupling/fit_pillarhall.py at the fitted parameters) ----------------
def pillarhall_data():
    from aldconf import Slit, Precursor, Process, simulate
    from aldconf.data import load_pillarhall, default_path
    df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
    sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0)
    prof = np.array([np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan) for g in (good[good.sheet == sh].sort_values("depth_um") for sh in sheets)])
    mean = np.nanmean(prof, 0); std = np.nanstd(prof, 0); plateau = np.nanmean(mean[(grid > 5) & (grid < 60)]); y = mean / plateau
    GAP, L, CYCLES, TPULSE = 500e-9, 250e-6, 500, 0.1; GPC = plateau * 1e-9 / CYCLES
    S_FIT, P_FIT = 8.5e-3, 75.0                      # figures/pillarhall_fit.png; exposure 7.5 Pa s per 0.1 s pulse
    def model(s, p):
        pr = Precursor(sticking=s, site_density=4e18, mass_amu=72.09, gpc=GPC); res = simulate(Process(Slit(GAP, L), pr, T=T_K, p_partial=p, t_dose=TPULSE, cycles=CYCLES, n_cells=400, n_substeps=60, narrowing=True))
        return np.interp(grid * 1e-6, res.x, res.thickness / res.thickness[:10].mean())
    fit = model(S_FIT, P_FIT); ok = ~np.isnan(y); rms = float(np.sqrt(np.mean((fit[ok] - y[ok]) ** 2)))
    # the same profile with the ideal-plate slit transport (coupling/fit_pillarhall_corrected.py): s 1.96e-2 at 32.3 Pa, same residual
    S_FIT_PL, P_FIT_PL = 1.96e-2, 32.3
    def model_plates(s, p):
        pr = Precursor(sticking=s, site_density=4e18, mass_amu=72.09, gpc=GPC); res = simulate(Process(Slit(GAP, L, knudsen_correction="plates"), pr, T=T_K, p_partial=p, t_dose=TPULSE, cycles=CYCLES, n_cells=400, n_substeps=60, narrowing=True))
        return np.interp(grid * 1e-6, res.x, res.thickness / res.thickness[:10].mean())
    fit_pl = model_plates(S_FIT_PL, P_FIT_PL); rms_pl = float(np.sqrt(np.mean((fit_pl[ok] - y[ok]) ** 2)))
    alt = {"1e-3": model(1e-3, P_FIT), "1e-2": model(1e-2, P_FIT)}
    print(f"pillarhall: {len(sheets)} sheets, plateau {plateau:.2f} nm, rms at fit {rms:.4f}")
    return {"AR": r(grid * 2, 1), "measured": r(np.where(ok, y, np.nan), 4), "std": r(std / plateau, 4), "fit": r(fit, 4), "fit_plates": r(fit_pl, 4), "s_fit_plates": S_FIT_PL, "p_fit_plates_Pa": P_FIT_PL, "rms_plates": r(rms_pl, 4), "alt": {k: r(v, 4) for k, v in alt.items()}, "s_fit": S_FIT, "p_fit_Pa": P_FIT, "rms": r(rms, 4), "plateau_nm": r(plateau, 2), "n_sheets": len(sheets), "source": "Zenodo 3974438, Puurunen and Ylivaara (CC BY 4.0)"}



# ---------------- night three: the seven follow-up packages, as scene data ----------------
def night3_data():
    import pandas as pd
    from scipy.spatial import Delaunay
    out = {}
    # 1. sticking sensitivity
    df = pd.read_csv(FIG / "sticking_sensitivity.csv")
    out["sticking"] = [{"kind": a.kind, "s": r(a.s, 5), "t90_centre": r(a.t90_centre, 3), "t90_worst": r(a.t90_worst, 3), "t90_flat": r(a.t90_flat_centre, 3) if a.t90_flat_centre == a.t90_flat_centre else None, "t_supply": r(a.t_supply, 4), "util": r(a.utilization, 3)} for a in df.itertuples()]
    out["sticking_refs"] = {"arts_band": [5e-4, 2e-3], "ylilammi": 5.72e-3, "vandalon": 4e-3, "pillarhall": 8.5e-3, "pillarhall_plates": 1.96e-2, "runs": 5e-3}
    # 2. step coverage and thickness after 50 cycles
    sc = pd.read_csv(FIG / "step_coverage_vs_pulse.csv"); curves = []
    for case, g in sc.groupby("case", sort=False):
        if "long" in case: continue
        g = g[g.t_s > 0.05]; idx = np.unique(np.geomspace(1, len(g), 70).astype(int) - 1)
        curves.append({"AR": int(g.AR.iloc[0]), "t": r(g.t_s.values[idx], 3), "sc": r(g.step_coverage_centre.values[idx], 4), "thick": r(g.thickness_bottom_50cyc_nm.values[idx], 3), "top": r(g.thickness_top_50cyc_nm.values[idx], 3)})
    tab = pd.read_csv(FIG / "step_coverage_table.csv") if (FIG / "step_coverage_table.csv").exists() else None
    out["stepcov"] = {"curves": curves, "table": tab.to_dict("records") if tab is not None else None, "gpc_nm": 0.103}
    # 3. cycle time and throughput
    cy = pd.read_csv(FIG / "cycle_time_wph.csv"); out["cycle"] = [{"recipe": a.recipe, "basis": a.purge_basis, "seg": r([a.t_TMA_s, a.t_purge1_s, a.t_H2O_s, a.t_purge2_s], 3), "cycle": r(a.cycle_s, 2), "nm": a.thickness_nm, "cycles": int(a.cycles), "wph": r(a.wph_no_overhead, 2), "wph60": r(a.wph_60s_overhead, 2)} for a in cy.itertuples()]
    w = pd.read_csv(FIG / "wph_vs_wall_ratio.csv"); out["wph"] = {"x": sorted(set(r(w.x.tolist(), 3))), "rows": [{"x": r(a.x, 3), "wall": r(a.wall_ratio, 2), "cycle": r(a.cycle_s, 2), "wph60": r(a.wph_60s_overhead, 3), "wph0": r(a.wph_no_overhead, 3)} for a in w.itertuples()][::2]}
    # 4. blanket saturation against literature
    bl = pd.read_csv(FIG / "blanket_saturation_literature.csv")
    out["exposure"] = [{"source": a.source, "T": r(a.T_C, 0), "kind": a.kind, "Pa_s": r(a.exposure_Pa_s, 4), "L": r(a.exposure_L, 0), "s": (r(a.s, 5) if a.s == a.s else None), "note": str(a.note)[:160]} for a in bl.itertuples()]
    # 5. Gakis reproduction: final coverage across the wafer, the surface pressure trace, the measured points
    from patterned_lib import _block, _field
    gk = {"cases": [], "meas_pos_mm": [-95.5, -76.4, -57.2, -38.2, -19.3, -0.7, 18.4, 37.4, 56.6, 75.5, 94.7], "meas_gpc": [0.9751, 0.9834, 0.9873, 0.9927, 0.9951, 0.9961, 0.9951, 0.9982, 0.9961, 0.9927, 0.9826], "meas_err": 0.003, "gpc_sat": 1.0, "source": "Gakis et al. 2018 (CERD 132, 795) and 2019 (CES 195, 399); measured points digitised from Fig. 10 of the 2019 manuscript"}
    for name, lab in (("pulse_25ms", "25 ms valve, s 3.44e-3"), ("pulse_60ms_x1.5", "60 ms, 1.5x amount"), ("pulse_60ms_x2.4", "60 ms, 2.4x amount"), ("pulse_25ms_s5e-3", "25 ms, s 5e-3")):
        case = ROOT / "reactor/validation_gakis" / name
        times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
        th_end = _field(_block((times[-1][1] / "T").read_text()), "theta", 100); tr, pc = [], []
        for tt, d in times:
            blk = _block((d / "T").read_text()); v = _field(blk, "value", 100); tr.append(tt); pc.append(float(v[0]) * 3.371e21 * KB * 573.15)
        gk["cases"].append({"name": lab, "r_mm": r(np.arange(100) + 0.5, 1), "theta_end": r(th_end, 5), "t": r(tr, 3), "p_centre_Pa": r(pc, 3)})
    out["gakis"] = gk
    # 6. 3D chamber: every second wafer face, triangulated, plus azimuthal bins
    c3 = pd.read_csv(FIG / "chamber3d_t90.csv"); c3 = c3[c3.r_mm <= 146.5].iloc[::3].reset_index(drop=True)   # the outermost face row (the edge drop) is left out so the triangulated rim stays clean
    tri = Delaunay(np.c_[c3.y_mm, c3.z_mm]).simplices; cen = np.c_[c3.y_mm, c3.z_mm][tri].mean(1); tri = tri[np.hypot(cen[:, 0], cen[:, 1]) <= 146.0]
    full = pd.read_csv(FIG / "chamber3d_t90.csv"); az = {}
    for R0 in (50, 100, 140):
        sel = full[(full.r_mm - R0).abs() <= 2.5]; bins = np.arange(-180, 181, 15); mid = []; val = []
        for a, b in zip(bins[:-1], bins[1:]):
            q = sel[(sel.azimuth_deg >= a) & (sel.azimuth_deg < b)]
            if len(q): mid.append((a + b) / 2); val.append(float(q.t90_bottom_s.mean()))
        az[str(R0)] = {"az": r(mid, 1), "t90": r(val, 4)}
    out["chamber3d"] = {"y": r(c3.y_mm, 1), "z": r(c3.z_mm, 1), "t90": r(c3.t90_bottom_s, 3), "dev": r(c3.dev_from_wedge_pct, 2), "tri": [int(v) for v in tri.ravel()], "n": int(len(c3)), "azimuthal": az, "port_deg": [-30, 30], "mean": r(full.t90_bottom_s.mean(), 3), "worst": r(full.t90_bottom_s.max(), 3), "wedge_worst": 14.05, "cells": 105984}
    # 7. tapered holes: coverage along the depth at 10 s and t90 per taper, chamber-driven
    from aldconf import TaperedHole, Precursor, Process, simulate_driven
    from patterned_lib import load_patterned
    t, f, params, pulse = load_patterned(ROOT / "reactor/showerhead_pulse_patterned_AR100_long"); n0 = f["value"][:, 0] * params["c0"]
    pr = Precursor(sticking=params["s"], site_density=params["N0"], mass_amu=72.09); depths = np.linspace(0, 1, 21); tp = []
    for ratio in (1.0, 0.9, 0.8, 0.6):
        proc = Process(TaperedHole(D_HOLE, AR * D_HOLE, bottom_ratio=ratio), pr, T=T_K, p_partial=1.0, n_cells=200)
        x, th, ht, hist = simulate_driven(proc, t, n0, substeps_per_sample=5, record_depths=depths)
        snaps = {str(q): r(hist[np.argmin(abs(ht - q))], 3) for q in (5, 10, 13)}
        tp.append({"ratio": ratio, "t90_bottom": r(t90_of(ht, hist[:, -1]), 3), "snap": snaps, "area_rel": r((1 + ratio) / 2, 3)})
    out["taper"] = {"depths": r(depths, 2), "holes": tp, "table": pd.read_csv(FIG / "tapered_hole.csv").to_dict("records")}
    # 8. word-line recesses: snapshots for the 100 nm slit with 365 nm recesses
    from aldconf import RecessedSlit, simulate_recessed
    g = RecessedSlit(100e-9, 5e-6, 29e-9, 365e-9, 51e-9)
    rr = simulate_recessed(g, pr, T_K, t, n0, substeps_per_sample=5, n_slit=200, n_recess=50, snapshot_times=(3.0, 6.0, 9.0, 12.0))
    si = np.linspace(0, 199, 40).astype(int); ri = np.linspace(0, 49, 10).astype(int)
    snaps = {str(k): {"slit": r(np.asarray(v[0])[si], 3), "recess": r(np.asarray(v[1])[np.ix_(si, ri)], 3)} for k, v in rr.snapshots.items()}
    out["recess"] = {"geom": {"w_nm": 100, "L_um": 5, "g_nm": 29, "L_r_nm": 365, "pitch_nm": 51}, "x_um": r(np.asarray(rr.x)[si] * 1e6, 3), "y_nm": r(np.asarray(rr.y)[ri] * 1e9, 1), "snaps": snaps, "table": pd.read_csv(FIG / "wordline_recess.csv").to_dict("records")}
    # 9. chemistry
    ch = pd.read_csv(FIG / "chemistry_rule_table.csv"); td = json.loads((ROOT / "reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah/results.json").read_text())
    out["chemistry"] = {"rows": [{"film": a.film, "precursor": a.precursor, "M": r(a.M_g_mol, 2), "T_C": a.T_C, "s": r(a.s, 5), "gpc": r(a.GPC_nm, 3), "N0": a.N0_per_m2, "vbar": r(a.vbar_m_s, 1), "rule_AR100": r(a.rule_AR100_wafer_s, 3), "rule_blanket": r(a.rule_blanket_s, 4), "note": str(a.note)[:200]} for a in ch.itertuples()], "tdmah_run": td}
    return out

def build_data():
    hero = hero_data()
    data = {"geometry": {"wafer_mm": 300, "chamber_mm": 400, "gap_mm": 15, "hole_nm": 100, "AR": AR, "pitch_nm": 180, "hole_density": HOLES, "open_area": r(HOLES * np.pi * (D_HOLE / 2) ** 2, 3), "wall_ratio": r(WALL, 2), "p_Pa": P_TOT, "T_K": T_K, "Q_sccm": 300, "U_m_s": U, "cells": 6000, "N0": N0, "c0_1pct": C0_1PCT},
            "hero": hero, "map": map_data(), "validation": validation_data(), "doe": doe_data(), "water": water_data(), "peald": peald_data(), "robust": robust_data(), "pillarhall": pillarhall_data(),
            "machball": {"AR": [20, 50, 100], "this": [1.20, 4.45, 16.2], "machball": [1.29, 4.46, 15.6], "source": "tests/test_machball_crosscheck.py"},
            "spatial": {"N0": N0, "U": U, "c0": C0_1PCT, "t_purge": 0.5, "f_zone": 0.25, "tau_rot": 1.0, "ratios": [10, 30, 100]}}
    # Python GP surrogate export (analysis/doe_gp.py) and its Pareto front, when present
    gpj = FIG / "doe_gp_python.json"
    if gpj.exists():
        G = json.loads(gpj.read_text())
        data["gp"] = {"software": G.get("software"), "coding": G.get("coding"), "factor_ranges": G.get("factor_ranges"),
                      "responses": [{"name": q["name"], "transform": q.get("transform"), "linear": q["linear"]["coef"], "constant": q["kernel"]["constant_value"], "length_scale": q["kernel"]["length_scale"], "noise": q["kernel"]["noise_level"],
                                     "X": r(q["X_train"], 4), "alpha": r(q["alpha"], 6), "K_inv": r(q["K_inv"], 5), "y_mean": q["y_normalisation"]["mean"], "y_std": q["y_normalisation"]["std"], "loo_r2": r(q["loo_r2"], 4), "loo_rmse": r(q["loo_rmse"], 5)} for q in G["responses"]]}
    lay = FIG / "pattern_layout.csv"
    if lay.exists():
        L = csv_rows(lay); rad = L.radius_mm.tolist()
        def col(c): return [None if (v != v) else r(v, 3) for v in L[c].tolist()]
        data["layout"] = {"radius_mm": r(rad, 1), "dense": {"AR": col("dense_centre_AR"), "density": col("dense_centre_holeDensity_per_m2"), "t90_flat": col("dense_centre_t90_flat_s"), "t90_bottom": col("dense_centre_t90_bottom_s"), "wall_avg": 19.2, "rule_avg": 5.10, "rule_local": 14.39, "t90_centre": 11.85, "t90_array_edge": 8.75},
                          "ramp": {"AR": col("ar_ramp_AR"), "density": col("ar_ramp_holeDensity_per_m2"), "t90_flat": col("ar_ramp_t90_flat_s"), "t90_bottom": col("ar_ramp_t90_bottom_s"), "wall_avg": 41.0, "rule_avg": 10.62, "rule_local": 14.36, "t90_centre": 12.50, "t90_edge": 7.60},
                          "source": "reactor/make_layout_cases.py, coupling/postprocess_layout.py, paper/assumptions.md night 3 part 1"}
    par = FIG / "doe_pareto.csv"
    if par.exists():
        pf = csv_rows(par); data["pareto"] = [{k: r(getattr(a, k), 4) for k in ("t90_s", "util", "precursor_over_minimum", "x", "Q_sccm", "p_Pa", "g_mm", "tma_sccm_s")} for a in pf.itertuples()]
    data["n3"] = night3_data()
    data["stats"] = {"blanket_t90": hero["blanket"]["t90_centre"], "blanket_t90_worst": hero["blanket"]["t90_worst"], "patterned_t90": hero["t90_bottom_centre"], "patterned_t90_worst": hero["t90_bottom_worst"], "wall_ratio": r(WALL, 1),
                     "underdose": r(hero["t90_bottom_centre"] / hero["blanket"]["t90_centre"], 1), "rule_band": data["map"]["patterned_ratio"], "blanket_band": data["map"]["blanket_ratio"], "fluidized_tau": data["validation"]["fluidized"]["tau"],
                     "util_range": data["doe"]["util_range"], "h2o_t90": data["water"]["patterned_h2o_t90_bottom"], "tma_t90": data["water"]["patterned_tma_t90_bottom"]}
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "data.js", "w", encoding="utf-8") as fh:
        fh.write("// Generated by docs/build_assets.py from the OpenFOAM cases, figures/*.csv and the aldconf feature solver. Do not edit by hand.\n")
        fh.write("window.ALD = " + json.dumps(data, separators=(",", ":")) + ";\n")
    print("data.js", (OUT / "data.js").stat().st_size, "bytes", data["stats"])


# ---------------- social preview ----------------
def build_preview():
    from PIL import Image, ImageDraw, ImageFont
    BG, INK, MUTED, ACC = "#0b1220", "#e6edf3", "#8b98a9", "#34d399"
    im = Image.new("RGB", (1280, 640), BG); d = ImageDraw.Draw(im)
    def font(sz, bold=True):
        for f in (["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/Library/Fonts/Arial Bold.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"] if bold else ["/System/Library/Fonts/Supplemental/Arial.ttf", "/Library/Fonts/Arial.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"]):
            try: return ImageFont.truetype(f, sz)
            except Exception: pass
        return ImageFont.load_default()
    d.text((70, 70), "SOHAM KAVATHEKAR  ·  ALD CHAMBER-TO-FEATURE MODEL", fill=MUTED, font=font(20))
    d.text((70, 120), "Every blanket-wafer ALD recipe", fill=INK, font=font(50)); d.text((70, 182), "hides a 40x under-dose.", fill=INK, font=font(50))
    d.text((70, 262), "OpenFOAM showerhead chamber, one hole model per point on the wafer, two-way coupled.", fill="#c9d4e2", font=font(21, False))
    d.text((70, 292), "The pulse that saturates a 3D-NAND wafer is 1.2 x (1 + wall area ratio) x one monolayer of supply.", fill="#c9d4e2", font=font(21, False))
    for i, (big, small) in enumerate([("0.35 s → 14 s", "blanket wafer vs AR-100 patterned wafer"), ("within 25 %", "rule across AR 20-100, 20x feed range"), ("τ = 1.12", "measured fluidized-bed ALD, 940 m², on the same axis")]):
        x = 70 + i * 390; d.rectangle([x, 380, x + 300, 384], fill=ACC); d.text((x, 400), big, fill=INK, font=font(30)); d.text((x, 444), small, fill=MUTED, font=font(15, False))
    # right: a schematic hole with the coverage front (no render dependency)
    hx, hy, hw, hh = 1040, 70, 110, 280
    d.rectangle([hx - 60, hy - 10, hx + hw + 60, hy + 12], fill="#1f2a3c"); d.rectangle([hx, hy, hx + hw, hy + hh], outline="#2a3649", width=3)
    for k in range(10):
        y0 = hy + k * hh / 10; c = ACC if k < 6 else ("#5b7fb5" if k < 8 else "#1f2a3c"); d.rectangle([hx + 4, y0 + 2, hx + hw - 4, y0 + hh / 10 - 2], fill=c)
    d.text((hx - 60, hy + hh + 14), "one hole, AR 100: the coverage", fill=MUTED, font=font(14, False)); d.text((hx - 60, hy + hh + 34), "front at 8 s of a 14 s pulse", fill=MUTED, font=font(14, False))
    d.text((70, 560), "github.com/Soham-ChemE/ald-conformality  ·  soham-cheme.github.io/ald-conformality", fill=MUTED, font=font(18, False))
    im.save(OUT / "social-preview.png", optimize=True); print("wrote social-preview.png")


def quantise(paths, colors=256, max_w=1400):
    from PIL import Image
    for p in paths:
        im = Image.open(p).convert("RGBA")
        if im.width > max_w: im = im.resize((max_w, round(im.height * max_w / im.width)), Image.LANCZOS)
        im.convert("RGB").quantize(colors=colors, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG).save(p, optimize=True); print("quantised", p.name, p.stat().st_size, "bytes")


def qa_readme():
    rep = json.loads((DOCS / "qa" / "qa.json").read_text()); lines = ["# QA for the live page", "", "Served from docs/ by docs/tools/page_qa.py (Playwright for Python) and loaded in headless Chromium. The screenshots use software WebGL with `window.__forceGL` so the live scenes render (frame timing is therefore not representative); the Lighthouse run has WebGL disabled, so it measures the page's static path, as the sibling pages' reports do. Screenshots in both themes: the hero, the first two body scenes, the 375 px hero.", ""]
    for page, pr in rep.get("pages", {}).items():
        lines += [f"## {page}", ""]
        if "lighthouse" in pr:
            lh = pr["lighthouse"]; lines += [f"Lighthouse (desktop preset): performance {lh['performance']}, accessibility {lh['accessibility']}, best practices {lh['bestPractices']}, SEO {lh['seo']}", "", "| Metric | Value |", "|---|---|", f"| First Contentful Paint | {lh['fcp']} |", f"| Largest Contentful Paint | {lh['lcp']} |", f"| Speed Index | {lh['si']} |", f"| Total Blocking Time | {lh['tbt']} |", f"| Cumulative Layout Shift | {lh['cls']} |", f"| Total byte weight | {lh['bytes']} |", ""]
        lines += ["| Theme | Console errors | Failed requests | Transfer (bytes) | Body scenes captured |", "|---|---|---|---|---|"]
        for theme in ("dark", "light"):
            t = pr.get(theme)
            if t: lines.append(f"| {theme} | {len(t['errors']) + len(t.get('mobileErrors', []))} | {len(t['failed'])} | {t['transferBytes']:,} | {', '.join(t['scenes'])} |")
        lines.append("")
    lines += ["Files: `light/` and `dark/` hold `index-hero.png`, `index-scene1-*.png`, `index-scene2-*.png` and `index-mobile-375-hero.png`; `lighthouse.report.html` and `.json` are the Lighthouse run; `qa.json` is the machine-readable summary."]
    (DOCS / "qa" / "README.md").write_text("\n".join(lines).replace("\u00a0", " ") + "\n", encoding="utf-8"); print("wrote qa/README.md")


FIGS_ON_PAGE = ["sticking_sensitivity", "step_coverage_vs_pulse", "cycle_time_wph", "blanket_saturation_literature", "validation_gakis", "chamber3d_t90_map", "tapered_hole", "wordline_recess"]


def copy_figures():
    """The night-three result figures shown as static cards on the page: copied from figures/ and quantised (Pages serves docs/ only)."""
    from PIL import Image
    for f in FIGS_ON_PAGE:
        src = FIG / f"{f}.png"
        if not src.exists(): continue
        im = Image.open(src).convert("RGB")
        if im.width > 1600: im = im.resize((1600, round(im.height * 1600 / im.width)), Image.LANCZOS)
        im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.FLOYDSTEINBERG).save(OUT / f"fig_{f}.png", optimize=True); print("copied", f)


def main(argv):
    only = [a for a in argv if a.startswith("--")]; do_all = not only
    if do_all or "--data" in only: build_data()
    if do_all or "--preview" in only: build_preview()
    if do_all or "--figs" in only: copy_figures()
    tool = [sys.executable, str(DOCS / "tools" / "page_qa.py"), "--docs", str(DOCS)]
    if do_all or "--hero-png" in only: subprocess.run(tool + ["--hero-png"], check=True); quantise([OUT / "hero.png"], max_w=1600)
    if do_all or "--snap" in only: subprocess.run(tool + ["--snap"], check=True); quantise(sorted(OUT.glob("viz_*_*.png")))
    flags = [f for f in ("--qa", "--lighthouse") if do_all or f in only]
    if flags: subprocess.run(tool + flags, check=True); qa_readme()


if __name__ == "__main__":
    main(sys.argv[1:])
