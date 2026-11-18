"""Figures 7 and 8 and the summary table of the LHAR re-analysis (draft_v2 Section 3.6), from figures/lhar_*.csv,
figures/ballistic_slope_calibration.csv, figures/pillarhall_fit_*.csv and data/lhar_profiles/. Writes
docs/paper/figures/fig7_slope_relation.png, fig8_lhar_refit.png, lhar_summary.csv, lhar_summary.md and lhar_audit.json.
  PYTHONPATH=src .venv/bin/python docs/paper/lhar_figures.py
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; OUT = ROOT / "docs/paper/figures"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5, "figure.dpi": 150, "savefig.dpi": 200})
audit = {}
cal = pd.read_csv(FIG / "ballistic_slope_calibration.csv"); b = cal[cal.model == "ballistic"]; h = cal[cal.model == "1d_hydraulic"]
st = pd.read_csv(FIG / "lhar_slope_table.csv"); rf = pd.read_csv(FIG / "lhar_refit.csv"); bf = pd.read_csv(FIG / "lhar_single_exposure_fit.csv")
cut = pd.read_csv(FIG / "pillarhall_fit_cutoff.csv"); corr = pd.read_csv(FIG / "pillarhall_fit_corrected.csv")
# ---- empirical ballistic slope relation ----
bb = b[(b.z50 >= 80) & (b.z50 <= 500) & (b.s >= 1e-3)]
X = np.c_[np.ones(len(bb)), np.log(bb.z50), np.log(bb.s)]; k, *_ = np.linalg.lstsq(X, np.log(bb.A), rcond=None)
err = np.exp(np.log(bb.A) - X @ k)
audit["A_continuum_this_code_median"] = round(float(h.A.median()), 2); audit["A_continuum_range"] = [round(float(h.A.min()), 2), round(float(h.A.max()), 2)]
audit["A_ballistic_range_80_500"] = [round(float(bb.A.min()), 1), round(float(bb.A.max()), 1)]
audit["A_ballistic_over_13.9_range"] = [round(float(bb.A.min() / 13.9), 2), round(float(bb.A.max() / 13.9), 2)]
audit["ballistic_relation_prefactor_expz_exps"] = [round(float(np.exp(k[0])), 2), round(float(k[1]), 3), round(float(k[2]), 3)]
audit["ballistic_relation_fit_error_range"] = [round(float(err.min()), 3), round(float(err.max()), 3)]; audit["ballistic_relation_rms_log"] = round(float(np.sqrt(np.mean(np.log(err) ** 2))), 3); audit["ballistic_relation_n"] = int(len(bb))
for s0 in (0.004, 0.008, 0.016):
    d = b[b.s == s0].sort_values("z50"); audit[f"A_ballistic_s{s0}_at_z100_200_300_400"] = [round(float(np.interp(z, d.z50, d.A)), 1) for z in (100, 200, 300, 400)]
for s0 in (0.0002, 0.0005):
    d = b[b.s == s0].sort_values("z50"); audit[f"A_ballistic_s{s0}_at_z500_1000"] = [round(float(np.interp(z, d.z50, d.A)), 1) for z in (500, 1000)]
allb = b[b.z50 >= 80]; audit["A_ballistic_full_range"] = [round(float(allb.A.min()), 1), round(float(allb.A.max()), 1)]
sc = pd.read_csv(FIG / "lhar_single_exposure_scans.csv")
fig, (a, c) = plt.subplots(1, 2, figsize=(10.5, 3.9))
cm = plt.cm.viridis
SS = sorted(b.s.unique())
for i, s0 in enumerate(SS):
    d = b[b.s == s0].sort_values("z50"); a.plot(d.z50, d.A, "-", color=cm(i / (len(SS) - 1)), lw=1.2, label=f"ballistic, s = {s0:g}")
z = np.linspace(60, 1150, 200)
a.axhline(13.9, color="k", ls="--", lw=1); a.text(620, 15.3, "continuum, Arts 2019 (13.9)", fontsize=7)
a.axhline(11.1, color="0.5", ls="--", lw=1); a.text(620, 7.2, "continuum, Gonsalves 2024 (11.1)", fontsize=7, color="0.4")
a.plot(h.z50, h.A, "kx", ms=4, label="1D hydraulic, this code")
a.plot(z, 13.9 * 0.375 * np.log(z), "r:", lw=1.3, label="13.9 x (3/8) ln(z$_{50}$)  (cut-off at the front)")
a.set_xlabel("half-coverage position z$_{50}$ = x$_{50}$/H"); a.set_ylabel("A = s / (d$\\theta$/d$\\tilde x$)$^2$ at $\\theta$ = 1/2"); a.set_ylim(0, 105); a.set_xscale("log"); a.set_xlim(60, 1200); a.set_title("(a) slope relation: ballistic (Machball trench) against continuum"); a.legend(ncol=3, fontsize=6, loc="upper left"); a.grid(alpha=0.3)
ref = st[(st.source == "Yim2020") & (st.H_nm == 500) & (st.magnification == 50) & (st.cycles == 500) & (st.pulse.astype(str) == "0.1")]
def g1(df, **kw):
    m = df
    for kk, vv in kw.items(): m = m[m[kk] == vv]
    return m.iloc[0]
methods = [("s_arts", "hydraulic\n(Arts, 13.9)"), ("s_gonsalves", "hydraulic\n(Gonsalves, 11.1)"), ("s_corr_front", "corrected,\ncut-off at front"), ("s_corr_pillar", "corrected, cut-off\nat pillars or 1 mm"), ("s_ballistic", "ballistic")]
for j_, (col, lab) in enumerate(methods):
    v = ref[col].values; c.plot(np.full(len(v), j_) + np.linspace(-0.14, 0.02, len(v)), v, "o", ms=3.5, color="tab:blue", alpha=0.7, label="Yim 2020, 300 C: 8 samples, tabulated slopes (bar: geometric mean)" if j_ == 0 else None)
    g = np.exp(np.log(v).mean()); c.plot([j_ - 0.25, j_ + 0.1], [g, g], "-", color="tab:blue", lw=2)
    audit[f"ref_group_{col}_gmean"] = float(f"{g:.3e}"); audit[f"ref_group_{col}_spread"] = round(float(np.exp(np.log(v).std(ddof=1))), 2); audit[f"ref_group_{col}_range"] = [float(f"{v.min():.2e}"), float(f"{v.max():.2e}")]
audit["ref_group_samples"] = sorted(int(x) for x in ref["sample"])
v1b = ref[ref.layout == "v1b"]; audit["ref_group_v1b_only_gmean_arts_ballistic"] = [float(f"{np.exp(np.log(v1b.s_arts).mean()):.3e}"), float(f"{np.exp(np.log(v1b.s_ballistic).mean()):.3e}"), int(len(v1b))]
for ds, mk, colr, lab in (("zenodo", "s", "tab:red", "Zenodo 3974438, 300 C"), ("ylilammi_al2o3", "D", "tab:purple", "Ylilammi 2018 Fig. 6, 300 C")):
    xs = [0.22, 3.22, 4.22]; vals = [g1(bf, dataset=ds, model="1d_hydraulic").s, g1(bf, dataset=ds, model="1d_corr_1mm").s, g1(bf, dataset=ds, model="ballistic").s]
    c.plot(xs, vals, mk, ms=6.5, mfc="none", color=colr, label=f"{lab}: full-profile fit, single exposure")
    mh = 8.5e-3 if ds == "zenodo" else g1(rf, dataset=ds, transport="hyd").s; mc = g1(cut, L_phys_m=0.001).s if ds == "zenodo" else g1(rf, dataset=ds, transport="plates:1e-3").s
    c.plot([0.36, 3.36], [mh, mc], mk, ms=6.5, color=colr, label=f"{lab}: full-profile fit, all cycles with narrowing")
ar0 = st[st["sample"] == 101].iloc[0]
c.errorbar([-0.36, 1.64, 3.64], [ar0.s_arts, ar0.s_corr_front, ar0.s_ballistic], yerr=[[ar0.s_arts - ar0.s_arts_lo, ar0.s_corr_front - ar0.s_corr_front_lo, ar0.s_ballistic - ar0.s_ballistic_lo], [ar0.s_arts_hi - ar0.s_arts, ar0.s_corr_front_hi - ar0.s_corr_front, ar0.s_ballistic_hi - ar0.s_ballistic]], fmt="^", ms=6, color="tab:green", capsize=2, label="Arts 2019 Fig. 2, about 220 C, digitised slope")
c.axhspan(3e-3, 5e-3, color="0.85", zorder=0); c.text(4.45, 4.35e-3, "SFG 4 $\\pm$ 1e-3", fontsize=6.5, va="center", ha="right")
c.axhline(3.44e-3, color="0.4", ls=":", lw=1); c.text(4.45, 3.15e-3, "reactor fit 3.4e-3", fontsize=6.5, color="0.3", ha="right", va="center")
c.set_yscale("log"); c.set_ylim(5e-4, 6e-2); c.set_xlim(-0.6, 4.5); c.set_xticks(range(5)); c.set_xticklabels([m[1] for m in methods], fontsize=7); c.set_ylabel("TMA sticking coefficient"); c.set_title("(b) TMA sticking coefficient by transport description and analysis"); c.legend(fontsize=5.8, loc="upper left"); c.grid(axis="y", which="both", alpha=0.3)
fig.tight_layout(); fig.savefig(OUT / "fig7_slope_relation.png"); plt.close(fig)
# ---- Fig. 8: profiles and single-exposure fits ----
from aldconf.data import load_pillarhall, default_path
df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
for sh in sheets:
    gg = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, gg.depth_um, gg.thickness_nm, left=np.nan, right=np.nan))
mean = np.nanmean(np.array(prof), 0); plateau = np.nanmean(mean[(grid > 5) & (grid < 60)])
d6 = pd.read_csv(ROOT / "data/lhar_profiles/ylilammi2018_fig6_al2o3.csv"); pl6 = d6.thickness_nm[(d6.x_um > 5) & (d6.x_um < 60)].mean()
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.7), sharey=True)
for ax, ds, xm, ym, ttl in ((axs[0], "zenodo", grid * 2, mean / plateau, "(a) Zenodo 3974438 (Yim 2020 sample 8), mean of 8 scans"), (axs[1], "ylilammi_al2o3", d6.x_um.values * 2, d6.thickness_nm.values / pl6, "(b) Ylilammi 2018 Fig. 6")):
    ax.plot(xm, ym, "k.", ms=3.5, alpha=0.7, label="measured")
    for mdl, sty, col, lab in (("1d_hydraulic", "--", "tab:blue", "1D, hydraulic D$_K$"), ("1d_corr_1mm", "-.", "tab:green", "1D, corrected (1 mm cut-off)"), ("ballistic", "-", "tab:red", "ballistic (Machball)")):
        cc = pd.read_csv(FIG / f"lhar_single_curve_{ds}_{mdl}.csv"); r = g1(bf, dataset=ds, model=mdl)
        ax.plot(cc.x_um * 2, cc.model_norm, sty, color=col, lw=1.2, label=f"{lab}: s = {r.s:.1e}, {r.exposure_Pa_s:.1f} Pa s, rms {r.rms_all:.3f} (front {r.rms_front:.3f})")
    if ds == "zenodo": ax.axvspan(300, 400, color="0.92", zorder=0); ax.text(350, 0.36, "not fitted\n(scans rise\nagain)", fontsize=6.5, ha="center", color="0.4")
    ax.set_xlim(0, 400); ax.set_ylim(-0.05, 1.45); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_xlabel("depth / gap"); ax.set_title(ttl); ax.legend(fontsize=6.3, loc="upper right", bbox_to_anchor=(1.0, 0.93)); ax.grid(alpha=0.3)
axs[0].set_ylabel("normalised thickness; coverage after one exposure")
fig.tight_layout(); fig.savefig(OUT / "fig8_lhar_refit.png"); plt.close(fig)
# ---- summary table ----
rows = []
def add(profile, basis, transport, s, rms=np.nan, expo=np.nan, note=""):
    rows.append(dict(profile=profile, basis=basis, transport=transport, s=s, rms=rms, exposure_Pa_s=expo, note=note))
for ds, name in (("zenodo", "Zenodo 3974438 (Al2O3 300 C; Yim 2020 sample 8)"), ("ylilammi_al2o3", "Ylilammi 2018 Fig. 6 (Al2O3 300 C)")):
    for mdl, tr in (("1d_hydraulic", "hydraulic"), ("1d_corr_1mm", "corrected, cut-off 1 mm (factor 2.85)"), ("ballistic", "ballistic (Machball trench AR 800)")):
        r = g1(bf, dataset=ds, model=mdl); add(name, "full profile, single exposure, no narrowing" + (", x <= 150 um" if ds == "zenodo" else ""), tr, r.s, r.rms_all, r.exposure_Pa_s, f"rms plateau / front / tail {r.rms_plateau:.3f} / {r.rms_front:.3f} / {r.rms_tail:.3f}")
r = g1(corr, fit="uncorrected_published"); add("Zenodo 3974438 (Al2O3 300 C; Yim 2020 sample 8)", "full profile, 500 cycles with narrowing", "hydraulic", r.s, r.rms, r.exposure_Pa_s)
r = g1(corr, fit="joint_corrected"); add("Zenodo 3974438 (Al2O3 300 C; Yim 2020 sample 8)", "full profile, 500 cycles with narrowing", "corrected, cut-off 250 um (factor 2.33)", r.s, r.rms, r.exposure_Pa_s)
for _, r in cut.sort_values("L_phys_m").iterrows(): add("Zenodo 3974438 (Al2O3 300 C; Yim 2020 sample 8)", "full profile, 500 cycles with narrowing", f"corrected, cut-off {r.L_phys_m*1e3:g} mm (factor {r.factor:.2f})", r.s, r.rms, r.exposure_Pa_s)
for ds, name, cyc in (("ylilammi_al2o3", "Ylilammi 2018 Fig. 6 (Al2O3 300 C)", 500), ("gao_al2o3", "Gao 2015 Fig. 5a (Al2O3 300 C, first-generation LHAR, 107 nm film in a 500 nm gap)", 1109), ("ylilammi_tio2", "Ylilammi 2018 Fig. 7 (TiO2 110 C)", 1000)):
    for tr, lab in (("hyd", "hydraulic"), ("plates:1e-3", "corrected, cut-off 1 mm (factor 2.85)")):
        r = g1(rf, dataset=ds, transport=tr); add(name, f"full profile, {cyc} cycles with narrowing", lab, r.s, r.rms, r.exposure_Pa_s, "poor fit: not Langmuir-like" if "tio2" in ds else "")
for col, lab in (("s_arts", "hydraulic (Arts slope relation, A = 13.9)"), ("s_gonsalves", "hydraulic (Gonsalves slope relation, A = 11.1)"), ("s_corr_front", "corrected, cut-off at the front"), ("s_corr_pillar", "corrected, cut-off at the pillar mean free path"), ("s_ballistic", "ballistic slope relation")):
    v = ref[col].values; add("Yim 2020 Table 1, 8 like-for-like samples (500 nm, 50x, 500 cycles, 0.1 s)", "tabulated front slope; geometric mean", lab, float(np.exp(np.log(v).mean())), note=f"x/÷ {np.exp(np.log(v).std(ddof=1)):.2f}; range {v.min():.1e} to {v.max():.1e}")
ar = st[st.source == "Arts2019"]
for _, r in ar.iterrows():
    for col, lab in (("s_arts", "hydraulic (Arts slope relation)"), ("s_corr_front", "corrected, cut-off at the front"), ("s_ballistic", "ballistic slope relation" + (" (outside calibration)" if r.ballistic_extrapolated else ""))):
        add(f"Arts 2019 Fig. 2, {r.species}-limited (about 220 C), digitised", "front slope from the digitised raster", lab, r[col], note=f"interval {r[col+'_lo']:.1e} to {r[col+'_hi']:.1e}; z50 {r.x50_tilde:.0f}")
add("Vandalon and Kessels 2016, SFG, 300 C", "optical, planar surface", "none", 4e-3, note="4 +/- 1e-3")
add("Gakis et al. 2019, reactor CFD fit, 300 C", "reactor scale", "none", 0.004 * 0.86, note="s0 0.004 x 86 % reaction probability")
S = pd.DataFrame(rows); S.to_csv(OUT / "lhar_summary.csv", index=False)
with open(OUT / "lhar_summary.md", "w") as fh:
    fh.write("| profile | basis | transport model | s | rms | exposure (Pa s) | note |\n|---|---|---|---|---|---|---|\n")
    for _, r in S.iterrows(): fh.write(f"| {r.profile} | {r.basis} | {r.transport} | {r.s:.2e} | {'' if np.isnan(r.rms) else f'{r.rms:.3f}'} | {'' if np.isnan(r.exposure_Pa_s) else f'{r.exposure_Pa_s:.2f}'} | {r.note} |\n")
for ds in ("zenodo", "ylilammi_al2o3"):
    hy = g1(bf, dataset=ds, model="1d_hydraulic"); ba = g1(bf, dataset=ds, model="ballistic"); c1 = g1(bf, dataset=ds, model="1d_corr_1mm")
    audit[f"{ds}_single_exposure_hyd_corr1mm_ball"] = [round(hy.exposure_Pa_s, 2), round(c1.exposure_Pa_s, 2), round(ba.exposure_Pa_s, 2)]; audit[f"{ds}_rms_plateau_front_tail_hyd"] = [round(hy.rms_plateau, 4), round(hy.rms_front, 4), round(hy.rms_tail, 4)]; audit[f"{ds}_rms_plateau_front_tail_ball"] = [round(ba.rms_plateau, 4), round(ba.rms_front, 4), round(ba.rms_tail, 4)]; audit[f"{ds}_msq_ratio_ball_over_hyd"] = round((ba.rms_all / hy.rms_all) ** 2, 2)
    audit[f"{ds}_single_s_hyd_corr1mm_ball"] = [float(f"{hy.s:.3e}"), float(f"{c1.s:.3e}"), float(f"{ba.s:.3e}")]; audit[f"{ds}_single_rms_hyd_corr_ball"] = [round(hy.rms_all, 4), round(c1.rms_all, 4), round(ba.rms_all, 4)]
    audit[f"{ds}_ballistic_over_hydraulic_s"] = round(ba.s / hy.s, 2); audit[f"{ds}_hydraulic_over_ballistic_exposure"] = round(hy.exposure_Pa_s / ba.exposure_Pa_s, 2); audit[f"{ds}_ballistic_over_corr1mm_s"] = round(ba.s / c1.s, 2)
for ds in ("ylilammi_al2o3", "gao_al2o3", "ylilammi_tio2"):
    a1, a2 = g1(rf, dataset=ds, transport="hyd"), g1(rf, dataset=ds, transport="plates:1e-3"); audit[f"{ds}_multicycle_s_hyd_corr"] = [float(f"{a1.s:.3e}"), float(f"{a2.s:.3e}")]; audit[f"{ds}_multicycle_rms"] = round(a1.rms, 4); audit[f"{ds}_multicycle_exposure_hyd_corr"] = [round(a1.exposure_Pa_s, 2), round(a2.exposure_Pa_s, 2)]
audit["arts_tma_z50_slope"] = [round(float(ar.iloc[0].x50_tilde)), round(float(ar.iloc[0].slope), 5)]; audit["arts_h2o_z50_slope"] = [round(float(ar.iloc[1].x50_tilde)), round(float(ar.iloc[1].slope), 5)]
audit["arts_tma_s_arts_corrfront_ball"] = [float(f"{ar.iloc[0][c_]:.2e}") for c_ in ("s_arts", "s_corr_front", "s_ballistic")]; audit["arts_h2o_s_arts_corrfront_ball"] = [float(f"{ar.iloc[1][c_]:.2e}") for c_ in ("s_arts", "s_corr_front", "s_ballistic")]
audit["zenodo_scans_s_hyd_range_gmean"] = [float(f"{sc.s_hydraulic.min():.2e}"), float(f"{sc.s_hydraulic.max():.2e}"), float(f"{np.exp(np.log(sc.s_hydraulic).mean()):.2e}")]; audit["zenodo_scans_s_ball_range_gmean"] = [float(f"{sc.s_ballistic.min():.2e}"), float(f"{sc.s_ballistic.max():.2e}"), float(f"{np.exp(np.log(sc.s_ballistic).mean()):.2e}")]; audit["zenodo_scans_ratio_range"] = [round(float(sc.ratio.min()), 2), round(float(sc.ratio.max()), 2)]; audit["zenodo_scans_prefer_hydraulic"] = [int((sc.rms_hydraulic < sc.rms_ballistic).sum()), int(len(sc))]
audit["set_ratio_hyd_ball"] = [round(audit["ref_group_s_arts_gmean"] / ar.iloc[0].s_arts, 2), round(audit["ref_group_s_ballistic_gmean"] / ar.iloc[0].s_ballistic, 2)]
audit["pillar_mfp_um_v1a_v1b_v2a"] = [round(np.sqrt(3) / 2 * p ** 2 / 4.0) for p in (28, 49, 98)]
json.dump(audit, open(OUT / "lhar_audit.json", "w"), indent=1); print(json.dumps(audit, indent=1)); print(S.to_string(float_format=lambda v: f"{v:.3g}"))
