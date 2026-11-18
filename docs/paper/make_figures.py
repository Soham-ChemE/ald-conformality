"""Figures of docs/paper/draft_v1.md, regenerated from the repository's result CSVs (figures/*.csv) and, for Fig. 6, from
the Zenodo 3974438 workbook through the aldconf solver. Every number plotted here is also listed in docs/paper/NUMBERS_AUDIT.md.

  PYTHONPATH=src .venv/bin/python docs/paper/make_figures.py        (about 1 min; Fig. 6 reruns the two PillarHall fits' curves)
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter

ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
FIG = ROOT / "figures"; OUT = ROOT / "docs" / "paper" / "figures"; OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5, "figure.dpi": 150, "savefig.dpi": 200})
audit = {}

# ---------- Fig. 1: pulse-time map and step coverage ----------
pm = pd.read_csv(FIG / "pulse_time_map.csv")
pm["tau"] = pm.t90_worst / pm.t_supply
fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.6))
for kind, mk, lab in (("blanket", "o", "blanket wafer"), ("patterned", "s", "patterned wafer, AR 20 / 50 / 100")):
    d = pm[pm.kind == kind]; a.loglog(d.t_supply, d.t90_worst, mk, ms=6, label=lab, mfc="none" if kind == "blanket" else None)
t = np.logspace(np.log10(0.015), np.log10(30), 50)
a.loglog(t, t, "k:", lw=1, label="supply limit  (1 + A$_{wall}$/A$_{wafer}$) N$_0$/(U c$_0$)")
a.loglog(t, 1.2 * t, "k-", lw=1, label="rule, 1.2 x supply limit")
a.set_xlabel("supply-limit time  (1 + A$_{wall}$/A$_{wafer}$) N$_0$ / (U c$_0$)  (s)"); a.set_ylabel("time to 90 % bottom coverage, worst face (s)")
a.set_title("(a) 16 chamber runs, 0.5-10 % TMA, 300 Pa, 300 sccm, 573 K"); a.legend(loc="upper left"); a.grid(True, which="both", alpha=0.3)
audit["fig1a_tau_patterned_min_max"] = [round(pm[pm.kind == "patterned"].tau.min(), 3), round(pm[pm.kind == "patterned"].tau.max(), 3)]
audit["fig1a_tau_blanket_min_max"] = [round(pm[pm.kind == "blanket"].tau.min(), 3), round(pm[pm.kind == "blanket"].tau.max(), 3)]
sc = pd.read_csv(FIG / "step_coverage_table.csv")
for lvl, mk in ((0.9, "o"), (0.95, "s"), (0.99, "^")):
    d = sc[sc.level == lvl]; b.plot(d.AR, d.t_bottom_worst_s, mk, ms=6, label=f"bottom reaches {int(lvl*100)} %, worst face")
d = sc[sc.level == 0.9]; b.plot(d.AR, d.rule_s, "k-", lw=1, label="rule 1.2 (1 + A$_{wall}$/A$_{wafer}$) N$_0$/(U c$_0$)")
b.plot(d.AR, d.t_flat_centre_s, "x", color="0.5", ms=6, label="flat top between the holes, 90 %")
b.set_xlabel("hole aspect ratio (100 nm holes, 180 nm pitch, half the wafer)"); b.set_ylabel("time (s)")
b.set_title("(b) 1 % TMA: 90, 95 and 99 % bottom coverage and the rule"); b.legend(); b.grid(alpha=0.3)
audit["fig1b_t99_over_rule"] = [round(v, 3) for v in (sc[sc.level == 0.99].t_bottom_worst_s.values / sc[sc.level == 0.99].rule_s.values)]
fig.tight_layout(); fig.savefig(OUT / "fig1_pulse_map.png"); plt.close(fig)

# ---------- Fig. 2: what does not move the patterned pulse ----------
ss = pd.read_csv(FIG / "sticking_sensitivity.csv"); doe = pd.read_csv(FIG / "doe_results.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.8))
for kind, mk, lab in (("blanket", "o", "blanket wafer"), ("patterned", "s", "patterned AR-100 wafer, hole bottom")):
    d = ss[ss.kind == kind].sort_values("s"); a.loglog(d.s, d.t90_worst, mk + "-", ms=6, label=lab, mfc="none" if kind == "blanket" else None)
d = ss[ss.kind == "patterned"].sort_values("s"); a.loglog(d.s, d.t90_flat_worst, "x--", color="0.5", ms=6, label="patterned wafer, flat top between holes")
a.axhline(ss[ss.kind == "patterned"].t_supply.iloc[0], color="k", ls=":", lw=1); a.text(5.5e-4, ss[ss.kind == "patterned"].t_supply.iloc[0] * 1.08, "supply limit, patterned", fontsize=7)
a.axhline(ss[ss.kind == "blanket"].t_supply.iloc[0], color="k", ls=":", lw=1); a.text(5.5e-4, ss[ss.kind == "blanket"].t_supply.iloc[0] * 1.08, "supply limit, blanket", fontsize=7)
for s0, lab in ((5e-4, "Arts 2019 (low)"), (2e-3, "Arts 2019 (high)"), (4e-3, "Vandalon 2016"), (5.72e-3, "Ylilammi 2018"), (8.5e-3, "this work, hydraulic"), (1.96e-2, "this work, corrected")):
    a.axvline(s0, color="0.85", lw=0.8, zorder=0)
a.set_xlabel("TMA sticking coefficient s"); a.set_ylabel("time to 90 % coverage, worst face (s)"); a.set_title("(a) sticking coefficient, 1 % TMA, 300 Pa, 300 sccm"); a.legend(loc="center left"); a.grid(True, which="both", alpha=0.3)
audit["fig2a_patterned_t90_worst_range"] = [float(d.t90_worst.min()), float(d.t90_worst.max())]
audit["fig2a_blanket_t90_worst_range"] = [float(ss[ss.kind == "blanket"].t90_worst.min()), float(ss[ss.kind == "blanket"].t90_worst.max())]
scat = b.scatter(doe.t_supply, doe.t90_worst, c=np.log10(doe.p), cmap="viridis", s=28, label="25 DOE runs (colour: log$_{10}$ pressure / Pa)")
t = np.logspace(np.log10(2), np.log10(60), 30); b.plot(t, 1.2 * t, "k-", lw=1, label="rule"); b.plot(t, t, "k:", lw=1, label="supply limit")
b.set_xscale("log"); b.set_yscale("log"); b.set_xlabel("supply-limit time (s)"); b.set_ylabel("time to 90 % bottom coverage, worst face (s)")
b.set_title("(b) DOE: 0.5-3 % TMA, 150-600 sccm, 100-1000 Pa, 10-25 mm gap"); plt.colorbar(scat, ax=b, pad=0.02).set_label("log$_{10}$ p (Pa)"); b.legend(loc="upper left"); b.grid(True, which="both", alpha=0.3)
for ax_ in (a, b):
    ax_.xaxis.set_minor_formatter(NullFormatter()); ax_.yaxis.set_minor_formatter(NullFormatter())
doe["tau"] = doe.t90_worst / doe.t_supply; audit["fig2b_doe_tau_min_max"] = [round(doe.tau.min(), 3), round(doe.tau.max(), 3)]
fig.tight_layout(); fig.savefig(OUT / "fig2_insensitivity.png"); plt.close(fig)

# ---------- Fig. 3: non-uniform layouts ----------
pl = pd.read_csv(FIG / "pattern_layout.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
a.plot(pl.radius_mm, pl.dense_centre_t90_bottom_s, "s-", ms=3, label="hole bottom"); a.plot(pl.radius_mm, pl.dense_centre_t90_flat_s, "x-", color="0.5", ms=3, label="flat surface")
a.axhline(14.05, color="k", ls="--", lw=1, label="uniform AR-100 wafer, bottom (14.05 s)"); a.axhline(5.1, color="tab:red", ls=":", lw=1, label="rule with the wafer-average ratio (5.1 s)"); a.axhline(14.4, color="tab:green", ls=":", lw=1, label="rule with the local ratio (14.4 s)")
a2 = a.twinx(); a2.plot(pl.radius_mm, pl.dense_centre_holeDensity_per_m2 / 1.78e13, color="tab:orange", lw=0.8, alpha=0.6); a2.set_ylabel("hole density / 1.78e13 m$^{-2}$", color="tab:orange"); a2.set_ylim(0, 1.6)
a.set_xlabel("radius (mm)"); a.set_ylabel("time to 90 % coverage (s)"); a.set_title("(a) dense AR-100 centre (r < 75 mm), bare edge"); a.legend(loc="upper right", bbox_to_anchor=(1.0, 0.93), framealpha=0.95); a.grid(alpha=0.3)
b.plot(pl.radius_mm, pl.ar_ramp_t90_bottom_s, "s-", ms=3, label="hole bottom"); b.plot(pl.radius_mm, pl.ar_ramp_t90_flat_s, "x-", color="0.5", ms=3, label="flat surface")
b.axhline(10.6, color="tab:red", ls=":", lw=1, label="rule, wafer-average ratio (10.6 s)"); b.plot(pl.radius_mm, 1.2 * (1 + pl.ar_ramp_AR * 0.559) * 11.99 / 56.92, color="tab:green", ls=":", lw=1, label="rule, local ratio")
b2 = b.twinx(); b2.plot(pl.radius_mm, pl.ar_ramp_AR, color="tab:orange", lw=0.8, alpha=0.6); b2.set_ylabel("local aspect ratio", color="tab:orange"); b2.set_ylim(0, 120)
b.set_xlabel("radius (mm)"); b.set_title("(b) aspect ratio 100 at the centre to 60 at the edge"); b.legend(loc="lower left", framealpha=0.95); b.grid(alpha=0.3)
audit["fig3_dense_centre_bottom_t90_centre_edge"] = [float(pl.dense_centre_t90_bottom_s.iloc[0]), float(pl.dense_centre_t90_bottom_s[pl.radius_mm == 74.5].iloc[0])]
audit["fig3_dense_centre_flat_t90_bare_edge_range"] = [float(pl.dense_centre_t90_flat_s[pl.radius_mm > 100].min()), float(pl.dense_centre_t90_flat_s[pl.radius_mm > 100].max())]
audit["fig3_ar_ramp_bottom_t90_centre_edge"] = [float(pl.ar_ramp_t90_bottom_s.iloc[0]), float(pl.ar_ramp_t90_bottom_s.iloc[-1])]
fig.tight_layout(); fig.savefig(OUT / "fig3_layouts.png"); plt.close(fig)

# ---------- Fig. 4: 3D chamber with a one-sided pump ----------
c3 = pd.read_csv(FIG / "chamber3d_t90.csv")
fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.6))
for r0, col in ((50, "tab:blue"), (100, "tab:green"), (140, "tab:red")):
    d = c3[(c3.r_mm - r0).abs() <= 2.5].copy(); d["bin"] = (d.azimuth_deg // 15) * 15 + 7.5
    g = d.groupby("bin").t90_bottom_s.mean(); a.plot(g.index, g.values, "o-", ms=3, color=col, label=f"r = {r0} mm ({len(d)} faces)")
    audit[f"fig4_p2p_pct_r{r0}"] = round(100 * (g.max() - g.min()) / g.mean(), 2)
a.set_xlabel("azimuth from the pump port (deg)"); a.set_ylabel("hole-bottom t90, 15-degree bin means (s)"); a.set_title("(a) azimuthal variation, 106k-cell chamber, 60-degree pump port"); a.legend(); a.grid(alpha=0.3)
sc_ = b.scatter(c3.y_mm, c3.z_mm, c=c3.t90_bottom_s, s=4, cmap="viridis", vmin=13.5, vmax=14.2); plt.colorbar(sc_, ax=b, pad=0.02).set_label("hole-bottom t90 (s)")
b.set_aspect("equal"); b.set_xlabel("y (mm), port at +x"); b.set_ylabel("z (mm)"); b.set_title("(b) wafer map, 6144 faces; worst 14.17 s vs 14.05 s wedge")
audit["fig4_worst_t90"] = float(c3.t90_bottom_s.max()); audit["fig4_min_t90"] = float(c3.t90_bottom_s.min())
fig.tight_layout(); fig.savefig(OUT / "fig4_chamber3d.png"); plt.close(fig)

# ---------- Fig. 5: three-method verification ----------
vm = pd.read_csv(FIG / "voxel_mc_crosscheck.csv")
holes = vm[(vm.kind == "hole") & (vm.update != "deposit") & (~vm.run.str.contains("plainfaces"))]
best = holes.sort_values("dt").groupby("AR").first().reset_index()   # smallest time step per AR
fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.4))
a.plot(best.AR, best.t90_machball, "o", ms=7, label="Machball (ballistic, Argonne), converged N$_z$")
a.plot(best.AR, best.t90, "s", ms=6, label="voxel Monte Carlo (this work), smallest dt")
a.plot(best.AR, best.t90_1d, "^", ms=6, label="1D Knudsen-Langmuir hole model")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("aspect ratio, 100 nm holes, 0.3 Pa TMA, 573 K, s = 5e-3"); a.set_ylabel("time to 90 % bottom coverage (s)"); a.set_title("(a) straight holes: three methods"); a.legend(); a.grid(True, which="both", alpha=0.3)
audit["fig5_voxel_over_machball_pct"] = [round(100 * (v / m - 1), 1) for v, m in zip(best.t90, best.t90_machball)]
audit["fig5_1d_over_machball_pct"] = [round(100 * (v / m - 1), 1) for v, m in zip(best.t90_1d, best.t90_machball)]
audit["fig5_1d_over_voxel_pct"] = [round(100 * (v / m - 1), 1) for v, m in zip(best.t90_1d, best.t90)]
sl = vm[vm.kind == "slit"].sort_values("dt").iloc[0]
# the corrected 1D slit model, run here: AR-50 slit, 100 nm gap, constant 0.3 Pa TMA, 573 K, s = 5e-3 (same conditions as the CSV rows)
from aldconf import Slit as _Slit, Precursor as _Pre, Process as _Proc, simulate_driven as _drv
def _slit_t90(correction):
    pr = _Pre(sticking=5e-3, site_density=4e18, mass_amu=72.09, gpc=0.0); n0 = 0.3 / (1.380649e-23 * 573.15)
    proc = _Proc(_Slit(100e-9, 5e-6, knudsen_correction=correction), pr, T=573.15, p_partial=0.3, t_dose=3.0, n_cells=200, narrowing=False)
    t = np.linspace(0, 3.0, 601); _, _, hist_t, hist = _drv(proc, t, np.full_like(t, n0), substeps_per_sample=5, record_depths=(1.0,))
    hist = np.asarray(hist).reshape(len(hist_t), -1)[:, 0]
    return float(np.interp(0.9, hist, hist_t))
t90_corr, t90_unc = _slit_t90("plates"), _slit_t90("none")
audit["fig5b_1d_slit_corrected_t90_run_here"] = round(t90_corr, 4); audit["fig5b_1d_slit_uncorrected_t90_run_here"] = round(t90_unc, 4)
vals = [sl.t90_machball, sl.t90, sl.t90_1d, t90_corr]
labels = ["Machball\nTrench", "voxel MC", "1D slit,\nD$_K$ = 2w$\\bar v$/3", "1D slit,\n(3/8) ln(L/w) factor"]
bars = b.bar(labels, vals, color=["tab:blue", "tab:orange", "tab:green", "tab:green"], alpha=[1, 1, 1, 0.5][0])
bars[3].set_alpha(0.5)
for x, v in zip(range(4), vals): b.text(x, v + 0.02, f"{v:.2f} s", ha="center", fontsize=8)
b.set_ylabel("time to 90 % bottom coverage (s)"); b.set_title("(b) slit, AR 50, 100 nm gap, same conditions"); b.set_ylim(0, 1.6); b.grid(axis="y", alpha=0.3)
audit["fig5b_slit_t90_machball_voxel_1d"] = [float(sl.t90_machball), float(sl.t90), float(sl.t90_1d)]
fig.tight_layout(); fig.savefig(OUT / "fig5_verification.png"); plt.close(fig)

# ---------- Fig. 6: PillarHall fit with the two slit transport models ----------
from aldconf import Slit, Precursor, Process, simulate, slit_knudsen_factor
from aldconf.data import load_pillarhall, default_path
df = load_pillarhall(default_path()); good = df[df.good & (df.depth_um >= 0) & (df.depth_um <= df.length_mm * 1000)]
sheets = [s for s in good.sheet.unique() if "L0.1mm" not in s]; grid = np.arange(0, 180, 2.0); prof = []
for sh in sheets:
    g = good[good.sheet == sh].sort_values("depth_um"); prof.append(np.interp(grid, g.depth_um, g.thickness_nm, left=np.nan, right=np.nan))
prof = np.array(prof); mean = np.nanmean(prof, 0); std = np.nanstd(prof, 0); plateau = np.nanmean(mean[(grid > 5) & (grid < 60)]); y = mean / plateau
GAP, L, CYCLES, TPULSE, T, KB = 500e-9, 250e-6, 500, 0.1, 573.15, 1.380649e-23; GPC = plateau * 1e-9 / CYCLES
fits = pd.read_csv(FIG / "pillarhall_fit_corrected.csv"); jc = fits[fits.fit == "joint_corrected"].iloc[0]; un = fits[fits.fit == "uncorrected_published"].iloc[0]
def model(s, p_Pa, correction):
    pr = Precursor(sticking=s, site_density=4e18, mass_amu=72.09, gpc=GPC)
    res = simulate(Process(Slit(GAP, L, knudsen_correction=correction), pr, T=T, p_partial=p_Pa, t_dose=TPULSE, cycles=CYCLES, n_cells=400, n_substeps=60, narrowing=True))
    th = res.thickness / res.thickness[:10].mean(); return np.interp(grid * 1e-6, res.x, th)
m_un = model(un.s, un.p_TMA_Pa, "none"); m_c = model(jc.s, jc.p_TMA_Pa, "plates")
ok = ~np.isnan(y); audit["fig6_rms_uncorrected"] = round(float(np.sqrt(np.mean((m_un[ok] - y[ok]) ** 2))), 4); audit["fig6_rms_corrected"] = round(float(np.sqrt(np.mean((m_c[ok] - y[ok]) ** 2))), 4)
audit["fig6_plateau_nm"] = round(float(plateau), 2); audit["fig6_n_sheets"] = len(sheets); audit["fig6_half_depth_AR"] = round(float(np.interp(0.5, y[::-1], grid[::-1]) * 2), 0)
audit["fig6_slit_factor_AR500"] = round(float(slit_knudsen_factor(GAP, L)), 3)
fig, a = plt.subplots(figsize=(6.2, 3.6))
a.errorbar(grid * 2, y, yerr=std / plateau, fmt=".", ms=3, color="k", alpha=0.6, label=f"Zenodo 3974438: Al$_2$O$_3$, 300 C, 500 x 0.1 s TMA, mean of {len(sheets)} profiles")
a.plot(grid * 2, m_un, "b--", lw=1.2, label=f"1D slit, D$_K$ = 2w$\\bar v$/3:  s = {un.s:.1e}, p$_{{TMA}}$ = {un.p_TMA_Pa:.0f} Pa, rms {audit['fig6_rms_uncorrected']:.3f}")
a.plot(grid * 2, m_c, "r-", lw=1.2, label=f"1D slit, (3/8) ln(L/w) correction:  s = {jc.s:.2e}, p$_{{TMA}}$ = {jc.p_TMA_Pa:.1f} Pa, rms {audit['fig6_rms_corrected']:.3f}")
a.set_xlabel("depth / gap  (aspect ratio; 500 nm gap)"); a.set_ylabel("thickness / plateau"); a.set_ylim(-0.05, 1.1); a.set_title("One measured profile, two slit transport models, the same residual"); a.legend(loc="lower left"); a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT / "fig6_pillarhall.png"); plt.close(fig)

json.dump(audit, open(OUT / "figure_audit.json", "w"), indent=1)
print(json.dumps(audit, indent=1)); print("saved", sorted(p.name for p in OUT.glob("fig*.png")))
