"""Figures 13 and 14 of draft_v4, drawn from result files only (nothing is recomputed here).
Figure 13: EDS line scans of the cycle series against the model with and without channel narrowing, from
docs/blindtest/v2/posthoc/eds_cycle_series_curves.json (written by eds_cycle_series.py; post hoc analysis).
Figure 14: diffusivity along a slit against kinetic theory, from figures/benchmark_slit_diffusion.csv and
figures/slit_diffusion_tma_n2.csv (written by benchmark_slit_diffusion.py).
Writes docs/paper/figures/fig13_eds_cycle_series.png and fig14_slit_benchmark.png.
  .venv/bin/python docs/paper/final_figures.py
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[2]; FIG = ROOT / "figures"; OUT = ROOT / "docs/paper/figures"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5, "figure.dpi": 150, "savefig.dpi": 200})

# ---- Figure 13: EDS cycle series ----
C = json.loads((ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").read_text())
fig, ax = plt.subplots(2, 3, figsize=(11, 5.6), sharey=True)
for c in C:
    j = (250, 500, 1000).index(c["cycles"]); i = ("left", "right").index(c["channel"]); a = ax[i, j]
    a.plot(c["z"], c["y"], ".", color="0.25", ms=3, label="EDS Al signal (Yim et al. 2020, Fig. S8)")
    a.plot(c["grid"], c["with_narrowing"], "C3-", lw=1.4, label="model with channel narrowing, s = 0.04")
    a.plot(c["grid"], c["without"], "C0--", lw=1.1, label="same s without narrowing (not refitted)")
    a.axvline(1.6 * c["z50"], color="0.7", lw=0.6, ls=":"); a.set_xlim(0, 1.7 * c["z50"]); a.set_ylim(-0.08, 1.25)
    a.set_title(f"({'abcdef'[3 * i + j]}) {c['cycles']} cycles, {c['channel']} channel, {100 * c['fill']:.0f} % of the gap filled")
    if i == 1: a.set_xlabel("distance from the entrance / gap")
    if j == 0: a.set_ylabel("relative thickness")
h, l = ax[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=8.5)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(OUT / "fig13_eds_cycle_series.png"); plt.close(fig)

# ---- Figure 14: slit diffusivity against kinetic theory ----
B = pd.read_csv(FIG / "benchmark_slit_diffusion.csv"); T = pd.read_csv(FIG / "slit_diffusion_tma_n2.csv"); N = B[B.pair == "Ne-Ar"]; X = B[B.pair == "He-Xe"]
fig, ax = plt.subplots(1, 3, figsize=(11, 3.5))
a = ax[0]
a.plot(N.width_2a, N.kinetic_theory, "k-", lw=1.2, label="kinetic theory (McCormack model)")
a.errorbar(N.width_2a, N.monte_carlo, yerr=N.monte_carlo * N.mc_rel_err, fmt="o", color="C3", ms=4, capsize=2, label="Monte Carlo, hard-sphere collisions")
R = N.dropna(subset=["monte_carlo_relaxation_collisions"]); a.plot(R.width_2a, R.monte_carlo_relaxation_collisions, "s", mfc="none", color="C2", ms=6, label="Monte Carlo, relaxation collisions")
a.plot(N.width_2a, N.hydraulic_bosanquet, "C0--", lw=1.1, label="hydraulic Knudsen + gas diffusion")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("slit width / mean free path"); a.set_ylabel("diffusivity along the slit (reference units)"); a.set_title("(a) Ne-Ar, heavy species"); a.legend(frameon=False)
a = ax[1]
a.axhline(1, color="k", lw=0.8); a.axhspan(0.95, 1.05, color="0.92", zorder=0)
a.errorbar(N.width_2a, N.mc_over_theory, yerr=N.mc_over_theory * N.mc_rel_err, fmt="o-", color="C3", ms=4, lw=0.8, capsize=2, label="Ne-Ar, hard spheres")
a.plot(R.width_2a, R.monte_carlo_relaxation_collisions / R.kinetic_theory, "s-", mfc="none", color="C2", ms=6, lw=0.8, label="Ne-Ar, relaxation collisions")
a.plot(N.width_2a, N.hydraulic_over_theory, "C0--", lw=1.1, label="Ne-Ar, hydraulic + gas diffusion")
a.plot(X.width_2a, X.mc_over_theory, "^:", color="0.5", ms=4, lw=0.8, label="He-Xe, hard spheres (not a tracer case)")
a.set_xscale("log"); a.set_ylim(0.4, 1.6); a.set_xlabel("slit width / mean free path"); a.set_ylabel("ratio to kinetic theory"); a.set_title("(b) ratio to the reference; band: 5 %"); a.legend(frameon=False, loc="upper left")
a = ax[2]
a.axhline(1, color="k", lw=0.8)
a.errorbar(T.pressure_Pa, T.mc_over_hydraulic, yerr=T.mc_over_hydraulic * T.mc_rel_err, fmt="o-", color="C3", ms=4, lw=0.8, capsize=2)
for p, r in zip(T.pressure_Pa, T.D_monte_carlo_m2_s / T.D_hydraulic_bosanquet_m2_s): a.annotate(f"{r:.2f}", (p, r), textcoords="offset points", xytext=(4, 5), fontsize=7)
a.set_xscale("log"); a.set_xlim(18, 90000); a.set_ylim(0.9, 2.6); a.set_xlabel("N$_2$ pressure / Pa"); a.set_ylabel("Monte Carlo / hydraulic + gas diffusion"); a.set_title("(c) TMA in N$_2$, 500 nm gap")
fig.tight_layout(); fig.savefig(OUT / "fig14_slit_benchmark.png"); plt.close(fig)
print("wrote fig13_eds_cycle_series.png, fig14_slit_benchmark.png")

# ---- Figure: sticking coefficient fitted profile by profile, without and with channel narrowing ----
P = pd.read_csv(FIG / "profile_fits_cut0.3.csv"); fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.7), sharey=True)
for a, fam, title in zip(ax, ("hydraulic", "gas", "ballistic"), ("(a) hydraulic model", "(b) carrier-gas model, N$_2$ at 300 Pa", "(c) gas-free model")):
    g = P[P.transport == fam]; top = g.s_grid_max.max()
    for nar, col, lab in ((False, "C0", "without narrowing"), (True, "C3", "with narrowing")):
        mk = "s" if nar else "o"
        for meth in ("EDS", "reflectometry"):
            q = g[(g.narrowing == nar) & (g.method == meth)]; x = q.fill + (0.006 if nar else -0.006)
            lo = np.clip(q.s_fit - q.s_p16, 0, None); hi = np.clip(q.s_p84 - q.s_fit, 0, None); open_top = q.s_p84 >= 0.99 * q.s_grid_max
            a.errorbar(x, q.s_fit, yerr=[lo, hi], fmt=mk, color=col, ms=4.5, mfc=col if meth == "EDS" else "none", lw=0.8, capsize=2, label=f"{meth}, {lab}")
            a.plot(x[open_top], np.full(open_top.sum(), 0.2 * 1.25), "^", color=col, ms=4)
        q = g[g.narrowing == nar]; c = np.polyfit(q.fill, np.log(q.s_fit), 1); ff = np.linspace(0.04, 0.41, 20); a.plot(ff, np.exp(np.polyval(c, ff)), "-" if not nar else "--", color=col, lw=0.9, alpha=0.7)
    a.set_yscale("log"); a.set_ylim(1e-3, 0.32); a.set_xlabel("fraction of the gap filled, f = 2t/H"); a.set_title(title); a.set_yticks([0.001, 0.003, 0.01, 0.03, 0.1, 0.3]); a.set_yticklabels(["0.001", "0.003", "0.01", "0.03", "0.1", "0.3"]); a.minorticks_off()
ax[0].set_ylabel("fitted sticking coefficient"); h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(OUT / "fig_profile_fits.png"); plt.close(fig)
print("wrote fig_profile_fits.png")

# ---- Figure: how strongly narrowing acts (the exponent as an empirical dial) ----
XS = pd.read_csv(FIG / "exponent_scan_summary.csv"); XS = XS[XS.fronts == "all 16 fronts"]; XB = pd.read_csv(FIG / "exponent_scan_bounds.csv"); XP = pd.read_csv(FIG / "exponent_scan_pair.csv")
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0)); a = ax[0]
for (fam, win), mk, col, dx in ((("hydraulic", 0.3), "o", "C0", -0.006), (("gas", 0.3), "s", "C3", 0.006), (("hydraulic", 0.1), "o", "C0", -0.014), (("gas", 0.1), "s", "C3", 0.014)):
    g = XS[(XS.transport == fam) & (XS.window == win)].sort_values("n")
    a.errorbar(g.n + dx, g.dlns_df, yerr=g.se, fmt=mk + ("-" if win == 0.3 else ":"), color=col, ms=4.5, mfc=col if win == 0.3 else "none", lw=0.9, capsize=2, label=f"{'hydraulic' if fam == 'hydraulic' else 'carrier-gas'} model, window to {win}")
a.axhline(0, color="k", lw=0.8); steep = XB[XB.n_upper_bound < 0.8].n_upper_bound; a.axvspan(-0.05, steep.max(), color="0.9", zorder=0)
for b in steep: a.plot([b, b], [-10.6, -9.6], color="k", lw=1.2)
pr = XP[XP.reading == "floor subtracted"]; ok = pr[pr.rms <= 1.5 * pr.rms.min()]; best = pr.loc[pr.rms.idxmin()]
a.plot([ok.n.min(), ok.n.max()], [-8.8, -8.8], color="C2", lw=2.5); a.plot(best.n, -8.8, "D", color="C2", ms=5)
a.text(0.02, -10.1, "upper bounds from the three steep fronts", fontsize=7.5, va="center"); a.text(0.02, -8.8, "second-reactor pair: best n and range", fontsize=7.5, va="center", color="C2")
for v, lab in ((0.92, "0.92"), (1.0, "1")): a.axvline(v, color="0.5", lw=0.7, ls="--")
a.text(0.96, 5.6, "theory", fontsize=7.5, ha="center", color="0.4")
a.set_xlim(-0.05, 1.06); a.set_ylim(-11, 6.5); a.set_xlabel("narrowing exponent n (0 = no narrowing)"); a.set_ylabel("slope of ln s against filled fraction, 16 fronts"); a.set_title("(a) main series refitted for each n"); a.legend(frameon=False, fontsize=7, loc="upper left")
a = ax[1]; ff = np.linspace(0.03, 0.5, 100)
for n_, ls, lab in ((1.0, "-", "steepest front, n = 1"), (0.92, "--", "n = 0.92"), (0.7, ":", "n = 0.7")): a.plot(ff, (1 - ff / 2) / (n_ * ff), "k" + ls, lw=1, label=lab)
T1 = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv"); T1 = T1[(T1.H_nm == 500) & (T1.magnification == 50)]; a.plot(2 * T1.S_IIb_nm / T1.H_nm, T1.slope_dtheta_dxtilde * T1.x50_tilde, "o", color="0.55", ms=4, mfc="none", label="reflectometry, tabulated slopes")
E_ = pd.read_csv(ROOT / "docs/blindtest/v2/posthoc/eds_slopes.csv"); a.plot(E_.fill, E_.slope_eds * E_.z50_gaps, "*", color="C0", ms=8, label="EDS line scans")
BQ = pd.read_csv(FIG / "profile_fits_beneq_slopes.csv"); BQ = BQ[BQ.n == 1.0]
for rd, mfc in (("floor subtracted", "C3"), ("as read", "none")):
    g = BQ[BQ.reading == rd]; a.plot(g.fill, g.slope_x_z50, "^", color="C3", ms=6, mfc=mfc, label=f"second reactor, {rd}")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(0.04, 0.5); a.set_ylim(1.5, 25); a.set_xticks([0.05, 0.1, 0.2, 0.4]); a.set_xticklabels(["0.05", "0.1", "0.2", "0.4"]); a.set_yticks([2, 5, 10, 20]); a.set_yticklabels(["2", "5", "10", "20"]); a.minorticks_off()
a.set_xlabel("fraction of the gap filled, f"); a.set_ylabel("front slope x half-thickness depth"); a.set_title("(b) measured fronts against the steepest front"); a.legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout(); fig.savefig(OUT / "fig_exponent_scan.png"); plt.close(fig); print("wrote fig_exponent_scan.png")
