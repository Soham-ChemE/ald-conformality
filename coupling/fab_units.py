"""Fab units from the existing chamber runs: step coverage and film thickness against pulse time, purge times
(hole, chamber, desorption) and the ALD cycle time / wafers per hour.

Reads only time directories that already exist (no OpenFOAM is run here) and the Python feature solver's
transport module. Writes:
    figures/step_coverage_vs_pulse.csv, figures/step_coverage_table.csv, figures/step_coverage_vs_pulse.png
    figures/purge_times.csv
    figures/cycle_time_wph.csv, figures/wph_vs_wall_ratio.csv, figures/cycle_time_wph.png
Run:  PYTHONPATH=src .venv/bin/python coupling/fab_units.py

Definitions used throughout
- step coverage SC(t) = thetaBottom / thetaFlat at a wafer face, evaluated at time t of a pulse that is still on.
  In the irreversible Langmuir model the coverage at time t of a long pulse is the coverage a pulse of length t
  leaves behind, up to the small amount adsorbed from the gas left in the chamber after the valve closes (the
  surface density over the patterned wafer is 0.019 c0 during the pulse, so this is negligible).
- thickness after N cycles of a given pulse = N x GPC x theta, GPC = 0.103 nm per cycle from the PillarHall record
  in this repository (plateau 51.5 nm after 500 cycles, coupling/fit_pillarhall.py; Puurunen 2005 reviews ~0.1 nm
  per cycle for TMA/H2O as the literature context).
- purge to a level = time after the pulse ends for the outlet precursor signal to fall to that fraction of its
  value at the instant the pulse ends (same metric as coupling/postprocess_desorption.py).
"""
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "coupling"), str(ROOT / "reactor"), str(ROOT / "src")]
from patterned_lib import load_patterned, _block, _field
from postprocess_langmuir import load as load_flat
from aldconf.transport import knudsen_diffusivity, mean_speed
FIG = ROOT / "figures"

# ---- constants shared with the rest of the repository ------------------------------------------------------
GPC_NM = 51.5 / 500                      # 0.103 nm per cycle, PillarHall plateau (paper/assumptions.md night 1 part 10)
N0, U, C0_1PCT = 4e18, 0.050, 3.79e20    # as in coupling/temporal_vs_spatial.py
Q_M3S, V_M3 = 3.54e-3, np.pi * 0.2 ** 2 * 0.015   # 300 sccm at 300 Pa, 573 K; chamber r 200 mm x 15 mm gap
A_WAFER = np.pi * 0.15 ** 2
T_K, M_TMA = 573.15, 72.09
COL = {20: "tab:blue", 50: "tab:orange", 100: "tab:green"}; LS = {20: "-", 50: "--", 100: "-."}

def first_time(arr, t, lvl):
    k = np.argmax(arr >= lvl); return float(t[k]) if arr.max() >= lvl else np.nan

def outlet(case):
    """Mean outlet-patch value against time (copied from coupling/postprocess_desorption.py; not imported because
    that script draws its figure on import)."""
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times]); o = np.zeros(len(t))
    for i, (_, d) in enumerate(times): o[i] = _field(_block((d / "T").read_text(), "outlet"), "value", 30).mean()
    return t, o

def pulse_of(case):
    return float(re.search(r"\(\(0 1\) \(([0-9.]+) 1\)", (case / "0/T").read_text()).group(1))

# ============================================================================================================
# 1. Step coverage and thickness against pulse time
# ============================================================================================================
print("== 1. step coverage and thickness")
cases1 = {20: ROOT / "reactor/map/patterned_AR20_x0.01", 50: ROOT / "reactor/map/patterned_AR50_x0.01",
          100: ROOT / "reactor/map/patterned_AR100_x0.01", "100 long": ROOT / "reactor/showerhead_pulse_patterned_AR100_long"}
pmap = pd.read_csv(FIG / "pulse_time_map.csv")
curves, rows1 = [], []
for key, case in cases1.items():
    AR = 100 if key == "100 long" else key
    t, f, p, pulse = load_patterned(case)
    on = t <= pulse + 1e-9                                 # only while the pulse is on
    t, fb, ff = t[on], f["thetaBottom"][on], f["thetaFlat"][on]
    sc = fb / np.maximum(ff, 1e-12)
    sc_worst = np.nanmin(sc[:, :149], axis=1)              # worst face (face 149 is the edge-spike face, excluded as in the map)
    row = pmap[(pmap.kind == "patterned") & (pmap.AR == AR) & (pmap.x == 0.01)].iloc[0]
    rule = 1.2 * row.t_supply                              # the map's t_supply column is already (1 + A_wall/A_wafer) N0/(U c0)
    for lvl in (0.90, 0.95, 0.99):
        rows1.append(dict(case=case.relative_to(ROOT).as_posix(), AR=AR, wall_ratio=round(row.wall_ratio, 2), pulse_s=pulse, write_dt_s=round(float(t[1] - t[0]), 4),
                          level=lvl, t_bottom_centre_s=first_time(fb[:, 0], t, lvl), t_bottom_worst_s=np.nanmax([first_time(fb[:, j], t, lvl) for j in range(149)]),
                          t_stepcov_centre_s=first_time(sc[:, 0], t, lvl), t_stepcov_worst_s=first_time(sc_worst, t, lvl),
                          t_flat_centre_s=first_time(ff[:, 0], t, lvl), rule_s=rule))
    curves.append(pd.DataFrame(dict(case=case.relative_to(ROOT).as_posix(), AR=AR, t_s=t, theta_bottom_centre=fb[:, 0], theta_flat_centre=ff[:, 0],
                                    step_coverage_centre=sc[:, 0], step_coverage_worst=sc_worst,
                                    thickness_bottom_50cyc_nm=50 * GPC_NM * fb[:, 0], thickness_top_50cyc_nm=50 * GPC_NM * ff[:, 0])))
    print(f"  {case.name}: pulse {pulse} s, write {t[1]-t[0]:.3f} s, rule {rule:.2f} s; SC90/95/99 centre "
          f"{first_time(sc[:,0],t,0.9):.3f} / {first_time(sc[:,0],t,0.95):.3f} / {first_time(sc[:,0],t,0.99):.3f} s; flat t90 {first_time(ff[:,0],t,0.9):.2f} s")
tab1 = pd.DataFrame(rows1); tab1.to_csv(FIG / "step_coverage_table.csv", index=False)
pd.concat(curves).to_csv(FIG / "step_coverage_vs_pulse.csv", index=False)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for cv in curves:
    AR = int(cv.AR.iloc[0]); long = "long" in cv.case.iloc[0]
    lab = f"AR {AR}, {11.2 if AR == 20 else 28.0 if AR == 50 else 55.9:.0f}x wall area" + (" (20 s run)" if long else "")
    ax[0].plot(cv.t_s, 100 * cv.step_coverage_centre, LS[AR], color=COL[AR], lw=2 if not long else 1, alpha=1 if not long else 0.6, label=lab)
    ax[1].plot(cv.t_s, cv.thickness_bottom_50cyc_nm, LS[AR], color=COL[AR], lw=2 if not long else 1, alpha=1 if not long else 0.6, label=lab)
    if not long:
        rule = tab1[(tab1.AR == AR)].rule_s.iloc[0]
        for a in ax: a.axvline(rule, color=COL[AR], lw=0.8, ls=":", alpha=0.9)
        ax[0].text(rule, 4, f" rule {rule:.1f} s", color=COL[AR], fontsize=7, rotation=90, va="bottom", ha="right")
ax[1].plot(curves[0].t_s, curves[0].thickness_top_50cyc_nm, "k:", lw=1, label="flat top of the array, any AR (50 cycles)")
for a in ax:
    a.set_xscale("log"); a.set_xlim(0.1, 25); a.set_xlabel("pulse time (s)"); a.grid(True, which="both", alpha=0.3)
ax[0].axhline(90, color="0.6", lw=0.8); ax[0].axhline(95, color="0.6", lw=0.8); ax[0].axhline(99, color="0.6", lw=0.8)
ax[0].text(0.105, 90.5, "90 %", fontsize=7, color="0.4"); ax[0].text(0.105, 95.5, "95 %", fontsize=7, color="0.4"); ax[0].text(0.105, 99.5, "99 %", fontsize=7, color="0.4")
ax[0].set_ylim(0, 105); ax[0].set_ylabel("step coverage, hole bottom / flat top (%)")
ax[0].set_title("Step coverage against TMA pulse (1 % TMA, 300 sccm, 300 Pa, wafer centre); dotted: rule", fontsize=9); ax[0].legend(fontsize=7, loc="upper left")
ax[1].set_ylabel("Al2O3 at the hole bottom after 50 cycles (nm)"); ax[1].set_ylim(0, 5.6)
ax[1].axhline(50 * GPC_NM, color="0.6", lw=0.8); ax[1].text(24, 50 * GPC_NM + 0.06, f"saturated: {50*GPC_NM:.2f} nm (0.103 nm per cycle)", fontsize=7, color="0.4", ha="right")
ax[1].set_title("Thickness at the hole bottom after 50 cycles of that pulse (GPC 0.103 nm, PillarHall record)", fontsize=9); ax[1].legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(FIG / "step_coverage_vs_pulse.png", dpi=150); print("  saved figures/step_coverage_vs_pulse.png")

# ============================================================================================================
# 2. Purge: (a) the hole, (b) the chamber, (c) desorption
# ============================================================================================================
print("== 2. purge")
purge_rows = []
# (a) the hole: Knudsen diffusion out of a 100 nm x 10 um hole with the mouth density set to zero, walls saturated
d, L = 100e-9, 10e-6
DK = float(knudsen_diffusivity(d, T_K, M_TMA)); vbar = float(mean_speed(T_K, M_TMA))
t_scale = L ** 2 / DK
tau1 = 4 * L ** 2 / (np.pi ** 2 * DK)             # slowest eigenmode, Dirichlet mouth / Neumann bottom
t_bottom_1pct_analytic = tau1 * np.log(4 / np.pi / 0.01)     # first mode: n_bottom = (4/pi) exp(-t/tau1)
# numerical: implicit (backward Euler) 1D transient, 200 cells, n = 1 at t = 0, n(0) = 0, dn/dx(L) = 0
from scipy.linalg import solve_banded
Nc = 200; dx = L / Nc; dt = 2e-8; n = np.ones(Nc); a = DK * dt / dx ** 2
diag = np.full(Nc, 1 + 2 * a); diag[-1] = 1 + a; diag[0] = 1 + 3 * a     # ghost at the mouth: n_ghost = -n_0
ab = np.zeros((3, Nc)); ab[0, 1:] = -a; ab[1] = diag; ab[2, :-1] = -a
tt, nb, nm = [0.0], [1.0], [1.0]
while tt[-1] < 60e-6:
    n = solve_banded((1, 1), ab, n); tt.append(tt[-1] + dt); nb.append(n[-1]); nm.append(n.mean())
tt, nb, nm = map(np.array, (tt, nb, nm))
t_bottom_e = tt[np.argmax(nb < np.exp(-1))]; t_bottom_1pct = tt[np.argmax(nb < 0.01)]; t_mean_1pct = tt[np.argmax(nm < 0.01)]
print(f"  (a) hole: vbar {vbar:.0f} m/s, D_K {DK:.3e} m2/s, L^2/D_K {t_scale*1e6:.2f} us; eigenmode tau1 {tau1*1e6:.2f} us, bottom to 1 % analytic {t_bottom_1pct_analytic*1e6:.2f} us;"
      f" numeric bottom to 1/e {t_bottom_e*1e6:.2f} us, to 1 % {t_bottom_1pct*1e6:.2f} us, hole mean to 1 % {t_mean_1pct*1e6:.2f} us")
purge_rows += [dict(part="hole", case="analytic", quantity="L^2/D_K, 100 nm x 10 um, TMA 573 K", value_s=t_scale, note=f"D_K = d vbar/3 = {DK:.3e} m2/s, vbar {vbar:.0f} m/s"),
               dict(part="hole", case="analytic", quantity="slowest eigenmode time constant 4 L^2/(pi^2 D_K)", value_s=tau1, note="Dirichlet mouth, Neumann bottom"),
               dict(part="hole", case="analytic", quantity="bottom density to 1 %, first mode", value_s=t_bottom_1pct_analytic, note="(4/pi) exp(-t/tau1) = 0.01"),
               dict(part="hole", case="numeric", quantity="bottom density to 1/e", value_s=t_bottom_e, note="implicit 1D transient, 200 cells, dt 20 ns"),
               dict(part="hole", case="numeric", quantity="bottom density to 1 %", value_s=t_bottom_1pct, note="same"),
               dict(part="hole", case="numeric", quantity="hole-mean density to 1 %", value_s=t_mean_1pct, note="same")]

# (b) the chamber: outlet signal after the pulse ends
def purge_times(case, label, extrap_window=(None, None)):
    pulse = pulse_of(case); to, o = outlet(case)
    i_end = np.argmin(abs(to - pulse)); o_end = o[i_end]; o_peak = o.max(); t_peak = to[np.argmax(o)]
    after = to > pulse + 1e-6
    # exponential decay constant of the tail: log-linear fit over the last fitted window (default: last 40 % of the purge)
    t_lo = extrap_window[0] if extrap_window[0] else pulse + 0.6 * (to[-1] - pulse); t_hi = extrap_window[1] if extrap_window[1] else to[-1]
    w = (to >= t_lo) & (to <= t_hi) & (o > 0)
    k_decay, lnA = np.polyfit(to[w] - pulse, np.log(o[w]), 1); k_decay = -k_decay
    out = dict(case=label, pulse_s=pulse, run_end_s=to[-1], outlet_at_pulse_end=o_end, outlet_peak=o_peak, t_peak_after_pulse_s=t_peak - pulse,
               tail_decay_rate_1_s=k_decay, tail_time_constant_s=1 / k_decay)
    for frac in (0.1, 0.01, 0.001):
        below = np.where(after & (o < frac * o_end))[0]
        if len(below):
            out[f"t_{frac:g}_s"] = to[below[0]] - pulse; out[f"t_{frac:g}_how"] = "run"
        else:
            out[f"t_{frac:g}_s"] = (lnA - np.log(frac * o_end)) / k_decay; out[f"t_{frac:g}_how"] = f"extrapolated (exp tail fitted {t_lo-pulse:.2f}-{t_hi-pulse:.2f} s after the pulse)"
    return out, (to, o, pulse)
chamber = {}
for label, case in (("TMA patterned AR-100, 20 s pulse", ROOT / "reactor/showerhead_pulse_patterned_AR100_long"),
                    ("TMA blanket, 0.5 s pulse", ROOT / "reactor/showerhead_pulse_langmuir_p0.5"),
                    ("TMA blanket, 1.26 s pulse (map)", ROOT / "reactor/map/blanket_x0.01"),
                    ("H2O blanket kdes 0, 5 s pulse", ROOT / "reactor/showerhead_h2o_blanket_kdes0")):
    r, tr = purge_times(case, label); chamber[label] = (r, tr)
    print(f"  (b) {label}: outlet at pulse end {r['outlet_at_pulse_end']:.3f} c0 (peak {r['outlet_peak']:.3f} at +{r['t_peak_after_pulse_s']:.2f} s); "
          f"to 10 % {r['t_0.1_s']:.2f} s [{r['t_0.1_how'][:4]}], 1 % {r['t_0.01_s']:.2f} s [{r['t_0.01_how'][:4]}], 0.1 % {r['t_0.001_s']:.2f} s [{r['t_0.001_how'][:4]}]; "
          f"tail decay {r['tail_decay_rate_1_s']:.2f} 1/s (tau {r['tail_time_constant_s']:.2f} s)")
    for frac in (0.1, 0.01, 0.001):
        purge_rows.append(dict(part="chamber", case=label, quantity=f"outlet to {100*frac:g} % of pulse-end value", value_s=r[f"t_{frac:g}_s"], note=r[f"t_{frac:g}_how"]))
    purge_rows.append(dict(part="chamber", case=label, quantity="outlet tail time constant", value_s=r["tail_time_constant_s"], note=f"decay rate {r['tail_decay_rate_1_s']:.2f} 1/s"))
purge_rows += [dict(part="chamber", case="reference", quantity="gap transit 0.015 m / 0.050 m/s", value_s=0.015 / 0.050, note="showerhead face speed"),
               dict(part="chamber", case="reference", quantity="chamber volume / volumetric flow", value_s=V_M3 / Q_M3S, note=f"V {V_M3:.3e} m3 (r 200 mm, 15 mm gap), Q {Q_M3S:.2e} m3/s")]
T_PURGE_TMA = chamber["TMA patterned AR-100, 20 s pulse"][0]["t_0.01_s"]

# (c) desorption: water outlet to 1 % for each kdes; where the 3 s run does not reach it, extrapolate with the
# simulated coverage decay (outlet in the desorption-limited tail = kdes theta N0 A_wafer / (Q c0))
print("  (c) desorption")
QC0 = Q_M3S * C0_1PCT
des = {}
for k in (0, 0.1, 1, 10):
    case = ROOT / f"reactor/showerhead_h2o_blanket_kdes{k:g}"
    t, th, face, s, vb, c0, N0c, pulse = load_flat(case); to, o = outlet(case)
    i_end = np.argmin(abs(t - pulse)); th_end = th[i_end, 0]; th_final = th[-1, 0]
    o_end = o[np.argmin(abs(to - pulse))]; after = to > pulse + 0.02
    below = np.where(after & (o < 0.01 * max(o_end, 1e-9)))[0]
    o_tail_model = k * th_final * N0c * A_WAFER / QC0 if k else 0.0
    if len(below):
        t_purge, how = to[below[0]] - pulse, "run"
    else:
        # coverage decay rate over the last 2 s of the purge (includes readsorption), then theta(t) = theta_end exp(-k_eff (t - t_end))
        w = t >= t[-1] - 2.0
        k_eff = -np.polyfit(t[w], np.log(th[w, 0]), 1)[0]
        th_target = 0.01 * o_end * QC0 / (k * N0c * A_WAFER)         # coverage at which the desorption-fed outlet equals 1 % of the pulse-end level
        t_purge = (t[-1] - pulse) + np.log(th_final / th_target) / k_eff
        how = f"extrapolated: desorption-limited tail kdes theta N0 A/(Q c0) with theta decaying at {k_eff:.3f} 1/s (fitted {t[-1]-pulse-2:.0f}-{t[-1]-pulse:.0f} s after the pulse)"
    des[k] = dict(kdes=k, t_purge=t_purge, how=how, th_end=th_end, th_final=th_final, o_end=o_end, o_last=o[-1], o_tail_model=o_tail_model)
    print(f"     kdes {k:4g} 1/s: coverage at pulse end {th_end:.3f}, 3 s after {th_final:.3f}; outlet 3 s after {o[-1]:.4f} c0 (desorption-tail model {o_tail_model:.4f}); purge to 1 %: {t_purge:.2f} s [{how[:12]}]")
    purge_rows.append(dict(part="desorption", case=f"H2O blanket kdes {k:g} 1/s", quantity="outlet to 1 % of pulse-end value", value_s=t_purge, note=how))
T_PURGE_H2O = des[0.1]["t_purge"]
pd.DataFrame(purge_rows).to_csv(FIG / "purge_times.csv", index=False); print("  saved figures/purge_times.csv")

# ============================================================================================================
# 3. Cycle time and throughput
# ============================================================================================================
print("== 3. cycle time and wafers per hour")
T_TMA_BLANKET, T_H2O_BLANKET = 0.345, 1.74           # reactor/showerhead_pulse_langmuir_p0.5 (0.35 s rounded in the paper), reactor/make_h2o_cases.py
T_TMA_PATT, T_H2O_PATT = 14.05, 13.4                 # showerhead_pulse_patterned_AR100_long, showerhead_h2o_patterned_AR100
T_TMA_OPT = 9.83                                     # figures/doe_gp_report.txt, same-supply optimum (2.51 %, 150 sccm, 631 Pa, 10 mm)
T_H2O_OPT = T_H2O_PATT * T_TMA_OPT / T_TMA_PATT      # water scaled by the same factor (not in the DOE)
recipes = [("blanket wafer, baseline", T_TMA_BLANKET, T_H2O_BLANKET),
           ("patterned AR-100, baseline", T_TMA_PATT, T_H2O_PATT),
           ("patterned AR-100, same-supply optimum", T_TMA_OPT, T_H2O_OPT)]
OVERHEAD = 60.0
rows3 = []
for name, ttma, th2o in recipes:
    for purge_label, tp_tma, tp_h2o in (("chamber 1 % (TMA), kdes 0.1 (H2O)", T_PURGE_TMA, T_PURGE_H2O), ("gas clearing only (kdes 0)", T_PURGE_TMA, des[0]["t_purge"])):
        cyc = ttma + tp_tma + th2o + tp_h2o
        for thick in (5.0, 10.0):
            n = int(np.ceil(thick / GPC_NM)); tw = n * cyc
            rows3.append(dict(recipe=name, purge_basis=purge_label, t_TMA_s=ttma, t_purge1_s=tp_tma, t_H2O_s=th2o, t_purge2_s=tp_h2o, cycle_s=cyc,
                              thickness_nm=thick, cycles=n, min_per_wafer=tw / 60, wph_no_overhead=3600 / tw, wph_60s_overhead=3600 / (tw + OVERHEAD)))
tab3 = pd.DataFrame(rows3); tab3.to_csv(FIG / "cycle_time_wph.csv", index=False)
for _, r in tab3[tab3.purge_basis.str.startswith("chamber")].iterrows():
    print(f"  {r.recipe:40s} {r.thickness_nm:4.0f} nm: cycle {r.cycle_s:5.2f} s ({r.t_TMA_s:.2f} + {r.t_purge1_s:.2f} + {r.t_H2O_s:.2f} + {r.t_purge2_s:.2f}), {r.cycles} cycles, "
          f"{r.min_per_wafer:5.1f} min/wafer, {r.wph_no_overhead:5.1f} WPH ({r.wph_60s_overhead:5.1f} with 60 s handling)")

# the rule as a throughput formula: WPH against wall-area ratio, 1 % and 3 % TMA, 5 nm film, purges from part 2
wall = np.logspace(0, 2.2, 60); phi_req = 1.2 * (1 + wall) * N0; J1 = U * C0_1PCT
N5 = int(np.ceil(5.0 / GPC_NM))
wrows = []
for x, J in ((0.01, J1), (0.03, 3 * J1)):
    t_half = phi_req / J                                 # same supply flux for both half-cycles, as in coupling/temporal_vs_spatial.py
    cyc = 2 * t_half + T_PURGE_TMA + T_PURGE_H2O
    cyc0 = 2 * t_half + T_PURGE_TMA + des[0]["t_purge"]
    for w, c, c0_ in zip(wall, cyc, cyc0):
        wrows.append(dict(x=x, wall_ratio=w, cycle_s=c, wph_no_overhead=3600 / (N5 * c), wph_60s_overhead=3600 / (N5 * c + OVERHEAD), cycle_kdes0_s=c0_, wph_kdes0_60s_overhead=3600 / (N5 * c0_ + OVERHEAD)))
wtab = pd.DataFrame(wrows); wtab.to_csv(FIG / "wph_vs_wall_ratio.csv", index=False)
for x in (0.01, 0.03):
    for w in (1, 11.2, 28, 55.9, 112):
        sub = wtab[wtab.x == x]; i = np.argmin(abs(sub.wall_ratio.values - w)); r = sub.iloc[i]
        print(f"  rule, {100*x:g} % TMA, wall {w:5.1f}x: cycle {r.cycle_s:5.1f} s, {r.wph_no_overhead:5.1f} WPH ({r.wph_60s_overhead:5.1f} with 60 s handling; kdes 0 purge: {r.wph_kdes0_60s_overhead:5.1f})")

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=(1, 1.15)))
seg_cols = ("tab:blue", "0.75", "tab:orange", "0.55"); seg_names = ("TMA pulse", "purge (chamber to 1 %)", "H2O pulse", "purge (kdes 0.1 1/s)")
base = tab3[(tab3.purge_basis.str.startswith("chamber")) & (tab3.thickness_nm == 5)].reset_index(drop=True)
xpos = np.arange(len(base)); bottom = np.zeros(len(base))
for col, nm, key in zip(seg_cols, seg_names, ("t_TMA_s", "t_purge1_s", "t_H2O_s", "t_purge2_s")):
    vals = base[key].values; ax[0].bar(xpos, vals, 0.6, bottom=bottom, color=col, edgecolor="white", linewidth=1.5, label=nm)
    for xx, b, v in zip(xpos, bottom, vals):
        if v > 1.2: ax[0].text(xx, b + v / 2, f"{v:.1f}", ha="center", va="center", fontsize=7, color="black" if col != "tab:blue" else "white")
    bottom += vals
for xx, r in base.iterrows():
    ax[0].text(xx, r.cycle_s + 0.5, f"{r.cycle_s:.1f} s per cycle\n{r.wph_60s_overhead:.1f} WPH at 5 nm", ha="center", va="bottom", fontsize=7)
ax[0].set_xticks(xpos); ax[0].set_xticklabels(["blanket\nbaseline", "AR 100\nbaseline", "AR 100\nsame-supply optimum"], fontsize=8)
ax[0].set_ylabel("time per ALD cycle (s)"); ax[0].set_ylim(0, 1.25 * base.cycle_s.max()); ax[0].legend(fontsize=7, loc="upper left")
ax[0].set_title("Al2O3 cycle: 1 % precursors, 300 sccm, 300 Pa (optimum 2.5 %, 150 sccm, 631 Pa)", fontsize=9); ax[0].grid(True, axis="y", alpha=0.3)
for x, c, lsx in ((0.01, "k", "-"), (0.03, "k", "--")):
    sub = wtab[wtab.x == x]
    ax[1].plot(sub.wall_ratio, sub.wph_60s_overhead, lsx, color=c, lw=2, label=f"{100*x:g} % TMA and H2O, 60 s handling per wafer")
    ax[1].plot(sub.wall_ratio, sub.wph_no_overhead, lsx, color="0.55", lw=1, label=f"{100*x:g} %, no handling overhead")
for w, lab in ((11.2, "AR 20"), (28.0, "AR 50"), (55.9, "AR 100")):
    ax[1].axvline(w, color="0.75", lw=0.8); ax[1].text(w, 10.5, lab, fontsize=7, ha="center", color="0.3")
for _, r in base.iterrows():
    if "patterned" in r.recipe and "baseline" in r.recipe: ax[1].plot([55.9], [r.wph_60s_overhead], "o", color="tab:green", ms=7, label="chamber run, AR-100 baseline (14.05 + 13.4 s)")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("patterned wall area / wafer area"); ax[1].set_ylabel("wafers per hour, 5 nm Al2O3 (49 cycles)")
ax[1].set_title(f"Rule as throughput: {N5} cycles x [2 x 1.2 (1 + A_wall/A_wafer) N0/(U c0) + {T_PURGE_TMA:.1f} s + {T_PURGE_H2O:.1f} s purges]", fontsize=8)
ax[1].legend(fontsize=7, loc="lower left"); ax[1].grid(True, which="both", alpha=0.3); ax[1].set_ylim(0.8, 12)
fig.tight_layout(); fig.savefig(FIG / "cycle_time_wph.png", dpi=150); print("  saved figures/cycle_time_wph.png")
