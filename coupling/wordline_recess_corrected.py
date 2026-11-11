"""Word-line recesses filled through a gate-line slit, restated with the corrected slit transport
(knudsen_correction="plates": D_slit = (3/8) ln(L/w) x 2 w vbar/3 for the slit and for every recess;
paper/wp/wp9b_slit_transport.md). Same geometries, drives and outputs as coupling/wordline_recess.py, whose numbers
(hydraulic-diameter diffusivity) stay in figures/wordline_recess.csv. The voxel Monte Carlo reference for w100_Lr365 at
0.3 Pa is slit floor 10.70 s, deepest recess 10.88 s (paper/wp/wp9_voxel_mc.md (d)).

Original description:

Geometry from the 3D-NAND literature (sources in paper/wp/wp4_geometry_chemistry.md): word-line (recess) height
29 nm and inter-word-line oxide 22 nm, i.e. a 51 nm vertical pitch (Lee and Shin 2022); 98 oxide/nitride pairs
for a 5 um stack (Reiter et al. 2024; 98 x 51 nm = 5.0 um); slit width 250 nm (Onto Innovation 2023, 2XX-layer
device); lateral recess length 365 nm (Samsung 32-layer V-NAND, 730 nm slit-to-slit space; TechInsights 2014) and
1 um (slit pitch above 2 um, Onto 2023). TMA parameters of the chamber runs (s 5e-3, N0 4e18).
Three drives: (A) the chamber's centre-face precursor density from reactor/showerhead_pulse_patterned_AR100_long
(1 % TMA, 20 s pulse), one-way, as for the straight hole (14.05 s); (B) a constant 0.3 Pa mouth density (feature-only
benchmark); and the limits (plain slit, area-equivalent slit). The supply rule is evaluated for a wafer tiled with
these slits at pitch w_s + 2 L_r (fully patterned and 50 % patterned as the reference wafer).
Output: figures/wordline_recess_corrected.png, figures/wordline_recess_corrected.csv."""
import sys, csv, json, time
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "coupling"))
from aldconf import (Slit, Precursor, Process, simulate, simulate_driven, RecessedSlit, AreaEquivalentSlit, simulate_recessed,
                     simulate_recessed_constant, time_to_coverage, mean_speed, slit_knudsen_factor)
from aldconf.kinetics import number_density
from patterned_lib import load_patterned

CASE = ROOT / "reactor/showerhead_pulse_patterned_AR100_long"
T, U, C0, N0, S = 573.15, 0.0501, 3.79e20, 4e18, 5e-3
WL, OX, L_S = 29e-9, 22e-9, 5e-6                         # word-line height, oxide, stack height
GEOMS = {
    "w250_Lr365": RecessedSlit(250e-9, L_S, WL, 365e-9, WL + OX),
    "w250_Lr1000": RecessedSlit(250e-9, L_S, WL, 1000e-9, WL + OX),
    "w100_Lr365": RecessedSlit(100e-9, L_S, WL, 365e-9, WL + OX),
    "w100_Lr1000": RecessedSlit(100e-9, L_S, WL, 1000e-9, WL + OX),
}
BASE = "w250_Lr1000"
P_CONST = 0.3
CORR = "plates"


def main():
    t, f, params, pulse = load_patterned(CASE)
    n0_series = f["value"][:, 0] * params["c0"]
    pr = Precursor(sticking=S, site_density=N0, mass_amu=72.09)
    vbar = mean_speed(T, 72.09)
    rows = []; keep = {}
    for name, g in GEOMS.items():
        t0 = time.time()
        pitch = g.width + 2 * g.recess_length
        ratio_full = g.wall_area_ratio(pitch); ratio_half = 0.5 * ratio_full
        rule_full = 1.2 * (1 + ratio_full) * N0 / (U * C0); rule_half = 1.2 * (1 + ratio_half) * N0 / (U * C0)
        # (A) chamber-driven
        snaps = (3.0, 6.0, 9.0) if name == BASE else ()
        rA = simulate_recessed(g, pr, T, t, n0_series, substeps_per_sample=5, n_slit=200, n_recess=50, snapshot_times=snaps, knudsen_correction=CORR)
        plain = simulate_driven(Process(Slit(g.width, L_S, knudsen_correction=CORR), pr, T=T, p_partial=1.0, n_cells=200), t, n0_series, substeps_per_sample=5, record_depths=(0.0, 1.0))
        areaeq = simulate_driven(Process(AreaEquivalentSlit(g.width, L_S, g.area_factor, knudsen_correction=CORR), pr, T=T, p_partial=1.0, n_cells=200), t, n0_series, substeps_per_sample=5, record_depths=(0.0, 1.0))
        tA = dict(slit_bottom=time_to_coverage(rA.t, rA.hist_slit_bottom), recess_top=time_to_coverage(rA.t, rA.hist_recess_bottom_top),
                  recess_mid=time_to_coverage(rA.t, rA.hist_recess_bottom_mid), recess_deep=time_to_coverage(rA.t, rA.hist_recess_bottom_deep),
                  recess_min=time_to_coverage(rA.t, rA.hist_recess_bottom_min), plain_slit=time_to_coverage(plain[2], plain[3][:, 1]),
                  area_eq_slit=time_to_coverage(areaeq[2], areaeq[3][:, 1]))
        # (B) constant 0.3 Pa
        rB = simulate_recessed_constant(g, pr, T, P_CONST, 60.0, n_substeps=1200, n_samples=300, n_slit=200, n_recess=50, knudsen_correction=CORR)
        tB = dict(slit_bottom=time_to_coverage(rB.t, rB.hist_slit_bottom), recess_deep=time_to_coverage(rB.t, rB.hist_recess_bottom_deep),
                  recess_min=time_to_coverage(rB.t, rB.hist_recess_bottom_min))
        n0c = number_density(P_CONST, T)
        # stoichiometric minimum at the effusion limit: molecules needed per mouth area / impingement rate
        walls_per_mouth = g.wall_area_per_width() / g.width
        t_effusion = walls_per_mouth * N0 / (n0c * vbar / 4.0)
        row = dict(geometry=name, knudsen_correction=CORR, slit_factor=float(slit_knudsen_factor(g.width, L_S)), recess_factor=float(slit_knudsen_factor(g.recess_gap, g.recess_length)), w_s_nm=g.width * 1e9, L_r_nm=g.recess_length * 1e9, g_r_nm=g.recess_gap * 1e9, p_r_nm=g.recess_pitch * 1e9,
                   L_s_um=L_S * 1e6, area_factor=g.area_factor, slit_pitch_um=pitch * 1e6, wall_ratio_full=ratio_full, wall_ratio_half=ratio_half,
                   rule_full_s=rule_full, rule_half_s=rule_half,
                   A_slit_bottom_t90=tA["slit_bottom"], A_recess_top_t90=tA["recess_top"], A_recess_mid_t90=tA["recess_mid"], A_recess_deep_t90=tA["recess_deep"],
                   A_recess_min_t90=tA["recess_min"], A_lag_s=tA["recess_min"] - tA["slit_bottom"], A_plain_slit_t90=tA["plain_slit"], A_area_eq_slit_t90=tA["area_eq_slit"],
                   B_slit_bottom_t90=tB["slit_bottom"], B_recess_deep_t90=tB["recess_deep"], B_recess_min_t90=tB["recess_min"],
                   B_exposure_Pa_s=tB["recess_min"] * P_CONST, B_effusion_limit_s=t_effusion, B_ratio_to_effusion=tB["recess_min"] / t_effusion)
        rows.append(row); keep[name] = (rA, rB)
        print(f"{name}: A_f {g.area_factor:.1f}, wall ratio (full tiling) {ratio_full:.0f}, rule {rule_full:.1f} s (50 %: {rule_half:.1f} s); "
              f"chamber-driven t90 slit bottom {tA['slit_bottom']:.2f} s, recess bottoms top/mid/deep {tA['recess_top']:.2f}/{tA['recess_mid']:.2f}/{tA['recess_deep']:.2f} s "
              f"(min {tA['recess_min']:.2f}), plain slit {tA['plain_slit']:.2f}, area-eq slit {tA['area_eq_slit']:.2f}; "
              f"constant 0.3 Pa: slit bottom {tB['slit_bottom']:.2f} s, recess min {tB['recess_min']:.2f} s = {tB['recess_min'] * P_CONST:.2f} Pa s, effusion limit {t_effusion:.3f} s  [{time.time() - t0:.0f} s]", flush=True)
    with open(ROOT / "figures/wordline_recess_corrected.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    # figure: base geometry, chamber-driven, three times
    rA, rB = keep[BASE]; g = GEOMS[BASE]
    times = sorted(rA.snapshots)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), gridspec_kw=dict(height_ratios=(1.3, 1)), constrained_layout=True)
    for ax, tt in zip(axes[0], times):
        ths, thr = rA.snapshots[tt]
        im = ax.imshow(thr.T, origin="lower", aspect="auto", extent=(0, L_S * 1e6, 0, g.recess_length * 1e9), vmin=0, vmax=1, cmap="viridis")
        ax.set_xlabel("depth along the slit (um)"); ax.set_ylabel("distance into the recess (nm)"); ax.set_title(f"t = {tt:g} s: coverage in the recesses", fontsize=9)
    fig.colorbar(im, ax=axes[0].tolist(), fraction=0.02, pad=0.01, label="coverage", location="right")
    ax = axes[1][0]
    for tt, c in zip(times, ("tab:green", "tab:purple", "tab:brown")):
        ths, thr = rA.snapshots[tt]
        ax.plot(rA.x * 1e6, ths, color=c, label=f"{tt:g} s, slit wall")
        ax.plot(rA.x * 1e6, thr[:, -1], color=c, ls="--", label=f"{tt:g} s, recess bottoms")
    ax.set_xlabel("depth along the slit (um)"); ax.set_ylabel("coverage"); ax.set_title("Slit wall (solid) and recess bottoms (dashed)", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ax = axes[1][1]
    ax.plot(rA.t, rA.hist_slit_bottom, "k", label="slit bottom (flat wall)")
    ax.plot(rA.t, rA.hist_recess_bottom_top, color="tab:blue", label="recess bottom at the slit mouth")
    ax.plot(rA.t, rA.hist_recess_bottom_deep, color="tab:red", label="recess bottom at the slit bottom")
    ax.plot(rA.t, rA.hist_slit_mouth, color="0.6", ls=":", label="slit mouth")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlim(0, 22); ax.set_xlabel("time (s)"); ax.set_ylabel("coverage")
    ax.set_title(f"Chamber-driven (1 % TMA), corrected slit transport: w_s {g.width * 1e9:.0f} nm, L_r {g.recess_length * 1e9:.0f} nm, L_s {L_S * 1e6:.0f} um", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ax = axes[1][2]
    ax.plot(rB.t, rB.hist_slit_bottom, "k", label="slit bottom (flat wall)")
    ax.plot(rB.t, rB.hist_recess_bottom_top, color="tab:blue", label="recess bottom at the slit mouth")
    ax.plot(rB.t, rB.hist_recess_bottom_deep, color="tab:red", label="recess bottom at the slit bottom")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlabel("time (s)"); ax.set_ylabel("coverage"); ax.set_xscale("log")
    ax.set_title(f"Constant {P_CONST:g} Pa mouth density (feature-only)", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both")
    fig.savefig(ROOT / "figures/wordline_recess_corrected.png", dpi=150); print("saved figures/wordline_recess_corrected.png")


if __name__ == "__main__":
    main()
