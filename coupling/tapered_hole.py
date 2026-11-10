"""Tapered (frustum) channel holes: time to 90 % bottom coverage for a 100 nm-mouth AR-100 hole whose bottom diameter
is 100, 90, 80 or 60 % of the mouth, driven by the chamber's centre-face precursor density from
reactor/showerhead_pulse_patterned_AR100_long (1 % TMA, 300 Pa, 573 K, 300 sccm, 20 s pulse), exactly as
docs/build_assets.py drives the straight hole; so the straight-hole result is the one-way counterpart of the chamber's
14.05 s. Also: the same holes at a constant 0.3 Pa mouth density (feature-only benchmark, as in the Machball check),
the frustum wall area relative to the straight hole, and the supply rule 1.2 (1 + 55.9 A_rel) N0/(U c0) with the
tapered wall area. The 90 % bottom diameter is the literature-anchored profile (cryogenic 3D-NAND channel-hole etch,
about 10 nm of top-to-bottom CD loss over 10 um; see paper/wp/wp4_geometry_chemistry.md); 80 and 60 % are a
parametric sweep. Output: figures/tapered_hole.png, figures/tapered_hole.csv."""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "coupling"))
from aldconf import Hole, TaperedHole, Precursor, Process, simulate, simulate_driven, dose_for_coverage
from patterned_lib import load_patterned

CASE = ROOT / "reactor/showerhead_pulse_patterned_AR100_long"
D, AR, T = 100e-9, 100, 573.15
U, C0, N0, S, WALL = 0.0501, 3.79e20, 4e18, 5e-3, 1.78e13 * np.pi * D * AR * D
RATIOS = (1.0, 0.9, 0.8, 0.6)


def t90_of(t, series, level=0.9):
    idx = np.where(series >= level)[0]
    if len(idx) == 0: return float("nan")
    i = idx[0]
    return float(t[i]) if i == 0 else float(t[i - 1] + (level - series[i - 1]) * (t[i] - t[i - 1]) / (series[i] - series[i - 1]))


def main():
    t, f, params, pulse = load_patterned(CASE)
    n0_series = f["value"][:, 0] * params["c0"]
    chamber_t90 = t90_of(t, f["thetaBottom"][:, 0])
    pr = Precursor(sticking=params["s"], site_density=params["N0"], mass_amu=72.09)
    depths = np.linspace(0, 1, 21)
    rows = []; profiles = {}; hists = {}
    for br in RATIOS:
        g = TaperedHole(D, AR * D, br)
        proc = Process(g, pr, T=T, p_partial=1.0, n_cells=200)
        x, theta, ht, hist = simulate_driven(proc, t, n0_series, substeps_per_sample=5, record_depths=depths)
        hists[br] = (ht, hist[:, -1])
        profiles[br] = {tt: hist[np.argmin(abs(ht - tt)), :] for tt in (4.0, 8.0, 12.0)}
        t90_driven = t90_of(ht, hist[:, -1])
        pc = Process(g, pr, T=T, p_partial=0.3, t_dose=1.0, n_cells=200, n_substeps=200, narrowing=False)
        t90_const = dose_for_coverage(pc, target=0.90, t_lo=1e-3, t_hi=1e3)
        a_rel = g.area_ratio_to_straight()
        rule = 1.2 * (1 + WALL * a_rel) * N0 / (U * C0)
        rows.append(dict(bottom_ratio=br, bottom_nm=br * D * 1e9, area_rel=a_rel, t90_driven_s=t90_driven,
                         t90_const_0p3Pa_s=t90_const, rule_s=rule, ratio_driven_to_rule=t90_driven / rule))
        print(f"bottom {br * 100:.0f} %: area {a_rel:.3f} x straight; driven bottom t90 {t90_driven:.2f} s; constant 0.3 Pa t90 {t90_const:.2f} s; rule {rule:.2f} s", flush=True)
    print(f"chamber (two-way) bottom t90 at the centre: {chamber_t90:.2f} s; straight one-way: {rows[0]['t90_driven_s']:.2f} s")
    import csv
    with open(ROOT / "figures/tapered_hole.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    (ROOT / "figures/tapered_hole.json").write_text(json.dumps(dict(chamber_bottom_t90_s=chamber_t90, rows=rows), indent=1))
    # figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    colors = {1.0: "k", 0.9: "tab:blue", 0.8: "tab:orange", 0.6: "tab:red"}
    for br in RATIOS:
        ht, hb = hists[br]
        axes[0].plot(ht, hb, color=colors[br], label=f"bottom {br * D * 1e9:.0f} nm (area {rows[RATIOS.index(br)]['area_rel']:.2f}x)")
    axes[0].axhline(0.9, color="0.6", ls=":"); axes[0].axvline(chamber_t90, color="0.4", ls="--", lw=0.8)
    axes[0].text(chamber_t90 + 0.2, 0.05, f"chamber two-way\n{chamber_t90:.2f} s", fontsize=7, color="0.3")
    axes[0].set_xlabel("time (s)"); axes[0].set_ylabel("bottom coverage"); axes[0].set_xlim(0, 22)
    axes[0].set_title("AR-100 hole bottoms driven by the chamber centre density", fontsize=9); axes[0].legend(fontsize=7); axes[0].grid(alpha=0.3)
    for br, ls in ((1.0, "-"), (0.6, "--")):
        for tt, c in zip((4.0, 8.0, 12.0), ("tab:green", "tab:purple", "tab:brown")):
            axes[1].plot(depths * AR, profiles[br][tt], ls=ls, color=c, label=f"{tt:g} s, bottom {br * D * 1e9:.0f} nm")
    axes[1].set_xlabel("depth / mouth diameter"); axes[1].set_ylabel("coverage"); axes[1].set_title("Coverage along the hole, straight (solid) vs 60 % bottom (dashed)", fontsize=9)
    axes[1].legend(fontsize=7, ncol=2); axes[1].grid(alpha=0.3)
    brs = [r["bottom_ratio"] for r in rows]
    axes[2].plot(brs, [r["t90_driven_s"] for r in rows], "o-", color="tab:blue", label="chamber-driven bottom t90")
    axes[2].plot(brs, [r["rule_s"] for r in rows], "s--", color="tab:gray", label="rule with frustum wall area")
    axes[2].plot(brs, [r["t90_const_0p3Pa_s"] for r in rows], "^-", color="tab:red", label="constant 0.3 Pa mouth, t90")
    axes[2].set_xlabel("bottom / mouth diameter"); axes[2].set_ylabel("time to 90 % bottom coverage (s)"); axes[2].invert_xaxis()
    axes[2].set_title("Taper: supply-limited chamber vs constant-density feature", fontsize=9); axes[2].legend(fontsize=7); axes[2].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(ROOT / "figures/tapered_hole.png", dpi=150); print("saved figures/tapered_hole.png")


if __name__ == "__main__":
    main()
