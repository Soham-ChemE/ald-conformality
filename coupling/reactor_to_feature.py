"""Reactor -> feature coupling (one way): the TMA number density at the wafer surface from the chamber model,
at several radii, drives the feature-scale Langmuir/Knudsen model of a hole of diameter d and aspect ratio AR.

Output per case, radius and AR: coverage profile down the hole at the end of the run, bottom-coverage history,
time for the hole bottom to reach 90 % coverage (nan if not reached). Compared with the flat-wafer t90.
"""
import sys, re
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "reactor"))
from aldconf import Hole, Precursor, Process, simulate_driven
from postprocess_langmuir import load as load_cov

OUT = ROOT / "figures"
D_HOLE = 100e-9
ARS = (20, 50, 100)
RADII = {0: "centre", 100: "r = 100 mm", 140: "r = 140 mm", 149: "r = 149 mm"}
CASES = {"base (300 sccm, 1 % TMA)": ROOT / "reactor/sweep/base/pulse",
         "3 % TMA": ROOT / "reactor/sweep/x0.03/pulse",
         "600 sccm": ROOT / "reactor/sweep/Q600/pulse"}

def run_case(case_dir):
    t, theta_w, face, s, vbar, c0, N0, pulse = load_cov(case_dir)
    res = {}
    for j, lab in RADII.items():
        n0_t = face[:, j] * c0                     # TMA number density at the wafer surface, 1/m3
        t90_flat = t[np.argmax(theta_w[:, j] >= 0.9)] if theta_w[:, j].max() >= 0.9 else np.nan
        for AR in ARS:
            pr = Precursor(sticking=s, site_density=N0, mass_amu=72.09)
            proc = Process(Hole(D_HOLE, AR * D_HOLE), pr, T=573.15, n_cells=200)
            x, th, ht, h = simulate_driven(proc, t, n0_t, substeps_per_sample=4, record_depths=(0.0, 0.5, 1.0))
            bottom = h[:, 2]
            t90_bottom = ht[np.argmax(bottom >= 0.9)] if bottom.max() >= 0.9 else np.nan
            res[(j, AR)] = dict(x=x, theta=th, ht=ht, hist=h, t90_bottom=t90_bottom, t90_flat=t90_flat, pulse=pulse)
    return res, pulse

if __name__ == "__main__":
    allres = {}
    for lab, cd in CASES.items():
        if not cd.exists(): print("missing", cd); continue
        allres[lab], pulse = run_case(cd)
        print(f"== {lab}: pulse {pulse:g} s")
        for (j, AR), r in allres[lab].items():
            print(f"   {RADII[j]:11s} AR {AR:3d}: flat t90 {r['t90_flat']:.2f} s | hole bottom coverage at end {r['theta'][-1]:.3f}, "
                  f"mid {np.interp(0.5, [0,1], [0,0]) if False else r['hist'][-1,1]:.3f}, t90 bottom {r['t90_bottom']:.2f} s")
    # Figure 1: coverage down the hole at the end of the run, base case, all radii, AR 100 and AR 50
    base = next(iter(allres.values()))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for AR, ls in ((50, "--"), (100, "-")):
        for j, lab in RADII.items():
            r = base[(j, AR)]
            axes[0].plot(r["x"] / D_HOLE, r["theta"], ls, label=f"AR {AR}, {lab}")
    axes[0].set_xlabel("depth / diameter"); axes[0].set_ylabel("coverage at end of pulse + purge"); axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Base case: coverage down a 100 nm hole", fontsize=9); axes[0].legend(fontsize=6, ncol=2)
    # Figure 2: bottom coverage vs time, centre, by AR, base case
    for AR in ARS:
        r = base[(0, AR)]; axes[1].plot(r["ht"], r["hist"][:, 2], label=f"hole bottom, AR {AR}")
    r = base[(0, 20)]; axes[1].plot(r["ht"], r["hist"][:, 0], "k:", label="hole mouth (= flat wafer)")
    axes[1].axvspan(0, base[(0, 20)]["pulse"], color="0.92"); axes[1].set_xlabel("time (s)"); axes[1].set_ylabel("coverage"); axes[1].set_ylim(0, 1.05)
    axes[1].set_title("Base case, wafer centre: bottom coverage vs time", fontsize=9); axes[1].legend(fontsize=7)
    # Figure 3: final bottom coverage vs AR for each case (centre and r=140)
    for lab, res in allres.items():
        axes[2].plot(ARS, [res[(0, AR)]["theta"][-1] for AR in ARS], "o-", label=f"{lab}, centre")
        axes[2].plot(ARS, [res[(140, AR)]["theta"][-1] for AR in ARS], "s--", alpha=0.6, label=f"{lab}, r = 140 mm")
    axes[2].set_xlabel("aspect ratio"); axes[2].set_ylabel("hole-bottom coverage after 1 s pulse"); axes[2].set_ylim(0, 1.05)
    axes[2].set_title("What a 1 s pulse leaves at the bottom of the hole", fontsize=9); axes[2].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(OUT / "coupling_hole_coverage.png", dpi=150); print("saved", OUT / "coupling_hole_coverage.png")
