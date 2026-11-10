"""Coverage across the 200 mm wafer of the axisymmetric Fiji F200 approximation, against the thickness profile measured by
Gakis et al. (Chem. Eng. Sci. 195, 399 (2019), Fig. 10, 300 C, 25 ms TMA valve, 5 s purges, experiment 11).

Measured values were digitised from Fig. 10 of the 2019 author manuscript (hal-02329135), 400 dpi crop, connected
components of the filled-square markers, axes calibrated on the frame (-100..100 mm, 0.85..1.1 A/cycle); the centroid
precision is 0.0003 A/cycle and the axis calibration about 0.002 A/cycle, so each value carries +/- 0.003 A/cycle. The
digitised non-uniformity (max - min)/min = 2.37 % against the 2.40 % stated in the paper is the check on the reading.
The 60 ms profile (experiment 12) is given in the paper only as the map of Fig. 9d (0.98 to 1.01 A/cycle) and the
non-uniformity 0.58 %; its wafer value "1 A/cycle at 300 C" (Section 4.1.1, Fig. 3) is used as the saturated GPC.

Model output per case: final coverage theta(r) (= predicted GPC / saturated GPC), non-uniformity over |r| <= 95 mm (the
measured span), time to 90 %, the impingement exposure int p_s dt of the wafer (their Fig. 12 gives 0.0118-0.0130
mol/m2 = 0.55-0.60 Pa s for the 25 ms valve) and the fraction of the delivered TMA captured by the wafer.
Output: figures/validation_gakis.png and .csv.   Run: .venv/bin/python reactor/validation_gakis/postprocess_gakis.py
"""
import re, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE.parent))
from postprocess_langmuir import patch_block, list_or_uniform, scalar_entry

KB = 1.380649e-23
CASES = ["pulse_25ms", "pulse_60ms_x1.5", "pulse_60ms_x2.4", "pulse_25ms_s5e-3"]
# Gakis 2019 Fig. 10, squares (300 C, 25 ms, 5 s purge): position (mm, loading door at +100) and GPC (A/cycle)
MEAS_POS = np.array([-95.5, -76.4, -57.2, -38.2, -19.3, -0.7, 18.4, 37.4, 56.6, 75.5, 94.7])
MEAS_GPC = np.array([0.9751, 0.9834, 0.9873, 0.9927, 0.9951, 0.9961, 0.9951, 0.9982, 0.9961, 0.9927, 0.9826])
MEAS_ERR = 0.003
GPC_SAT = 1.00          # 60 ms value at 300 C (their text and Fig. 3); Fig. 9d map spans 0.98-1.01
NONUNI_MEAS = {"25 ms": 2.40, "60 ms": 0.58}; NONUNI_THEIR_MODEL = {"25 ms": 2.45, "60 ms": 0.83}

def read_list(txt, key):
    m = re.search(r"\n%s\s*\n\(\n(.*?)\n\)" % key, txt, re.S)
    return m.group(1)

def wafer_face_radii(case):
    pm = case / "constant/polyMesh"
    pts = np.array([list(map(float, l.strip("()").split())) for l in re.search(r"\n\(\n(.*?)\n\)", (pm / "points").read_text(), re.S).group(1).splitlines()])
    faces_txt = re.search(r"\n\(\n(.*?)\n\)", (pm / "faces").read_text(), re.S).group(1).splitlines()
    bnd = (pm / "boundary").read_text()
    m = re.search(r"wafer\s*\{[^}]*nFaces\s+(\d+);[^}]*startFace\s+(\d+);", bnd, re.S)
    n, start = int(m.group(1)), int(m.group(2))
    r = np.empty(n)
    for i in range(n):
        ids = list(map(int, faces_txt[start + i].split("(")[1].rstrip(")").split()))
        c = pts[ids].mean(0); r[i] = np.hypot(c[1], c[2])
    return r

def load(case):
    r = wafer_face_radii(case)
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times]); n = len(r)
    theta = np.zeros((len(t), n)); face = np.zeros((len(t), n)); s = vbar = c0 = N0 = None
    for i, (_, d) in enumerate(times):
        blk = patch_block((d / "T").read_text())
        if s is None and "langmuirSticking" in blk:
            s, vbar, c0, N0 = (scalar_entry(blk, k) for k in ("s", "vbar", "c0", "N0"))
        th = list_or_uniform(blk, "theta", n); va = list_or_uniform(blk, "value", n)
        theta[i] = 0 if th is None else th; face[i] = 0 if va is None else va
    prm = dict(l.split() for l in (case / "params.txt").read_text().splitlines())
    return r, t, theta, face, s, vbar, c0, N0, {k: float(v) for k, v in prm.items()}

if __name__ == "__main__":
    rows, curves = [], {}
    for name in CASES:
        case = HERE / name
        if not (case / "params.txt").exists() or not any(re.fullmatch(r"0\.\d+|[1-9][0-9.]*", d.name) for d in case.iterdir()):
            print("skip", name); continue
        r, t, theta, face, s, vbar, c0, N0, prm = load(case)
        T = prm["T"]; p_s = face * c0 * KB * T                       # TMA partial pressure at the wafer faces, Pa
        th_end = theta[-1]; sel = r <= 0.095
        nonuni = (th_end[sel].max() - th_end[sel].min()) / th_end[sel].min() * 100
        t90 = np.array([t[np.argmax(theta[:, j] >= 0.9)] if theta[:, j].max() >= 0.9 else np.nan for j in range(len(r))])
        expo = np.trapezoid(p_s, t, axis=0)                           # Pa s impingement exposure per face
        uptake = (vbar / 4) * s * (1 - theta) * face * c0            # molecules/m2/s
        dt = np.diff(t, prepend=0)
        # face areas: 1 mm radial rings (wafer blocks have 20 + 80 uniform radial cells)
        dr = 1e-3; area = 2 * np.pi * r * dr
        captured = np.sum((uptake * dt[:, None]).sum(0) * area)
        delivered = prm["n_mol"] * 6.02214e23
        j_c, j_e = np.argmin(r), np.argmax(r)
        rows.append(dict(case=name, amount_rel_25ms=prm["amount"], s=s, N0=N0, TMA_molecules=delivered,
                         theta_end_centre=th_end[j_c], theta_end_edge=th_end[j_e], theta_end_min=th_end[sel].min(), theta_end_max=th_end[sel].max(),
                         nonuniformity_pct=nonuni, t90_centre=t90[j_c], t90_worst=np.nanmax(t90[sel]) if np.isfinite(t90[sel]).any() else np.nan,
                         exposure_centre_Pa_s=expo[j_c], exposure_min_Pa_s=expo.min(), exposure_max_Pa_s=expo.max(),
                         exposure_nonuni_pct=(expo.max() - expo.min()) / expo.min() * 100, captured_fraction=captured / delivered,
                         p_s_peak_centre_Pa=p_s[:, j_c].max(), theta_at_pulse_end_0p3s=theta[np.argmin(abs(t - 0.3)), j_c]))
        curves[name] = (r, t, theta, p_s)
        print(f"{name}: theta_end centre {th_end[j_c]:.4f} edge {th_end[j_e]:.4f}, nonuni {nonuni:.2f} %, t90 centre {t90[j_c]}, "
              f"exposure centre {expo[j_c]:.3f} Pa s ({expo.min():.3f}-{expo.max():.3f}), captured {captured/delivered*100:.1f} %, peak p_s {p_s[:, j_c].max():.3f} Pa")
    df = pd.DataFrame(rows)
    meas = pd.DataFrame([dict(case="measured, Gakis 2019 Fig. 10, 25 ms (digitised)", theta_end_min=MEAS_GPC.min() / GPC_SAT, theta_end_max=MEAS_GPC.max() / GPC_SAT,
                              theta_end_centre=MEAS_GPC[5] / GPC_SAT, nonuniformity_pct=(MEAS_GPC.max() - MEAS_GPC.min()) / MEAS_GPC.min() * 100),
                         dict(case="measured, Gakis 2019, 25 ms (stated non-uniformity)", nonuniformity_pct=NONUNI_MEAS["25 ms"]),
                         dict(case="measured, Gakis 2019, 60 ms (stated)", nonuniformity_pct=NONUNI_MEAS["60 ms"], theta_end_centre=1.0),
                         dict(case="their CFD + kinetics, 25 ms", nonuniformity_pct=NONUNI_THEIR_MODEL["25 ms"]),
                         dict(case="their CFD + kinetics, 60 ms", nonuniformity_pct=NONUNI_THEIR_MODEL["60 ms"]),
                         dict(case="their CFD, 25 ms impingement exposure (Fig. 12, 2018)", exposure_min_Pa_s=0.0118 * 46.45, exposure_max_Pa_s=0.0130 * 46.45, exposure_nonuni_pct=10.7)])
    out = pd.concat([df, meas], ignore_index=True)
    out.to_csv(ROOT / "figures/validation_gakis.csv", index=False)
    print(out.to_string())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
    col = {"pulse_25ms": BLUE, "pulse_60ms_x1.5": ORANGE, "pulse_60ms_x2.4": AQUA, "pulse_25ms_s5e-3": YELLOW}
    lab = {"pulse_25ms": "25 ms valve (0.63 mg TMA), s = 3.44e-3", "pulse_60ms_x1.5": "60 ms, 1.5 x amount", "pulse_60ms_x2.4": "60 ms, 2.4 x amount", "pulse_25ms_s5e-3": "25 ms, s = 5e-3"}
    ax = axes[0]
    for name, (r, t, theta, p_s) in curves.items():
        ax.plot(r * 1e3, theta[-1], color=col[name], lw=2, label=lab[name])
    ax.errorbar(np.abs(MEAS_POS[MEAS_POS < 0]), MEAS_GPC[MEAS_POS < 0] / GPC_SAT, yerr=MEAS_ERR, fmt="s", color="#1a1a19", ms=5, label="measured 25 ms, side away from the loading door (Fig. 10)")
    ax.errorbar(np.abs(MEAS_POS[MEAS_POS > 0]), MEAS_GPC[MEAS_POS > 0] / GPC_SAT, yerr=MEAS_ERR, fmt="o", mfc="none", color="#1a1a19", ms=5, label="measured 25 ms, loading-door side")
    ax.set_xlabel("radius (mm)"); ax.set_ylabel("coverage after the cycle = GPC / GPC(60 ms)"); ax.set_ylim(0.9, 1.02)
    ax.set_title("Final TMA coverage across the 200 mm wafer", fontsize=9); ax.legend(fontsize=6.5, loc="lower left"); ax.grid(alpha=0.25)
    ax = axes[1]
    for name, (r, t, theta, p_s) in curves.items():
        ax.plot(t, theta[:, np.argmin(r)], color=col[name], lw=2, label=lab[name] + ", centre")
        ax.plot(t, theta[:, np.argmax(r)], color=col[name], lw=1, ls="--")
    ax.set_xlim(0, 0.6); ax.set_xlabel("time from the valve opening (s)"); ax.set_ylabel("coverage theta"); ax.set_title("Coverage vs time (dashed: wafer edge)", fontsize=9)
    ax.legend(fontsize=6.5, loc="lower right"); ax.grid(alpha=0.25)
    ax = axes[2]
    for name, (r, t, theta, p_s) in curves.items():
        ax.plot(t, p_s[:, np.argmin(r)], color=col[name], lw=2, label=lab[name])
    ax.axhspan(0, 9.6, color="#f0efec", zorder=0); ax.text(0.45, 9.0, "Ar at 9.6 Pa", fontsize=7, color="#555555", va="top")
    ax.set_xlim(0, 0.6); ax.set_xlabel("time (s)"); ax.set_ylabel("TMA partial pressure at the wafer centre (Pa)")
    ax.set_title("Surface TMA pressure (their Fig. 7b: 0.005-0.028 Torr = 0.7-3.7 Pa at 25 ms)", fontsize=9); ax.legend(fontsize=6.5); ax.grid(alpha=0.25)
    fig.suptitle("Axisymmetric Fiji F200 approximation vs Gakis et al. 2018/2019: 300 C, 9.6 Pa Ar, 180 sccm, 200 mm wafer", fontsize=10)
    fig.tight_layout(); fig.savefig(ROOT / "figures/validation_gakis.png", dpi=150); print("saved", ROOT / "figures/validation_gakis.png")
