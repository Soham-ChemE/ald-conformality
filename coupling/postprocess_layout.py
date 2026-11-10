"""Spatially varying pattern layouts (reactor/showerhead_pulse_layout_*): per-radius saturation times of the flat
surface and the hole bottoms, the area-weighted wall-area ratio from the per-face ARField / holeDensityField, and the
supply-limit rule t = 1.2 (1 + A_wall/A_wafer) N0/(U c0) against the worst and the deepest-hole t90 on the wafer.
Writes figures/pattern_layout.csv and figures/pattern_layout.png.
  .venv/bin/python coupling/postprocess_layout.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "coupling")]
from patterned_lib import load_patterned

OUT = ROOT / "figures"
U, N0_RULE, RULE_CONST = 0.0501, 4e18, 1.2        # m/s, sites/m2, rule constant (paper, pulse-time map)
T90_UNIFORM_AR100 = 14.05                         # s, hole bottom, uniform AR-100 wafer at 1 % TMA (figures/pulse_time_map.csv)
T90_BLANKET = 0.35                                # s, blanket wafer centre
r_mm = np.arange(150) + 0.5
area_w = 2 * np.pi * r_mm * 1.0                   # face area weights (annuli of 1 mm width)

CASES = {
    "dense_centre": ("showerhead_pulse_layout_dense_centre", "dense centre: 1.78e13 /m2 for r < 75 mm, ramp to 0 at 100 mm, bare beyond; AR 100"),
    "ar_ramp": ("showerhead_pulse_layout_ar_ramp", "AR ramp: AR 100 at the centre to 60 at 150 mm; 1.78e13 /m2 everywhere"),
}
UNIFORM = ("showerhead_pulse_patterned_AR100_long", "uniform AR 100, 1.78e13 /m2 (reference run)")


def t90_profile(arr, t, mask=None):
    """First write time at which each face reaches 90 % coverage; nan if never, or where mask is False."""
    out = np.full(arr.shape[1], np.nan)
    for j in range(arr.shape[1]):
        if mask is not None and not mask[j]: continue
        if arr[:, j].max() >= 0.9: out[j] = t[np.argmax(arr[:, j] >= 0.9)]
    return out


def analyse(casedir):
    t, f, prm, pulse = load_patterned(casedir)
    AR, hd = f["ARField"][-1], f["holeDensityField"][-1]        # constant in time; 0/T carries only the scalars
    holes = hd > 0
    wall_local = hd * np.pi * prm["d"] * AR * prm["d"]              # local wall area per unit wafer area
    wall_ratio = float(np.sum(wall_local * area_w) / np.sum(area_w))  # area-weighted over the wafer
    rule = RULE_CONST * (1 + wall_ratio) * prm["N0"] / (U * prm["c0"])
    t90_b = t90_profile(f["thetaBottom"], t, holes)
    t90_f = t90_profile(f["thetaFlat"], t)
    j_deep = int(np.argmax(np.where(holes, AR, -1)))                  # deepest holes; ties -> innermost face
    return dict(t=t, f=f, prm=prm, pulse=pulse, AR=AR, hd=hd, holes=holes, wall_local=wall_local,
                wall_ratio=wall_ratio, rule=rule, t90_bottom=t90_b, t90_flat=t90_f, j_deep=j_deep,
                open_area=float(np.sum(hd * np.pi * prm["d"] ** 2 / 4 * area_w) / np.sum(area_w)))


res = {k: analyse(ROOT / "reactor" / v[0]) for k, v in CASES.items()}
if (ROOT / "reactor" / UNIFORM[0]).exists():
    res["uniform"] = analyse(ROOT / "reactor" / UNIFORM[0])

# ---- table
rows = []
for k, R in res.items():
    tb, tf = R["t90_bottom"], R["t90_flat"]
    jw = int(np.nanargmax(tb))
    # the rule with the wafer-average wall-area ratio (rule_s) and with the LOCAL ratio of the face in question
    local_rule = RULE_CONST * (1 + R["wall_local"]) * R["prm"]["N0"] / (U * R["prm"]["c0"])
    rows.append(dict(layout=k, pulse_s=R["pulse"], open_area_pct=100 * R["open_area"], wall_ratio=R["wall_ratio"],
                     rule_s=R["rule"], t90_bottom_centre=tb[0], t90_bottom_worst=tb[jw], r_worst_mm=r_mm[jw],
                     t90_bottom_deepest=tb[R["j_deep"]], r_deepest_mm=r_mm[R["j_deep"]],
                     rule_over_worst=R["rule"] / tb[jw], rule_over_deepest=R["rule"] / tb[R["j_deep"]],
                     local_rule_worst_s=local_rule[jw], local_rule_over_worst=local_rule[jw] / tb[jw],
                     t90_bottom_edge=tb[149], local_rule_edge_s=local_rule[149],
                     t90_flat_centre=tf[0], t90_flat_edge=tf[149], run_end_s=R["t"][-1]))
tab = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.float_format", lambda x: f"{x:.3f}")
print(tab.to_string(index=False))

# ---- the open-field question for the dense-centre layout
R = res["dense_centre"]
bare = ~R["holes"]; dense = R["hd"] >= 0.999 * R["hd"].max(); ramp = R["holes"] & ~dense
tf, tb, t, f = R["t90_flat"], R["t90_bottom"], R["t"], R["f"]
k03 = np.argmin(abs(t - 0.3))
print("\nDense-centre layout, flat-surface t90 by region (blanket wafer reference 0.35 s centre / 0.38 s edge;"
      " uniform AR-100 wafer flat regions 1.95 s):")
print(f"  dense array r < 75 mm:   t90_flat {np.nanmin(tf[dense]):.2f} to {np.nanmax(tf[dense]):.2f} s;"
      f" hole bottom t90 {np.nanmin(tb[dense]):.2f} to {np.nanmax(tb[dense]):.2f} s")
print(f"  ramp 75-100 mm:          t90_flat {np.nanmin(tf[ramp]):.2f} to {np.nanmax(tf[ramp]):.2f} s;"
      f" hole bottom t90 {np.nanmin(tb[ramp]):.2f} to {np.nanmax(tb[ramp]):.2f} s")
print(f"  bare wafer r > 100 mm:   t90_flat {np.nanmin(tf[bare]):.2f} (r = {r_mm[bare][np.nanargmin(tf[bare])]:.1f} mm)"
      f" to {np.nanmax(tf[bare]):.2f} s (r = {r_mm[bare][np.nanargmax(tf[bare])]:.1f} mm)")
print(f"  surface TMA density / c0 at 0.3 s: r = 0.5 mm {f['value'][k03, 0]:.3f}, r = 100.5 mm {f['value'][k03, 100]:.3f},"
      f" r = 125.5 mm {f['value'][k03, 125]:.3f}, r = 149.5 mm {f['value'][k03, 149]:.3f}")
print(f"  flat coverage at 0.35 s: r = 100.5 mm {f['thetaFlat'][np.argmin(abs(t-0.35)), 100]:.3f},"
      f" r = 125.5 mm {f['thetaFlat'][np.argmin(abs(t-0.35)), 125]:.3f}, r = 149.5 mm {f['thetaFlat'][np.argmin(abs(t-0.35)), 149]:.3f}")
if "uniform" in res:
    print(f"  hole bottom t90 at the centre: dense-centre layout {tb[0]:.2f} s vs uniform wafer {res['uniform']['t90_bottom'][0]:.2f} s")

# ---- csv
csv = pd.DataFrame({"radius_mm": r_mm})
for k in CASES:
    R = res[k]
    csv[f"{k}_AR"] = R["AR"]; csv[f"{k}_holeDensity_per_m2"] = R["hd"]
    csv[f"{k}_t90_flat_s"] = R["t90_flat"]; csv[f"{k}_t90_bottom_s"] = R["t90_bottom"]
csv.to_csv(OUT / "pattern_layout.csv", index=False, float_format="%.6g")

# ---- figure
COL = {"dense_centre": "#2a78d6", "ar_ramp": "#eb6834"}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
ax = axes[0]
for k, (cd, lab) in CASES.items():
    R = res[k]
    ax.plot(r_mm, R["t90_bottom"], "-", color=COL[k], lw=2, label=f"{k.replace('_', ' ')}: hole bottom")
    ax.plot(r_mm, R["t90_flat"], "--", color=COL[k], lw=1.6, label=f"{k.replace('_', ' ')}: flat surface")
ax.axhline(T90_UNIFORM_AR100, color="0.45", ls=":", lw=1.2)
ax.text(2, T90_UNIFORM_AR100 + 0.3, f"uniform AR-100 wafer, hole bottom: {T90_UNIFORM_AR100} s", fontsize=8, color="0.3")
ax.axhline(T90_BLANKET, color="0.45", ls=":", lw=1.2)
ax.text(2, T90_BLANKET + 0.3, f"blanket wafer: {T90_BLANKET} s", fontsize=8, color="0.3")
ax.set_xlabel("radius (mm)"); ax.set_ylabel("time to 90 % coverage (s)"); ax.set_xlim(0, 150)
ax.set_ylim(0, max(18.5, np.nanmax([res[k]["t90_bottom"].max() for k in CASES]) + 4))
ax.set_title("Saturation time across the wafer, 1 % TMA, 16 s pulse", fontsize=10); ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.25)
ax = axes[1]
for k, (cd, lab) in CASES.items():
    R = res[k]
    ax.plot(r_mm, R["wall_local"], "-", color=COL[k], lw=2, label=lab)
ax.axhline(res["uniform"]["wall_ratio"] if "uniform" in res else 1.78e13 * np.pi * 100e-9 * 100 * 100e-9,
           color="0.45", ls=":", lw=1.2)
ax.text(2, 57, "uniform AR-100 wafer: 55.9", fontsize=8, color="0.3")
for k in CASES:
    ax.text(103, 26 - 4 * list(CASES).index(k), f"{k.replace('_', ' ')}: wafer average {res[k]['wall_ratio']:.1f}, rule {res[k]['rule']:.2f} s",
            fontsize=8, color=COL[k])
ax.set_xlabel("radius (mm)"); ax.set_ylabel("local wall area / wafer area  (holeDensity x pi d x AR d)"); ax.set_xlim(0, 150); ax.set_ylim(0, 62)
ax.set_title("The two patterns as the supply rule sees them", fontsize=10); ax.legend(fontsize=7.5, loc="lower left"); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(OUT / "pattern_layout.png", dpi=150)
print("saved", OUT / "pattern_layout.png", "and", OUT / "pattern_layout.csv")
