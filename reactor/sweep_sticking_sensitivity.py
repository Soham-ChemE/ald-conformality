"""Sticking-coefficient sensitivity of the saturating pulse at the baseline operating point (300 Pa, 573 K, 300 sccm,
1 % TMA, 15 mm gap): blanket wafer and the 3D-NAND-like AR-100 patterned wafer (wall-area ratio 55.9) at
s = 5e-4, 1e-3, 2e-3, 5e-3, 8.5e-3. The range spans the Arts et al. 2019 lumped-model estimate (0.5-2e-3) and the
PillarHall fit of coupling/fit_pillarhall.py (8.5e-3); 5e-3 is the value used in every other chamber run.

Cases are copies of the two templates with only the wafer BC line, the pulse table and endTime/writeInterval changed:
  patterned: reactor/showerhead_pulse_patterned_AR100_long  -> 24 s pulse, endTime 26, writes every 0.05 s
             (the existing long case, 20 s pulse / endTime 22, is reused as the s = 5e-3 point)
  blanket:   reactor/showerhead_pulse_langmuir_p0.5          -> 8 s pulse,  endTime 9,  writes every 0.01 s
Each run is one scalarTransportFoam call through the OpenFOAM launcher (no shell spawned inside the environment).
Re-runnable: a case whose final time directory exists is not rebuilt or rerun.
  .venv/bin/python reactor/sweep_sticking_sensitivity.py
Outputs: reactor/sticking/<case>/ (fields, log, runtime.txt), figures/sticking_sensitivity.csv, figures/sticking_sensitivity.png.

Post-processing per case: t90 = first write time at which coverage >= 0.9 (resolution one write interval);
hole bottoms and flat surface of the patterned wafer at the centre face and the worst face over faces 0-148
(face 149 is the outermost half-face next to the bare floor and is excluded, as in the rest of the repository);
t_supply = (1 + A_wall/A_wafer) N0/(U c0) with U = 0.0501 m/s read from the frozen flow (11.99 s patterned,
0.2107 s blanket); tau_sat = t90_worst / t_supply; utilization up to saturation = 1/tau_sat
(site-equivalents reacted / site-equivalents delivered, as in reactor/doe_patterned.py and coupling/validate_supply_rule.py).
"""
import re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT), str(ROOT.parent / "coupling")]
from patterned_lib import load_patterned
from postprocess_langmuir import load as load_flat

OUTDIR = ROOT / "sticking"; OUTDIR.mkdir(exist_ok=True)
FIG = ROOT.parent / "figures"
LAUNCHER = Path.home() / "Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam"
PAT_SRC = ROOT / "showerhead_pulse_patterned_AR100_long"
BLK_SRC = ROOT / "showerhead_pulse_langmuir_p0.5"
S_PAT = (5e-4, 1e-3, 2e-3, 8.5e-3)              # 5e-3 is PAT_SRC itself
S_BLK = (5e-4, 1e-3, 2e-3, 5e-3, 8.5e-3)
PAT_PULSE, PAT_END, PAT_WI = 24.0, 26.0, 0.05
BLK_PULSE, BLK_END, BLK_WI = 8.0, 9.0, 0.01
U, N0, C0, D_HOLE, AR, HOLES = 0.0501, 4e18, 3.79e20, 100e-9, 100, 1.78e13
WALL = HOLES * np.pi * D_HOLE * AR * D_HOLE                 # 55.9
T_MONO = N0 / (U * C0)                                     # 0.2107 s
T_SUPPLY_PAT = (1 + WALL) * T_MONO                         # 11.99 s
S_ARTS = (0.5e-3, 2e-3); S_PILLARHALL = 8.5e-3
COMMON = "r 0; vbar 410; N0 4e18; c0 3.79e20;"
NFACE = 149


def build(name, src, kind, s, pulse, end, wi):
    c = OUTDIR / name
    if (c / f"{end:g}" / "T").exists():
        return c, False
    if c.exists(): shutil.rmtree(c)
    c.mkdir()
    for sub in ("0", "constant", "system"): shutil.copytree(src / sub, c / sub)
    if kind == "patterned":
        pat = r"    wafer      \{ type patternedSurface;[^\n]*\n"
        bc = f"type patternedSurface; s {s:g}; {COMMON} d {D_HOLE:g}; AR {AR}; holeDensity {HOLES:g}; nCells 200; value uniform 0;"
    else:
        pat = r"    wafer      \{ type langmuirSticking;[^\n]*\n"
        bc = f"type langmuirSticking; s {s:g}; {COMMON} value uniform 0;"
    t = c / "0/T"; txt = t.read_text()
    txt, n1 = re.subn(pat, "    wafer      { " + bc + " }\n", txt)
    txt, n2 = re.subn(r"table \(\(0 1\) \([0-9.]+ 1\) \([0-9.e-]+ 0\) \([0-9.]+ 0\)\)",
                      f"table ((0 1) ({pulse:g} 1) ({pulse + 1e-4:g} 0) (100 0))", txt)
    txt, n3 = re.subn(r"// 0\.1 s pulse at unit concentration, then purge",
                      f"// {pulse:g} s pulse at unit concentration, then purge to {end:g} s; sticking sensitivity s = {s:g}", txt)
    assert n1 == n2 == n3 == 1, (name, n1, n2, n3); t.write_text(txt)
    cd = c / "system/controlDict"; d = cd.read_text()
    d, m1 = re.subn(r"endTime\s+[0-9.]+;[^\n]*", f"endTime         {end:g};        // {pulse:g} s pulse + {end - pulse:g} s purge", d)
    d, m2 = re.subn(r"writeInterval\s+[0-9.]+;", f"writeInterval   {wi:g};", d)
    assert m1 == m2 == 1, (name, m1, m2); cd.write_text(d)
    return c, True


def run(c):
    # one command through the launcher; the bundled bash segfaults when a shell is spawned inside the environment
    t0 = time.time()
    rc = subprocess.run([str(LAUNCHER), "-c", f"cd '{c}' && scalarTransportFoam > log.scalarTransportFoam 2>&1"]).returncode
    dt = time.time() - t0
    (c / "runtime.txt").write_text(f"scalarTransportFoam exit {rc}, wall time {dt:.0f} s\n")
    return rc, dt


def t90(arr, j, t):
    col = arr[:, j]
    return float(t[np.argmax(col >= 0.9)]) if col.max() >= 0.9 else float("nan")


def t90_profile(arr, t):
    return np.array([t90(arr, j, t) for j in range(NFACE)])


def runtime_of(c):
    f = c / "runtime.txt"
    if f.exists():
        m = re.search(r"wall time (\d+) s", f.read_text()); return float(m.group(1)) if m else np.nan
    f = c / "time.txt"                                       # the reused long case was timed with /usr/bin/time
    if f.exists():
        m = re.search(r"real\s+([0-9.]+)", f.read_text()); return float(m.group(1)) if m else np.nan
    return np.nan


rows = []
# ---- blanket wafer ----
for s in S_BLK:
    c, fresh = build(f"blanket_s{s:g}", BLK_SRC, "blanket", s, BLK_PULSE, BLK_END, BLK_WI)
    if fresh:
        print(f"running {c.name} ...", flush=True); rc, dt = run(c); print(f"{c.name}: exit {rc}, {dt:.0f} s", flush=True)
    t, th, face, s_read, vbar, c0, n0, pulse = load_flat(c)
    assert abs(s_read - s) < 1e-12 * max(1, s) and abs(c0 - C0) < 1e12 and n0 == N0, (c.name, s_read, c0, n0)
    prof = t90_profile(th, t)
    rows.append(dict(kind="blanket", s=s, case=str(c.relative_to(ROOT.parent)), pulse=pulse, end_time=float(t[-1]),
                     t_supply=T_MONO, wall_ratio=0.0,
                     t90_centre=prof[0], t90_worst=np.nanmax(prof), t90_worst_face_mm=float(np.nanargmax(prof) + 0.5),
                     t90_flat_centre=np.nan, t90_flat_worst=np.nan,
                     theta_end_centre=float(th[-1, 0]), theta_min_end=float(th[-1, :NFACE].min()),
                     wall_time_s=runtime_of(c)))
    r = rows[-1]; r["tau_sat"] = r["t90_worst"] / r["t_supply"]; r["utilization"] = 1 / r["tau_sat"]
    r["saturated_in_pulse"] = bool(r["t90_worst"] <= pulse)
    print(f"blanket s={s:g}: t90 centre {r['t90_centre']:.3f} s, worst {r['t90_worst']:.3f} s (r = {r['t90_worst_face_mm']:.1f} mm), "
          f"tau_sat {r['tau_sat']:.2f}, utilization {100*r['utilization']:.1f} %, theta(end) centre {r['theta_end_centre']:.4f}", flush=True)

# ---- patterned AR-100 wafer ----
pat_cases = [(PAT_SRC, 5e-3)] + [(build(f"patterned_AR100_s{s:g}", PAT_SRC, "patterned", s, PAT_PULSE, PAT_END, PAT_WI), s) for s in S_PAT]
for c, s in pat_cases:
    if isinstance(c, tuple):
        c, fresh = c
        if fresh:
            print(f"running {c.name} ...", flush=True); rc, dt = run(c); print(f"{c.name}: exit {rc}, {dt/60:.1f} min", flush=True)
    t, f, prm, pulse = load_patterned(c)
    assert abs(prm["s"] - s) < 1e-12 * max(1, s) and prm["AR"] == AR and abs(prm["holeDensity"] - HOLES) < 1e9, (c.name, prm)
    pb, pf = t90_profile(f["thetaBottom"], t), t90_profile(f["thetaFlat"], t)
    rows.append(dict(kind="patterned", s=s, case=str(c.relative_to(ROOT.parent)), pulse=pulse, end_time=float(t[-1]),
                     t_supply=T_SUPPLY_PAT, wall_ratio=WALL,
                     t90_centre=pb[0], t90_worst=np.nanmax(pb), t90_worst_face_mm=float(np.nanargmax(pb) + 0.5),
                     t90_flat_centre=pf[0], t90_flat_worst=np.nanmax(pf),
                     theta_end_centre=float(f["thetaBottom"][-1, 0]), theta_min_end=float(f["thetaBottom"][-1, :NFACE].min()),
                     wall_time_s=runtime_of(c)))
    r = rows[-1]; r["tau_sat"] = r["t90_worst"] / r["t_supply"]; r["utilization"] = 1 / r["tau_sat"]
    r["saturated_in_pulse"] = bool(r["t90_worst"] <= pulse)
    print(f"patterned s={s:g}: bottom t90 centre {r['t90_centre']:.2f} s, worst {r['t90_worst']:.2f} s (r = {r['t90_worst_face_mm']:.1f} mm), "
          f"flat {r['t90_flat_centre']:.2f} s (worst {r['t90_flat_worst']:.2f}), tau_sat {r['tau_sat']:.3f}, "
          f"utilization {100*r['utilization']:.1f} %, bottom theta(end) centre {r['theta_end_centre']:.4f}", flush=True)

df = pd.DataFrame(rows).sort_values(["kind", "s"]).reset_index(drop=True)
df["t90_centre_over_supply"] = df.t90_centre / df.t_supply
df["t90_worst_over_supply"] = df.t90_worst / df.t_supply
df["t90_flat_centre_over_supply"] = df.t90_flat_centre / df.t_supply
df.to_csv(FIG / "sticking_sensitivity.csv", index=False)
print(df.to_string())

# ---- figure ----
b, p = df[df.kind == "blanket"], df[df.kind == "patterned"]
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.6))
for a in (ax, ax2):
    a.axvspan(*S_ARTS, color="tab:green", alpha=0.12, lw=0, label="Arts et al. 2019 (JVST A 37, 030908), LHAR front slope, Tset 275 C: 0.5-2e-3")
    a.axvline(5.72e-3, color="tab:purple", ls="-.", lw=0.8, label="Ylilammi et al. 2018 (JAP 123, 205301), LHAR fit, 300 C: 5.72e-3")
    a.axvline(4e-3, color="tab:olive", ls="-.", lw=0.8, label="Vandalon and Kessels 2016 (APL 108, 011607), SFG, 300 C: 4e-3 (+/- 1e-3)")
    a.axvline(S_PILLARHALL, color="tab:red", ls="--", lw=1, label="PillarHall fit (this work, coupling/fit_pillarhall.py), 300 C: 8.5e-3")
    a.axvline(5e-3, color="0.5", ls=":", lw=1, label="value used in the chamber runs: 5e-3")
ax.plot(p.s, p.t90_worst, "s-", color="tab:blue", ms=7, label="patterned AR-100 wafer, hole bottoms (worst face)")
ax.plot(p.s, p.t90_centre, "s", mfc="none", color="tab:blue", ms=7, label="patterned AR-100 wafer, hole bottoms (centre)")
ax.plot(p.s, p.t90_flat_centre, "^-", color="tab:orange", ms=7, label="patterned wafer, flat surface between holes")
ax.plot(b.s, b.t90_worst, "o-", color="tab:gray", ms=7, label="blanket wafer (worst face)")
ax.axhline(T_SUPPLY_PAT, color="tab:blue", lw=0.8, ls="-.", label=f"supply-limit time, patterned: (1 + {WALL:.1f}) N0/(U c0) = {T_SUPPLY_PAT:.2f} s")
ax.axhline(T_MONO, color="tab:gray", lw=0.8, ls="-.", label=f"supply-limit time, blanket: N0/(U c0) = {T_MONO:.3f} s")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("TMA sticking coefficient s"); ax.set_ylabel("time to 90 % coverage (s)")
ax.set_title("Saturating pulse vs sticking coefficient, 1 % TMA, 300 sccm, 300 Pa, 573 K", fontsize=9)
ax.set_ylim(0.02, 40); ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=6, loc="lower left", ncol=2)
ax2.axhspan(1.05, 1.27, color="tab:blue", alpha=0.10, lw=0, label="pulse-time map, 12 patterned cases at s = 5e-3: 1.05-1.27")
ax2.plot(p.s, p.t90_worst_over_supply, "s-", color="tab:blue", ms=7, label="patterned AR-100 hole bottoms (worst) / 11.99 s")
ax2.plot(p.s, p.t90_centre_over_supply, "s", mfc="none", color="tab:blue", ms=7, label="patterned AR-100 hole bottoms (centre) / 11.99 s")
ax2.plot(b.s, b.t90_worst_over_supply, "o-", color="tab:gray", ms=7, label=f"blanket wafer (worst) / {T_MONO:.3f} s")
ax2.axhline(1, color="0.6", lw=0.8)
ax2.set_xscale("log"); ax2.set_yscale("log"); ax2.set_xlabel("TMA sticking coefficient s"); ax2.set_ylabel("t90 / supply-limit time")
ax2.set_title("Saturation time relative to the supply limit", fontsize=9)
ax2.set_ylim(0.9, 12); ax2.grid(True, which="both", alpha=0.3); ax2.legend(fontsize=6, loc="upper right")
for a in (ax, ax2):                                        # tick at the five run values only (after the log scale is set)
    a.set_xticks(list(S_BLK)); a.set_xticklabels([f"{v:g}" for v in S_BLK]); a.set_xticks([], minor=True); a.set_xlim(4e-4, 1.05e-2)
fig.tight_layout(); fig.savefig(FIG / "sticking_sensitivity.png", dpi=150); print("saved", FIG / "sticking_sensitivity.png")
