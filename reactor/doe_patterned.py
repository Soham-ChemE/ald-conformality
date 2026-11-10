"""Face-centred central composite design over the chamber operating point for a 3D-NAND-like AR-100 wafer, with the
two-way multiscale wall. Factors (coded -1..+1): log10 TMA fraction x (0.5 % .. 3 %), log10 N2 flow Q (150 .. 600 sccm),
log10 pressure p (100 .. 1000 Pa), gap g (10 .. 25 mm). 16 factorial + 8 axial + 1 centre = 25 runs.
Responses per run: t90 of the hole bottoms at the centre and the worst point (saturating pulse), edge lag
(t90 at r = 140 mm vs centre), flat-region t90, precursor utilization over the saturating pulse (site-equivalents
needed / molecules supplied up to t90_worst), and the supply-limit time. Resumable: finished cases are skipped.
Run inside the OpenFOAM environment:
  ~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam -c '.venv/bin/python reactor/doe_patterned.py'
"""
import itertools, re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT), str(ROOT.parent / "coupling")]
import sweep_operating as so
from patterned_lib import load_patterned
DOE = ROOT / "doe"; DOE.mkdir(exist_ok=True)
OUT = ROOT.parent / "figures"
TPL_PULSE = ROOT / "showerhead_pulse_langmuir_p0.5"
N0, S, HOLES, D_HOLE, AR = 4e18, 5e-3, 1.78e13, 100e-9, 100
WALL = HOLES * np.pi * D_HOLE * AR * D_HOLE          # 55.9
KB, T = 1.380649e-23, 573.15

LO = dict(x=np.log10(0.005), Q=np.log10(150.0), p=np.log10(100.0), g=10.0)
HI = dict(x=np.log10(0.03), Q=np.log10(600.0), p=np.log10(1000.0), g=25.0)
def decode(c):
    out = {}
    for k in ("x", "Q", "p", "g"):
        v = LO[k] + (c[k] + 1) / 2 * (HI[k] - LO[k])
        out[k] = v if k == "g" else 10 ** v
    return out
runs = []
for f in itertools.product((-1, 1), repeat=4): runs.append(dict(zip(("x", "Q", "p", "g"), f)))
for k in ("x", "Q", "p", "g"):
    for a in (-1, 1):
        c = dict(x=0, Q=0, p=0, g=0); c[k] = a; runs.append(c)
runs.append(dict(x=0, Q=0, p=0, g=0))
design = pd.DataFrame([{**{f"c_{k}": v for k, v in c.items()}, **decode(c)} for c in runs])
design.index.name = "run"; design.to_csv(OUT / "doe_design.csv")

def build_pulse(case, c, d):
    pu = case / "pulse"
    if pu.exists(): shutil.rmtree(pu)
    pu.mkdir(); (pu / "0").mkdir(); (pu / "constant").mkdir()
    shutil.copytree(TPL_PULSE / "system", pu / "system")
    shutil.copy(TPL_PULSE / "0/T", pu / "0/T")
    shutil.copy(TPL_PULSE / "constant/transportProperties", pu / "constant/transportProperties")
    so.sub(pu / "constant/transportProperties", [("DT              8.4e-3;", f"DT              {d['D']:.5g};")])
    t_supply = (1 + WALL) * N0 / (d["U"] * d["c0"])
    pulse = round(1.6 * t_supply + 0.5, 2)
    t = pu / "0/T"; s = t.read_text()
    s, n1 = re.subn(r"    wafer      \{ type langmuirSticking;[^\n]*\n",
        f"    wafer      {{ type patternedSurface; s {S:g}; r 0; vbar 410; N0 {N0:g}; c0 {d['c0']:.4g}; d {D_HOLE:g}; AR {AR}; holeDensity {HOLES:g}; nCells 150; value uniform 0; }}\n", s)
    s, n2 = re.subn(r"uniformValueFraction constant [0-9.eE+-]+;", f"uniformValueFraction constant {d['f_in']:.6g};", s)
    s, n3 = re.subn(r"table \(\(0 1\) \([0-9.]+ 1\) \([0-9.e-]+ 0\) \([0-9.]+ 0\)\)", f"table ((0 1) ({pulse} 1) ({pulse + 1e-4} 0) (1000 0))", s)
    assert n1 == n2 == n3 == 1, (case, n1, n2, n3); t.write_text(s)
    cd = pu / "system/controlDict"; dd = cd.read_text()
    dd, m1 = re.subn(r"endTime\s+[0-9.]+;", f"endTime         {pulse + 1.0};", dd)
    dd, m2 = re.subn(r"writeInterval\s+[0-9.]+;", f"writeInterval   {max(0.005, round(pulse / 150, 3))};", dd)
    dd = re.sub(r"functions\s*\{.*\}\s*$", "functions {}\n", dd, flags=re.S)
    assert m1 == m2 == 1; cd.write_text(dd)
    return pulse, t_supply

def t90(arr, j, t):
    k = np.argmax(arr[:, j] >= 0.9); return float(t[k]) if arr[:, j].max() >= 0.9 else float("nan")

rows = []
for i, c in enumerate(runs):
    dec = decode(c); d = so.derived(dec)
    case = DOE / f"run{i:02d}"
    done = (case / "pulse" / "log.scalarTransportFoam").exists() and "End" in (case / "pulse" / "log.scalarTransportFoam").read_text()[-200:]
    t0 = time.time()
    if not done:
        # steady flow (reuse sweep_operating's builder; it keeps a finished steady run)
        so.SWEEP = DOE
        so.build(f"run{i:02d}", dec, d)                 # builds steady (if needed) and a Langmuir pulse we then replace
        pulse, t_supply = build_pulse(case, c, d)
        rc, dt, err = so.run(case)
        if rc: print(f"run{i:02d} FAILED: {err}", flush=True); continue
    else:
        t_supply = (1 + WALL) * N0 / (d["U"] * d["c0"]); pulse = round(1.6 * t_supply + 0.5, 2)
    t, f, prm, _ = load_patterned(case / "pulse")
    tb = np.array([t90(f["thetaBottom"], j, t) for j in range(150)])
    t90c, t90w, t90e = tb[0], np.nanmax(tb[:149]), tb[140]
    util = (1 + WALL) * N0 / (d["U"] * d["c0"] * t90w) if t90w == t90w else np.nan
    rows.append({**{f"c_{k}": v for k, v in c.items()}, **dec, "U": d["U"], "D": d["D"], "Pe": d["Pe"], "t_supply": t_supply, "pulse": pulse,
                 "t90_centre": t90c, "t90_worst": t90w, "edge_lag_pct": (t90e / t90c - 1) * 100, "t90_flat": t90(f["thetaFlat"], 0, t),
                 "utilization": util, "tau_sat": t90w / t_supply, "run": i})
    print(f"run{i:02d} x={dec['x']*100:.2f}% Q={dec['Q']:.0f} p={dec['p']:.0f} g={dec['g']:.1f} | t_supply {t_supply:.2f} s, t90 worst {t90w:.2f} s (tau {t90w/t_supply:.2f}), edge lag {(t90e/t90c-1)*100:+.1f} %, util {util*100:.0f} % | {time.time()-t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "doe_results.csv", index=False)
print(pd.DataFrame(rows).round(3).to_string())
