"""Pulse-time map for the patterned wafer: time for AR-20/50/100 hole bottoms to reach 90 % coverage versus TMA
fraction (0.5, 1, 3, 10 %), two-way multiscale wall, baseline chamber (300 Pa, 300 sccm, 15 mm). Plus blanket wafers.
Run inside the OpenFOAM environment:
  ~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam -c '.venv/bin/python reactor/sweep_pattern_map.py'
"""
import re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT), str(ROOT.parent / "coupling")]
from patterned_lib import load_patterned
from postprocess_langmuir import load as load_flat
BASE = ROOT / "sweep/base/pulse"
MAP = ROOT / "map"; MAP.mkdir(exist_ok=True)
OUT = ROOT.parent / "figures"
KB, T = 1.380649e-23, 573.15
P_TOT, U, N0, HOLES, D_HOLE = 300.0, 0.0501, 4e18, 1.78e13, 100e-9
XS = (0.005, 0.01, 0.03, 0.10); ARS = (20, 50, 100)

def build(name, bc, pulse, end, wi, c0):
    c = MAP / name
    if c.exists(): shutil.rmtree(c)
    c.mkdir()
    for sub in ("0", "constant", "system"): shutil.copytree(BASE / sub, c / sub)
    t = c / "0/T"; s = t.read_text()
    s, n1 = re.subn(r"    wafer      \{ type langmuirSticking;[^\n]*\n", "    wafer      { " + bc + " }\n", s)
    s, n2 = re.subn(r"table \(\(0 1\) \([0-9.]+ 1\) \([0-9.e-]+ 0\) \([0-9.]+ 0\)\)", f"table ((0 1) ({pulse} 1) ({pulse + 1e-4} 0) (1000 0))", s)
    assert n1 == n2 == 1, (name, n1, n2); t.write_text(s)
    cd = c / "system/controlDict"; d = cd.read_text()
    d, m1 = re.subn(r"endTime\s+[0-9.]+;", f"endTime         {end};", d)
    d, m2 = re.subn(r"writeInterval\s+[0-9.]+;", f"writeInterval   {wi};", d)
    assert m1 == m2 == 1; cd.write_text(d)
    return c

def run(c):
    with open(c / "log.scalarTransportFoam", "w") as log:
        t0 = time.time(); rc = subprocess.run(["scalarTransportFoam"], cwd=c, stdout=log, stderr=subprocess.STDOUT).returncode
    return rc, time.time() - t0

def t90(arr, j, t):
    k = np.argmax(arr[:, j] >= 0.9); return float(t[k]) if arr[:, j].max() >= 0.9 else float("nan")

rows = []
for x in XS:
    c0 = x * P_TOT / (KB * T)
    t_mono = N0 / (U * c0)
    # blanket wafer
    pulse = max(0.3, 6 * t_mono); end = pulse + 1.0; wi = max(0.002, pulse / 300)
    c = build(f"blanket_x{x:g}", f"type langmuirSticking; s 5e-3; r 0; vbar 410; N0 {N0:g}; c0 {c0:.4g}; value uniform 0;", pulse, end, wi, c0)
    rc, dt = run(c)
    t, th, face, *_ = load_flat(c)
    rows.append(dict(kind="blanket", x=x, AR=0, wall_ratio=0, t_supply=t_mono, pulse=pulse, t90_centre=t90(th, 0, t), t90_worst=np.nanmax([t90(th, j, t) for j in range(149)]), rc=rc))
    print(f"blanket x={x:g}: t_supply {t_mono:.3f} s, t90 {rows[-1]['t90_centre']:.3f} s (worst {rows[-1]['t90_worst']:.3f}), rc={rc}, {dt:.0f}s", flush=True)
    for AR in ARS:
        wall = HOLES * np.pi * D_HOLE * AR * D_HOLE
        t_supply = (1 + wall) * t_mono
        pulse = round(1.6 * t_supply + 0.5, 2); end = pulse + 1.0; wi = max(0.005, round(pulse / 200, 3))
        bc = f"type patternedSurface; s 5e-3; r 0; vbar 410; N0 {N0:g}; c0 {c0:.4g}; d {D_HOLE:g}; AR {AR}; holeDensity {HOLES:g}; nCells 200; value uniform 0;"
        c = build(f"patterned_AR{AR}_x{x:g}", bc, pulse, end, wi, c0)
        rc, dt = run(c)
        t, f, prm, _ = load_patterned(c)
        tb = [t90(f["thetaBottom"], j, t) for j in range(149)]
        rows.append(dict(kind="patterned", x=x, AR=AR, wall_ratio=wall, t_supply=t_supply, pulse=pulse,
                         t90_centre=tb[0], t90_worst=np.nanmax(tb), t90_flat=t90(f["thetaFlat"], 0, t), rc=rc))
        print(f"patterned AR={AR} x={x:g}: wall {wall:.0f}x, t_supply {t_supply:.2f} s, pulse {pulse} s -> bottom t90 {tb[0]:.2f} s (worst {rows[-1]['t90_worst']:.2f}), flat {rows[-1]['t90_flat']:.2f}, rc={rc}, {dt:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "pulse_time_map.csv", index=False)
print(pd.DataFrame(rows).to_string())
