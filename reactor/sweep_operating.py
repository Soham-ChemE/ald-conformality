"""One-factor-at-a-time sweep of the reactor operating point around the baseline
(300 Pa, 300 sccm N2, 15 mm gap, 1 % TMA = 3 Pa), with the self-limiting wafer (s = 5e-3, N0 = 4e18).

For each case: mesh (gap sets the axial cell count at 0.5 mm), steady laminar N2 flow (nu and inlet speed from p, Q),
then a 1.0 s TMA pulse + 1.0 s purge with the langmuirSticking wall (D and c0 from p, x) and a Danckwerts inlet whose
value fraction is recomputed from U, D and the wall distance. Metrics: time to 90 % coverage at centre / r = 140 mm /
worst on wafer, TMA capture fraction, final coverage. Output: figures/reactor_sweep_operating.csv and .png.

Run inside the OpenFOAM environment:
  ~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam -c '.venv/bin/python reactor/sweep_operating.py'
"""
import json, re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SWEEP = ROOT / "sweep"
OUT = ROOT.parent / "figures"
TPL_STEADY = ROOT / "showerhead_steady"
TPL_PULSE = ROOT / "showerhead_pulse_langmuir_p0.5"

KB, R, T = 1.380649e-23, 8.314, 573.15
M_N2, MU_N2 = 0.028, 2.9e-5
R_WAFER = 0.150
A_SH = np.pi * R_WAFER ** 2
D_REF, P_REF = 8.4e-3, 300.0           # TMA in N2 at 300 Pa, 573 K
S, N0, VBAR = 5e-3, 4e18, 410.0
PULSE, PURGE = 1.0, 1.0

BASE = dict(p=300.0, Q=300.0, g=15.0, x=0.01)
CASES = {"base": BASE}
for p in (100.0, 1000.0): CASES[f"p{int(p)}"] = {**BASE, "p": p}
for Q in (150.0, 600.0): CASES[f"Q{int(Q)}"] = {**BASE, "Q": Q}
for g in (10.0, 25.0): CASES[f"g{int(g)}"] = {**BASE, "g": g}
for x in (0.005, 0.03): CASES[f"x{x:g}"] = {**BASE, "x": x}

def derived(c):
    p, Q, g, x = c["p"], c["Q"], c["g"], c["x"]
    rho = p * M_N2 / (R * T); nu = MU_N2 / rho
    Qc = Q * 1.6667e-8 * (101325.0 / p) * (T / 273.15)      # m3/s at chamber conditions
    U = Qc / A_SH
    D = D_REF * P_REF / p
    c0 = x * p / (KB * T)
    nx = int(round(g / 0.5)); delta = (g * 1e-3 / nx) / 2
    f_in = U / (U + D / delta)
    return dict(rho=rho, nu=nu, U=U, D=D, c0=c0, nx=nx, delta=delta, f_in=f_in, Pe=U * g * 1e-3 / D,
                supply=U * c0, t_mono=N0 / (U * c0))

def sub(path, pairs):
    s = path.read_text()
    for old, new in pairs:
        assert old in s, (path, old)
        s = s.replace(old, new)
    path.write_text(s)

def build(name, c, d):
    case = SWEEP / name
    st, pu = case / "steady", case / "pulse"
    steady_done = st.exists() and any(re.fullmatch(r"[1-9][0-9]*", q.name) for q in st.iterdir())
    if not steady_done:
        if st.exists(): shutil.rmtree(st)
        shutil.copytree(TPL_STEADY, st, ignore=shutil.ignore_patterns("[0-9]*", "log.*", "postProcessing", "dynamicCode"))
        (st / "0").mkdir(exist_ok=True)
        for f in ("U", "p"): shutil.copy(TPL_STEADY / "0" / f, st / "0" / f)
        sub(st / "system/blockMeshDict", [("G  15;", f"G  {c['g']:g};"), ("(30 150 1)", f"({d['nx']} 150 1)"), ("(30 50 1)", f"({d['nx']} 50 1)")])
        sub(st / "constant/transportProperties", [("nu              1.65e-2;", f"nu              {d['nu']:.5g};")])
        sub(st / "0/U", [("uniform (-0.050 0 0)", f"uniform (-{d['U']:.6g} 0 0)")])
    if pu.exists(): shutil.rmtree(pu)
    pu.mkdir(); (pu / "0").mkdir(); (pu / "constant").mkdir()
    shutil.copytree(TPL_PULSE / "system", pu / "system")
    shutil.copy(TPL_PULSE / "0/T", pu / "0/T")
    shutil.copy(TPL_PULSE / "constant/transportProperties", pu / "constant/transportProperties")
    sub(pu / "constant/transportProperties", [("DT              8.4e-3;", f"DT              {d['D']:.5g};")])
    sub(pu / "0/T", [("uniformValueFraction constant 1.486e-3;", f"uniformValueFraction constant {d['f_in']:.6g};"),
                     ("c0 3.79e20;", f"c0 {d['c0']:.4g};"),
                     ("table ((0 1) (0.5 1) (0.5001 0) (10 0))", f"table ((0 1) ({PULSE} 1) ({PULSE + 1e-4} 0) (10 0))")])
    sub(pu / "system/controlDict", [("endTime         1.5;", f"endTime         {PULSE + PURGE};"), ("writeInterval   0.005;", "writeInterval   0.01;")])
    # drop the useless patch sampler
    cd = (pu / "system/controlDict").read_text()
    cd = re.sub(r"functions\s*\{.*\}\s*$", "functions {}\n", cd, flags=re.S)
    (pu / "system/controlDict").write_text(cd)
    return case

def _foam(cmd, cwd):
    """Run one OpenFOAM utility directly (no shell: the app's bundled bash segfaults when spawned from here)."""
    with open(cwd / f"log.{cmd[0]}", "w") as log:
        return subprocess.run(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT).returncode

def run(case):
    t0 = time.time()
    st, pu = case / "steady", case / "pulse"
    if not any(re.fullmatch(r"[1-9][0-9]*", q.name) for q in st.iterdir()):
        for cmd in (["blockMesh"], ["simpleFoam"]):
            rc = _foam(cmd, st)
            if rc: return rc, time.time() - t0, f"{cmd[0]} failed, see {st / ('log.' + cmd[0])}"
    last = max((q for q in st.iterdir() if re.fullmatch(r"[1-9][0-9]*", q.name)), key=lambda q: int(q.name))
    if (pu / "constant/polyMesh").exists(): shutil.rmtree(pu / "constant/polyMesh")
    shutil.copytree(st / "constant/polyMesh", pu / "constant/polyMesh")
    for f in ("U", "phi"): shutil.copy(last / f, pu / "0" / f)
    rc = _foam(["scalarTransportFoam"], pu)
    return rc, time.time() - t0, ("" if rc == 0 else f"scalarTransportFoam failed, see {pu / 'log.scalarTransportFoam'}")

sys.path.insert(0, str(ROOT))
from postprocess_langmuir import load as load_cov

def metrics(case, d):
    t, theta, face, s, vbar, c0, N0_, pulse = load_cov(case / "pulse")
    def t90(j): 
        k = np.argmax(theta[:, j] >= 0.9); return t[k] if theta[:, j].max() >= 0.9 else np.nan
    t90r = np.array([t90(j) for j in range(150)])
    uptake = (vbar / 4) * s * (1 - theta) * face * c0
    dt = np.diff(t, prepend=0)
    r = (np.arange(150) + 0.5) * 1e-3
    absorbed = 2 * np.pi * np.sum((uptake * dt[:, None]).sum(0) * r) * 1e-3
    injected = d["U"] * c0 * A_SH * PULSE
    return dict(t90_centre=float(t90r[0]), t90_r140=float(t90r[140]), t90_worst=float(np.nanmax(t90r[:149])),
                t90_spread_pct=float((np.nanmax(t90r[:149]) - t90r[0]) / t90r[0] * 100) if t90r[0] == t90r[0] else np.nan,
                theta_final_centre=float(theta[-1, 0]), captured_pct=float(absorbed / injected * 100),
                theta_at_0p1s=float(theta[np.argmin(abs(t - 0.1)), 0]))

if __name__ == "__main__":
    SWEEP.mkdir(exist_ok=True)
    rows = []
    for name, c in CASES.items():
        d = derived(c)
        case = build(name, c, d)
        rc, dt_, err = run(case)
        print(f"{name:8s} p={c['p']:g} Q={c['Q']:g} g={c['g']:g} x={c['x']:g} | U={d['U']:.4f} D={d['D']:.2e} Pe={d['Pe']:.3f} "
              f"t_mono={d['t_mono']:.3f}s | rc={rc} {dt_:.0f}s", flush=True)
        if rc != 0:
            print(err, flush=True); continue
        m = metrics(case, d)
        rows.append({**c, **{k: d[k] for k in ("U", "D", "Pe", "t_mono", "f_in")}, **m, "name": name})
        print("   ", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in m.items()}, flush=True)
    import pandas as pd
    df = pd.DataFrame(rows).set_index("name")
    df.to_csv(OUT / "reactor_sweep_operating.csv")
    print(df[["p", "Q", "g", "x", "Pe", "t_mono", "t90_centre", "t90_worst", "t90_spread_pct", "captured_pct", "theta_at_0p1s"]].to_string())
