"""Temperature sweep (item 3c): 200, 300, 400 C at fixed 300 Pa, 300 sccm, 1 % TMA, 15 mm gap, blanket and
patterned AR-100 wafers. Temperature enters through: carrier density/viscosity (nu ~ T^1.75/p), inlet speed
(U ~ T at fixed sccm), diffusivity (D ~ T^1.75/p), precursor density (c0 ~ 1/T at fixed fraction), thermal speed
(vbar ~ sqrt(T)) and the Danckwerts fraction. Supply U c0 is T-independent at fixed sccm and fraction, so the rule
predicts the same patterned pulse time at every T; the blanket wafer (kinetics-limited) should speed up with vbar.
Sticking is kept at 5e-3 (its own T dependence is a chemistry input we do not have for this range).
Builds steady + pulse cases under reactor/temperature/; run with run_temperature_cases.sh."""
import re, shutil, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); import sweep_operating as so
OUTDIR = ROOT / "temperature"; OUTDIR.mkdir(exist_ok=True)
KB, R, M_N2, MU_REF, T_REF = 1.380649e-23, 8.314, 0.028, 2.9e-5, 573.15
P, Q, X, G = 300.0, 300.0, 0.01, 15.0
N0, S, HOLES, D_HOLE, AR = 4e18, 5e-3, 1.78e13, 100e-9, 100
WALL = HOLES * np.pi * D_HOLE * AR * D_HOLE
for TC in (200, 300, 400):
    T = TC + 273.15
    rho = P * M_N2 / (R * T); mu = MU_REF * (T / T_REF) ** 0.7; nu = mu / rho
    Qc = Q * 1.6667e-8 * (101325.0 / P) * (T / 273.15); U = Qc / so.A_SH
    D = so.D_REF * (so.P_REF / P) * (T / T_REF) ** 1.75
    c0 = X * P / (KB * T); vbar = np.sqrt(8 * R * T / (np.pi * 0.07209))
    nx = 30; delta = (G * 1e-3 / nx) / 2; f_in = U / (U + D / delta)
    d = dict(rho=rho, nu=nu, U=U, D=D, c0=c0, nx=nx, delta=delta, f_in=f_in, Pe=U * G * 1e-3 / D, supply=U * c0, t_mono=N0 / (U * c0))
    for kind in ("blanket", "patterned"):
        name = f"T{TC}_{kind}"; case = OUTDIR / name
        if case.exists(): shutil.rmtree(case)
        so.SWEEP = OUTDIR; so.build(name, dict(p=P, Q=Q, g=G, x=X), d)   # steady case + a template pulse
        pu = case / "pulse"; t = pu / "0/T"; s = t.read_text()
        if kind == "blanket":
            bc = f"type langmuirSticking; s {S:g}; r 0; vbar {vbar:.0f}; N0 {N0:g}; c0 {c0:.4g}; value uniform 0;"; pulse = 1.5
        else:
            bc = f"type patternedSurface; s {S:g}; r 0; vbar {vbar:.0f}; N0 {N0:g}; c0 {c0:.4g}; d {D_HOLE:g}; AR {AR}; holeDensity {HOLES:g}; nCells 150; value uniform 0;"
            pulse = round(1.6 * (1 + WALL) * d["t_mono"] + 0.5, 2)
        s, n1 = re.subn(r"    wafer      \{ type langmuirSticking;[^\n]*\n", "    wafer      { " + bc + " }\n", s)
        s, n2 = re.subn(r"table \(\(0 1\) \([0-9.]+ 1\) \([0-9.e-]+ 0\) \([0-9.]+ 0\)\)", f"table ((0 1) ({pulse} 1) ({pulse + 1e-4} 0) (1000 0))", s)
        assert n1 == n2 == 1; t.write_text(s)
        cd = pu / "system/controlDict"; dd = cd.read_text()
        dd = re.sub(r"endTime\s+[0-9.]+;", f"endTime         {pulse + 1.0};", dd)
        dd = re.sub(r"writeInterval\s+[0-9.]+;", f"writeInterval   {max(0.005, round(pulse / 150, 3))};", dd)
        dd = re.sub(r"functions\s*\{.*\}\s*$", "functions {}\n", dd, flags=re.S); cd.write_text(dd)
        # steady transportProperties and U were written by so.build from d (nu, U); confirm
        print(f"{name}: T={T:.0f} K nu={nu:.3e} U={U:.4f} D={D:.3e} c0={c0:.3e} vbar={vbar:.0f} supply={U*c0:.3e} t_mono={d['t_mono']:.3f} pulse={pulse}")
