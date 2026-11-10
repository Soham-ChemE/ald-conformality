"""Axisymmetric approximation of the Ultratech/Cambridge NanoTech Fiji F200 chamber of Gakis et al. (Chem. Eng. Res.
Des. 132, 795 (2018), doi 10.1016/j.cherd.2018.02.031; Chem. Eng. Sci. 195, 399 (2019), doi 10.1016/j.ces.2018.09.037):
steady Ar flow (simpleFoam) and TMA pulses on the frozen flow (scalarTransportFoam) with the langmuirSticking wafer.

Geometry. The papers give no chamber dimensions; everything below was measured on Fig. 5a of the 2018 paper (flow-field
cross-section, rendered at 300 dpi, red-arrow extents per pixel row) with the 200 mm substrate (dome top, 330 px) as the
scale, 1.65 px/mm, reading uncertainty about +/- 3 mm on every length:
  top inlet tube        radius 20 mm, 82 mm long above the cone apex plane
  conical ceiling       from r = 20 mm at 199 mm above the wafer plane to r = 128 mm at 35 mm above it
  wide chamber          radius 128 mm, from 35 mm above the wafer plane to 50 mm below it (the real taper, 33-69 mm
                        below, is replaced by a step at 50 mm)
  substrate chuck       cylinder r = 100 mm, 50 mm deep (the real chuck is a dome 124 mm deep: r 95 mm at 47 mm below
                        the wafer, 59 mm at 96 mm; its lower half is omitted)
  lower pumping duct    radius 106 mm, from 50 mm to 205 mm below the wafer plane; outlet on its floor
Approximations forced by axisymmetry (all stated in the deliverable): the three inlets (top 100 sccm, loading door 50,
side 30 with the reactants) become one axial inlet carrying 180 sccm and the TMA; the loading-door box, the off-axis
side-inlet tube and the recirculation driven by the pulse momentum are not represented; the carrier flow is frozen at
its steady value (the real TMA pulse is a slug of nearly pure TMA at up to 450 sccm for 40 ms, which raises the outlet
pressure by 18 % in their Fig. 9); isothermal gas at 573 K (their walls are at 543 K, substrate 573 K).

Conditions (their Table 1): Ar, 72 mTorr = 9.6 Pa (Pirani 50 mTorr corrected), substrate 300 C, 180 sccm at 24 C and
1 bar -> Q = 0.0603 m3/s at 9.6 Pa and 573 K (their pump speed 57.1 L/s at the 543 K outlet is the same number);
rho = 8.05e-5 kg/m3, mu_Ar(573 K) = 3.6e-5 Pa s, nu = 0.45 m2/s; inlet speed 48 m/s (Re 4.3 on the tube diameter,
their Re 2.05 at the side inlet; Mach 0.11); D_TMA = 8.4e-3 * 300/9.6 = 0.263 m2/s (the N2 value scaled 1/p, taken for
Ar as well). Mean free path 1.4 mm; first cell at the wafer 2 mm.
TMA pulse: their feeding-system model delivers 0.63 mg = 8.74e-6 mol per 25 ms valve opening with the flow-rate and
mole-fraction history of their Fig. 8 (digitised here as the SHAPE table, normalised); the 60 ms valve is not quantified
in the paper and is run at 1.5x and 2.4x the 25 ms amount (bounds: headspace depletion vs linear in valve time).
Wafer: langmuirSticking with s = 3.44e-3 (their fitted s0 = 0.004 times the 86 % reaction probability of the adsorbed
TMA at 300 C, Chem. Eng. Sci. 2019 Section 4.1.2), N0 = 4.1e18 /m2 (one TMA per Al deposited: their GPC 1.0 A/cycle at
their density 3500 kg/m3), vbar = 410 m/s. A sensitivity case uses this repository's s = 5e-3.

  .venv/bin/python reactor/validation_gakis/make_gakis_cases.py          # build steady + pulse cases
  .venv/bin/python reactor/validation_gakis/make_gakis_cases.py --run    # build, mesh, simpleFoam, then the pulses
Every OpenFOAM command goes through the launcher as one command (the app's bash segfaults if a shell is spawned inside
the environment). Post-process with reactor/validation_gakis/postprocess_gakis.py.
"""
import re, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TPL_STEADY = ROOT / "showerhead_steady"
TPL_PULSE = ROOT / "showerhead_pulse_langmuir_p0.5"
LAUNCHER = Path.home() / "Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam"

KB, R, NA = 1.380649e-23, 8.314, 6.02214e23
T = 573.15
P = 9.6                                   # Pa (72 mTorr)
M_AR, MU_AR = 0.039948, 3.6e-5
SCCM = 180.0                              # 100 top + 50 loading door + 30 side, MFC standard 24 C and 1 bar
Q = SCCM * 1e-6 / 60 * (1e5 / P) * (T / 297.15)   # m3/s at chamber conditions
RHO = P * M_AR / (R * T); NU = MU_AR / RHO
D_TMA = 8.4e-3 * 300.0 / P
S_EFF, S_REPO, N0, VBAR = 3.44e-3, 5e-3, 4.1e18, 410.0
N_TMA_25MS = 0.63e-3 / 72.09             # mol per 25 ms valve opening (their feeding-system model)

# geometry, mm (x axial, 0 at the wafer plane, positive toward the inlet)
R_TUBE, R_WAFER, R_DUCT, R_CH = 20.0, 100.0, 106.0, 128.0
X_APEX, X_EDGE, X_TUBE, X_STEP, X_OUT = 199.0, 35.0, 281.0, -50.0, -205.0
def ceiling(r): return X_APEX - (X_APEX - X_EDGE) * (r - R_TUBE) / (R_CH - R_TUBE)
X_C100, X_C106 = ceiling(R_WAFER), ceiling(R_DUCT)
THETA = 2.5                                # half wedge angle, deg

# normalised TMA molar-flow history read off Fig. 8 of the 2018 paper (flow rate x TMA mole fraction, peak = 1);
# reading uncertainty about 20 % on each ordinate, the integral is renormalised to the stated 0.63 mg so only the
# shape matters
SHAPE = [(0, 0), (0.005, 0.3), (0.012, 1.0), (0.020, 0.9), (0.030, 0.45), (0.040, 0.15), (0.050, 0.05),
         (0.060, 0.02), (0.10, 0.01), (0.30, 0.003), (0.3001, 0.0), (10, 0.0)]
SHAPE_INT = np.trapezoid([v for _, v in SHAPE], [t for t, _ in SHAPE])   # s

CASES = {   # name: (TMA amount relative to the 25 ms valve, s)
    "pulse_25ms": (1.0, S_EFF),
    "pulse_60ms_x1.5": (1.5, S_EFF),
    "pulse_60ms_x2.4": (2.4, S_EFF),
    "pulse_25ms_s5e-3": (1.0, S_REPO),
}
PULSE_END, DT, WRITE = 1.0, 4e-4, 0.01   # three residence times after the 0.3 s pulse tail; Courant 0.25 at the wafer

def wedge_pts(x, r):
    cy, sz = np.cos(np.radians(THETA)), np.sin(np.radians(THETA))
    return f"({x:g} {r*cy:.6g} {-r*sz:.6g})", f"({x:g} {r*cy:.6g} {r*sz:.6g})"

def blockmesh():
    # vertex table: (x, r); axis points (r = 0) are single vertices, others come as back/front pairs
    P2 = [(0, 0), (0, R_TUBE), (0, R_WAFER), (0, R_DUCT), (0, R_CH),                      # 0-4 wafer plane
          (X_APEX, 0), (X_APEX, R_TUBE), (X_C100, R_WAFER), (X_C106, R_DUCT), (X_EDGE, R_CH),  # 5-9 ceiling
          (X_TUBE, 0), (X_TUBE, R_TUBE),                                                   # 10-11 inlet
          (X_STEP, 0), (X_STEP, R_WAFER), (X_STEP, R_DUCT), (X_STEP, R_CH),                # 12-15 step plane
          (X_OUT, 0), (X_OUT, R_WAFER), (X_OUT, R_DUCT)]                                   # 16-18 outlet plane
    verts, idx = [], {}
    for i, (x, r) in enumerate(P2):
        if r == 0:
            idx[i] = (len(verts), len(verts)); verts.append(f"({x:g} 0 0)")
        else:
            b, f = wedge_pts(x, r); idx[i] = (len(verts), len(verts) + 1); verts += [b, f]
    def hexblock(p0, p1, p2, p3, n, grading="(1 1 1)"):
        # p0..p3 counter-clockwise in the (x, r) plane: p0 = (x low, r low), p1 = (x high, r low), p2 = (x high, r high), p3 = (x low, r high)
        b = [idx[p][0] for p in (p0, p1, p2, p3)]; f = [idx[p][1] for p in (p0, p1, p2, p3)]
        return f"    hex ({b[0]} {b[1]} {b[2]} {b[3]} {f[0]} {f[1]} {f[2]} {f[3]}) ({n[0]} {n[1]} 1) simpleGrading {grading}"
    NX_CONE, NX_TUBE, NX_STEP, NX_DUCT = 60, 20, 25, 40
    blocks = [
        hexblock(0, 5, 6, 1, (NX_CONE, 20), "(3 1 1)"),       # A: axis column under the inlet tube
        hexblock(5, 10, 11, 6, (NX_TUBE, 20)),                 # T: inlet tube
        hexblock(1, 6, 7, 2, (NX_CONE, 80), "(3 1 1)"),       # B: wafer region under the cone
        hexblock(2, 7, 8, 3, (NX_CONE, 6), "(3 1 1)"),        # C1: r 100..106 above the annulus
        hexblock(3, 8, 9, 4, (NX_CONE, 22), "(3 1 1)"),       # C2: r 106..128
        hexblock(13, 2, 3, 14, (NX_STEP, 6)),                  # D1: annulus beside the chuck, r 100..106
        hexblock(14, 3, 4, 15, (NX_STEP, 22)),                 # D2: r 106..128
        hexblock(16, 12, 13, 17, (NX_DUCT, 100)),              # E1: lower duct under the chuck
        hexblock(17, 13, 14, 18, (NX_DUCT, 6)),                # E2: lower duct r 100..106
    ]
    def face(p, q, r_, s_):   # quad from four (x,r) point ids on the back/front wedge faces, mixed as needed
        return "(" + " ".join(str(v) for v in (p, q, r_, s_)) + ")"
    B, F = (lambda i: idx[i][0]), (lambda i: idx[i][1])
    boundary = f"""
boundary
(
    inlet      {{ type patch; faces (({B(10)} {B(11)} {F(11)} {F(10)})); }}
    wafer      {{ type wall;  faces (({B(0)} {B(1)} {F(1)} {F(0)}) ({B(1)} {B(2)} {F(2)} {F(1)})); }}
    chuckSide  {{ type wall;  faces (({B(13)} {B(2)} {F(2)} {F(13)})); }}
    chuckBottom {{ type wall; faces (({B(12)} {B(13)} {F(13)} {F(12)})); }}
    ceiling    {{ type wall;  faces (({B(6)} {B(7)} {F(7)} {F(6)}) ({B(7)} {B(8)} {F(8)} {F(7)}) ({B(8)} {B(9)} {F(9)} {F(8)})); }}
    tubeWall   {{ type wall;  faces (({B(6)} {B(11)} {F(11)} {F(6)})); }}
    sideWall   {{ type wall;  faces (({B(4)} {B(9)} {F(9)} {F(4)}) ({B(15)} {B(4)} {F(4)} {F(15)})); }}
    stepFloor  {{ type wall;  faces (({B(14)} {B(15)} {F(15)} {F(14)})); }}
    ductWall   {{ type wall;  faces (({B(18)} {B(14)} {F(14)} {F(18)})); }}
    outlet     {{ type patch; faces (({B(16)} {B(17)} {F(17)} {F(16)}) ({B(17)} {B(18)} {F(18)} {F(17)})); }}
    axis       {{ type empty; faces (({B(0)} {B(5)} {B(5)} {B(0)}) ({B(5)} {B(10)} {B(10)} {B(5)}) ({B(16)} {B(12)} {B(12)} {B(16)})); }}
    front      {{ type wedge; faces (
        ({F(0)} {F(5)} {F(6)} {F(1)}) ({F(5)} {F(10)} {F(11)} {F(6)}) ({F(1)} {F(6)} {F(7)} {F(2)}) ({F(2)} {F(7)} {F(8)} {F(3)}) ({F(3)} {F(8)} {F(9)} {F(4)})
        ({F(13)} {F(2)} {F(3)} {F(14)}) ({F(14)} {F(3)} {F(4)} {F(15)}) ({F(16)} {F(12)} {F(13)} {F(17)}) ({F(17)} {F(13)} {F(14)} {F(18)})); }}
    back       {{ type wedge; faces (
        ({B(0)} {B(1)} {B(6)} {B(5)}) ({B(5)} {B(6)} {B(11)} {B(10)}) ({B(1)} {B(2)} {B(7)} {B(6)}) ({B(2)} {B(3)} {B(8)} {B(7)}) ({B(3)} {B(4)} {B(9)} {B(8)})
        ({B(13)} {B(14)} {B(3)} {B(2)}) ({B(14)} {B(15)} {B(4)} {B(3)}) ({B(16)} {B(17)} {B(13)} {B(12)}) ({B(17)} {B(18)} {B(14)} {B(13)})); }}
);
"""
    return ("/* Axisymmetric Fiji F200 approximation, 5 degree wedge; see make_gakis_cases.py for the measured dimensions. */\n"
            "FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }\n"
            "scale 0.001;\nvertices\n(\n" + "\n".join(f"    {v}   // {i}" for i, v in enumerate(verts)) + "\n);\nblocks\n(\n"
            + "\n".join(blocks) + "\n);\nedges ();\n" + boundary + "mergePatchPairs ();\n")

STEADY = HERE / "steady"
U_IN = Q / (np.pi * (R_TUBE * 1e-3) ** 2)
WALLS = ["wafer", "chuckSide", "chuckBottom", "ceiling", "tubeWall", "sideWall", "stepFloor", "ductWall"]

def build_steady():
    if STEADY.exists(): shutil.rmtree(STEADY)
    for sub in ("0", "constant", "system"): (STEADY / sub).mkdir(parents=True)
    (STEADY / "system/blockMeshDict").write_text(blockmesh())
    for f in ("controlDict", "fvSchemes", "fvSolution"): shutil.copy(TPL_STEADY / "system" / f, STEADY / "system" / f)
    for f in ("momentumTransport", "turbulenceProperties"): shutil.copy(TPL_STEADY / "constant" / f, STEADY / "constant" / f)
    (STEADY / "constant/transportProperties").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object transportProperties; }\n"
        f"// Ar at {P} Pa, {T:.0f} K: rho = {RHO:.3e} kg/m3, mu = {MU_AR:.2e} Pa s -> nu = {NU:.4g} m2/s\n"
        f"transportModel  Newtonian;\nnu              {NU:.5g};\n")
    wall_u = "\n".join(f"    {w:<12}{{ type noSlip; }}" for w in WALLS)
    (STEADY / "0/U").write_text(
        "FoamFile { version 2.0; format ascii; class volVectorField; object U; }\ndimensions [0 1 -1 0 0 0 0];\n"
        f"// {SCCM:g} sccm Ar at {P} Pa, {T:.0f} K: Q = {Q:.4g} m3/s through the r = {R_TUBE:g} mm top tube -> {U_IN:.3g} m/s downward (-x)\n"
        "internalField   uniform (0 0 0);\nboundaryField\n{\n"
        f"    inlet       {{ type fixedValue; value uniform (-{U_IN:.5g} 0 0); }}\n{wall_u}\n"
        "    outlet      { type pressureInletOutletVelocity; value uniform (0 0 0); }\n"
        "    axis        { type empty; }\n    front       { type wedge; }\n    back        { type wedge; }\n}\n")
    wall_p = "\n".join(f"    {w:<12}{{ type zeroGradient; }}" for w in WALLS)
    (STEADY / "0/p").write_text(
        "FoamFile { version 2.0; format ascii; class volScalarField; object p; }\ndimensions [0 2 -2 0 0 0 0];\n"
        "internalField   uniform 0;\nboundaryField\n{\n    inlet       { type zeroGradient; }\n" + wall_p + "\n"
        "    outlet      { type totalPressure; p0 uniform 0; }\n    axis        { type empty; }\n    front       { type wedge; }\n    back        { type wedge; }\n}\n")
    # the template controlDict runs 2000 iterations with residual control; keep it, write every 500
    return STEADY

def build_pulse(name, amount, s):
    c = HERE / name
    if c.exists(): shutil.rmtree(c)
    for sub in ("0", "constant", "system"): (c / sub).mkdir(parents=True)
    for f in ("fvSchemes", "fvSolution"): shutil.copy(TPL_PULSE / "system" / f, c / "system" / f)
    cd = (TPL_PULSE / "system/controlDict").read_text()
    cd = re.sub(r"endTime\s+[0-9.]+;[^\n]*", f"endTime         {PULSE_END:g};   // pulse tail ends at 0.3 s, residence time {0.0136/Q:.2f} s, three residence times later the gas has cleared", cd)
    cd = re.sub(r"deltaT\s+[0-9.e-]+;", f"deltaT          {DT:g};", cd)
    cd = re.sub(r"writeInterval\s+[0-9.]+;", f"writeInterval   {WRITE:g};", cd)
    cd = re.sub(r"functions\s*\{.*\}\s*$", "functions {}\n", cd, flags=re.S)
    (c / "system/controlDict").write_text(cd)
    (c / "constant/transportProperties").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary; object transportProperties; }\n"
        f"// TMA at {P} Pa, {T:.0f} K: 8.4e-3 m2/s (300 Pa, N2) scaled as 1/p; the Ar value is taken equal\nDT              {D_TMA:.4g};\n")
    n_mol = amount * N_TMA_25MS
    c0 = n_mol * NA / (Q * SHAPE_INT)                   # number density of TMA in the inlet stream at shape = 1
    dx_tube = (X_TUBE - X_APEX) * 1e-3 / 20; delta = dx_tube / 2
    f_in = U_IN / (U_IN + D_TMA / delta)
    table = " ".join(f"({t:g} {v:g})" for t, v in SHAPE)
    walls = "\n".join(f"    {w:<12}{{ type zeroGradient; }}" for w in WALLS if w != "wafer")
    (c / "0/T").write_text(
        "FoamFile { version 2.0; format ascii; class volScalarField; object T; }\n"
        "// T = TMA number density / c0, c0 the inlet value at the peak of the pulse shape.\n"
        f"// {amount:g} x the 25 ms valve amount (0.63 mg = {N_TMA_25MS:.3e} mol): {n_mol:.3e} mol = {n_mol*NA:.3e} molecules,\n"
        f"// = {n_mol*NA/(np.pi*(R_WAFER*1e-3)**2)/N0:.1f} monolayer-equivalents of the {2*R_WAFER:g} mm wafer at N0 = {N0:.2g}.\n"
        f"// c0 = {c0:.4g} /m3 (TMA partial pressure {c0*KB*T:.3g} Pa at shape = 1 against {P} Pa of Ar: the real pulse is a\n"
        "// slug of nearly pure TMA; the passive-scalar, frozen-flow treatment is an approximation stated in the deliverable).\n"
        f"// Danckwerts inlet U (c_face - c0 shape(t)) = D (c_face - c_cell)/delta: f = U/(U + D/delta) = {U_IN:.3g}/({U_IN:.3g} + {D_TMA/delta:.3g}) = {f_in:.4g}\n"
        "dimensions [0 0 0 1 0 0 0];\ninternalField   uniform 0;\nboundaryField\n{\n"
        f"    inlet\n    {{\n        type                 uniformMixed;\n        uniformValue         table ({table});\n"
        f"        uniformGradient      constant 0;\n        uniformValueFraction constant {f_in:.5g};\n    }}\n"
        f"    wafer       {{ type langmuirSticking; s {s:g}; r 0; vbar {VBAR:g}; N0 {N0:.3g}; c0 {c0:.4g}; value uniform 0; }}\n"
        + walls + "\n    outlet      { type inletOutlet; inletValue uniform 0; value uniform 0; }\n"
        "    axis        { type empty; }\n    front       { type wedge; }\n    back        { type wedge; }\n}\n")
    (c / "params.txt").write_text(f"amount {amount}\ns {s}\nn_mol {n_mol:.6e}\nc0 {c0:.6e}\nQ {Q:.6e}\nU_in {U_IN:.6e}\nD {D_TMA:.6e}\n"
                                  f"f_in {f_in:.6e}\nN0 {N0:.3e}\nvbar {VBAR}\nP {P}\nT {T}\nshape_int {SHAPE_INT:.6e}\n")
    return c

def foam(case, cmd):
    t0 = time.time()
    rc = subprocess.run([str(LAUNCHER), "-c", f"cd '{case}' && {cmd} > log.{cmd.split()[0]} 2>&1"]).returncode
    return rc, time.time() - t0

def run_steady():
    log = []
    for cmd in ("blockMesh", "checkMesh", "simpleFoam"):
        rc, dt = foam(STEADY, cmd); log.append(f"{cmd}: exit {rc}, {dt:.0f} s"); print(log[-1], flush=True)
        if rc: break
    (STEADY / "runtime.txt").write_text("\n".join(log) + "\n")

def run_pulse(c):
    last = max((q for q in STEADY.iterdir() if re.fullmatch(r"[1-9][0-9]*", q.name)), key=lambda q: int(q.name))
    if (c / "constant/polyMesh").exists(): shutil.rmtree(c / "constant/polyMesh")
    shutil.copytree(STEADY / "constant/polyMesh", c / "constant/polyMesh")
    for f in ("U", "phi"): shutil.copy(last / f, c / "0" / f)
    rc, dt = foam(c, "scalarTransportFoam")
    (c / "runtime.txt").write_text(f"scalarTransportFoam exit {rc}, wall time {dt:.0f} s, steady iteration {last.name}\n")
    print(f"{c.name}: exit {rc}, {dt/60:.1f} min", flush=True)
    return rc

if __name__ == "__main__":
    print(f"Q = {Q:.4g} m3/s, U_in = {U_IN:.3g} m/s, nu = {NU:.3g}, Re_tube = {U_IN*2*R_TUBE*1e-3/NU:.2f}, D = {D_TMA:.3g} m2/s, "
          f"residence V/Q ~ {0.0136/Q:.2f} s, shape integral {SHAPE_INT:.4f} s", flush=True)
    build_steady()
    cases = [build_pulse(n, a, s) for n, (a, s) in CASES.items()]
    print("built: steady", *[c.name for c in cases], flush=True)
    if "--run" in sys.argv:
        if not any(re.fullmatch(r"[1-9][0-9]*", q.name) for q in STEADY.iterdir()): run_steady()
        for c in cases: run_pulse(c)
