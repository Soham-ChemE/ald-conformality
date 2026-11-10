"""Full 3D showerhead chamber with a one-sided pumping port (work package 6): build the O-grid mesh, run the steady
carrier flow, build the TMA pulse case on the frozen flow and run it. Every step is skipped when its output exists, so
the script resumes after an interruption.

  .venv/bin/python reactor/chamber3d/make_and_run.py              # all steps (mesh, steady, pulse build, pulse run)
  .venv/bin/python reactor/chamber3d/make_and_run.py --timing     # mesh, steady, then a 0.2 s copy of the pulse case
                                                                  # (reactor/chamber3d/pulse_timing) to estimate the wall time
  .venv/bin/python reactor/chamber3d/make_and_run.py --no-pulse   # stop after building the pulse case

Geometry and settings (see reactor/chamber3d/steady/system/blockMeshDict): axis along x, wafer at x = 0, showerhead at
x = 15 mm, radius r = sqrt(y^2 + z^2) as in the wedge, so the patternedSurface BC profiles carry over. Patches:
showerhead (r < 150 mm, x = 15 mm), topWall, wafer (r < 150 mm, x = 0, one patch), floor, sideWall (r = 200 mm) and
outlet, a 60 degree sector of the side wall centred on +y (azimuth 0) carved with topoSet + createPatch.
Flow: simpleFoam laminar, same BCs and settings as reactor/showerhead_steady. Pulse: scalarTransportFoam with the
0/T of reactor/showerhead_pulse_patterned_AR100_long (Danckwerts uniformMixed inlet with the value fraction
recomputed for this mesh's first-cell height, patternedSurface wafer, inletOutlet port), 1 % TMA, 16 s pulse, run to
17 s, written every 0.1 s (the wedge wrote every 0.05 s; 0.1 s halves the 3D output, about 15 MB per write).
Every OpenFOAM command goes through the app launcher in one subprocess call, as reactor/make_layout_cases.py does.
"""
import re, shutil, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REACTOR = HERE.parent
STEADY = HERE / "steady"
PULSE = HERE / "pulse"
WEDGE_PULSE = REACTOR / "showerhead_pulse_patterned_AR100_long"
LAUNCHER = Path.home() / "Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam"
PULSE_S, END_S, WRITE_INTERVAL = 16.0, 17.0, 0.1
U_IN, D_TMA = 0.050, 0.0084          # showerhead speed (0/U) and TMA diffusivity (constant/transportProperties)
RUNLOG = HERE / "runtimes.txt"


def foam(case, cmd, logname=None):
    """Run one OpenFOAM command in `case` through the launcher; returns (exit code, wall seconds)."""
    logname = logname or f"log.{cmd.split()[0]}"
    t0 = time.time()
    rc = subprocess.run([str(LAUNCHER), "-c", f"cd '{case}' && {cmd} > {logname} 2>&1"]).returncode
    dt = time.time() - t0
    with RUNLOG.open("a") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {case.name:14s} {cmd:60s} exit {rc}  {dt:8.1f} s\n")
    print(f"  {case.name}: {cmd}  ->  exit {rc}, {dt:.1f} s", flush=True)
    if rc != 0:
        sys.exit(f"{cmd} failed in {case}; see {case / logname}")
    return rc, dt


def time_dirs(case):
    out = []
    for d in case.iterdir():
        if re.fullmatch(r"[0-9.]+(e-?[0-9]+)?", d.name) and d.is_dir():
            out.append((float(d.name), d))
    return sorted(out)


def mesh():
    bnd = STEADY / "constant/polyMesh/boundary"
    if bnd.exists() and "outlet" in bnd.read_text():
        print("mesh: exists (outlet patch present), skipping", flush=True); return
    print("mesh: blockMesh, topoSet, createPatch, checkMesh", flush=True)
    foam(STEADY, "blockMesh")
    foam(STEADY, "topoSet")
    foam(STEADY, "createPatch -overwrite", "log.createPatch")
    foam(STEADY, "checkMesh")
    log = (STEADY / "log.checkMesh").read_text()
    for key in ("cells:", "Mesh non-orthogonality", "Max skewness", "Max aspect ratio", "Mesh OK", "FAILED"):
        for line in log.splitlines():
            if key in line: print("   ", line.strip())


def steady():
    tds = [t for t in time_dirs(STEADY) if t[0] > 0 and (t[1] / "phi").exists()]
    if tds:
        print(f"steady: converged field at time {tds[-1][1].name} exists, skipping", flush=True); return tds[-1][1]
    print("steady: simpleFoam", flush=True)
    foam(STEADY, "simpleFoam")
    foam(STEADY, "postProcess -func 'flowRatePatch(name=outlet)' -latestTime", "log.flowRate.outlet")
    foam(STEADY, "postProcess -func 'flowRatePatch(name=showerhead)' -latestTime", "log.flowRate.showerhead")
    tds = [t for t in time_dirs(STEADY) if t[0] > 0 and (t[1] / "phi").exists()]
    return tds[-1][1]


def first_cell_half_height():
    """Showerhead face to first cell centre, from the blockMeshDict gap G and cell count nG (uniform grading)."""
    s = (STEADY / "system/blockMeshDict").read_text()
    G = float(re.search(r"^G\s+([0-9.]+);", s, re.M).group(1))
    nG = int(re.search(r"^nG\s+([0-9]+);", s, re.M).group(1))
    return 1e-3 * G / nG / 2


def write_T(case, delta):
    f = U_IN / (U_IN + D_TMA / delta)
    src = (WEDGE_PULSE / "0/T").read_text()
    wafer = re.search(r"    wafer      \{ type patternedSurface;[^\n]*\n", src).group(0)
    txt = f"""FoamFile {{ version 2.0; format ascii; class volScalarField; object T; }}
// T is the TMA mole fraction normalised to its showerhead value during the pulse (dimensionless carrier).
// Same BCs as reactor/showerhead_pulse_patterned_AR100_long on the full 3D chamber with a one-sided port.
dimensions [0 0 0 1 0 0 0];
internalField   uniform 0;
boundaryField
{{
    // {PULSE_S:g} s pulse at unit concentration, then purge to {END_S:g} s
    // Danckwerts total-flux inlet: U (c_face - c0(t)) = D (c_face - c_cell)/delta, so the TMA entering the chamber is
    // exactly U*c0*A during the pulse. As a mixed BC: value fraction f = U/(U + D/delta) with delta = {delta*1e3:.4g} mm (first cell
    // centre of this mesh, 12 cells over the 15 mm gap): f = {U_IN:g}/({U_IN:g} + {D_TMA/delta:.4g}) = {f:.6g}.
    showerhead
    {{
        type                 uniformMixed;
        uniformValue         table ((0 1) ({PULSE_S:g} 1) ({PULSE_S + 1e-4:g} 0) (100 0));
        uniformGradient      constant 0;
        uniformValueFraction constant {f:.6g};
    }}
    topWall    {{ type zeroGradient; }}
    // self-limiting ALD surface: s on open sites, coverage theta integrated per face (reactor/lib/aldWallKinetics)
    // TMA: vbar 410 m/s at 573 K; N0 4e18 sites/m2; c0 = TMA number density in the feed at T = 1 (3 Pa at 573 K -> 3.79e20 1/m3).
{wafer}    floor      {{ type zeroGradient; }}
    sideWall   {{ type zeroGradient; }}
    outlet     {{ type inletOutlet; inletValue uniform 0; value uniform 0; }}
}}
"""
    (case / "0/T").write_text(txt)
    return f


def build_pulse(case, steady_time, end=END_S):
    if (case / "0/T").exists() and (case / "0/phi").exists():
        print(f"{case.name}: exists, skipping build", flush=True); return
    print(f"{case.name}: building from {STEADY.name}/{steady_time.name} and {WEDGE_PULSE.name}", flush=True)
    if case.exists(): shutil.rmtree(case)
    (case / "0").mkdir(parents=True); (case / "constant").mkdir(); (case / "system").mkdir()
    shutil.copytree(STEADY / "constant/polyMesh", case / "constant/polyMesh")
    shutil.copy(WEDGE_PULSE / "constant/transportProperties", case / "constant/transportProperties")
    for f in ("fvSchemes", "fvSolution"): shutil.copy(WEDGE_PULSE / "system" / f, case / "system" / f)
    cd = (WEDGE_PULSE / "system/controlDict").read_text()
    cd, m1 = re.subn(r"endTime\s+[0-9.]+;[^\n]*", f"endTime         {end:g};        // {PULSE_S:g} s pulse + purge", cd)
    cd, m2 = re.subn(r"writeInterval\s+[0-9.]+;", f"writeInterval   {WRITE_INTERVAL:g};", cd)
    assert m1 == m2 == 1
    (case / "system/controlDict").write_text(cd)
    shutil.copy(steady_time / "U", case / "0/U")
    shutil.copy(steady_time / "phi", case / "0/phi")
    delta = first_cell_half_height()
    f = write_T(case, delta)
    print(f"  Danckwerts value fraction {f:.6g} (delta {delta*1e3:.4g} mm)", flush=True)
    # face centres of every patch (boundaryField of C) for the post-processing
    foam(case, "postProcess -func writeCellCentres -time 0", "log.writeCellCentres")


def run_pulse(case, end):
    tds = time_dirs(case)
    if tds and abs(tds[-1][0] - end) < 1e-6 and (tds[-1][1] / "T").exists():
        print(f"{case.name}: already ran to {end:g} s, skipping", flush=True); return
    print(f"{case.name}: scalarTransportFoam to {end:g} s", flush=True)
    rc, dt = foam(case, "scalarTransportFoam")
    (case / "runtime.txt").write_text(f"scalarTransportFoam to {end:g} s: exit {rc}, wall time {dt:.0f} s\n")
    return dt


if __name__ == "__main__":
    mesh()
    st = steady()
    if "--timing" in sys.argv:
        tc = HERE / "pulse_timing"
        build_pulse(tc, st, end=0.2)
        dt = run_pulse(tc, 0.2)
        if dt:
            print(f"timing: 0.2 s of simulation took {dt:.0f} s wall; the {END_S:g} s run would take about "
                  f"{dt * END_S / 0.2 / 3600:.2f} h serial", flush=True)
        sys.exit(0)
    build_pulse(PULSE, st)
    if "--no-pulse" in sys.argv: sys.exit(0)
    run_pulse(PULSE, END_S)
