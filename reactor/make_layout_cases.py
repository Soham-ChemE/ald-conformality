"""Spatially varying hole pattern across the wafer: build (and with --run, run) two chamber cases that use the
ARProfile / holeDensityProfile tables of the patternedSurface boundary condition.

  dense_centre : AR 100 everywhere; 1.78e13 holes/m2 for r < 75 mm, ramping linearly to 0 between 75 and 100 mm,
                 bare wafer (flat Langmuir surface) beyond.
  ar_ramp      : 1.78e13 holes/m2 everywhere; AR 100 at the centre ramping linearly to 60 at r = 150 mm.

Both are copies of reactor/showerhead_pulse_patterned_AR100_long (1 % TMA, 300 sccm, 300 Pa, 15 mm gap, frozen
carrier flow in 0/U and 0/phi) with the pulse shortened to 16 s and the run to 17 s; s, N0, c0, d and nCells
unchanged. Each run is one scalarTransportFoam call made through the OpenFOAM launcher:
  .venv/bin/python reactor/make_layout_cases.py          # build only
  .venv/bin/python reactor/make_layout_cases.py --run    # build and run sequentially (10-20 min each)
Post-process with coupling/postprocess_layout.py.
"""
import re, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "showerhead_pulse_patterned_AR100_long"
LAUNCHER = Path.home() / "Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam"
PULSE, END, WRITE_INTERVAL = 16.0, 17.0, 0.05
COMMON = "s 5e-3; r 0; vbar 410; N0 4e18; c0 3.79e20; d 100e-9; nCells 200;"
LAYOUTS = {
    "showerhead_pulse_layout_dense_centre":
        f"type patternedSurface; {COMMON} AR 100; holeDensity 1.78e13; "
        "holeDensityProfile table ((0 1.78e13) (0.075 1.78e13) (0.1 0) (0.15 0)); value uniform 0;",
    "showerhead_pulse_layout_ar_ramp":
        f"type patternedSurface; {COMMON} AR 100; holeDensity 1.78e13; "
        "ARProfile table ((0 100) (0.15 60)); value uniform 0;",
}


def build(name, bc):
    c = ROOT / name
    if c.exists(): shutil.rmtree(c)
    c.mkdir()
    for sub in ("0", "constant", "system"): shutil.copytree(SRC / sub, c / sub)
    t = c / "0/T"; s = t.read_text()
    s, n1 = re.subn(r"    wafer      \{ type patternedSurface;[^\n]*\n", "    wafer      { " + bc + " }\n", s)
    s, n2 = re.subn(r"table \(\(0 1\) \([0-9.]+ 1\) \([0-9.e-]+ 0\) \([0-9.]+ 0\)\)",
                    f"table ((0 1) ({PULSE:g} 1) ({PULSE + 1e-4:g} 0) (100 0))", s)
    s, n3 = re.subn(r"// 0\.1 s pulse at unit concentration, then purge", f"// {PULSE:g} s pulse at unit concentration, then purge to {END:g} s", s)
    assert n1 == n2 == n3 == 1, (name, n1, n2, n3); t.write_text(s)
    cd = c / "system/controlDict"; d = cd.read_text()
    d, m1 = re.subn(r"endTime\s+[0-9.]+;[^\n]*", f"endTime         {END:g};        // {PULSE:g} s pulse + {END - PULSE:g} s purge", d)
    d, m2 = re.subn(r"writeInterval\s+[0-9.]+;", f"writeInterval   {WRITE_INTERVAL:g};", d)
    assert m1 == m2 == 1; cd.write_text(d)
    return c


def run(c):
    # one command through the launcher; the bundled bash segfaults when a shell is spawned inside the environment
    t0 = time.time()
    rc = subprocess.run([str(LAUNCHER), "-c", f"cd '{c}' && scalarTransportFoam > log.scalarTransportFoam 2>&1"]).returncode
    dt = time.time() - t0
    (c / "runtime.txt").write_text(f"scalarTransportFoam exit {rc}, wall time {dt:.0f} s\n")
    return rc, dt


if __name__ == "__main__":
    cases = [build(n, bc) for n, bc in LAYOUTS.items()]
    print("built:", *[c.name for c in cases], flush=True)
    if "--run" in sys.argv:
        for c in cases:
            print(f"running {c.name} ...", flush=True)
            rc, dt = run(c)
            print(f"{c.name}: exit {rc}, {dt/60:.1f} min", flush=True)
