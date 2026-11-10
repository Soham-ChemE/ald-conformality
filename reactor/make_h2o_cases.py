"""Build the H2O half-cycle cases from the sweep base pulse case (TMA settings replaced)."""
import re, shutil, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BASE = ROOT / "sweep/base/pulse"
# H2O at 573 K: vbar 821 m/s; D(H2O-N2) 2.5e-2 m2/s at 300 Pa (estimate); s ~ 1e-4 (Arts 2019 JVSTA); 1 % H2O -> c0 3.79e20
# Danckwerts f = U/(U + D/delta) = 0.05/(0.05 + 2.5e-2/2.5e-4) = 4.998e-4
CASES = {
    "showerhead_h2o_blanket": ("type langmuirSticking; s 1e-4; r 0; vbar 821; N0 4e18; c0 3.79e20; value uniform 0;", 5.0, 6.0, 0.02),
    "showerhead_h2o_patterned_AR100": ("type patternedSurface; s 1e-4; r 0; vbar 821; N0 4e18; c0 3.79e20; d 100e-9; AR 100; holeDensity 1.78e13; nCells 200; value uniform 0;", 60.0, 62.0, 0.2),
}
for name, (bc, P, end, wi) in CASES.items():
    c = ROOT / name
    if c.exists(): shutil.rmtree(c)
    c.mkdir(); 
    for sub in ("0", "constant", "system"): shutil.copytree(BASE / sub, c / sub)
    t = c / "0/T"; s = t.read_text()
    s, n1 = re.subn(r"    wafer      \{ type langmuirSticking;[^\n]*\n", "    wafer      { " + bc + " }\n", s)
    s, n2 = re.subn(r"table \(\(0 1\) \([0-9.]+ 1\) \([0-9.e-]+ 0\) \([0-9.]+ 0\)\)", f"table ((0 1) ({P} 1) ({P + 1e-4} 0) (1000 0))", s)
    s, n3 = re.subn(r"uniformValueFraction constant [0-9.eE+-]+;", "uniformValueFraction constant 4.998e-4;", s)
    assert n1 == n2 == n3 == 1, (name, n1, n2, n3)
    t.write_text(s)
    tp = c / "constant/transportProperties"
    tp.write_text(re.sub(r"DT\s+[0-9.eE+-]+;", "DT              2.5e-2;", tp.read_text()))
    cd = c / "system/controlDict"; d = cd.read_text()
    d, m1 = re.subn(r"endTime\s+[0-9.]+;", f"endTime         {end};", d)
    d, m2 = re.subn(r"writeInterval\s+[0-9.]+;", f"writeInterval   {wi};", d)
    assert m1 == m2 == 1
    cd.write_text(d)
    print(name, "built:", re.search(r"wafer\s+\{[^}]*\}", s).group(0)[:80], "| pulse", P, "| endTime", end)
