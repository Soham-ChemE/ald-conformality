"""Water half-cycle on the blanket wafer with reversible adsorption: kdes = 0 (irreversible), 0.1, 1 and 10 1/s.
Shows the gas-phase water that reappears after the pulse (the real purge problem) at the outlet, and how long the
surface takes to settle. Uses the kdes parameter of langmuirSticking (developed in a side library, since merged into aldWallKinetics)."""
import re, shutil
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BASE = ROOT / "showerhead_h2o_blanket"
for kdes in (0.0, 0.1, 1.0, 10.0):
    c = ROOT / f"showerhead_h2o_blanket_kdes{kdes:g}"
    if c.exists(): shutil.rmtree(c)
    c.mkdir()
    for sub in ("0", "constant", "system"): shutil.copytree(BASE / sub, c / sub)
    t = c / "0/T"; s = t.read_text()
    s, n = re.subn(r"type langmuirSticking; s 1e-4; r 0;", f"type langmuirSticking; s 1e-4; r 0; kdes {kdes:g};", s)
    assert n == 1; t.write_text(s)
    cd = c / "system/controlDict"; d = cd.read_text()
    d = d.replace('libs            ("libaldWallKinetics.dylib");', 'libs            ("libaldWallKinetics.dylib");')
    d, m = re.subn(r"endTime\s+[0-9.]+;", "endTime         8.0;", d)   # 5 s pulse + 3 s purge
    assert 'libaldWallKineticsDev' in d and m == 1; cd.write_text(d)
    print("built", c.name)
