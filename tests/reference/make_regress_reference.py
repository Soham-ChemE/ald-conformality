"""Write the regression reference for reactor/showerhead_pulse_layout_regress from a local run of the case.

    PYTHONPATH=src .venv/bin/python tests/reference/make_regress_reference.py

Reads every written time with coupling/patterned_lib.load_patterned and stores the wafer patch at t = 2 s
(regress_t2.csv: face, thetaFlat, thetaBottom, value) and the first 20 faces' thetaBottom at t = 1 s, mid-pulse
(regress_t1_thetaBottom_first20.csv). tests/test_regression_reference.py compares a fresh run against these files.
"""
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "coupling"))
from patterned_lib import load_patterned  # noqa: E402

CASE = ROOT / "reactor" / "showerhead_pulse_layout_regress"
OUT = Path(__file__).resolve().parent
N_MID = 20


def main():
    t0 = time.time()
    t, f, params, pulse = load_patterned(CASE)
    print(f"parsed {len(t)} times in {time.time() - t0:.1f} s; pulse {pulse} s; params {params}")
    i2 = int(np.argmin(np.abs(t - 2.0)))
    i1 = int(np.argmin(np.abs(t - 1.0)))
    assert abs(t[i2] - 2.0) < 1e-9 and abs(t[i1] - 1.0) < 1e-9, (t[i2], t[i1])
    face = np.arange(f["thetaFlat"].shape[1])
    np.savetxt(OUT / "regress_t2.csv",
               np.column_stack([face, f["thetaFlat"][i2], f["thetaBottom"][i2], f["value"][i2]]),
               fmt=["%d", "%.8g", "%.8g", "%.8g"], delimiter=",", header="face,thetaFlat,thetaBottom,value",
               comments="")
    np.savetxt(OUT / "regress_t1_thetaBottom_first20.csv",
               np.column_stack([face[:N_MID], f["thetaBottom"][i1, :N_MID]]),
               fmt=["%d", "%.8g"], delimiter=",", header="face,thetaBottom", comments="")
    print("thetaFlat t=2: min %.6g max %.6g" % (f["thetaFlat"][i2].min(), f["thetaFlat"][i2].max()))
    print("thetaBottom t=2: min %.6g max %.6g" % (f["thetaBottom"][i2].min(), f["thetaBottom"][i2].max()))
    print("value t=2: min %.6g max %.6g" % (f["value"][i2].min(), f["value"][i2].max()))
    print("thetaBottom t=1 first 20:", np.array2string(f["thetaBottom"][i1, :N_MID], precision=5))


if __name__ == "__main__":
    main()
