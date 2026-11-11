"""Export the centre-face gas-density history of reactor/showerhead_pulse_patterned_AR100_long as a CSV.

    PYTHONPATH=src .venv/bin/python lib/aldkin/examples/export_centre_face_drive.py

Reads every written time of the long AR-100 case with coupling/patterned_lib.load_patterned and writes
lib/aldkin/examples/centre_face_drive.csv with the columns

    t_s          time (s) of the written step
    T_face       face value of the transported scalar T = c/c0 on the innermost wafer face (face 0, r = 0.67 mm)
    thetaBottom  the boundary condition's own hole-bottom coverage on that face (for comparison)
    thetaFlat    the boundary condition's flat coverage on that face

and a header line with the case's kinetic parameters. The CSV drives examples/plain_c_example.c and
tests/test_aldkin.py; the parameters in the header are the ones the library must be created with to reproduce
the case (s, r, vbar, N0, c0, d, AR, holeDensity, nCells).
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "coupling"))
from patterned_lib import load_patterned  # noqa: E402

CASE = ROOT / "reactor" / "showerhead_pulse_patterned_AR100_long"
OUT = Path(__file__).resolve().parent / "centre_face_drive.csv"
FACE = 0


def main():
    t, f, params, pulse = load_patterned(CASE)
    keys = ("s", "vbar", "N0", "c0", "d", "AR", "holeDensity", "nCells")
    header = "# " + " ".join(f"{k}={params[k]:.8g}" for k in keys) + f" r=0 pulse_s={pulse:g} face={FACE}\n"
    header += "t_s,T_face,thetaBottom,thetaFlat"
    data = np.column_stack([t, f["value"][:, FACE], f["thetaBottom"][:, FACE], f["thetaFlat"][:, FACE]])
    np.savetxt(OUT, data, fmt=["%.6g", "%.8g", "%.8g", "%.8g"], delimiter=",", header=header, comments="")
    i90 = int(np.argmax(f["thetaBottom"][:, FACE] >= 0.9))
    print(f"wrote {OUT.relative_to(ROOT)}: {len(t)} times, t = {t[0]:g}..{t[-1]:g} s, pulse {pulse:g} s")
    print(f"centre face: T_face max {f['value'][:, FACE].max():.5g}, thetaBottom reaches 0.9 at t = {t[i90]:g} s")


if __name__ == "__main__":
    main()
