"""Loader for the Zenodo 3974438 PillarHall workbook (Al2O3, TMA/H2O, 300 C, 500 cycles, LHAR3 500 nm gap).

Sheet layout (as received): rows 0-7 are a header block; from row 8, column 0 is X (mm) along the
scan, column 2 is thickness (nm), column 3 is the reflectometry RMSE, column 5 is X re-zeroed by the
authors at the cavity entrance (mm), column 6 repeats the thickness. Sheet names encode the channel
length: L0.1mm, L0.5mm, L1mm, L5mm. Nominal gap 500 nm.
"""
from pathlib import Path
import re
import pandas as pd

GAP_NM = 500.0


def load_pillarhall(path):
    """Return one tidy DataFrame: sheet, length_mm, x_raw_mm, depth_um, thickness_nm, rmse, good."""
    xl = pd.ExcelFile(path)
    frames = []
    for sh in xl.sheet_names:
        m = re.search(r"L([\d.]+)mm", sh)
        if not m:
            continue
        df = xl.parse(sh, header=None)
        body = df.iloc[8:, [0, 2, 3, 5]].apply(pd.to_numeric, errors="coerce")
        body.columns = ["x_raw_mm", "thickness_nm", "rmse", "x_shift_mm"]
        body = body.dropna(subset=["x_raw_mm", "thickness_nm"])
        body["sheet"] = sh.strip()
        body["length_mm"] = float(m.group(1))
        body["depth_um"] = body["x_shift_mm"] * 1000.0
        # Reflectometry fit failures show up as thicknesses of thousands of nm (500 cycles of Al2O3 at
        # 300 C is about 50 nm) or as large RMSE; the record says to reject poor fits.
        body["good"] = (body["thickness_nm"] < 100.0) & (body["rmse"].fillna(0) <= 2.0)
        frames.append(body)
    out = pd.concat(frames, ignore_index=True)
    return out[["sheet", "length_mm", "x_raw_mm", "depth_um", "thickness_nm", "rmse", "good"]]


def default_path():
    return Path(__file__).resolve().parents[2] / "data" / "zenodo_3974438" / "V0001_LHAR3_Al2O3_300C_500cyc.xlsx"
