"""Regression test for the OpenFOAM chain: reactor/showerhead_pulse_layout_regress (steady N2 flow frozen into a
1 s TMA pulse over a patternedSurface wafer, AR 50, 2 s run) against tests/reference/regress_t2.csv, written from a
local run by tests/reference/make_regress_reference.py. Skipped when the case has not been run; the CI openfoam job
runs the steady case, copies U and phi, runs the pulse and then this test. Tolerance: ALDCONF_REGRESS_ATOL (1e-6)."""
import os
from pathlib import Path

import numpy as np
import pytest

from patterned_lib import load_patterned

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "reactor" / "showerhead_pulse_layout_regress"
REF = Path(__file__).resolve().parent / "reference"
ATOL = float(os.environ.get("ALDCONF_REGRESS_ATOL", "1e-6"))

pytestmark = pytest.mark.skipif(
    not (CASE / "2" / "T").exists(),
    reason=f"{CASE.relative_to(ROOT)}/2/T not found: run the steady case and the regression case first "
           "(reproduce.sh --regress or the CI openfoam job)")


@pytest.fixture(scope="module")
def run():
    t, f, params, pulse = load_patterned(CASE)
    return dict(t=t, f=f, params=params, pulse=pulse)


def test_case_setup_is_the_reference_setup(run):
    assert run["pulse"] == 1.0
    assert run["params"]["AR"] == 50 and run["params"]["holeDensity"] == 1 and run["params"]["nCells"] == 200
    assert run["t"][-1] == pytest.approx(2.0)


def test_wafer_fields_at_t2_match_reference(run):
    ref = np.genfromtxt(REF / "regress_t2.csv", delimiter=",", names=True)
    i = int(np.argmin(np.abs(run["t"] - 2.0)))
    assert abs(run["t"][i] - 2.0) < 1e-9
    for col in ("thetaFlat", "thetaBottom", "value"):
        got = run["f"][col][i]
        assert got.shape == ref[col].shape
        err = np.max(np.abs(got - ref[col]))
        assert err <= ATOL, f"{col}: max |diff| {err:.3e} > {ATOL:g} (face {int(np.argmax(np.abs(got - ref[col])))})"


def test_mid_pulse_bottom_coverage_matches_reference(run):
    ref = np.genfromtxt(REF / "regress_t1_thetaBottom_first20.csv", delimiter=",", names=True)
    i = int(np.argmin(np.abs(run["t"] - 1.0)))
    assert abs(run["t"][i] - 1.0) < 1e-9
    got = run["f"]["thetaBottom"][i, :len(ref)]
    err = np.max(np.abs(got - ref["thetaBottom"]))
    assert err <= ATOL, f"thetaBottom at t = 1 s: max |diff| {err:.3e} > {ATOL:g}"
    assert np.all(got < 1.0) and np.all(got > 0.99)       # mid-pulse: nearly, not fully, saturated
