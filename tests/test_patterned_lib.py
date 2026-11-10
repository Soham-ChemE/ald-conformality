"""coupling/patterned_lib.load_patterned on a static fixture written in the three text forms the patternedSurface
boundary condition produces: the case author's one-line dictionary (t = 0), uniform scalars without the per-face
pattern fields (t = 0.5) and nonuniform lists with ARField / holeDensityField present (t = 1)."""
from pathlib import Path

import numpy as np
import pytest

from patterned_lib import NR, _block, _field, _scalar, load_patterned

CASE = Path(__file__).resolve().parent / "fixtures" / "patterned_case"
J = np.arange(NR)


@pytest.fixture(scope="module")
def loaded():
    return load_patterned(CASE)


def test_times_are_sorted_and_only_dirs_with_T_count(loaded):
    t, f, params, pulse = loaded
    assert t.tolist() == [0.0, 0.5, 1.0]


def test_pulse_length_read_from_the_showerhead_table(loaded):
    assert loaded[3] == 1.5


def test_every_field_has_ntime_by_150_shape(loaded):
    _, f, _, _ = loaded
    assert set(f) == {"thetaFlat", "thetaBottom", "thetaMid", "Keff", "value", "ARField", "holeDensityField"}
    for k, v in f.items():
        assert v.shape == (3, NR), k


def test_scalar_parameters_come_from_the_first_patterned_block(loaded):
    _, _, params, _ = loaded
    assert params == pytest.approx(dict(s=5e-3, vbar=410, N0=4e18, c0=3.79e20, d=100e-9, AR=100, holeDensity=1.78e13,
                                        nCells=200))


def test_initial_time_has_no_coverage_fields_so_zeros(loaded):
    _, f, _, _ = loaded
    for k in ("thetaFlat", "thetaBottom", "thetaMid", "Keff", "value"):
        assert np.all(f[k][0] == 0.0), k


def test_uniform_scalar_form_is_broadcast_to_every_face(loaded):
    _, f, _, _ = loaded
    assert np.all(f["thetaFlat"][1] == 0.25)
    assert np.all(f["thetaBottom"][1] == 0.0625)
    assert np.all(f["thetaMid"][1] == 0.125)
    assert np.all(f["Keff"][1] == 0.01)
    assert np.all(f["value"][1] == 0.5)


def test_absent_pattern_fields_fall_back_to_the_scalar_AR_and_holeDensity(loaded):
    _, f, _, _ = loaded
    assert np.all(f["ARField"][:2] == 100.0)
    assert np.all(f["holeDensityField"][:2] == 1.78e13)


def test_nonuniform_lists_are_read_face_by_face(loaded):
    _, f, _, _ = loaded
    assert np.allclose(f["thetaFlat"][2], 0.5 + 0.001 * J)
    assert np.allclose(f["thetaBottom"][2], 0.2 + 0.002 * J)
    assert np.allclose(f["thetaMid"][2], 0.3 + 0.001 * J)
    assert np.allclose(f["value"][2], 0.9 - 0.001 * J)
    assert np.all(f["Keff"][2] == 0.01)


def test_hole_density_profile_ramp_is_read_from_the_per_face_field(loaded):
    _, f, _, _ = loaded
    hd = f["holeDensityField"][2]
    assert np.all(hd[:75] == 1.78e13)
    assert np.allclose(hd[75:100], 1.78e13 * (100 - J[75:100]) / 25.0)
    assert np.all(hd[100:] == 0.0)
    assert np.all(f["ARField"][2] == 100.0)


def test_value_field_is_the_wafer_value_not_the_outlet_or_showerhead_value(loaded):
    """The other patches also carry 'value' entries (0.3 and 0.1 in the fixture); the wafer block must be isolated."""
    _, f, _, _ = loaded
    assert np.all(f["value"][1] == 0.5)
    assert f["value"][2][0] == pytest.approx(0.9)


def test_block_extracts_one_line_and_multi_line_wafer_blocks():
    one = _block((CASE / "0/T").read_text())
    multi = _block((CASE / "1/T").read_text())
    assert "patternedSurface" in one and "holeDensityProfile" in one and "zeroGradient" not in one
    assert "patternedSurface" in multi and "thetaHole" in multi and "inletOutlet" not in multi


def test_field_and_scalar_helpers_on_a_minimal_block():
    blk = "\n        s 0.005;\n        foo uniform 2.5;\n        bar nonuniform List<scalar> \n3\n(\n1\n2\n3\n)\n;\n"
    assert _scalar(blk, "s") == 0.005
    assert np.isnan(_scalar(blk, "missing"))
    assert np.all(_field(blk, "foo", n=4) == 2.5)
    assert _field(blk, "bar").tolist() == [1.0, 2.0, 3.0]
    assert np.all(_field(blk, "absent", n=7, default=-1.0) == -1.0)
