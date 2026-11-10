"""The pure functions of coupling/postprocess_layout.py (area-weighted wall-area ratio, supply rule, open area,
t90 profiles) on synthetic per-face profiles."""
import numpy as np
import pytest

import postprocess_layout as pl

D, HD, AR = 100e-9, 1.78e13, 100.0
NR = 150


def test_uniform_reference_wafer_wall_ratio_is_55_9():
    wall = pl.local_wall_ratio(np.full(NR, HD), np.full(NR, AR), D)
    assert wall.shape == (NR,)
    assert pl.area_weighted_mean(wall) == pytest.approx(55.9, abs=0.05)
    assert pl.area_weighted_mean(wall) == pytest.approx(HD * np.pi * D * AR * D, rel=1e-12)


def test_scalar_inputs_give_the_same_number():
    assert pl.local_wall_ratio(HD, AR, D) == pytest.approx(55.92, abs=0.01)


def test_area_weighting_uses_annulus_area_not_face_count():
    """Holes only inside r < 75 mm cover a quarter of the wafer area (half the faces): the average is 55.9/4."""
    hd = np.where(pl.r_mm < 75, HD, 0.0)
    wall = pl.local_wall_ratio(hd, AR, D)
    assert pl.area_weighted_mean(wall) == pytest.approx(55.92 * 0.25, rel=1e-3)
    assert pl.area_weighted_mean(wall) != pytest.approx(wall.mean(), rel=1e-2)


def test_area_weights_sum_to_the_wafer_area():
    assert np.sum(pl.area_w) == pytest.approx(np.pi * 150.0 ** 2, rel=1e-12)


def test_supply_rule_matches_the_reference_wafer_pulse():
    """1.2 (1 + 55.92) 4e18 / (0.0501 x 3.79e20): 14.4 s for the uniform AR-100 wafer at 1 % TMA."""
    t = pl.supply_rule(55.92, 4e18, 3.79e20)
    assert t == pytest.approx(14.39, abs=0.02)
    assert pl.supply_rule(0.0, 4e18, 3.79e20) == pytest.approx(14.39 / 56.92, rel=1e-2)   # blanket wafer, 0.25 s


def test_supply_rule_is_elementwise_on_a_profile():
    wall = np.array([0.0, 10.0, 55.92])
    t = pl.supply_rule(wall, 4e18, 3.79e20, U=0.0501, const=1.0)
    assert t.shape == (3,)
    assert np.allclose(t / t[0], 1 + wall)


def test_open_area_fraction_of_the_reference_pattern():
    """1.78e13 holes/m2 of 100 nm diameter: 14 % open area."""
    assert pl.open_area_fraction(np.full(NR, HD), D) == pytest.approx(0.1398, abs=5e-4)
    assert pl.open_area_fraction(np.where(pl.r_mm < 75, HD, 0.0), D) == pytest.approx(0.1398 / 4, rel=2e-3)


def test_t90_profile_first_crossing_and_mask():
    t = np.array([0.0, 1.0, 2.0, 3.0])
    arr = np.array([[0.0, 0.0, 0.0], [0.95, 0.5, 0.0], [1.0, 0.9, 0.0], [1.0, 1.0, 0.5]])
    out = pl.t90_profile(arr, t)
    assert out[0] == 1.0 and out[1] == 2.0 and np.isnan(out[2])
    masked = pl.t90_profile(arr, t, mask=np.array([True, False, True]))
    assert masked[0] == 1.0 and np.isnan(masked[1]) and np.isnan(masked[2])
