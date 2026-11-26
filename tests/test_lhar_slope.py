"""aldconf.lhar: the slope-to-sticking relation of the complete lateral-channel model (table lhar_slope_table.csv)."""
import numpy as np
import pytest
from aldconf.lhar import sticking_from_slope, slope_curve, narrowing_limit, relative_sensitivity, max_fill_for, TRANSPORTS


def test_reduces_to_the_slope_method_without_narrowing_and_spot():
    s, sl = slope_curve(200, 0.0, "hydraulic", 0.0)
    assert np.all(np.abs(13.9 * sl ** 2 / s - 1) < 0.16) and np.median(np.abs(13.9 * sl ** 2 / s - 1)) < 0.05   # the hydraulic kernel is the model behind s = 13.9 slope^2


def test_round_trip_on_table_nodes():
    for tr in TRANSPORTS:
        s, sl = slope_curve(200, 0.1, tr, 11.0); k = len(s) // 2
        assert abs(sticking_from_slope(sl[k], 200, 0.1, tr)["s"] / s[k] - 1) < 1e-6


def test_slope_loses_its_grip_on_s_as_the_gap_fills():
    sens = [sticking_from_slope(slope_curve(200, f, "gas300")[1][5], 200, f, "gas300")["sensitivity"] for f in (0.0, 0.1, 0.2, 0.3)]
    assert sens[0] > 0.4 and sens[0] > sens[1] > sens[2] > sens[3] and sens[3] < 0.15
    s, sl = slope_curve(200, 0.4, "gas300"); assert sl[-1] / sl[0] < 1.1                      # 40 % filled: the slope is the same for every s
    r = sticking_from_slope(sl[-1] * 1.2, 200, 0.4, "gas300"); assert r["s"] is None and "ceiling" in r["note"]
    assert np.all(sl * 200 <= narrowing_limit(0.4) * 1.05)                                    # the model respects the narrowing limit


def test_published_slope_gives_more_than_the_slope_method():
    r = sticking_from_slope(0.019, 208, 0.19, "gas300")                                       # Yim 2020 sample 4 (nominal gap)
    assert r["slope_method_s"] == pytest.approx(13.9 * 0.019 ** 2) and r["s"] > 2 * r["slope_method_s"] and r["sensitivity"] < 0.3


def test_bad_input():
    with pytest.raises(ValueError): slope_curve(500, 0.1)
    with pytest.raises(ValueError): slope_curve(200, 0.1, "nonsense")
    with pytest.raises(ValueError): narrowing_limit(0)


def test_one_number_rule_follows_the_table():
    for tr in ("hydraulic", "gas300"):
        s0, sl0 = slope_curve(200, 0.0, tr); k = 6; base = np.log(sl0[k + 1] / sl0[k - 1]) / np.log(s0[k + 1] / s0[k - 1])
        for f in (0.1, 0.2, 0.3):
            s, sl = slope_curve(200, f, tr); ratio = np.log(sl[k + 1] / sl[k - 1]) / np.log(s[k + 1] / s[k - 1]) / base
            assert abs(ratio - relative_sensitivity(f, 200, sl0[k])) < 0.2
    f = max_fill_for(200, 0.025, 0.8); assert 0.08 < f < 0.14 and abs(relative_sensitivity(f, 200, 0.025) - 0.8) < 1e-6
