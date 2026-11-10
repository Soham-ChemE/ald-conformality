"""aldconf.rule_pulse against the live page's calculator (docs/index.html, function rule()) and against the
pulse-time map, figures/pulse_time_map.csv (U = 0.0501 m/s, 300 Pa, 573.15 K, 300 sccm, 1.78e13 holes/m2, 100 nm)."""
import math
from pathlib import Path

import pandas as pd
import pytest

from aldconf import rule_pulse
from aldconf.surrogate import showerhead_velocity

ROOT = Path(__file__).resolve().parents[1]
MAP = pd.read_csv(ROOT / "figures" / "pulse_time_map.csv")
P_PA, T_K, Q_SCCM, HD, D = 300.0, 573.15, 300.0, 1.78e13, 100e-9
KB = 1.380649e-23


def docs_calculator(d, AR, den, x, Q, p, T, N0):
    """Line-for-line port of rule() in docs/index.html."""
    wall = den * math.pi * d * AR * d
    c0 = x * p / (KB * T)
    U = Q * 1.6667e-8 * (101325 / p) * (T / 273.15) / (math.pi * 0.15 * 0.15)
    ts = (1 + wall) * N0 / (U * c0)
    return dict(wall=wall, c0=c0, U=U, ts=ts, tp=1.2 * ts)


@pytest.mark.parametrize("x,Q,p,T,AR,d,den,N0", [
    (0.01, 300, 300, 573.15, 100, 100e-9, 1.78e13, 4e18),
    (0.0122, 300, 316, 573.15, 100, 100e-9, 1.78e13, 4e18),
    (0.03, 150, 1000, 473.15, 20, 50e-9, 5e13, 2e18),
    (0.005, 600, 100, 673.15, 60, 200e-9, 1e12, 4e18),
    (0.1, 300, 300, 573.15, 0, 100e-9, 0.0, 4e18),
])
def test_rule_pulse_equals_the_docs_calculator(x, Q, p, T, AR, d, den, N0):
    ours = rule_pulse(x, Q, p, T, AR, d, den, N0=N0)
    web = docs_calculator(d, AR, den, x, Q, p, T, N0)
    assert ours["wall_ratio"] == pytest.approx(web["wall"], rel=1e-12)
    assert ours["c0"] == pytest.approx(web["c0"], rel=1e-12)
    assert ours["U_m_s"] == pytest.approx(web["U"], rel=1e-12)
    assert ours["t_supply_s"] == pytest.approx(web["ts"], rel=1e-12)
    assert ours["t_pulse_s"] == pytest.approx(web["tp"], rel=1e-12)


def test_showerhead_velocity_is_the_map_velocity():
    assert showerhead_velocity(Q_SCCM, P_PA, T_K) == pytest.approx(0.0501, rel=1e-3)


def test_map_has_the_16_runs():
    assert len(MAP) == 16
    assert sorted(MAP.x.unique().tolist()) == [0.005, 0.01, 0.03, 0.1]
    assert sorted(MAP.AR.unique().tolist()) == [0, 20, 50, 100]


@pytest.mark.parametrize("row", MAP.to_dict("records"), ids=[f"{r.kind}_x{r.x}_AR{r.AR}" for r in MAP.itertuples()])
def test_rule_matches_pulse_time_map_supply_time(row):
    """1.2 x the csv's t_supply must equal rule_pulse's t_pulse within 0.5 % for the csv's x and AR."""
    den = HD if row["kind"] == "patterned" else 0.0
    r = rule_pulse(row["x"], Q_SCCM, P_PA, T_K, row["AR"], D, den)
    assert r["wall_ratio"] == pytest.approx(row["wall_ratio"], abs=1e-9)
    assert r["t_pulse_s"] == pytest.approx(1.2 * row["t_supply"], rel=5e-3)


def test_map_simulations_sit_above_the_supply_limit_within_the_stated_band():
    """The rule's claim: t90_worst / t_supply between 1.05 and 1.27 for every patterned run."""
    pat = MAP[MAP.kind == "patterned"]
    ratio = pat.t90_worst / pat.t_supply
    assert ratio.min() > 1.0 and ratio.max() < 1.3, ratio.describe()
