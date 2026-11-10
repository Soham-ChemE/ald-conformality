"""Supply-limit rule and the numpy-only GP surrogate read from figures/doe_gp_python.json."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from aldconf import load_surrogate, rule_pulse
from aldconf.surrogate import DEFAULT_JSON, code_factors

ROOT = Path(__file__).resolve().parents[1]
MATLAB_CENTRE_T90 = 11.62   # figures/doe_gp_report.txt, fitrgp centre point


def test_rule_reproduces_reference_wafer_pulse():
    """Reference wafer (1.78e13 holes/m2, 100 nm, AR 100) at 1 % TMA, 300 sccm, 300 Pa, 573.15 K: 14.4 s."""
    r = rule_pulse(0.01, 300, 300, 573.15, 100, 100e-9, 1.78e13)
    assert abs(r["wall_ratio"] - 55.9) < 0.1
    assert abs(r["t_pulse_s"] - 14.4) < 0.1, r


def test_rule_velocity_matches_doe_table():
    """U below the showerhead must match the U column of the DOE table within 1 %."""
    T = pd.read_csv(ROOT / "figures/doe_results.csv")
    for _, row in T.iterrows():
        U = rule_pulse(row.x, row.Q, row.p, 573.15, 100, 100e-9, 1.78e13)["U_m_s"]
        assert abs(U / row.U - 1) < 0.01


def test_json_predictor_matches_sklearn_at_training_points():
    """The numpy predictor from the JSON reproduces the scikit-learn model's t90 at the 25 runs within 0.5 %."""
    sklearn = pytest.importorskip("sklearn")  # noqa: F841
    sys.path.insert(0, str(ROOT / "analysis"))
    import doe_gp
    T = pd.read_csv(ROOT / "figures/doe_results.csv")
    X, ys, fits = doe_gp.fit_all(T)
    sur = load_surrogate(DEFAULT_JSON)
    t_sk = 10 ** doe_gp.predict(*fits["log10_t90_worst"], X)
    t_js = sur.predict_coded(X)["t90_s"]
    assert np.all(np.abs(t_js / t_sk - 1) < 0.005), np.abs(t_js / t_sk - 1).max()
    # and the JSON hyperparameters are the ones the refit produced
    spec = json.loads(Path(DEFAULT_JSON).read_text())
    ls_json = spec["responses"][0]["kernel"]["length_scale"]
    ls_fit = fits["log10_t90_worst"][1].kernel_.k1.k2.length_scale
    assert np.allclose(ls_json, ls_fit, rtol=1e-3)


def test_centre_point_matches_matlab():
    """Centre of the coded cube (x = 1.22 %, 300 sccm, 316 Pa, 17.5 mm): MATLAB fitrgp gives 11.62 s."""
    sur = load_surrogate()
    pr = sur.predict(0.005 * 6 ** 0.5, 300.0, 100 * 10 ** 0.5, 17.5)
    assert abs(pr["t90_s"] / MATLAB_CENTRE_T90 - 1) < 0.05, pr
    assert pr["t90_sd_factor"] > 1.0
    assert 0.8 < pr["utilization"] < 0.9
    assert abs(pr["edge_lag_pct"]) < 1.0


def test_coding_round_trip_and_training_points_recovered():
    """The surrogate interpolates the runs: predicted t90 within 3 % of the simulated value at every run."""
    T = pd.read_csv(ROOT / "figures/doe_results.csv")
    C = code_factors(T.x.to_numpy(), T.Q.to_numpy(), T.p.to_numpy(), T.g.to_numpy())
    assert np.allclose(C, T[["c_x", "c_Q", "c_p", "c_g"]].to_numpy(), atol=1e-9)
    sur = load_surrogate()
    pr = sur.predict(T.x.to_numpy(), T.Q.to_numpy(), T.p.to_numpy(), T.g.to_numpy())
    assert np.all(np.abs(pr["t90_s"] / T.t90_worst.to_numpy() - 1) < 0.03)
