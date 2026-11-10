"""The aldconf command line, run exactly as the README shows it: python -m aldconf with PYTHONPATH=src."""
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def run(*args):
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    return subprocess.run([sys.executable, "-m", "aldconf", *args], cwd=ROOT, env=env, capture_output=True,
                          text=True, timeout=300)


def number_after(text, label):
    m = re.search(re.escape(label) + r"\s*=?\s*(-?[0-9.]+(?:e[+-]?\d+)?)", text)
    assert m, (label, text)
    return float(m.group(1))


@pytest.fixture(scope="module")
def rule_out():
    r = run("rule", "--x", "0.01", "--Q", "300", "--p", "300", "--T", "573.15", "--AR", "100", "--d", "100e-9",
            "--density", "1.78e13")
    assert r.returncode == 0, r.stderr
    return r.stdout


def test_rule_prints_the_reference_wafer_pulse(rule_out):
    """1 % TMA, 300 sccm, 300 Pa, 573.15 K, AR 100, 100 nm holes at 1.78e13 /m2: t_pulse = 14.38 s."""
    t_pulse = number_after(rule_out, "t_pulse = 1.2 x t_supply =")
    assert abs(t_pulse - 14.38) < 0.05, rule_out


def test_rule_prints_wall_ratio_velocity_and_supply_time(rule_out):
    assert abs(number_after(rule_out, "A_wall/A_wafer =") - 55.92) < 0.01
    assert abs(number_after(rule_out, "U =") - 0.0501) < 0.0002
    t_supply = number_after(rule_out, "N0 / (U c0) =")
    assert abs(1.2 * t_supply - number_after(rule_out, "t_pulse = 1.2 x t_supply =")) < 2e-3


def test_rule_factor_option_scales_the_pulse():
    r = run("rule", "--x", "0.01", "--Q", "300", "--p", "300", "--AR", "100", "--d", "100e-9",
            "--density", "1.78e13", "--factor", "1.0")
    assert r.returncode == 0, r.stderr
    t_pulse = number_after(r.stdout, "t_pulse = 1 x t_supply =")
    assert abs(t_pulse - 14.38 / 1.2) < 0.05


def test_surrogate_readme_point_is_11_86_s():
    """The README example, x = 1.2 %, 300 sccm, 316 Pa, 17.5 mm (the design centre to three figures): 11.86 s."""
    r = run("surrogate", "--x", "0.012", "--Q", "300", "--p", "316", "--g", "17.5")
    assert r.returncode == 0, r.stderr
    t90 = number_after(r.stdout, "t90 worst point =")
    assert abs(t90 / 11.86 - 1) < 0.02, r.stdout


def test_surrogate_exact_centre_matches_matlab():
    """Exact centre of the coded cube, x = 0.005*sqrt(6) = 1.2247 %, p = 316.23 Pa: MATLAB fitrgp gives 11.62 s."""
    r = run("surrogate", "--x", "0.012247", "--Q", "300", "--p", "316.23", "--g", "17.5")
    assert r.returncode == 0, r.stderr
    t90 = number_after(r.stdout, "t90 worst point =")
    assert abs(t90 / 11.62 - 1) < 0.02, r.stdout
    lag = number_after(r.stdout, "edge lag =")
    assert abs(lag) < 1.0
    util = number_after(r.stdout, "utilization =")
    assert 80 < util < 90
    assert "fit LOO R2:" in r.stdout


def test_surrogate_readme_example_runs():
    r = run("surrogate", "--x", "0.012", "--Q", "300", "--p", "316", "--g", "17.5")
    assert r.returncode == 0, r.stderr
    lo_hi = re.search(r"sd factor ([0-9.]+): ([0-9.]+) to ([0-9.]+) s", r.stdout)
    assert lo_hi and float(lo_hi.group(1)) > 1.0 and float(lo_hi.group(2)) < float(lo_hi.group(3))


def test_pareto_front_is_monotone_in_both_objectives():
    r = run("pareto", "--max-lag", "5")
    assert r.returncode == 0, r.stderr
    lines = r.stdout.strip().splitlines()
    assert lines[0].startswith("Pareto front, edge lag <= 5 %, 21^4 grid,")
    n = int(re.search(r"grid, (\d+) points", lines[0]).group(1))
    rows = [list(map(float, ln.split())) for ln in lines[2:]]
    assert len(rows) == n and n >= 3
    t90 = [row[0] for row in rows]
    prec = [row[2] for row in rows]
    assert all(a < b for a, b in zip(t90, t90[1:]))          # sorted by pulse time
    assert all(a >= b for a, b in zip(prec, prec[1:]))       # slower points never use more precursor (4 printed digits)
    assert prec[0] > prec[-1]
    assert all(abs(row[1] * row[2] - 1) < 1e-3 for row in rows)   # precursor/min = 1/util


def test_pareto_grid_option_changes_the_grid():
    r = run("pareto", "--max-lag", "5", "--n", "5")
    assert r.returncode == 0, r.stderr
    assert "5^4 grid" in r.stdout.splitlines()[0]


def test_no_subcommand_is_a_usage_error_and_help_works():
    assert run().returncode == 2
    h = run("--help")
    assert h.returncode == 0 and "rule" in h.stdout and "surrogate" in h.stdout and "pareto" in h.stdout
