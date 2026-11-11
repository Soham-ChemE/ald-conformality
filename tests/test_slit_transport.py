"""Parallel-plate correction of the slit Knudsen diffusivity (Slit(knudsen_correction="plates")): formula values,
backward compatibility of the default, and the 1D slit model against Machball's ballistic Trench."""
import numpy as np
import pytest
from aldconf import Slit, Precursor, Process, simulate, dose_for_coverage, slit_knudsen_factor, RecessedSlit, simulate_recessed_constant

S, N0, T, P_PA, W = 5e-3, 4e18, 573.15, 0.3, 100e-9
PR = Precursor(sticking=S, site_density=N0, mass_amu=72.09)


def test_factor_is_three_eighths_log_aspect_ratio():
    assert slit_knudsen_factor(1.0, 1.0) == 1.0 and slit_knudsen_factor(2.0, 1.0) == 1.0
    assert abs(slit_knudsen_factor(100e-9, 5e-6) - 0.375 * np.log(50)) < 1e-12
    assert abs(Slit(500e-9, 250e-6, knudsen_correction="plates").knudsen_factor() - 0.375 * np.log(500)) < 1e-12
    assert Slit(500e-9, 250e-6).knudsen_factor() == 1.0
    with pytest.raises(ValueError):
        Slit(1e-7, 1e-6, knudsen_correction="log")


def test_default_reproduces_the_uncorrected_slit():
    """The WP9 reference value (figures/voxel_mc_crosscheck_runs.json, refs/1d_slit_AR50): t90 1.3813 s."""
    proc = Process(Slit(W, 50 * W), PR, T=T, p_partial=P_PA, t_dose=1.0, n_cells=200, n_substeps=200, narrowing=False)
    assert abs(dose_for_coverage(proc, 0.90, 1e-3, 1e3) / 1.3813 - 1) < 0.01


def test_plates_correction_is_the_knudsen_factor():
    a = simulate(Process(Slit(W, 50 * W, knudsen_correction="plates"), PR, T=T, p_partial=P_PA, t_dose=1.0, n_cells=100, n_substeps=100, narrowing=False))
    b = simulate(Process(Slit(W, 50 * W), PR, T=T, p_partial=P_PA, t_dose=1.0, n_cells=100, n_substeps=100, narrowing=False, knudsen_factor=0.375 * np.log(50)))
    assert np.abs(a.theta - b.theta).max() < 1e-12


def test_recess_correction_without_recesses_is_the_corrected_slit():
    geom = RecessedSlit(200e-9, 5e-6, 0.0, 0.0, 50e-9)
    res = simulate_recessed_constant(geom, PR, T, 0.3, 2.0, n_substeps=400, n_samples=100, n_slit=200, knudsen_correction="plates")
    ref = simulate(Process(Slit(200e-9, 5e-6, knudsen_correction="plates"), PR, T=T, p_partial=0.3, t_dose=2.0, n_cells=200, n_substeps=400, narrowing=False))
    assert np.abs(res.theta_slit - ref.theta).max() < 1e-9


def _t_bottom(geom, level=0.9):
    """Time at which the bottom coverage itself (not the step coverage) reaches `level`: Machball's end criterion."""
    from scipy.optimize import brentq
    def f(log_t):
        proc = Process(geom, PR, T=T, p_partial=P_PA, t_dose=float(np.exp(log_t)), n_cells=200, n_substeps=200, narrowing=False)
        return simulate(proc).theta[-1] - level
    return float(np.exp(brentq(f, np.log(1e-3), np.log(1e3), xtol=1e-4)))


@pytest.mark.parametrize("AR", [20, 50])
def test_corrected_slit_matches_machball_trench_within_5_percent(AR):
    machball = pytest.importorskip("machball")
    from machball import ALDIdeal
    from machball.ballistic import Trench
    ald = ALDIdeal(beta0=S, MM=72.09, T=T, p0=P_PA, s0=1.0 / N0)
    times, cov = ald.saturation_ballistic(Trench(AR, 8 * AR), endcov=0.90, verbose=False)
    t_mb = float(times[-1])
    t_corr = _t_bottom(Slit(W, AR * W, knudsen_correction="plates"))
    assert abs(t_corr / t_mb - 1) < 0.05, (AR, t_corr, t_mb)
    if AR == 50:  # the uncorrected model is the 24 % of WP9 (f)
        t_none = _t_bottom(Slit(W, AR * W))
        assert t_none / t_mb > 1.2
