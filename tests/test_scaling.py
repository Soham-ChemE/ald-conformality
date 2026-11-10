"""Scaling laws the model must reproduce before it is trusted on data."""
import numpy as np
import pytest
from aldconf import Slit, Precursor, Process, simulate, penetration_depth, step_coverage, dose_for_coverage, fit_sticking


def make(gap=500e-9, L=200e-6, s=0.05, t_dose=0.05, p=1.0, cycles=1, n_cells=300):
    pr = Precursor(sticking=s)
    return Process(Slit(gap, L), pr, T=573.15, p_partial=p, t_dose=t_dose, cycles=cycles,
                   n_cells=n_cells, n_substeps=150, narrowing=False)


def test_profile_is_monotone_and_saturates_at_mouth():
    res = simulate(make())
    assert res.theta[0] > 0.99
    assert np.all(np.diff(res.thickness) <= 1e-15)


def test_penetration_scales_as_sqrt_dose():
    """Diffusion-limited front: PD50 ~ sqrt(dose) (Gordon et al. 2003)."""
    pd1 = penetration_depth(simulate(make(t_dose=0.02)))
    pd4 = penetration_depth(simulate(make(t_dose=0.08)))
    ratio = pd4 / pd1
    assert 1.8 < ratio < 2.2, ratio


def test_dose_scales_as_aspect_ratio_squared():
    """Dose for 95% step coverage ~ AR^2 in the diffusion-limited regime (Gordon et al. 2003).

    At small s the front is partly reaction-limited and the exponent drops below 2 (measured 3.2x
    for s = 0.05 at AR 20 -> 40), so the law is tested where it should hold, s = 0.5."""
    t20 = dose_for_coverage(make(L=20 * 500e-9, s=0.5), target=0.95)
    t40 = dose_for_coverage(make(L=40 * 500e-9, s=0.5), target=0.95)
    ratio = t40 / t20
    assert 3.5 < ratio < 4.6, ratio


def test_sticking_is_recovered_from_synthetic_profile():
    truth = 2e-3
    proc = make(s=truth, t_dose=0.5)
    res = simulate(proc)
    s_fit, rms = fit_sticking(res.x, res.thickness, proc, s_bounds=(1e-5, 1e-1))
    assert abs(np.log10(s_fit / truth)) < 0.05, (s_fit, rms)


def test_recombination_caps_coverage():
    """With radical loss r, more dose stops helping: coverage at the bottom plateaus below 1."""
    base = make(L=100 * 500e-9, s=0.05)
    pr = Precursor(sticking=0.05, recombination=5e-3)
    p1 = Process(**{**base.__dict__, "precursor": pr, "t_dose": 0.5})
    p2 = Process(**{**base.__dict__, "precursor": pr, "t_dose": 5.0})
    sc1, sc2 = step_coverage(simulate(p1)), step_coverage(simulate(p2))
    assert sc2 < 0.9
    assert sc2 - sc1 < 0.1
