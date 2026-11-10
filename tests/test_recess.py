"""Nested slit + lateral-recess model (aldconf.recess): limits and consistency."""
import numpy as np
from aldconf import (Slit, Precursor, Process, simulate, RecessedSlit, AreaEquivalentSlit, simulate_recessed,
                     simulate_recessed_constant, time_to_coverage)

T = 573.15
PR = Precursor(sticking=5e-3, site_density=4e18, mass_amu=72.09)


def test_no_recess_is_the_plain_slit():
    geom = RecessedSlit(200e-9, 5e-6, 0.0, 0.0, 50e-9)
    res = simulate_recessed_constant(geom, PR, T, 0.3, 2.0, n_substeps=400, n_samples=100, n_slit=200)
    ref = simulate(Process(Slit(200e-9, 5e-6), PR, T=T, p_partial=0.3, t_dose=2.0, n_cells=200, n_substeps=400, narrowing=False))
    assert np.abs(res.theta_slit - ref.theta).max() < 1e-9


def test_short_recesses_are_an_area_equivalent_slit():
    """Recess AR 1 (30 nm gap, 30 nm long, 60 nm pitch; area factor 1.5): the interior sits at the mouth density,
    so the nested model must match a slit with 1.5 x the wall sink, at an unsaturated time and at the bottom."""
    geom = RecessedSlit(200e-9, 5e-6, 30e-9, 30e-9, 60e-9)
    assert abs(geom.area_factor - 1.5) < 1e-12
    res = simulate_recessed_constant(geom, PR, T, 0.3, 0.5, n_substeps=400, n_samples=100, n_slit=200, n_recess=20)
    ref = simulate(Process(AreaEquivalentSlit(200e-9, 5e-6, geom.area_factor), PR, T=T, p_partial=0.3, t_dose=0.5,
                           n_cells=200, n_substeps=400, narrowing=False))
    assert 0.05 < ref.theta[-1] < 0.95          # compare on the front, not on a saturated profile
    assert np.abs(res.theta_slit - ref.theta).max() < 0.03
    assert abs(res.theta_slit[-1] / ref.theta[-1] - 1) < 0.03
    # and the recess bottoms follow their local wall within a few per cent
    assert np.abs(res.theta_recess[:, -1] - res.theta_slit).max() < 0.03


def test_recess_bottoms_lag_the_slit_wall_and_order_with_depth():
    geom = RecessedSlit(250e-9, 5e-6, 29e-9, 1000e-9, 51e-9)
    res = simulate_recessed_constant(geom, PR, T, 0.3, 20.0, n_substeps=800, n_samples=200, n_slit=100, n_recess=40)
    t_slit = time_to_coverage(res.t, res.hist_slit_bottom)
    t_top, t_deep = time_to_coverage(res.t, res.hist_recess_bottom_top), time_to_coverage(res.t, res.hist_recess_bottom_deep)
    assert t_top < t_slit < t_deep
    assert time_to_coverage(res.t, res.hist_recess_bottom_min) == t_deep
    # recesses add wall area: the slit bottom saturates later than in the plain slit
    plain = simulate_recessed_constant(RecessedSlit(250e-9, 5e-6, 0.0, 0.0, 51e-9), PR, T, 0.3, 20.0, n_substeps=800, n_samples=200, n_slit=100)
    assert time_to_coverage(plain.t, plain.hist_slit_bottom) < t_slit


def test_driven_zero_density_leaves_surface_bare():
    geom = RecessedSlit(250e-9, 5e-6, 29e-9, 365e-9, 51e-9)
    res = simulate_recessed(geom, PR, T, [0.0, 1.0], [0.0, 0.0], substeps_per_sample=10, n_slit=50, n_recess=10)
    assert res.theta_slit.max() == 0.0 and res.theta_recess.max() == 0.0
