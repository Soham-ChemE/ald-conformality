"""Tapered (frustum) holes: backward compatibility of Hole, the frustum wall area, and the area-weighted gas solve."""
import numpy as np
from aldconf import Hole, TaperedHole, Precursor, Process, simulate, knudsen_diffusivity
from aldconf.solver import _solve_gas, _area_terms

T = 573.15
PR = Precursor(sticking=5e-3, site_density=4e18, mass_amu=72.09)


def _proc(g, t_dose):
    return Process(g, PR, T=T, p_partial=0.3, t_dose=t_dose, n_cells=200, n_substeps=100, narrowing=False)


def test_straight_tapered_hole_is_the_hole():
    """bottom_ratio = 1 takes the uniform-area path and reproduces Hole exactly (no change to earlier results)."""
    a = simulate(_proc(Hole(100e-9, 10e-6), 5.0)).theta
    b = simulate(_proc(TaperedHole(100e-9, 10e-6, 1.0), 5.0)).theta
    assert np.array_equal(a, b)


def test_small_taper_continuous_with_straight():
    """The area-weighted form at bottom_ratio -> 1 agrees with the uniform form (same discretisation limit)."""
    a = simulate(_proc(Hole(100e-9, 10e-6), 5.0)).theta
    b = simulate(_proc(TaperedHole(100e-9, 10e-6, 0.999), 5.0)).theta
    assert np.abs(a - b).max() < 0.01, np.abs(a - b).max()


def test_frustum_wall_area():
    g = TaperedHole(100e-9, 10e-6, 0.6)
    x = np.linspace(0, g.length, 20001); r = g.diameter_at(x) / 2
    numeric = np.trapezoid(2 * np.pi * r * np.sqrt(1 + np.gradient(r, x) ** 2), x)
    assert abs(g.wall_area() / numeric - 1) < 1e-4
    assert abs(g.area_ratio_to_straight() - 0.8) < 1e-5      # (1 + 0.6)/2 to within the slant term (1.6e-6)
    assert TaperedHole(100e-9, 10e-6, 1.0).wall_area() == Hole(100e-9, 10e-6).wall_area()


def test_area_weighted_solve_conserves_flux():
    """Zero wall sink: the steady density in a frustum must be uniform (the conservative form is exact there)."""
    g = TaperedHole(100e-9, 10e-6, 0.5)
    x = (np.arange(200) + 0.5) * g.length / 200
    op = g.initial_opening(x); A_face, A_cell = _area_terms(g, op)
    D = knudsen_diffusivity(g.hydraulic_diameter(op), T, 72.09)
    D_face = np.empty(201); D_face[1:-1] = 2 * D[:-1] * D[1:] / (D[:-1] + D[1:]); D_face[0] = D[0]; D_face[-1] = D[-1]
    n = _solve_gas(x, D_face, np.zeros(200), 1.0, A_face, A_cell)
    assert np.abs(n - 1.0).max() < 1e-9


def test_taper_slows_saturation_at_constant_mouth_density():
    """At fixed mouth density a narrower bottom loses more conductance (D ~ d, A ~ d^2) than wall area (~d):
    the bottom saturates later. Checked at 10 s, 0.3 Pa, AR 100."""
    th_straight = simulate(_proc(Hole(100e-9, 10e-6), 10.0)).theta[-1]
    th_taper = simulate(_proc(TaperedHole(100e-9, 10e-6, 0.6), 10.0)).theta[-1]
    assert th_taper < th_straight
