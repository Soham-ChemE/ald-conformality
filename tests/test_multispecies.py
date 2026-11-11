"""Two-species feature model (aldconf.multispecies): limits, byproduct mass balance, and the sign of re-adsorption."""
import numpy as np
import pytest
from aldconf import Hole, Precursor, Process, simulate, simulate_driven
from aldconf.multispecies import (Byproduct, simulate_two_species, simulate_two_species_driven, advance_coverages,
                                  knudsen_number, HCL_TICL4_HALFSTEP)

T = 573.15
TMA = Precursor(sticking=5e-3, site_density=4e18, mass_amu=72.09)
TICL4 = Precursor(name="TiCl4", sticking=0.10, site_density=1.24e18, mass_amu=189.68, sigma=7.04e-10)
HOLE = Hole(100e-9, 10e-6)


def _proc(pr=TMA, t_dose=5.0, n_substeps=100):
    return Process(HOLE, pr, T=T, p_partial=0.3, t_dose=t_dose, n_cells=200, n_substeps=n_substeps, narrowing=False)


def test_no_byproduct_action_is_the_single_species_solver():
    """s_B = 0 (any nu) and nu = 0 (any s_B, no mouth byproduct) reproduce aldconf.simulate to 1e-9."""
    ref = simulate(_proc()).theta
    for bp in (Byproduct("CH4", 16.04, nu=3.0, sticking=0.0), Byproduct("HCl", 36.46, nu=0.0, sticking=0.3, k_des=10.0),
               Byproduct("HCl", 36.46, nu=4.0, sticking=0.0, k_des=0.0)):
        res = simulate_two_species(_proc(), bp)
        assert np.abs(res.theta_P - ref).max() < 1e-9, bp
        assert np.all(res.theta_B == 0.0)


def test_fast_desorption_is_the_single_species_solver():
    """k_B -> infinity: blocked sites are released at once; the precursor profile is the single-species one."""
    ref = simulate(_proc()).theta
    res = simulate_two_species(_proc(), Byproduct("HCl", 36.46, nu=4.0, sticking=0.5, k_des=1e12))
    assert np.abs(res.theta_P - ref).max() < 1e-9
    assert res.theta_B.max() < 1e-9


def test_driven_matches_simulate_driven_without_byproduct_action():
    t = np.linspace(0, 3.0, 61)
    n0 = 7.9e19 * np.clip(t / 0.5, 0, 1) * (t < 2.5)
    x, th, ht, hist = simulate_driven(_proc(), t, n0, substeps_per_sample=5, record_depths=(0.0, 0.5, 1.0))
    res = simulate_two_species_driven(_proc(), Byproduct("CH4", nu=3.0), t, n0, substeps_per_sample=5,
                                      record_depths=(0.0, 0.5, 1.0))
    assert np.abs(res.theta_P - th).max() < 1e-12
    assert np.abs(res.hist_P - hist).max() < 1e-12
    assert np.allclose(res.t, ht)


def test_byproduct_mass_balance_closes():
    """produced + desorbed = escaped + re-adsorbed + etched, to 1e-6 (quasi-steady finite volumes conserve it), with
    and without etching, with and without a byproduct density at the mouth."""
    for bp in (Byproduct("HCl", 36.46, nu=4.0, sticking=0.05, k_des=10.0),
               Byproduct("HCl", 36.46, nu=4.0, sticking=0.05, k_des=1.0, k_etch=1e-3),
               Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=100.0, mouth_fraction=0.03)):
        res = simulate_two_species(_proc(TICL4, t_dose=2.0, n_substeps=200), bp)
        b = res.balance
        assert b["produced"] > 0 and b["escaped"] > 0 and b["adsorbed"] > 0
        assert abs(b["closure"]) < 1e-6, b
        if bp.k_etch == 0:
            assert b["etched"] == 0.0
        else:
            assert b["etched"] > 0
    # the surface inventory agrees with the gas-side bookkeeping up to the frozen-gas splitting error, which is first
    # order in the substep (the byproduct rates a_B ~ 50 1/s set the substep here, not the precursor's 0.3 1/s)
    errs = []
    for n_sub in (1000, 4000):
        res = simulate_two_species(_proc(TICL4, t_dose=1.0, n_substeps=n_sub), Byproduct("HCl", 36.46, nu=4.0, sticking=0.05, k_des=10.0))
        hole = res.process.geometry
        wall_per_mouth = hole.wall_area() / hole.cross_section(hole.diameter)
        inventory = TICL4.site_density * res.theta_B.mean() * wall_per_mouth
        gas_side = res.balance["adsorbed"] - res.balance["desorbed"]
        errs.append(abs(inventory - gas_side) / gas_side)
    assert errs[1] < 0.05 and errs[1] < 0.6 * errs[0], errs


def test_coverage_advance_is_exact():
    """The closed-form 2 x 2 advance equals a fine explicit integration and the single-species exponential."""
    a_P, a_B, e, kB, dt = 0.7, 2.3, 0.1, 1.5, 0.8
    thP, thB = 0.2, 0.1
    # reference by many small Euler steps
    p, b = thP, thB
    n = 200000
    h = dt / n
    for _ in range(n):
        f = 1 - p - b
        p, b = p + h * (a_P * f - e * p), b + h * (a_B * f - kB * b)
    P, B = advance_coverages(np.array([thP]), np.array([thB]), a_P, a_B, e, kB, dt)
    assert abs(P[0] - p) < 1e-5 and abs(B[0] - b) < 1e-5
    # degenerate cases: singular matrix (k_B = e = 0), no gas at all, equal eigenvalues
    P, B = advance_coverages(np.array([0.3]), np.array([0.2]), 0.0, 0.0, 0.0, 0.0, 1.0)
    assert P[0] == 0.3 and B[0] == 0.2
    P, B = advance_coverages(np.array([0.3]), np.array([0.0]), 0.5, 0.0, 0.0, 0.0, 2.0)
    assert abs(P[0] - (1 - 0.7 * np.exp(-1.0))) < 1e-14 and B[0] == 0.0
    P, B = advance_coverages(np.array([0.1]), np.array([0.1]), 0.5, 0.5, 0.0, 0.5, 1.0)   # a_P + e = a_B + k_B
    assert np.isfinite(P[0]) and np.isfinite(B[0]) and 0 < P[0] < 1 and 0 < B[0] < 1


def test_hcl_readsorption_lengthens_the_ticl4_pulse_and_etching_more():
    """TiCl4 on the AR-100 hole: HCl blocking free sites delays the bottom; etching delays it further."""
    proc = _proc(TICL4, t_dose=14.0, n_substeps=280)
    base = simulate_two_species(proc, Byproduct("HCl", 36.46, nu=0.5, sticking=0.0))
    block = simulate_two_species(proc, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=1.0))
    etch = simulate_two_species(proc, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=1.0, k_etch=1e-3))
    t_base, t_block, t_etch = (r.time_to_bottom_coverage(0.9) for r in (base, block, etch))
    assert np.isfinite(t_base) and np.isfinite(t_block)
    assert t_block > t_base * 1.01
    assert np.isnan(t_etch) or t_etch > t_block          # etching can hold the bottom below 90 % for good
    assert etch.theta_P[-1] < block.theta_P[-1] < base.theta_P[-1] + 1e-12
    assert block.hist_B.max() > 0.01
    # the byproduct diffuses faster than the precursor: its gas profile is flatter than the precursor's deficit
    assert block.n_B.max() > 0


def test_mouth_byproduct_blocks_even_without_production():
    """A byproduct density at the mouth (chamber accumulation) with nu = 0 still blocks sites."""
    proc = _proc(TICL4, t_dose=14.0, n_substeps=280)
    ref = simulate_two_species(proc, Byproduct("HCl", 36.46, nu=0.0, sticking=0.1, k_des=1.0, mouth_fraction=0.0))
    acc = simulate_two_species(proc, Byproduct("HCl", 36.46, nu=0.0, sticking=0.1, k_des=1.0, mouth_fraction=0.1))
    assert acc.hist_B.max() > 0.01 and ref.hist_B.max() == 0.0
    assert acc.time_to_bottom_coverage(0.9) > ref.time_to_bottom_coverage(0.9)


def test_knudsen_number_is_large():
    """100 nm feature, 100 to 1,000 Pa of N2 at 573 K: Kn between 100 and 2,000."""
    assert knudsen_number(1000.0, T, 100e-9) > 100
    assert knudsen_number(100.0, T, 100e-9) < 2000
    assert HCL_TICL4_HALFSTEP.nu == 0.5
