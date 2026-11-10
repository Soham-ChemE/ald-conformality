"""Independent cross-check: Machball's ballistic view-factor model (Yanguas-Gil & Elam, Argonne) versus the
aldconf Knudsen-continuum solver for a cylindrical hole at constant mouth density. Same s, N0, T, p, no fitting."""
import numpy as np
import pytest
from aldconf import Hole, Precursor, Process, dose_for_coverage

machball = pytest.importorskip("machball")
from machball import ALDIdeal
from machball.ballistic import Via

S, N0, T, P_PA, D = 5e-3, 4e18, 573.15, 0.3, 100e-9


@pytest.mark.parametrize("AR,Nz,tol", [(20, 40, 0.12), (50, 100, 0.08)])
def test_saturation_time_matches_machball(AR, Nz, tol):
    pr = Precursor(sticking=S, site_density=N0, mass_amu=72.09)
    proc = Process(Hole(D, AR * D), pr, T=T, p_partial=P_PA, t_dose=1.0, n_cells=200, n_substeps=200, narrowing=False)
    t_ours = dose_for_coverage(proc, target=0.95, t_lo=1e-3, t_hi=1e3)
    ald = ALDIdeal(beta0=S, MM=72.09, T=T, p0=P_PA, s0=1.0 / N0)
    times, cov = ald.saturation_ballistic(Via(AR, Nz), endcov=0.95, verbose=False)
    t_mb = times[-1]
    assert abs(t_ours - t_mb) / t_mb < tol, (AR, t_ours, t_mb)
