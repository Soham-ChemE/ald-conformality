"""Tests of the standalone kinetics library lib/aldkin through the ctypes binding aldconf.aldkin.

Skipped with a clear reason when the library has not been built (cmake -S lib/aldkin -B lib/aldkin/build &&
cmake --build lib/aldkin/build; ALDKIN_LIBRARY or ALDKIN_BUILD_DIR point elsewhere).

Checks: the flat Langmuir face against 1 - exp(-s vbar/4 n t/N0); the radial profile interpolation; the hole model
against aldconf.solver.simulate_driven on the long AR-100 case's centre-face drive (lib/aldkin/examples/
centre_face_drive.csv): bottom t90 within 1e-3 s and the whole coverage history within 1e-6 with the same
(exponential) update; the boundary condition's own forward-Euler scheme against the coverages the OpenFOAM run
recorded; the mixed-BC mapping; zero density and zero hole density; error handling.
"""
import numpy as np
import pytest

from aldconf import aldkin
from aldconf.transport import mean_speed, K_B, AMU

if not aldkin.is_available():
    pytest.skip("libaldkin not built: " + aldkin.BUILD_COMMAND + " (or set ALDKIN_LIBRARY / ALDKIN_BUILD_DIR)",
                allow_module_level=True)

from aldconf.aldkin import AldKin, Params, AldKinError, drive  # noqa: E402
from aldconf.solver import simulate_driven, Process  # noqa: E402
from aldconf.geometry import Hole  # noqa: E402
from aldconf.kinetics import Precursor  # noqa: E402

DRIVE = aldkin.ROOT / "lib" / "aldkin" / "examples" / "centre_face_drive.csv"

# the long case's parameters (header line of the CSV)
S, VBAR, N0, C0, D, AR, DENS, NCELLS = 5e-3, 410.0, 4e18, 3.79e20, 1e-7, 100.0, 1.78e13, 200
MASS = 72.09
T_GAS = VBAR ** 2 * np.pi * MASS * AMU / (8 * K_B)      # temperature at which mean_speed(T, MASS) = VBAR


def tma(**kw) -> Params:
    base = dict(s=S, r=0.0, vbar=VBAR, N0=N0, c0=C0, d=D, AR=AR, holeDensity=DENS, nCells=NCELLS)
    base.update(kw)
    return Params(**base)


@pytest.fixture(scope="module")
def drive_csv():
    d = np.genfromtxt(DRIVE, delimiter=",", names=True, skip_header=1)
    return d["t_s"], d["T_face"], d["thetaBottom"], d["thetaFlat"]


def test_version_and_defaults():
    assert aldkin.version() == "1.1.0"
    p = aldkin.default_params()
    assert p == Params()
    assert mean_speed(T_GAS, MASS) == pytest.approx(VBAR, rel=1e-14)


def test_flat_langmuir_analytic():
    """Exponential scheme at constant density equals the closed form to round-off (checked to 1e-9)."""
    n = 2e18
    k = AldKin(tma(flatOnly=True, c0=1.0, scheme="exponential"))
    k.set_faces(n=4)
    a = S * (VBAR / 4) * n / N0
    dt = 1e-3
    worst = 0.0
    for i in range(1, 5001):
        flux, keff = k.step(np.full(4, n), dt)
        worst = max(worst, np.abs(k.theta("flat") - (1 - np.exp(-a * i * dt))).max())
    assert worst <= 1e-9
    th = k.theta("flat")
    assert keff == pytest.approx((VBAR / 4) * S * (1 - th), rel=1e-14)
    assert flux == pytest.approx(keff * n, rel=1e-14)


def test_flat_langmuir_euler_is_first_order():
    n = 2e18
    a = S * (VBAR / 4) * n / N0
    errs = []
    for dt in (4e-3, 2e-3, 1e-3):
        k = AldKin(tma(flatOnly=True, c0=1.0, scheme="euler"))
        k.set_faces(n=1)
        for _ in range(int(round(2.0 / dt))):
            k.step([n], dt)
        errs.append(abs(k.theta("flat")[0] - (1 - np.exp(-a * 2.0))))
    assert errs[0] / errs[1] == pytest.approx(2.0, abs=0.05) and errs[1] / errs[2] == pytest.approx(2.0, abs=0.05)
    # one Euler step is the boundary condition's line
    k = AldKin(tma(flatOnly=True, kdes=0.3, c0=C0))
    k.set_faces(n=1)
    k.set_theta("flat", [0.25])
    k.step([0.4], 1e-3)
    lit = 0.25 + 1e-3 * ((VBAR / 4.0) * S * (1.0 - 0.25) * C0 * 0.4 / N0 - 0.3 * 0.25)
    assert k.theta("flat")[0] == lit


def test_profile_interpolation():
    """ARProfile table ((0 100) (0.075 100) (0.1 60) (0.15 60)): linear in between, clamped outside; independent
    density table; NULL keeps the scalar; negative values clamp to 0; density 0 makes a flat face."""
    k = AldKin(tma())
    radius = np.array([-1.0, 0.0, 0.05, 0.0875, 0.1, 0.12, 0.2])
    k.set_faces(radius)
    k.set_profile("AR", [0, 0.075, 0.1, 0.15], [100, 100, 60, 60])
    ar, dens = k.face_pattern()
    assert ar == pytest.approx([100, 100, 100, 80, 60, 60, 60], abs=1e-12)
    assert np.all(dens == DENS)
    k.set_profile("density", [0, 0.1, 0.1, 0.15], [DENS, DENS, 0.0, -5.0])
    ar, dens = k.face_pattern()
    assert ar == pytest.approx([100, 100, 100, 80, 60, 60, 60], abs=1e-12)
    assert dens.tolist() == [DENS, DENS, DENS, DENS, DENS, 0.0, 0.0]
    for _ in range(50):
        flux, keff = k.step(np.full(7, 0.05), 1e-3)
    bottom, flat = k.theta("bottom"), k.theta("flat")
    assert bottom[5] == 0.0 and bottom[6] == 0.0 and bottom[2] > 0.0
    assert keff[5] == (VBAR / 4.0) * (S * (1.0 - flat[5]) + 0.0)
    assert keff[2] > keff[5]
    assert bottom[4] > bottom[2]                                   # AR 60 fills faster than AR 100
    # the two-table convenience on one radius vector
    k2 = AldKin(tma())
    k2.set_faces(radius)
    k2.set_profiles([0, 0.075, 0.1, 0.15], AR=[100, 100, 60, 60], density=None)
    assert k2.face_pattern()[0] == pytest.approx(ar, abs=1e-12)
    k2.set_profile("AR", [], [])
    assert np.all(k2.face_pattern()[0] == AR)
    with pytest.raises(AldKinError, match="ascending"):
        k.set_profile("AR", [0.1, 0.0], [1, 2])


def test_hole_matches_simulate_driven_on_centre_face_drive(drive_csv):
    """Same discretisation and update as aldconf.solver: the whole coverage history agrees to 1e-6 and the bottom
    t90 to the sample (1e-3 s tolerance, the histories are on a 50 ms grid)."""
    t, T, _, _ = drive_csv
    proc = Process(geometry=Hole(D, AR * D), precursor=Precursor(sticking=S, site_density=N0, mass_amu=MASS),
                   T=T_GAS, n_cells=NCELLS)
    x, th_py, ht_py, hist_py = simulate_driven(proc, t, C0 * T, substeps_per_sample=5)
    k = AldKin(tma(scheme="exponential"))
    k.set_faces(n=1)
    ht, hist = drive(k, t, T, substeps_per_sample=5, record=("flat", "bottom", "hole"))
    assert np.allclose(ht, ht_py)
    mouth = hist["hole"][:, 0, 0]
    bottom = hist["bottom"][:, 0]
    assert np.abs(bottom - hist_py[:, 1]).max() <= 1e-6
    assert np.abs(mouth - hist_py[:, 0]).max() <= 1e-6
    assert np.abs(k.theta("hole")[0] - th_py).max() <= 1e-6
    t90_py = ht_py[np.argmax(hist_py[:, 1] >= 0.9)]
    t90 = ht[np.argmax(bottom >= 0.9)]
    assert abs(t90 - t90_py) <= 1e-3
    assert t90_py == pytest.approx(14.05, abs=1e-9)                # paper: 14.05 s at the centre


def test_euler_scheme_matches_the_recorded_boundary_condition(drive_csv):
    """The boundary condition's forward-Euler scheme at its 1 ms step, fed the 50 ms samples of its own face value
    (linearly interpolated), lands on the recorded hole-bottom t90 within a write interval and on the recorded
    coverages within the sampling error."""
    t, T, thetaB_rec, thetaF_rec = drive_csv
    k = AldKin(tma(scheme="euler"))
    k.set_faces([0.00067])
    dt = 1e-3
    nsteps = int(round(t[-1] / dt))
    bottom_at_writes, flat_at_writes = [0.0], [0.0]
    t90 = None
    for i in range(1, nsteps + 1):
        tnow = i * dt
        k.step([float(np.interp(tnow - dt, t, T))], dt)
        b = k.theta("bottom")[0]
        if t90 is None and b >= 0.9:
            t90 = tnow
        if abs(tnow / 0.05 - round(tnow / 0.05)) < 1e-6:
            bottom_at_writes.append(b)
            flat_at_writes.append(k.theta("flat")[0])
    t90_rec = t[np.argmax(thetaB_rec >= 0.9)]
    assert t90_rec == 14.05
    assert abs(t90 - t90_rec) <= 0.1
    bottom_at_writes = np.array(bottom_at_writes)
    assert len(bottom_at_writes) == len(t)
    assert np.abs(bottom_at_writes - thetaB_rec).max() < 0.05
    assert np.abs(np.array(flat_at_writes) - thetaF_rec).max() < 0.02


def test_mixed_bc_mapping_and_evaluate():
    k = AldKin(tma())
    k.set_faces(n=2)
    k.step([0.02, 0.05], 1e-3)
    Dd = np.array([0.0084 * 4000.0, 0.0084 * 4000.0])
    rv, rg, vf = k.mixed_bc(Dd)
    keff = k.keff()
    assert np.all(rv == 0.0) and np.all(rg == 0.0)
    assert np.all(vf == keff / (keff + Dd))
    # evaluate re-forms K_eff with new face values but leaves the coverages alone
    hole_before = k.theta("hole").copy()
    flux, keff2 = k.evaluate([0.0, 0.0])
    assert np.all(k.theta("hole") == hole_before)
    assert np.all(flux == 0.0) and np.all(keff2 > 0.0)
    # flat-only with desorption: refValue = J_des / K in T units
    kw = AldKin(tma(flatOnly=True, kdes=0.5, s=1e-4))
    kw.set_faces(n=1)
    for _ in range(100):
        kw.step([0.3], 1e-3)
    th = kw.theta("flat")[0]
    rv, rg, vf = kw.mixed_bc([33.6])
    K = (VBAR / 4) * (1e-4 * (1 - th))
    assert rv[0] == pytest.approx((0.5 * th * N0 / C0) / K, rel=1e-12)
    assert kw.desorption_flux()[0] == pytest.approx(0.5 * th * N0, rel=1e-12)


def test_zero_density_zero_holes_and_clone():
    k = AldKin(tma())
    k.set_faces(n=2)
    for _ in range(5):
        flux, keff = k.step([0.0, -1.0], 1e-3)
    assert np.all(k.theta("flat") == 0.0) and np.all(k.theta("bottom") == 0.0)
    assert np.all(flux == 0.0) and np.isfinite(keff).all() and keff[0] == keff[1] > 0
    for p in (tma(holeDensity=0.0), tma(AR=0.0)):
        kf = AldKin(p)
        kf.set_faces(n=1)
        for _ in range(20):
            _, keff = kf.step([0.2], 1e-3)
        flat = kf.theta("flat")[0]
        assert kf.theta("bottom")[0] == 0.0 and kf.theta("mid")[0] == 0.0 and flat > 0
        assert keff[0] == (VBAR / 4.0) * (S * (1.0 - flat) + 0.0)
        assert np.all(kf.hole_gas(0) == 0.0)
    # clone continues identically; restart from the written state too
    k = AldKin(tma(nCells=40))
    k.set_faces([0.0, 0.1])
    for _ in range(200):
        k.step([0.1, 0.05], 1e-3)
    c = k.clone()
    r = AldKin(tma(nCells=40))
    r.set_faces([0.0, 0.1])
    r.set_theta("hole", k.theta("hole"))
    r.set_theta("flat", k.theta("flat"))
    for _ in range(50):
        _, a = k.step([0.1, 0.05], 1e-3)
        _, b = c.step([0.1, 0.05], 1e-3)
        _, d = r.step([0.1, 0.05], 1e-3)
    assert np.all(a == b) and np.all(a == d)
    assert np.all(k.theta("hole") == c.theta("hole")) and np.all(k.theta("hole") == r.theta("hole"))


def test_errors_are_exceptions_not_crashes():
    with pytest.raises(AldKinError, match="vbar"):
        AldKin(tma(vbar=-1))
    with pytest.raises(AldKinError, match="kdes"):
        AldKin(tma(kdes=0.1))                     # desorption only in flat-only mode
    with pytest.raises(ValueError):
        Params(scheme="rk4").to_struct()
    k = AldKin(tma())
    with pytest.raises(AldKinError, match="set_faces"):
        k.step([0.1], 1e-3)
    k.set_faces(n=1)
    with pytest.raises(ValueError):
        k.step([0.1, 0.2], 1e-3)
    with pytest.raises(AldKinError):
        k.step([0.1], -1.0)
    with pytest.raises(AldKinError):
        k.step([float("nan")], 1e-3)
    with pytest.raises(AldKinError, match="out of range"):
        k.hole_gas(3)
    k.step([0.1], 1e-3)                            # still usable
    assert k.n_faces == 1 and k.n_cells == NCELLS
    assert k.params == tma()


# ---- two-species option (library 1.1) --------------------------------------------------------------------------------
TICL4_M, HCL_M = 189.68, 36.46


def ticl4(**kw) -> Params:
    """TiCl4-like hole parameters at the temperature T_GAS (vbar from the same mean_speed as the Python twin)."""
    base = dict(s=0.10, r=0.0, vbar=mean_speed(T_GAS, TICL4_M), N0=1.24e18, c0=C0, d=D, AR=AR, holeDensity=DENS, nCells=NCELLS,
                nu=0.5, vbarB=mean_speed(T_GAS, HCL_M))
    base.update(kw)
    return Params(**base)


def test_two_species_fields_round_trip_and_defaults_inactive():
    p = ticl4(sB=0.1, kB=1.0, kEtch=1e-3, cB0=1e18)
    k = AldKin(p)
    assert k.params == p
    assert aldkin.default_params().sB == 0.0 and aldkin.default_params().nu == 0.0
    # inactive fields (sB = kEtch = 0) leave the single-species path bit-identical
    a, b = AldKin(tma()), AldKin(tma(nu=4.0, vbarB=600.0, kB=5.0, cB0=1e18))
    a.set_faces(n=1); b.set_faces(n=1)
    for _ in range(100):
        fa, ka = a.step([0.1], 1e-3)
        fb, kb = b.step([0.1], 1e-3)
    assert ka[0] == kb[0] and fa[0] == fb[0]
    assert np.all(a.theta("hole") == b.theta("hole")) and np.all(b.theta("blocked_hole") == 0.0)
    assert b.byproduct_flux()[0] > 0.0 and a.byproduct_flux()[0] == 0.0
    with pytest.raises(AldKinError, match="vbarB"):
        AldKin(ticl4(sB=0.1, vbarB=0.0))


def test_two_species_hole_matches_the_python_twin():
    """Library (exponential scheme, 1 ms steps) against aldconf.multispecies.simulate_two_species on the same AR-100
    hole at a constant mouth density: precursor and blocked coverages within 1e-9 over 2 s, byproduct gas within 1e-9
    relative, and the bookkeeping closes."""
    from aldconf.multispecies import Byproduct, simulate_two_species
    p = ticl4(sB=0.1, kB=1.0, kEtch=1e-3, scheme="exponential")
    k = AldKin(p)
    k.set_faces(n=1)
    Tf, dt, nsteps = 0.1, 1e-3, 2000
    for _ in range(nsteps):
        k.step([Tf], dt)
    pr = Precursor(name="TiCl4", sticking=0.10, site_density=1.24e18, mass_amu=TICL4_M)
    bp = Byproduct("HCl", HCL_M, nu=0.5, sticking=0.1, k_des=1.0, k_etch=1e-3)
    proc = Process(geometry=Hole(D, AR * D), precursor=pr, T=T_GAS, p_partial=Tf * C0 * K_B * T_GAS, t_dose=nsteps * dt,
                   n_cells=NCELLS, n_substeps=nsteps, narrowing=False)
    res = simulate_two_species(proc, bp)
    assert np.abs(k.theta("hole")[0] - res.theta_P).max() < 1e-9
    assert np.abs(k.theta("blocked_hole")[0] - res.theta_B).max() < 1e-9
    assert res.theta_B.max() > 1e-3                      # the blocking is not trivially zero
    # the last step solved the byproduct gas with the coverages before its update, as the Python twin's n_B
    nB = k.hole_byproduct_gas(0)
    assert np.abs(nB - res.n_B).max() <= 1e-9 * res.n_B.max()
    b = k.byproduct_balance(0)
    assert abs((b["produced"] + b["desorbed"]) - (b["escaped"] + b["adsorbed"] + b["etched"])) <= 1e-9 * (b["produced"] + b["desorbed"])
    assert b["etched"] > 0.0


def test_two_species_mouth_density_blocks_the_flat_face():
    k = AldKin(ticl4(sB=0.1, kB=1.0, scheme="exponential"))
    k.set_faces(n=2)
    k.set_byproduct_mouth([0.0, 3e18])
    for _ in range(50):
        k.step([0.1, 0.1], 1e-3)
    bf = k.theta("blocked_flat")
    assert bf[0] == 0.0 and bf[1] > 0.0
    assert k.theta("blocked_bottom").shape == (2,)
    k.set_byproduct_mouth(None)
