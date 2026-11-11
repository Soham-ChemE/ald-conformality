"""Tests of the 3D voxel Monte Carlo feature solver lib/aldmc through the ctypes binding aldconf.voxel.

Skipped with a clear reason when the library has not been built (cmake -S lib/aldmc -B lib/aldmc/build &&
cmake --build lib/aldmc/build; ALDMC_LIBRARY or ALDMC_BUILD_DIR point elsewhere).

Checks: a straight AR-20 hole at 0.3 Pa against the 1D Knudsen solver's t90 (within 15 %: a ballistic solver sits
about 10 % above the continuum one at AR 20, Machball's own t90 is 8.5 % above it, and the small particle budget of
the test scatters by 1-2 % between seeds) and, when Machball is installed, against Machball's t95 (within 5 %, the
tight check between the two ballistic codes); a flat surface against the Langmuir law; an open box with no walls loses every particle; a periodic slit
conserves the particle count; the same seed gives bitwise identical results single-threaded and with four threads;
the geometry builders and error handling.
"""
import numpy as np
import pytest

from aldconf import voxel

if not voxel.is_available():
    pytest.skip("libaldmc not built: " + voxel.BUILD_COMMAND + " (or set ALDMC_LIBRARY / ALDMC_BUILD_DIR)",
                allow_module_level=True)

from aldconf import Hole, Precursor, Process, dose_for_coverage, mean_speed  # noqa: E402
from aldconf.kinetics import number_density  # noqa: E402

S, N0, T, P_PA, D = 5e-3, 4e18, 573.15, 0.3, 100e-9
MASS = 72.09
VBAR = float(mean_speed(T, MASS))
N_MOUTH = float(number_density(P_PA, T))


def params(**kw):
    base = dict(s=S, N0=N0, vbar=VBAR, seed=3, threads=0)      # 0: min(8, cpu count), or ALDMC_THREADS
    base.update(kw)
    return voxel.Params(**base)


def test_version_and_library():
    assert voxel.version() == "1.0.0"
    assert voxel.library_path() is not None


@pytest.fixture(scope="module")
def ar20_voxel():
    """AR-20 hole, 16 voxels across, 5000 particles per 5 ms step (about 10 s on four threads)."""
    g = voxel.hole(D, 20 * D, D / 16)
    res = voxel.simulate_constant_voxel(g, params(), N_MOUTH, 2.0, 0.005, M=5000, stop_at_bottom=0.97)
    return res


def test_ar20_hole_t90_against_1d(ar20_voxel):
    pr = Precursor(sticking=S, site_density=N0, mass_amu=MASS)
    proc = Process(Hole(D, 20 * D), pr, T=T, p_partial=P_PA, t_dose=1.0, n_cells=200, n_substeps=200, narrowing=False)
    t90_1d = dose_for_coverage(proc, target=0.90, t_lo=1e-3, t_hi=1e3)
    t90_mc = voxel.time_to_level(ar20_voxel.t, ar20_voxel.bottom, 0.9)
    assert abs(t90_mc - t90_1d) / t90_1d < 0.15, (t90_mc, t90_1d)
    assert t90_mc > t90_1d                                  # ballistic transport is the slower one at low AR
    # a Monte Carlo run is bookkept: every particle escaped or stuck (no recombination, closed geometry)
    c = ar20_voxel.counts
    assert c["lost"] == 0 and c["left"] == 0 and c["escaped"] + c["stuck"] == ar20_voxel.rate["particles"]
    # the coverage falls with depth at the end of the run and the floor is the last to saturate
    prof = ar20_voxel.theta[ar20_voxel.slice_weight > 0]
    assert prof[0] > prof[len(prof) // 2] > 0.9


def test_ar20_hole_t95_against_machball(ar20_voxel):
    machball = pytest.importorskip("machball")
    from machball import ALDIdeal
    from machball.ballistic import Via
    ald = ALDIdeal(beta0=S, MM=MASS, T=T, p0=P_PA, s0=1.0 / N0)
    times, _ = ald.saturation_ballistic(Via(20, 160), endcov=0.95, verbose=False)   # Nz 160: Machball converged to 0.5 %
    t95_mc = voxel.time_to_level(ar20_voxel.t, ar20_voxel.bottom, 0.95)
    assert abs(t95_mc - times[-1]) / times[-1] < 0.05, (t95_mc, times[-1])


def test_flat_surface_follows_langmuir():
    """theta(t) = 1 - exp(-s n vbar/4 t / N0): every particle hits the floor exactly once and leaves."""
    g = voxel.flat_surface(8, 5e-9)
    p = params(s=0.5, seed=11, threads=2)
    rate = p.s * N_MOUTH * p.vbar / 4.0 / p.N0
    mc = voxel.VoxelMC(g, p)
    assert mc.n_faces == 64 and mc.mouth_area == pytest.approx(64 * 25e-18)
    assert mc.site_area == pytest.approx(64 * 25e-18)
    dt = 0.1 / rate
    for _ in range(10):
        mc.step(N_MOUTH, dt, 64000)
    theta = mc.face_coverage().mean()
    assert theta == pytest.approx(1.0 - np.exp(-1.0), abs=5e-3)
    c = mc.counts()
    assert c["left"] == 0 and c["escaped"] + c["stuck"] == 640000
    assert mc.bottom_coverage() == pytest.approx(theta)


def test_open_box_loses_every_particle():
    g = voxel.open_box(8, 8, 8, 5e-9)
    mc = voxel.VoxelMC(g, params())
    assert mc.n_faces == 0
    mc.step(1e20, 1e-3, 5000)
    c = mc.counts()
    assert c["stuck"] == 0 and c["lost"] == 0 and c["escaped"] + c["left"] == 5000
    assert mc.bottom_coverage() == 0.0


def test_periodic_slit_conserves_particles():
    g = voxel.slit(50e-9, 1e-6, 5e-9, nx=4)
    assert g.periodic == (True, False) and g.shape == (4, 12, 201)
    mc = voxel.VoxelMC(g, params())
    for _ in range(5):
        mc.step(N_MOUTH, 0.05, 20000)
    c = mc.counts()
    assert c["escaped"] + c["stuck"] + c["lost"] + c["left"] == 100000
    assert c["left"] == 0 and c["lost"] == 0 and c["stuck"] > 0 and c["escaped"] > 0
    z, prof, wt = mc.slice_coverage()
    assert len(z) == 201 and wt[200] == 0.0 and prof[0] > prof[150]


def test_recombination_only_loses_everything():
    g = voxel.slit(50e-9, 1e-6, 5e-9)
    mc = voxel.VoxelMC(g, params(s=0.0, r=1.0, area_correction=False))
    mc.step(N_MOUTH, 0.05, 5000)
    assert mc.counts() == dict(escaped=0, stuck=0, lost=5000, left=0)


def test_same_seed_is_deterministic_and_thread_invariant():
    g = voxel.hole(D, 4 * D, D / 10)

    def run(threads):
        mc = voxel.VoxelMC(g, params(seed=12345, threads=threads))
        for _ in range(6):
            mc.step(N_MOUTH, 0.05, 30000)
        return mc.face_coverage(), mc.counts()

    a, ca = run(1)
    b, cb = run(1)
    c, cc = run(4)
    assert np.array_equal(a, b) and ca == cb
    assert np.array_equal(a, c) and ca == cc
    mc = voxel.VoxelMC(g, params(seed=999, threads=1))
    for _ in range(6):
        mc.step(N_MOUTH, 0.05, 30000)
    assert not np.array_equal(a, mc.face_coverage())


def test_area_correction_recovers_the_smooth_cylinder():
    g = voxel.hole(D, 20 * D, D / 24)
    smooth = np.pi * D * 20 * D + np.pi * D ** 2 / 4
    corrected = voxel.VoxelMC(g, params()).site_area
    plain = voxel.VoxelMC(g, params(area_correction=False)).site_area
    assert abs(corrected / smooth - 1.0) < 0.03
    assert 1.15 < plain / smooth < 1.30            # staircase: 4/pi on the wall


def test_geometry_builders():
    g = voxel.hole(D, 10 * D, D / 16)
    assert g.shape == (18, 18, 161) and g.meta["voxels_across"] == 16
    assert g.occupancy[:, :, 160].all() and g.occupancy[9, 9, 0] == 0 and g.occupancy[0, 0, 0] == 1
    tp = voxel.tapered_hole(D, 10 * D, 0.5, D / 16)
    assert tp.occupancy[:, :, 0].sum() < tp.occupancy[:, :, 159].sum()      # narrower at the bottom
    with pytest.raises(ValueError):
        voxel.tapered_hole(D, 10 * D, 0.0, D / 16)
    r = voxel.recessed_slit(100e-9, 1e-6, 29e-9, 365e-9, 51e-9, 5e-9)
    assert r.periodic == (True, False) and r.meta["n_recesses"] == 20 and r.meta["g_r_actual"] == pytest.approx(30e-9)
    k0, k1 = r.meta["recess_k"][0]
    y0, y1 = r.meta["y_slit"]
    assert r.occupancy[0, 1, k0] == 0 and r.occupancy[0, 1, k0 - 1] == 1 and r.occupancy[0, y0, 0] == 0
    fa = voxel.from_array(np.zeros((2, 3, 4)), 1e-9, periodic=(True, True))
    assert fa.occupancy.dtype == np.uint8 and fa.shape == (2, 3, 4)
    with pytest.raises(ValueError):
        voxel.from_array(np.zeros((2, 3)), 1e-9)


def test_error_handling():
    g = voxel.hole(D, 4 * D, D / 8)
    with pytest.raises(voxel.AldMCError):
        voxel.VoxelMC(g, params(s=2.0))
    with pytest.raises(ValueError):
        voxel.VoxelMC(g, params(update="nonsense"))
    mc = voxel.VoxelMC(g, params())
    with pytest.raises(voxel.AldMCError):
        mc.step(-1.0, 0.1, 10)
    with pytest.raises(voxel.AldMCError):
        mc.set_face_coverage(np.full(mc.n_faces, 2.0))
    mc.set_face_coverage(np.full(mc.n_faces, 0.25))
    assert mc.bottom_coverage() == pytest.approx(0.25)
    mc.reset()
    assert mc.bottom_coverage() == 0.0 and mc.time == 0.0
    mc.step(0.0, 0.5, 100)                      # zero density: time advances, nothing launched
    assert mc.time == pytest.approx(0.5) and mc.counts()["escaped"] == 0


def test_driven_result_unpacks_like_simulate_driven():
    g = voxel.hole(D, 4 * D, D / 8)
    t = np.linspace(0, 0.5, 6)
    res = voxel.simulate_driven_voxel(g, params(), t, np.full(6, N_MOUTH), M=5000, record_depths=(0.0, 0.5, 1.0),
                                      snapshot_times=(0.3,))
    x, theta, ht, hist = res
    assert len(x) == g.shape[2] and theta.shape == x.shape and hist.shape == (5, 3) and len(ht) == 5
    assert 0.3 in res.snapshots and res.bottom.shape == (5,)
    assert res.rate["particles"] == 25000
