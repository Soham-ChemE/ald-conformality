"""Tests of the pillar-cavity Monte Carlo lib/pillarmc/pillarmc.cpp (a single-file program, compiled here into a temporary
directory; skipped when no C++ compiler is on the path).

Checks: the gas model's tracer diffusion coefficient against first-order Chapman-Enskog theory and its thermalisation;
the ideal slit against the front position and tail that the Machball cross-check established (docs/paper/run_pillar_mc.py);
that the optional arguments in positional and key=value form give identical output and that the defaults are the cosine
law without gas; that carrier-gas collisions shorten the penetration and remove the tail monotonically with pressure and
approach the Bosanquet (hydraulic Knudsen plus gas diffusion) penetration at high pressure; that an unknown option fails.
"""
import shutil, subprocess
from pathlib import Path
import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[1] / "lib/pillarmc/pillarmc.cpp"
CXX = shutil.which("clang++") or shutil.which("g++") or shutil.which("c++")
if CXX is None:
    pytest.skip("no C++ compiler found for lib/pillarmc", allow_module_level=True)


@pytest.fixture(scope="module")
def exe(tmp_path_factory):
    out = tmp_path_factory.mktemp("pillarmc") / "pillarmc"
    subprocess.check_call([CXX, "-O2", "-std=c++17", "-pthread", "-o", str(out), str(SRC)])
    return out


def run(exe, tmp, name, *extra, s=0.02, E=3.0, N=40, M=40000, L=200):
    f = tmp / f"{name}.csv"
    subprocess.check_call([str(exe), str(s), str(E), str(N), str(M), "0.5", str(L), "0", "2", "0", "1", "2", str(f), *extra],
                          stderr=subprocess.DEVNULL)
    d = np.loadtxt(f, delimiter=",", comments="#")
    th = d[-1, 1:]; z = np.arange(len(th)) + 0.5
    i = int(np.argmax(th < 0.5)); z50 = z[i - 1] + (th[i - 1] - 0.5) / (th[i - 1] - th[i])
    return z50, float(np.interp(z50 + 60, z, th)), f.read_text()


def test_gas_diffusion_matches_chapman_enskog_and_thermalises(exe):
    out = subprocess.check_output([str(exe), "diffusion", "300", "6000"], text=True)
    ratio = float(out.split("ratio")[1].split()[0]); energy = float(out.split("c1^2)")[1].split()[0])
    assert abs(ratio - 1.0) < 0.05                            # TMA in N2: first-order Chapman-Enskog within statistics
    assert abs(energy - 1.0) < 0.03                           # the tracer keeps the Maxwell energy: the collision rule is in detailed balance


def test_ideal_slit_front_and_tail(exe, tmp_path):
    z50, tail, _ = run(exe, tmp_path, "ideal")
    assert 189 < z50 < 201                                    # 195 heights at s = 0.02, 3.0 Pa s (Machball: within 2.4 %)
    assert 0.06 < tail < 0.12                                 # the free-molecular tail 60 heights ahead of the front


def test_option_forms_agree_and_defaults_are_cosine_no_gas(exe, tmp_path):
    a = run(exe, tmp_path, "a", "2", "0.1", "0", N=6, M=8000)[2]
    b = run(exe, tmp_path, "b", "ncos=2", "fspec=0.1", N=6, M=8000)[2]
    c = run(exe, tmp_path, "c", N=6, M=8000)[2]
    d = run(exe, tmp_path, "d", "1", "0", "0", "0", "p=0", N=6, M=8000)[2]
    body = lambda t: t.split("\n", 1)[1]
    assert body(a) == body(b) and body(c) == body(d) and body(a) != body(c)
    assert subprocess.run([str(exe), "0.02", "1", "2", "100", "0.5", "50", "0", "2", "0", "1", "1", str(tmp_path / "x.csv"), "bogus=1"],
                          stderr=subprocess.DEVNULL).returncode != 0


def test_carrier_gas_cuts_tail_and_penetration(exe, tmp_path):
    r = [run(exe, tmp_path, f"p{p}", f"p={p}", M=40000 if p < 5000 else 12000, N=40 if p < 5000 else 25) for p in (0, 300, 3000, 10000)]
    z = [v[0] for v in r]; tail = [v[1] for v in r]
    assert z[0] > z[1] > z[2] > z[3] and tail[0] > 2 * tail[1] and tail[2] < 0.005
    # Bosanquet at 10000 Pa: gas diffusion (Chapman-Enskog, hard spheres) in series with the hydraulic Knudsen diffusivity
    kB, amu, T = 1.380649e-23, 1.66053906660e-27, 573.15
    vbar = np.sqrt(8 * kB * T / (np.pi * 72.09 * amu)); mu = 72.09 * 28.0134 / (72.09 + 28.0134) * amu
    Dg = 3 / (8 * 10000 / (kB * T) * (4.825e-10) ** 2) * np.sqrt(kB * T / (2 * np.pi * mu)); Dk = 2 * 0.5e-6 * vbar / 3
    D = 1 / (1 / Dg + 1 / Dk); J1 = 1 / np.sqrt(2 * np.pi * 72.09 * amu * kB * T)
    # diffusion-limited front of an irreversible Langmuir film: x50^2 = 2 D H n t / (2 N0), with n t = 4 J1 E / vbar
    x50 = np.sqrt(2 * D * 0.5e-6 * 4 * J1 * 3.0 / vbar / (2 * 4e18)) / 0.5e-6
    assert abs(z[3] / x50 - 1) < 0.15


def test_equilibrium_flux_is_uniform_with_and_without_gas(exe, tmp_path):
    # a short closed cavity with a negligible sticking probability is in equilibrium with the reservoir: every wall element
    # receives the reservoir flux, so the coverage is uniform and equal to 1 - exp(-s J E / N0), with or without carrier gas
    kB, amu, T = 1.380649e-23, 1.66053906660e-27, 573.15
    expect = 1 - np.exp(-1e-5 / np.sqrt(2 * np.pi * 72.09 * amu * kB * T) * 10.0 / 4e18)
    for name, extra in (("eq_free", ()), ("eq_gas", ("p=3000",))):
        f = tmp_path / f"{name}.csv"
        subprocess.check_call([str(exe), "1e-5", "10", "4", "60000", "0.5", "20", "0", "2", "0", "1", "2", str(f), *extra], stderr=subprocess.DEVNULL)
        th = np.loadtxt(f, delimiter=",", comments="#")[-1, 1:]
        assert abs(th.mean() / expect - 1) < 0.03 and th.std() / th.mean() < 0.05


def test_film_on_the_plates(exe, tmp_path):
    # (1) a uniform 50 nm film in a 500 nm gap is a 400 nm gap: same half-coverage depth in micrometres at the same exposure
    # (2) over a tapered film the equilibrium flux stays uniform (tilted-normal cosine law and the profile intersection are consistent)
    def x50(name, H, *extra):
        f = tmp_path / f"{name}.csv"
        subprocess.check_call([str(exe), "0.02", "3.0", "30", "40000", str(H), "200", "0", "2", "0", "1", "2", str(f), *extra], stderr=subprocess.DEVNULL)
        th = np.loadtxt(f, delimiter=",", comments="#")[-1, 1:]; i = int(np.argmax(th < 0.5)); return ((i - 0.5) + (th[i - 1] - 0.5) / (th[i - 1] - th[i])) * H
    t = tmp_path / "t05.txt"; np.savetxt(t, np.full(400, 0.05))
    assert abs(x50("film", 0.5, f"tfile={t}") / x50("bare400", 0.4) - 1) < 0.03
    tp = tmp_path / "taper.txt"; np.savetxt(tp, 0.12 * np.exp(-(np.arange(40) + 0.5) / 15)); f = tmp_path / "eq.csv"
    subprocess.check_call([str(exe), "1e-5", "10", "4", "60000", "0.5", "20", "0", "2", "0", "1", "2", str(f), f"tfile={tp}"], stderr=subprocess.DEVNULL)
    th = np.loadtxt(f, delimiter=",", comments="#")[-1, 1:]
    kB, amu, T = 1.380649e-23, 1.66053906660e-27, 573.15; expect = 1 - np.exp(-1e-5 / np.sqrt(2 * np.pi * 72.09 * amu * kB * T) * 10.0 / 4e18)
    assert abs(th.mean() / expect - 1) < 0.03 and th.std() / th.mean() < 0.05
    bad = tmp_path / "bad.txt"; np.savetxt(bad, np.zeros(7))
    assert subprocess.run([str(exe), "0.02", "1", "2", "100", "0.5", "50", "0", "2", "0", "1", "1", str(tmp_path / "x.csv"), f"tfile={bad}"], stderr=subprocess.DEVNULL).returncode != 0
