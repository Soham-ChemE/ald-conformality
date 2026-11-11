"""Cross-checks of the 3D voxel Monte Carlo feature solver (lib/aldmc, aldconf.voxel) against this repository's 1D
Knudsen/Langmuir solver and Argonne's Machball (ballistic view factors).

Parts (each well under 20 minutes on 8 threads; results are cached in figures/voxel_mc_crosscheck_runs.json so a part
can be rerun alone and the CSV and figure are assembled from the cache):
  holes        straight AR 20 / 50 / 100 holes at the Machball-check conditions (s 5e-3, N0 4e18, 573 K, 0.3 Pa TMA,
               100 nm): time to 90 % and 95 % floor coverage, voxel vs 1D vs Machball
  convergence  AR 100: particles per step M, voxel resolution (16 / 24 / 32 across), time step, update scheme,
               seed, and the uncorrected staircase (area_correction off)
  tapered      frustum, bottom 60 % of the mouth, constant 0.3 Pa, vs the 1D frustum model (figures/tapered_hole.csv)
  recess       word-line slit w100_Lr365 at constant 0.3 Pa vs the nested recess model (figures/wordline_recess.csv)
  driven       AR-100 hole driven by the chamber's centre-face density history (the long case, 50 ms samples
               subsampled to 250 ms), vs aldconf.simulate_driven's 14.05 s and its depth profiles at 5, 10 and 13 s
  slit         AR-50 slit (no staircase) vs Machball's Trench and the 1D slit model
  machball     Machball's own wall discretisation (Nz) and the converged reference values

    PYTHONPATH=src .venv/bin/python coupling/voxel_mc_crosscheck.py --part holes,convergence,tapered,recess,driven
    PYTHONPATH=src .venv/bin/python coupling/voxel_mc_crosscheck.py --figure      # CSV + PNG from the cache

Output: figures/voxel_mc_crosscheck.csv (one row per voxel run with its references), figures/voxel_mc_crosscheck.png,
figures/voxel_mc_crosscheck_runs.json (cache, with the histories), figures/voxel_mc_profiles.csv (driven profiles).
"""
import sys, json, time, argparse, csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "coupling"))
from aldconf import (Hole, TaperedHole, Slit, Precursor, Process, dose_for_coverage, simulate_driven, mean_speed,  # noqa: E402
                     RecessedSlit, simulate_recessed_constant, time_to_coverage)
from aldconf.kinetics import number_density  # noqa: E402
from aldconf import voxel  # noqa: E402

FIG = ROOT / "figures"
CACHE = FIG / "voxel_mc_crosscheck_runs.json"
DRIVE = ROOT / "lib/aldkin/examples/centre_face_drive.csv"

S, N0, T, P_PA, D = 5e-3, 4e18, 573.15, 0.3, 100e-9
MASS = 72.09
VBAR = float(mean_speed(T, MASS))
N0_GAS = float(number_density(P_PA, T))
C0 = 3.79e20                                   # 1 % TMA at 300 Pa, 573 K (the chamber runs' c0)
PR = Precursor(sticking=S, site_density=N0, mass_amu=MASS)
THREADS = voxel.default_threads()


def load_cache():
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def save_cache(c):
    CACHE.write_text(json.dumps(c, indent=1))


def params(seed=1, **kw):
    base = dict(s=S, N0=N0, vbar=VBAR, seed=seed, threads=THREADS)
    base.update(kw)
    return voxel.Params(**base)


def smooth_area(kind, **g):
    if kind == "hole":
        return np.pi * g["d"] * g["L"] + np.pi * g["d"] ** 2 / 4
    if kind == "tapered":
        return TaperedHole(g["d"], g["L"], g["br"]).wall_area() + np.pi * (g["br"] * g["d"]) ** 2 / 4
    return float("nan")


def record(res, mc, extra):
    r = dict(extra)
    r.update(t90=voxel.time_to_level(res.t, res.bottom, 0.9), t95=voxel.time_to_level(res.t, res.bottom, 0.95),
             wall_s=res.rate["wall_s"], particles=res.rate["particles"], threads=res.rate["threads"],
             pps=res.rate["particles_per_s"], pps_thread=res.rate["particles_per_s_per_thread"],
             seg_per_particle=res.rate["segments_per_particle"], n_faces=mc.n_faces, mouth_area=mc.mouth_area,
             site_area=mc.site_area, counts=res.counts, t_hist=[float(x) for x in res.t], bottom_hist=[float(x) for x in res.bottom])
    return r


def run_hole_constant(AR, across, dt, M, seed=1, br=1.0, t_end=None, **pkw):
    g = voxel.tapered_hole(D, AR * D, br, D / across) if br < 1 else voxel.hole(D, AR * D, D / across)
    p = params(seed=seed, **pkw)
    t_end = t_end or 4.0 * (AR / 20) ** 2 * 1.2
    print(f"  voxel hole AR {AR} br {br} {across} across ({g.shape}), dt {dt}, M {M}, seed {seed}, {pkw or ''}", flush=True)
    res = voxel.simulate_constant_voxel(g, p, N0_GAS, t_end, dt, M=M, stop_at_bottom=0.96)
    mc = res.solver
    kind = "tapered" if br < 1 else "hole"
    rec = record(res, mc, dict(case=f"AR{AR}" + (f"_br{br}" if br < 1 else ""), kind=kind, AR=AR, bottom_ratio=br, across=across,
                              voxel_nm=D / across * 1e9, dt=dt, M=M, seed=seed, mouth_density=N0_GAS,
                              inventory_ratio=mc.site_area / smooth_area(kind, d=D, L=AR * D, br=br),
                              area_correction=p.area_correction, update=p.update))
    print(f"    t90 {rec['t90']:.3f} s, t95 {rec['t95']:.3f} s; {rec['wall_s']:.0f} s wall, {rec['pps_thread']:.0f} p/s/thread, "
          f"{rec['seg_per_particle']:.0f} segments/particle, inventory {rec['inventory_ratio']:.4f}", flush=True)
    return rec


def one_d_hole(AR, br=1.0):
    g = TaperedHole(D, AR * D, br) if br < 1 else Hole(D, AR * D)
    proc = Process(g, PR, T=T, p_partial=P_PA, t_dose=1.0, n_cells=200, n_substeps=200, narrowing=False)
    return dict(t90=dose_for_coverage(proc, 0.90, 1e-3, 1e3), t95=dose_for_coverage(proc, 0.95, 1e-3, 1e3))


def machball_hole(AR, Nz=None, structure="via"):
    """Machball's t90 / t95 for a via (or trench) of aspect ratio AR with Nz wall sections. Nz = 2 AR is what
    tests/test_machball_crosscheck.py uses; the times still move by 3-6 % between Nz = 2 AR and 8 AR (see part
    `machball`), so the reference here is Nz = 8 AR (16 AR for the AR-100 via)."""
    if Nz is None:
        Nz = 16 * AR if AR >= 100 else 8 * AR
    try:
        from machball import ALDIdeal
        from machball.ballistic import Via, Trench
    except ImportError:
        return dict(t90=float("nan"), t95=float("nan"), Nz=Nz)
    ald = ALDIdeal(beta0=S, MM=MASS, T=T, p0=P_PA, s0=1.0 / N0)
    out = dict(Nz=Nz)
    t0 = time.time()
    for level in (0.9, 0.95):
        st = Via(AR, Nz) if structure == "via" else Trench(AR, Nz)
        times, cov = ald.saturation_ballistic(st, endcov=level, verbose=False)
        out[f"t{int(level * 100)}"] = float(times[-1])
    out["wall_s"] = time.time() - t0
    return out


def part_machball(cache):
    """Machball's own discretisation: t90 / t95 of the AR 20 / 50 / 100 vias and the AR-50 trench against Nz, and the
    converged values stored as the references mb_AR*, mb_slit_AR50 (replacing the Nz = 2 AR values)."""
    refs = cache.setdefault("refs", {})
    table = refs.setdefault("mb_nz", {})
    for AR, factors in ((20, (2, 8, 32)), (50, (2, 8, 16)), (100, (2, 4, 8, 16))):
        for f in factors:
            k = f"via_AR{AR}_Nz{f * AR}"
            if k not in table:
                table[k] = machball_hole(AR, f * AR); save_cache(cache)
            print(f"Machball via AR {AR} Nz {f * AR}: t90 {table[k]['t90']:.3f} t95 {table[k]['t95']:.3f} [{table[k].get('wall_s', 0):.0f} s]", flush=True)
        refs[f"mb_AR{AR}"] = dict(table[f"via_AR{AR}_Nz{factors[-1] * AR}"])
    for f in (2, 8):
        k = f"trench_AR50_Nz{f * 50}"
        if k not in table:
            table[k] = machball_hole(50, f * 50, "trench"); save_cache(cache)
        print(f"Machball trench AR 50 Nz {f * 50}: t90 {table[k]['t90']:.3f} t95 {table[k]['t95']:.3f}", flush=True)
    refs["mb_slit_AR50"] = dict(table["trench_AR50_Nz400"])
    save_cache(cache)


# ------------------------------------------------------------------------------------------------------------------
def part_holes(cache):
    refs = cache.setdefault("refs", {})
    # two time steps per hole at the same total particle budget: dt = t90/50 and t90/200 (the frozen-coverage bias is
    # first order in dt; the fine step is the headline, the coarse one shows the bias)
    for AR, across, dt, M in ((20, 16, 0.02, 100000), (50, 16, 0.05, 200000), (100, 16, 0.1, 200000),
                              (20, 16, 0.005, 25000), (50, 16, 0.0125, 50000), (100, 16, 0.025, 50000)):
        key = f"hole_AR{AR}_x{across}_dt{dt}_M{M}_s1"
        if key not in cache:
            cache[key] = run_hole_constant(AR, across, dt, M); save_cache(cache)
        if f"1d_AR{AR}" not in refs:
            refs[f"1d_AR{AR}"] = one_d_hole(AR); save_cache(cache)
        if f"mb_AR{AR}" not in refs:
            t0 = time.time(); refs[f"mb_AR{AR}"] = machball_hole(AR); refs[f"mb_AR{AR}"]["wall_s"] = time.time() - t0; save_cache(cache)
        r, d1, mb = cache[key], refs[f"1d_AR{AR}"], refs[f"mb_AR{AR}"]
        print(f"AR {AR} (dt {dt}): t90 voxel {r['t90']:.3f} / 1D {d1['t90']:.3f} / Machball {mb['t90']:.3f}; "
              f"t95 {r['t95']:.3f} / {d1['t95']:.3f} / {mb['t95']:.3f}", flush=True)


def part_convergence(cache):
    AR = 100
    runs = []
    for dt, M in ((0.025, 50000), (0.05, 100000), (0.2, 400000)):           # time step at equal total particles
        runs.append((f"hole_AR100_x16_dt{dt}_M{M}_s1", dict(across=16, dt=dt, M=M)))
    for across in (24, 32):                                                 # resolution at M 200000
        runs.append((f"hole_AR100_x{across}_dt0.1_M200000_s1", dict(across=across, dt=0.1, M=200000)))
    for M in (25000, 50000, 100000, 400000):                                # particles per step (200000 is the holes run)
        runs.append((f"hole_AR100_x16_dt0.1_M{M}_s1", dict(across=16, dt=0.1, M=M)))
    runs.append(("hole_AR100_x16_dt0.1_M200000_s2", dict(across=16, dt=0.1, M=200000, seed=2)))
    runs.append(("hole_AR100_x16_dt0.1_M200000_s1_deposit", dict(across=16, dt=0.1, M=200000, update="deposit")))
    runs.append(("hole_AR100_x16_dt0.1_M200000_s1_plainfaces", dict(across=16, dt=0.1, M=200000, area_correction=False)))
    for key, kw in runs:
        if key not in cache:
            cache[key] = run_hole_constant(AR, **kw); save_cache(cache)
        print(f"{key}: t90 {cache[key]['t90']:.3f} t95 {cache[key]['t95']:.3f} ({cache[key]['wall_s']:.0f} s)", flush=True)


def part_tapered(cache):
    refs = cache.setdefault("refs", {})
    key = "tapered_AR100_br0.6_x16_dt0.1_M200000_s1"
    if key not in cache:
        cache[key] = run_hole_constant(100, 16, 0.1, 200000, br=0.6, t_end=30.0); save_cache(cache)
    if "1d_AR100_br0.6" not in refs:
        refs["1d_AR100_br0.6"] = one_d_hole(100, 0.6); save_cache(cache)
    r, d1 = cache[key], refs["1d_AR100_br0.6"]
    print(f"tapered 0.6: t90 voxel {r['t90']:.3f} / 1D {d1['t90']:.3f}; t95 {r['t95']:.3f} / {d1['t95']:.3f}", flush=True)


def recess_faces(mc, g):
    """Face selections for the word-line slit: floor, deepest slit-wall slice, recess end walls per recess."""
    f = mc.faces()
    h = g.voxel
    y_end_lo, y_end_hi = g.meta["recess_end_y"]
    k_of = (f["z"] / h - 0.5).round().astype(int)
    k_of[f["direction"] == 4] -= 0                            # floor faces sit at the far face of the last gas voxel
    end_left = (f["direction"] == 3) & (np.abs(f["y"] - y_end_lo * h) < 0.5 * h)
    end_right = (f["direction"] == 2) & (np.abs(f["y"] - y_end_hi * h) < 0.5 * h)
    per_recess = []
    for (k0, k1) in g.meta["recess_k"]:
        inz = (f["z"] >= k0 * h) & (f["z"] <= k1 * h)
        per_recess.append(np.where((end_left | end_right) & inz)[0])
    return f, per_recess


def part_recess(cache):
    refs = cache.setdefault("refs", {})
    WL, OX, L_S = 29e-9, 22e-9, 5e-6
    geom = RecessedSlit(100e-9, L_S, WL, 365e-9, WL + OX)
    h = 5e-9
    g = voxel.recessed_slit(geom.width, geom.depth, geom.recess_gap, geom.recess_length, geom.recess_pitch, h, nx=4)
    key = "recess_w100_Lr365_h5_dt0.1_M100000_s1"
    if key not in cache:
        p = params()
        mc = voxel.VoxelMC(g, p)
        f, per_recess = recess_faces(mc, g)
        w = f["weight"]
        print(f"  voxel recessed slit {g.shape}, {mc.n_faces} faces, {g.meta['n_recesses']} recesses per wall, g_r {g.meta['g_r_actual']*1e9:.0f} nm, "
              f"p_r {g.meta['p_r_actual']*1e9:.0f} nm, L_r {g.meta['L_r_actual']*1e9:.0f} nm; end-wall faces per recess {len(per_recess[0])}", flush=True)
        dt, M = 0.1, 100000
        n = int(round(24.0 / dt)); t_series = np.linspace(0, n * dt, n + 1)
        hist_t, floor, wall_deep, rec_top, rec_mid, rec_deep, rec_min = [], [], [], [], [], [], []
        mid = len(per_recess) // 2
        for i in range(1, len(t_series)):
            mc.step(N0_GAS, dt, M)
            th = mc.face_coverage()
            z, prof, wt = mc.slice_coverage()
            kw = np.where(wt > 0)[0]
            ends = np.array([np.average(th[ix], weights=w[ix]) for ix in per_recess])
            hist_t.append(t_series[i]); floor.append(mc.bottom_coverage()); wall_deep.append(prof[kw.max()])
            rec_top.append(ends[0]); rec_mid.append(ends[mid]); rec_deep.append(ends[-1]); rec_min.append(ends.min())
            if i % 20 == 0:
                r = mc.rate()
                print(f"    t {t_series[i]:5.1f} s floor {floor[-1]:.3f} recess deep {rec_deep[-1]:.3f} min {rec_min[-1]:.3f}  "
                      f"{r['particles_per_s']:.0f} p/s wall {r['wall_s']:.0f} s", flush=True)
            if rec_min[-1] > 0.96 and floor[-1] > 0.96:
                break
        rate = mc.rate()
        rec = dict(case="w100_Lr365", kind="recess", across=g.meta["voxels_across"], voxel_nm=h * 1e9, dt=dt, M=M, seed=1,
                   mouth_density=N0_GAS, shape=list(g.shape), n_recesses=g.meta["n_recesses"], g_r_nm=g.meta["g_r_actual"] * 1e9,
                   p_r_nm=g.meta["p_r_actual"] * 1e9, L_r_nm=g.meta["L_r_actual"] * 1e9,
                   t90=voxel.time_to_level(hist_t, floor, 0.9), t95=voxel.time_to_level(hist_t, floor, 0.95),
                   t90_wall_deep=voxel.time_to_level(hist_t, wall_deep, 0.9), t90_recess_top=voxel.time_to_level(hist_t, rec_top, 0.9),
                   t90_recess_mid=voxel.time_to_level(hist_t, rec_mid, 0.9), t90_recess_deep=voxel.time_to_level(hist_t, rec_deep, 0.9),
                   t90_recess_min=voxel.time_to_level(hist_t, rec_min, 0.9), wall_s=rate["wall_s"], particles=rate["particles"],
                   threads=rate["threads"], pps=rate["particles_per_s"], pps_thread=rate["particles_per_s_per_thread"],
                   seg_per_particle=rate["segments_per_particle"], n_faces=mc.n_faces, mouth_area=mc.mouth_area, site_area=mc.site_area,
                   inventory_ratio=mc.site_area / (g.shape[0] * h * (2 * g.meta["L_actual"] * (1 - g.meta["g_r_actual"] / g.meta["p_r_actual"])
                                                                   + g.meta["n_recesses"] * 2 * (2 * g.meta["L_r_actual"] + g.meta["g_r_actual"]) + g.meta["w_actual"])),
                   counts=mc.counts(), t_hist=[float(x) for x in hist_t], bottom_hist=[float(x) for x in floor],
                   recess_deep_hist=[float(x) for x in rec_deep], recess_min_hist=[float(x) for x in rec_min], wall_deep_hist=[float(x) for x in wall_deep])
        cache[key] = rec; save_cache(cache)
    if "nested_w100_Lr365" not in refs:
        out = {}
        for label, gg in (("nominal", geom), ("rounded", RecessedSlit(g.meta["w_actual"], g.meta["L_actual"], g.meta["g_r_actual"], g.meta["L_r_actual"], g.meta["p_r_actual"]))):
            rB = simulate_recessed_constant(gg, PR, T, P_PA, 60.0, n_substeps=1200, n_samples=300, n_slit=200, n_recess=50)
            out[label] = dict(slit_bottom=time_to_coverage(rB.t, rB.hist_slit_bottom), recess_deep=time_to_coverage(rB.t, rB.hist_recess_bottom_deep),
                              recess_min=time_to_coverage(rB.t, rB.hist_recess_bottom_min), recess_top=time_to_coverage(rB.t, rB.hist_recess_bottom_top),
                              recess_mid=time_to_coverage(rB.t, rB.hist_recess_bottom_mid), area_factor=gg.area_factor,
                              t=[float(x) for x in rB.t], slit_bottom_hist=[float(x) for x in rB.hist_slit_bottom],
                              recess_deep_hist=[float(x) for x in rB.hist_recess_bottom_deep])
        refs["nested_w100_Lr365"] = out; save_cache(cache)
    r, nb = cache[key], refs["nested_w100_Lr365"]
    print(f"recess w100_Lr365 at 0.3 Pa: voxel floor t90 {r['t90']:.2f} s, deepest wall slice {r['t90_wall_deep']:.2f}, recess ends top/mid/deep "
          f"{r['t90_recess_top']:.2f}/{r['t90_recess_mid']:.2f}/{r['t90_recess_deep']:.2f} (min {r['t90_recess_min']:.2f}); nested model nominal "
          f"slit bottom {nb['nominal']['slit_bottom']:.2f}, recess deep {nb['nominal']['recess_deep']:.2f}; rounded dims {nb['rounded']['slit_bottom']:.2f} / "
          f"{nb['rounded']['recess_deep']:.2f}", flush=True)


def load_drive():
    rows = np.loadtxt(DRIVE, delimiter=",", skiprows=2)
    return rows[:, 0], rows[:, 1] * C0, rows[:, 2]


def part_driven(cache):
    refs = cache.setdefault("refs", {})
    t, n0, theta_bc = load_drive()
    depths = np.linspace(0, 1, 41)
    if "1d_driven" not in refs:
        proc = Process(Hole(D, 100 * D), PR, T=T, p_partial=1.0, n_cells=200)
        x, th, ht, hist = simulate_driven(proc, t, n0, substeps_per_sample=5, record_depths=depths)
        sel = {q: int(np.argmin(np.abs(ht - q))) for q in (5.0, 10.0, 13.0)}
        refs["1d_driven"] = dict(t90=voxel.time_to_level(ht, hist[:, -1], 0.9), t95=voxel.time_to_level(ht, hist[:, -1], 0.95),
                                 t90_bc_recorded=voxel.time_to_level(t, theta_bc, 0.9), depths=[float(d) for d in depths],
                                 profiles={str(q): [float(v) for v in hist[i]] for q, i in sel.items()},
                                 t=[float(v) for v in ht], bottom_hist=[float(v) for v in hist[:, -1]])
        save_cache(cache)
    for key, across, dtr, M in (("driven_AR100_x16_dt0.25_M400000_s1", 16, 0.25, 400000),
                                ("driven_AR100_x16_dt0.05_M80000_s1", 16, 0.05, 80000),
                                ("driven_AR100_x24_dt0.25_M400000_s1", 24, 0.25, 400000)):
        if key in cache:
            continue
        g = voxel.hole(D, 100 * D, D / across)
        print(f"  voxel driven hole {across} across {g.shape}, history resampled to {dtr} s, M {M}", flush=True)
        res = voxel.simulate_driven_voxel(g, params(), t, n0, M=M, resample_dt=dtr, record_depths=depths,
                                          snapshot_times=(5.0, 10.0, 13.0), t_stop=16.0, verbose=False)
        mc = res.solver
        z, prof, wt = mc.slice_coverage()
        rec = record(res, mc, dict(case="driven_AR100", kind="driven", AR=100, across=across, voxel_nm=D / across * 1e9, dt=dtr, M=M, seed=1,
                                  inventory_ratio=mc.site_area / smooth_area("hole", d=D, L=100 * D), area_correction=True, update="exponential"))
        rec["depths"] = [float(d) for d in depths]
        rec["profiles"] = {str(q): [float(v) for v in np.interp(depths * 100 * D, z[wt > 0], p_[wt > 0])] for q, p_ in res.snapshots.items()}
        rec["hist_depths"] = [[float(v) for v in row] for row in res.hist]
        cache[key] = rec; save_cache(cache)
        print(f"    t90 {rec['t90']:.3f} s (1D {refs['1d_driven']['t90']:.3f}); {rec['wall_s']:.0f} s wall, {rec['pps_thread']:.0f} p/s/thread", flush=True)



def part_slit(cache):
    """Slit (trench) of gap 100 nm, AR 50: axis-aligned walls have no staircase, so voxel vs Machball's Trench isolates
    the time-stepping and sampling from the area correction; the 1D slit model (D_K = 2 w vbar/3) is the third line."""
    refs = cache.setdefault("refs", {})
    W, AR = 100e-9, 50
    for dt, M in ((0.05, 100000), (0.0125, 25000)):
        key = f"slit_AR50_x16_dt{dt}_M{M}_s1"
        if key in cache:
            continue
        g = voxel.slit(W, AR * W, W / 16, nx=4)
        print(f"  voxel slit AR {AR} 16 across {g.shape}, dt {dt}, M {M}", flush=True)
        res = voxel.simulate_constant_voxel(g, params(), N0_GAS, 60.0, dt, M=M, stop_at_bottom=0.96)
        mc = res.solver
        smooth = g.shape[0] * (W / 16) * (2 * AR * W + W)
        cache[key] = record(res, mc, dict(case="slit_AR50", kind="slit", AR=AR, bottom_ratio=1.0, across=16, voxel_nm=W / 16 * 1e9, dt=dt, M=M, seed=1,
                                         mouth_density=N0_GAS, inventory_ratio=mc.site_area / smooth, area_correction=True, update="exponential"))
        save_cache(cache)
        print(f"    t90 {cache[key]['t90']:.3f} t95 {cache[key]['t95']:.3f} ({cache[key]['wall_s']:.0f} s, {cache[key]['pps_thread']:.0f} p/s/thread)", flush=True)
    if "1d_slit_AR50" not in refs:
        proc = Process(Slit(W, AR * W), PR, T=T, p_partial=P_PA, t_dose=1.0, n_cells=200, n_substeps=200, narrowing=False)
        refs["1d_slit_AR50"] = dict(t90=dose_for_coverage(proc, 0.90, 1e-3, 1e3), t95=dose_for_coverage(proc, 0.95, 1e-3, 1e3)); save_cache(cache)
    if "mb_slit_AR50" not in refs:
        refs["mb_slit_AR50"] = machball_hole(50, 400, "trench"); save_cache(cache)
    for dt, M in ((0.05, 100000), (0.0125, 25000)):
        r = cache[f"slit_AR50_x16_dt{dt}_M{M}_s1"]
        print(f"slit AR 50 (dt {dt}): t90 voxel {r['t90']:.3f} / 1D {refs['1d_slit_AR50']['t90']:.3f} / Machball {refs['mb_slit_AR50']['t90']:.3f}; "
              f"t95 {r['t95']:.3f} / {refs['1d_slit_AR50']['t95']:.3f} / {refs['mb_slit_AR50']['t95']:.3f}", flush=True)


# ------------------------------------------------------------------------------------------------------------------
def write_outputs(cache):
    refs = cache.get("refs", {})
    rows = []
    for key, r in cache.items():
        if key == "refs":
            continue
        row = {k: v for k, v in r.items() if not isinstance(v, (list, dict))}
        row["run"] = key
        row["escaped"], row["stuck"], row["lost"], row["left"] = (r["counts"][k] for k in ("escaped", "stuck", "lost", "left"))
        if r["kind"] in ("hole", "tapered"):
            ref = refs.get(f"1d_AR{r['AR']}" + (f"_br{r['bottom_ratio']}" if r["bottom_ratio"] < 1 else ""), {})
            row["t90_1d"] = ref.get("t90"); row["t95_1d"] = ref.get("t95")
            mb = refs.get(f"mb_AR{r['AR']}", {}) if r["bottom_ratio"] == 1 else {}
            row["t90_machball"] = mb.get("t90"); row["t95_machball"] = mb.get("t95")
        elif r["kind"] == "slit":
            row["t90_1d"] = refs.get("1d_slit_AR50", {}).get("t90"); row["t95_1d"] = refs.get("1d_slit_AR50", {}).get("t95")
            row["t90_machball"] = refs.get("mb_slit_AR50", {}).get("t90"); row["t95_machball"] = refs.get("mb_slit_AR50", {}).get("t95")
        elif r["kind"] == "driven":
            row["t90_1d"] = refs.get("1d_driven", {}).get("t90"); row["t95_1d"] = refs.get("1d_driven", {}).get("t95")
            row["t90_bc_recorded"] = refs.get("1d_driven", {}).get("t90_bc_recorded")
        elif r["kind"] == "recess":
            nb = refs.get("nested_w100_Lr365", {})
            for lab in ("nominal", "rounded"):
                if lab in nb:
                    row[f"nested_{lab}_slit_bottom_t90"] = nb[lab]["slit_bottom"]; row[f"nested_{lab}_recess_deep_t90"] = nb[lab]["recess_deep"]
                    row[f"nested_{lab}_recess_min_t90"] = nb[lab]["recess_min"]
        rows.append(row)
    cols = ["run", "case", "kind", "AR", "bottom_ratio", "across", "voxel_nm", "dt", "M", "seed", "update", "area_correction", "t90", "t95",
            "t90_1d", "t95_1d", "t90_machball", "t95_machball", "t90_bc_recorded", "t90_wall_deep", "t90_recess_top", "t90_recess_mid",
            "t90_recess_deep", "t90_recess_min", "nested_nominal_slit_bottom_t90", "nested_nominal_recess_deep_t90", "nested_nominal_recess_min_t90",
            "nested_rounded_slit_bottom_t90", "nested_rounded_recess_deep_t90", "nested_rounded_recess_min_t90", "inventory_ratio", "n_faces",
            "mouth_area", "site_area", "particles", "threads", "wall_s", "pps", "pps_thread", "seg_per_particle", "escaped", "stuck", "lost", "left",
            "n_recesses", "g_r_nm", "p_r_nm", "L_r_nm", "mouth_density"]
    for r in rows:
        for c in cols:
            r.setdefault(c, "")
    rows.sort(key=lambda r: (str(r["kind"]), str(r["case"]), r["run"]))
    with open(FIG / "voxel_mc_crosscheck.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    # driven profiles
    d1 = refs.get("1d_driven")
    drv = [k for k in cache if k.startswith("driven_")]
    if d1 and drv:
        with open(FIG / "voxel_mc_profiles.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            header = ["rel_depth"] + [f"theta_1d_{q}s" for q in ("5.0", "10.0", "13.0")]
            for k in drv:
                header += [f"theta_{k}_{q}s" for q in ("5.0", "10.0", "13.0")]
            w.writerow(header)
            for i, d in enumerate(d1["depths"]):
                row = [f"{d:.3f}"] + [f"{d1['profiles'][q][i]:.4f}" for q in ("5.0", "10.0", "13.0")]
                for k in drv:
                    row += [f"{cache[k]['profiles'][q][i]:.4f}" if q in cache[k]["profiles"] else "" for q in ("5.0", "10.0", "13.0")]
                w.writerow(row)
    plot(cache)
    print(f"wrote {FIG / 'voxel_mc_crosscheck.csv'} ({len(rows)} runs), voxel_mc_crosscheck.png" + (", voxel_mc_profiles.csv" if d1 and drv else ""))


def plot(cache):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    refs = cache.get("refs", {})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    # (1) straight holes: t90 and t95 voxel / 1D / Machball
    ax = axes[0, 0]
    ARs = [20, 50, 100]; xs = np.arange(3); wbar = 0.25
    for off, (lab, getter, col) in enumerate((("voxel MC (16 across, dt = t90/200)", lambda AR, lv: cache.get({20: "hole_AR20_x16_dt0.005_M25000_s1", 50: "hole_AR50_x16_dt0.0125_M50000_s1", 100: "hole_AR100_x16_dt0.025_M50000_s1"}[AR], {}).get(lv, np.nan), "tab:blue"),
                                                ("1D Knudsen (this repo)", lambda AR, lv: refs.get(f"1d_AR{AR}", {}).get(lv, np.nan), "tab:orange"),
                                                ("Machball (ballistic)", lambda AR, lv: refs.get(f"mb_AR{AR}", {}).get(lv, np.nan), "tab:green"))):
        v90 = [getter(AR, "t90") for AR in ARs]; v95 = [getter(AR, "t95") for AR in ARs]
        ax.bar(xs + (off - 1) * wbar, v95, wbar, color=col, alpha=0.35)
        ax.bar(xs + (off - 1) * wbar, v90, wbar, color=col, label=lab)
    ax.set_xticks(xs); ax.set_xticklabels([f"AR {a}" for a in ARs]); ax.set_yscale("log"); ax.set_ylabel("time to 90 % (solid) / 95 % (pale) floor coverage, s")
    ax.set_title("Straight 100 nm holes at 0.3 Pa TMA, 573 K, s = 5e-3", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3, which="both")
    # (2) convergence with M, resolution, dt
    ax = axes[0, 1]
    ref = refs.get("1d_AR100", {}).get("t90", np.nan); mb = refs.get("mb_AR100", {}).get("t90", np.nan)
    Ms, tM = [], []
    for M in (25000, 50000, 100000, 200000, 400000):
        k = f"hole_AR100_x16_dt0.1_M{M}_s1"
        if k in cache: Ms.append(M); tM.append(cache[k]["t90"])
    ax.semilogx(Ms, tM, "o-", color="tab:blue", label="16 across, dt 0.1 s, vs particles per step")
    for across, mk in ((24, "s"), (32, "^")):
        k = f"hole_AR100_x{across}_dt0.1_M200000_s1"
        if k in cache: ax.plot([200000], [cache[k]["t90"]], mk, color="tab:red", ms=8, label=f"{across} across, M 200000")
    for dt, M, mk in ((0.025, 50000, "<"), (0.05, 100000, "v"), (0.2, 400000, "P")):
        k = f"hole_AR100_x16_dt{dt}_M{M}_s1"
        if k in cache: ax.plot([M], [cache[k]["t90"]], mk, color="tab:purple", ms=8, label=f"dt {dt} s, M {M}")
    for key, lab, mk, col in (("hole_AR100_x16_dt0.1_M200000_s2", "seed 2", "D", "tab:cyan"), ("hole_AR100_x16_dt0.1_M200000_s1_deposit", "deposit update", "x", "k"),
                              ("hole_AR100_x16_dt0.1_M200000_s1_plainfaces", "no area correction", "*", "tab:brown")):
        if key in cache: ax.plot([200000], [cache[key]["t90"]], mk, color=col, ms=9, label=lab)
    ax.axhline(ref, color="tab:orange", ls="--", label=f"1D Knudsen {ref:.2f} s"); ax.axhline(mb, color="tab:green", ls=":", label=f"Machball {mb:.2f} s")
    ax.set_xlabel("particles per step M"); ax.set_ylabel("AR-100 floor t90, s"); ax.set_title("Convergence (AR 100, constant 0.3 Pa)", fontsize=9)
    ax.legend(fontsize=6, ncol=2); ax.grid(alpha=0.3)
    # (3) histories: AR 100 voxel vs 1D (constant); tapered
    ax = axes[0, 2]
    for key, lab, col in (("hole_AR100_x16_dt0.1_M200000_s1", "AR 100 straight, voxel", "tab:blue"), ("tapered_AR100_br0.6_x16_dt0.1_M200000_s1", "AR 100 frustum 60 %, voxel", "tab:red")):
        if key in cache: ax.plot(cache[key]["t_hist"], cache[key]["bottom_hist"], color=col, label=lab)
    for rk, lab, col in (("1d_AR100", "1D straight t90", "tab:orange"), ("1d_AR100_br0.6", "1D frustum t90", "tab:brown")):
        if rk in refs: ax.axvline(refs[rk]["t90"], color=col, ls="--", lw=0.9, label=f"{lab} {refs[rk]['t90']:.2f} s")
    if "mb_AR100" in refs: ax.axvline(refs["mb_AR100"]["t90"], color="tab:green", ls=":", label=f"Machball t90 {refs['mb_AR100']['t90']:.2f} s")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlabel("dose time, s"); ax.set_ylabel("floor coverage"); ax.set_title("Floor coverage at constant 0.3 Pa", fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # (4) recess
    ax = axes[1, 0]
    rk = "recess_w100_Lr365_h5_dt0.1_M100000_s1"
    if rk in cache:
        r = cache[rk]
        ax.plot(r["t_hist"], r["bottom_hist"], color="tab:blue", label=f"voxel slit floor (t90 {r['t90']:.2f} s)")
        ax.plot(r["t_hist"], r["recess_deep_hist"], color="tab:red", label=f"voxel deepest recess end (t90 {r['t90_recess_deep']:.2f} s)")
        ax.plot(r["t_hist"], r["recess_min_hist"], color="tab:red", ls=":", label=f"voxel slowest recess end (t90 {r['t90_recess_min']:.2f} s)")
    nb = refs.get("nested_w100_Lr365", {})
    if "nominal" in nb:
        ax.plot(nb["nominal"]["t"], nb["nominal"]["slit_bottom_hist"], color="tab:orange", ls="--", label=f"nested model slit bottom ({nb['nominal']['slit_bottom']:.2f} s)")
        ax.plot(nb["nominal"]["t"], nb["nominal"]["recess_deep_hist"], color="tab:brown", ls="--", label=f"nested model deepest recess ({nb['nominal']['recess_deep']:.2f} s)")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlim(0, 20); ax.set_xlabel("dose time, s"); ax.set_ylabel("coverage")
    ax.set_title("Word-line slit w100_Lr365 at constant 0.3 Pa (voxel 5 nm: g_r 30, p_r 50 nm)", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # (5) driven histories
    ax = axes[1, 1]
    d1 = refs.get("1d_driven")
    if d1: ax.plot(d1["t"], d1["bottom_hist"], color="tab:orange", label=f"1D simulate_driven, t90 {d1['t90']:.2f} s")
    for key, lab, col in (("driven_AR100_x16_dt0.25_M400000_s1", "voxel 16 across, 0.25 s, M 4e5", "tab:blue"), ("driven_AR100_x16_dt0.05_M80000_s1", "voxel 16 across, 0.05 s, M 8e4", "tab:cyan"),
                          ("driven_AR100_x24_dt0.25_M400000_s1", "voxel 24 across, 0.25 s, M 4e5", "tab:purple")):
        if key in cache: ax.plot(cache[key]["t_hist"], cache[key]["bottom_hist"], color=col, label=f"{lab}: t90 {cache[key]['t90']:.2f} s")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlim(0, 16); ax.set_xlabel("time, s"); ax.set_ylabel("hole-bottom coverage")
    ax.set_title("AR-100 hole driven by the chamber's centre-face density", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # (6) driven profiles
    ax = axes[1, 2]
    if d1:
        for q, col in (("5.0", "tab:blue"), ("10.0", "tab:orange"), ("13.0", "tab:green")):
            ax.plot(d1["depths"], d1["profiles"][q], color=col, ls="--", label=f"1D at {float(q):.0f} s")
            k = "driven_AR100_x16_dt0.25_M400000_s1"
            if k in cache and q in cache[k]["profiles"]:
                ax.plot(cache[k]["depths"], cache[k]["profiles"][q], color=col, label=f"voxel at {float(q):.0f} s")
    ax.set_xlabel("relative depth"); ax.set_ylabel("coverage"); ax.set_title("Coverage along the depth, chamber-driven", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.suptitle("3D voxel Monte Carlo feature solver (lib/aldmc) against the 1D Knudsen solver and Machball", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "voxel_mc_crosscheck.png", dpi=150)


PARTS = dict(holes=part_holes, convergence=part_convergence, tapered=part_tapered, recess=part_recess, driven=part_driven, slit=part_slit, machball=part_machball)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="", help="comma-separated: " + ",".join(PARTS))
    ap.add_argument("--figure", action="store_true", help="only assemble the CSV and figure from the cache")
    ap.add_argument("--threads", type=int, default=0)
    a = ap.parse_args()
    global THREADS
    if a.threads > 0:
        THREADS = a.threads
    cache = load_cache()
    if not a.figure:
        for name in [p for p in a.part.split(",") if p]:
            t0 = time.time()
            print(f"== {name} ({THREADS} threads)", flush=True)
            PARTS[name](cache)
            print(f"== {name} done in {time.time() - t0:.0f} s", flush=True)
    write_outputs(cache)


if __name__ == "__main__":
    main()
