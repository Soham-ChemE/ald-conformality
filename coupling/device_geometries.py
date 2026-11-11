"""Work package 10: the 3D voxel Monte Carlo feature solver (lib/aldmc, aldconf.voxel) on two device geometries from
the literature: gate-all-around (GAA) nanosheet stacks at the replacement-metal-gate step and a 3D-NAND word-line
cell (gate-line slit, word-line recesses and the channel-hole pillars standing in them), with the chamber's TMA,
HfO2/TDMAH and TiN/TiCl4 parameters, at a constant 0.3 Pa mouth density and under the chamber's centre-face density
history (the long AR-100 case, as in coupling/voxel_mc_crosscheck.py part `driven`).

Geometries (occupancy arrays through voxel.from_array; every dimension is rounded to whole voxels and the actual
values are recorded in `meta`):

  gaa_stack   unit cell of a nanosheet stack after channel release: N sheets of thickness t_s at vertical spacing g_s
              and width W, suspended g_bot above the substrate, under a gate trench of length Lg (along y, closed by
              the gate spacers / inner spacers, which are solid) and height H_top above the top sheet (the mouth is
              the top of the trench); periodic along x (the sheet-width direction) with period W + gap_x, gap_x
              being the gas gap to the neighbouring stack. Gas enters at the trench top, descends beside the stack
              and fills the inter-sheet slots from their two open edges (the sheet-width edges); the slot ends along
              y are the inner spacers.
  nand_cell   the word-line structure of aldconf.recess / voxel.recessed_slit (slit of width w_s and depth L_s with
              lateral recesses of height g_r, length L_r, pitch p_r on both walls) plus the channel-hole pillars:
              solid cylinders of diameter d_h running the full depth, in `rows` staggered rows per side at a pitch
              `pitch_h` along the slit (the x period). The pillars are what the nested 1D model cannot carry.

Runs (cached in figures/device_geometries_runs.json so a part can be rerun alone; `--figure` assembles the CSV and
figures from the cache):

    PYTHONPATH=src .venv/bin/python coupling/device_geometries.py --part gaa,nand5,nand13,refs
    PYTHONPATH=src .venv/bin/python coupling/device_geometries.py --figure

Outputs: figures/device_geometries.png (the geometries), figures/gaa_coverage.png (coverage maps and histories),
figures/device_geometries.csv (every number), figures/device_geometries_runs.json (cache with histories).
Threads: ALDMC_THREADS (default min(8, cpu)); every run here is designed for well under 20 minutes on 6 threads.
"""
import sys, json, time, argparse, csv, os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from aldconf import voxel, mean_speed, Precursor  # noqa: E402
from aldconf.kinetics import number_density  # noqa: E402
from aldconf.recess import RecessedSlit, simulate_recessed, simulate_recessed_constant, time_to_coverage  # noqa: E402

FIG = ROOT / "figures"
CACHE = FIG / "device_geometries_runs.json"
DRIVE = ROOT / "lib/aldkin/examples/centre_face_drive.csv"

T = 573.15                                      # K, the chamber temperature of every run in the paper
P_PA = 0.3                                      # Pa, the constant-mouth benchmark (Machball check conditions)
N0_GAS = float(number_density(P_PA, T))         # 3.79e19 m^-3
C0 = 3.79e20                                    # 1 % precursor at 300 Pa, 573 K (the chamber runs' c0)
U_C0 = 1.90e19                                  # molecules m^-2 s^-1 at 1 %: U = 0.050 m/s times c0 (wp4, wp2)

# chemistries: s, N0 [sites/m2], mass [amu], growth per cycle [nm] (wp4 Table 3; wp2 for Al2O3)
CHEM = {
    "TMA":   dict(s=5e-3, N0=4.0e18,  mass=72.09,  gpc_nm=0.103, film="Al2O3"),
    "TDMAH": dict(s=7e-3, N0=2.77e18, mass=354.7,  gpc_nm=0.10,  film="HfO2"),
    "TiCl4": dict(s=0.10, N0=1.24e18, mass=189.7,  gpc_nm=0.030, film="TiN"),
}
THREADS = voxel.default_threads()


def vbar_of(chem):
    return float(mean_speed(T, CHEM[chem]["mass"]))


def params(chem, seed=1, **kw):
    c = CHEM[chem]
    base = dict(s=c["s"], N0=c["N0"], vbar=vbar_of(chem), seed=seed, threads=THREADS)
    base.update(kw)
    return voxel.Params(**base)


def rule_pulse(chem, wall_ratio):
    """The paper's design rule 1.2 (1 + A_wall/A_wafer) N0/(U c0) at 1 % precursor, in seconds."""
    return 1.2 * (1.0 + wall_ratio) * CHEM[chem]["N0"] / U_C0


def effusion_limit(chem, wall_ratio_cell):
    """Stoichiometric minimum exposure time at the mouth at 0.3 Pa: (A_wall/A_mouth) N0 / (n vbar/4)."""
    return wall_ratio_cell * CHEM[chem]["N0"] / (N0_GAS * vbar_of(chem) / 4.0)


# -------------------------------------------------------------------------------------------------------------------
# geometry builders
def _nvox(length, voxel, minimum=1):
    return max(minimum, int(round(length / voxel)))


def gaa_stack(n_sheets, t_s, g_s, W, Lg, gap_x, H_top, g_bot, h, pad=1, **extra):
    """GAA nanosheet unit cell (see the module docstring). Returns a VoxelGeometry periodic along x with meta
    recording the actual (rounded) dimensions and the index ranges used by gaa_regions."""
    nt, ng, nW, nLg = _nvox(t_s, h), _nvox(g_s, h), _nvox(W, h), _nvox(Lg, h)
    ngap, nH, nb = _nvox(gap_x, h), _nvox(H_top, h), _nvox(g_bot, h)
    nx = nW + ngap
    ny = nLg + 2 * pad
    nz = nH + n_sheets * nt + (n_sheets - 1) * ng + nb + 1
    occ = np.ones((nx, ny, nz), dtype=np.uint8)
    occ[:, pad:pad + nLg, :nz - 1] = 0                                   # the gate trench, floor at the last slice
    x0 = ngap // 2
    sheets = []
    k = nH
    for i in range(n_sheets):
        occ[x0:x0 + nW, pad:pad + nLg, k:k + nt] = 1                      # a sheet (solid Si)
        sheets.append((k, k + nt))
        k += nt + ng
    gaps = [(sheets[i][1], sheets[i + 1][0]) for i in range(n_sheets - 1)] + [(sheets[-1][1], nz - 1)]
    meta = dict(kind="gaa", n_sheets=n_sheets, t_s=t_s, g_s=g_s, W=W, Lg=Lg, gap_x=gap_x, H_top=H_top, g_bot=g_bot,
                t_s_actual=nt * h, g_s_actual=ng * h, W_actual=nW * h, Lg_actual=nLg * h,
                gap_x_actual=ngap * h, H_top_actual=nH * h, g_bot_actual=nb * h, pitch_x=nx * h,
                sheet_x=(x0, x0 + nW), trench_y=(pad, pad + nLg), sheet_k=sheets, gap_k=gaps, pad=pad, **extra)
    return voxel.from_array(occ, h, periodic=(True, False), **meta)


def gaa_smooth_area(meta):
    """Exposed area of the smooth cell (m2): trench spacer walls, substrate floor, sheet faces and edges."""
    m = meta
    depth = m["H_top_actual"] + m["n_sheets"] * m["t_s_actual"] + (m["n_sheets"] - 1) * m["g_s_actual"] + m["g_bot_actual"]
    spacers = 2 * m["pitch_x"] * depth - 2 * m["n_sheets"] * m["W_actual"] * m["t_s_actual"]   # the sheets abut the spacers
    floor = m["pitch_x"] * m["Lg_actual"]
    sheets = m["n_sheets"] * (2 * m["W_actual"] * m["Lg_actual"] + 2 * m["t_s_actual"] * m["Lg_actual"])
    return spacers + floor + sheets


def _disc_mask(nx, ny, cx, cy, d_vox, periodic_x=True):
    """Boolean (nx, ny) mask of a digital disc of diameter d_vox voxels centred at (cx, cy) (voxel units, centre
    coordinates), with exactly round(pi r^2) voxels (as voxel.hole), wrapping in x."""
    x = np.arange(nx) + 0.5
    y = np.arange(ny) + 0.5
    dx = x[:, None] - cx
    if periodic_x:
        dx = dx - nx * np.round(dx / nx)
    dy = y[None, :] - cy
    r2 = dx ** 2 + dy ** 2
    count = int(round(np.pi * (d_vox / 2.0) ** 2))
    order = np.argsort(r2.ravel(), kind="stable")[:count]
    mask = np.zeros(nx * ny, dtype=bool)
    mask[order] = True
    return mask.reshape(nx, ny)


def nand_cell(w_s, L_s, g_r, L_r, p_r, d_h, pitch_h, rows, h, pad=1, oxide_first=True, **extra):
    """Word-line slit with recesses on both walls (voxel.recessed_slit's construction) and channel-hole pillars:
    `rows` staggered rows of solid cylinders (diameter d_h) per side, the rows at y = L_r (j + 1/2)/rows into the
    recess, row j shifted by half the period along the slit; the x period is pitch_h."""
    nw, nL, ng, nLr = _nvox(w_s, h), _nvox(L_s, h), _nvox(g_r, h), _nvox(L_r, h)
    npitch = _nvox(p_r, h, ng + 1)
    nx = _nvox(pitch_h, h)
    ny = nw + 2 * nLr + 2 * pad
    y0 = pad + nLr
    occ = np.ones((nx, ny, nL + 1), dtype=np.uint8)
    occ[:, y0:y0 + nw, :nL] = 0
    ox = npitch - ng
    k = 0
    recesses = []
    while k < nL:
        start = k + (ox if oxide_first else 0)
        end = min(start + ng, nL)
        if start >= nL:
            break
        if end - start == ng:
            occ[:, pad:y0, start:end] = 0
            occ[:, y0 + nw:y0 + nw + nLr, start:end] = 0
            recesses.append((start, end))
        k += npitch
    # pillars
    dvox = d_h / h
    pillars = []
    pmask = np.zeros((nx, ny), dtype=bool)
    for j in range(rows):
        cy_rel = nLr * (j + 0.5) / rows                                   # distance into the recess (voxels)
        cx = (nx / 2.0) if j % 2 == 0 else 0.0
        for side in (0, 1):
            cy = (y0 - cy_rel) if side == 0 else (y0 + nw + cy_rel)
            m = _disc_mask(nx, ny, cx, cy, dvox)
            pmask |= m
            pillars.append((cx * h, cy * h))
    occ[pmask, :] = 1
    meta = dict(kind="nand", w_s=w_s, L_s=L_s, g_r=g_r, L_r=L_r, p_r=p_r, d_h=d_h, pitch_h=pitch_h, rows=rows,
                w_actual=nw * h, L_actual=nL * h, g_r_actual=ng * h, p_r_actual=npitch * h,
                L_r_actual=nLr * h, d_h_voxels=dvox, pitch_h_actual=nx * h, n_recesses=len(recesses),
                recess_k=recesses, pad=pad, y_slit=(y0, y0 + nw), recess_end_y=(pad, y0 + nw + nLr),
                pillar_xy=pillars, pillar_voxels_per_slice=int(pmask.sum()), **extra)
    return voxel.from_array(occ, h, periodic=(True, False), **meta)


def nand_smooth_area(meta):
    """Smooth exposed area of the cell (m2): slit walls between recesses, slit floor, recess floors and ceilings
    (minus the pillar footprints), recess end walls and the pillar walls inside the recesses."""
    m = meta
    h_cell = m["pitch_h_actual"]
    n = m["n_recesses"]
    slit_wall = 2 * h_cell * (m["L_actual"] - n * m["g_r_actual"])
    floor = h_cell * m["w_actual"]
    disc = np.pi * m["d_h"] ** 2 / 4
    recess_plates = n * 2 * 2 * (h_cell * m["L_r_actual"] - m["rows"] * disc)        # two sides, floor + ceiling, minus the pillar footprints
    recess_ends = n * 2 * h_cell * m["g_r_actual"]
    pillar_wall = n * m["g_r_actual"] * len(m["pillar_xy"]) * np.pi * m["d_h"]       # pillar walls inside every recess
    return slit_wall + floor + recess_plates + recess_ends + pillar_wall


# -------------------------------------------------------------------------------------------------------------------
# face regions (index arrays into the solver's face list, area-weighted means are taken over them)
def gaa_regions(mc, g):
    f = mc.faces()
    h = g.voxel
    m = g.meta
    x, y, z, d, w = f["x"], f["y"], f["z"], f["direction"], f["weight"]
    xs0, xs1 = (v * h for v in m["sheet_x"])
    xmid = 0.5 * (xs0 + xs1)
    W = xs1 - xs0
    regions = {}
    # top of the top sheet
    k0, k1 = m["sheet_k"][0]
    regions["sheet0_top"] = np.where((d == 4) & (np.abs(z - k0 * h) < 0.5 * h) & (x > xs0) & (x < xs1))[0]
    for i, (ka, kb) in enumerate(m["gap_k"]):
        name = f"gap{i}" if i < len(m["gap_k"]) - 1 else "gap_bottom"
        inz = (z > ka * h - 0.5 * h) & (z < kb * h + 0.5 * h) & (x > xs0) & (x < xs1)
        plates = inz & ((d == 4) | (d == 5))                             # the sheet undersides / tops facing the slot
        regions[name + "_all"] = np.where(plates)[0]
        regions[name + "_mid"] = np.where(plates & (np.abs(x - xmid) < 0.1 * W))[0]
        regions[name + "_edge"] = np.where(plates & (np.abs(x - xmid) > 0.4 * W))[0]
    # stack-to-stack gas gap: spacer walls there are flat-like; the substrate floor under the gap
    kf = g.shape[2] - 1
    regions["floor_under_stack"] = np.where((d == 4) & (np.abs(z - kf * h) < 0.5 * h) & (x > xs0) & (x < xs1))[0]
    regions["floor_between_stacks"] = np.where((d == 4) & (np.abs(z - kf * h) < 0.5 * h) & ((x < xs0) | (x > xs1)))[0]
    regions = {k: v for k, v in regions.items() if len(v) > 0}
    return f, regions


def nand_regions(mc, g):
    f = mc.faces()
    h = g.voxel
    m = g.meta
    x, y, z, d, w = f["x"], f["y"], f["z"], f["direction"], f["weight"]
    ya, yb = (v * h for v in m["y_slit"])
    y_lo, y_hi = (v * h for v in m["recess_end_y"])
    Lr = m["L_r_actual"]
    kf = g.shape[2] - 1
    regions = {}
    regions["slit_floor"] = np.where((d == 4) & (np.abs(z - kf * h) < 0.5 * h) & (y > ya) & (y < yb))[0]
    in_slit = (y > ya) & (y < yb)
    deep_left = (y < ya) & (y < y_lo + 0.5 * Lr)
    deep_right = (y > yb) & (y > y_hi - 0.5 * Lr)
    deep = deep_left | deep_right
    n = len(m["recess_k"])
    picks = {"top": 0, "mid": n // 2, "deep": n - 1}
    ends = np.zeros(len(x), dtype=bool)
    ends |= (d == 3) & (np.abs(y - y_lo) < 0.5 * h)
    ends |= (d == 2) & (np.abs(y - y_hi) < 0.5 * h)
    rec_all = []
    for i, (ka, kb) in enumerate(m["recess_k"]):
        inz = (z > ka * h - 0.6 * h) & (z < kb * h + 0.6 * h) & ~in_slit
        rec_all.append(np.where(inz & deep)[0])
        for lab, idx in picks.items():
            if i == idx:
                regions[f"recess_{lab}_deephalf"] = np.where(inz & deep)[0]
                regions[f"recess_{lab}_end"] = np.where(inz & ends)[0]
    regions["recess_all_deephalf"] = rec_all
    return f, regions


# -------------------------------------------------------------------------------------------------------------------
def load_cache():
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def save_cache(c):
    CACHE.write_text(json.dumps(c, indent=1))


def load_drive():
    rows = np.loadtxt(DRIVE, delimiter=",", skiprows=2)
    return rows[:, 0], rows[:, 1] * C0


def drive_series(drive, dt, t_end):
    """(t, n0) sampled every dt: constant 0.3 Pa, or the chamber's centre-face history interpolated (its last
    sample, 22 s, is held constant beyond)."""
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    if drive == "constant":
        return t, np.full(n + 1, N0_GAS)
    td, nd = load_drive()
    return t, np.interp(t, td, nd, right=nd[-1])


def wmean(theta, w, idx):
    return float(np.average(theta[idx], weights=w[idx])) if len(idx) else float("nan")


def run_cell(g, chem, drive, dt, M, t_end, regions_fn, seed=1, stop_level=0.995, snapshot_every=0, min_steps=5,
             verbose_every=50, label="", max_wall_s=1080.0, **pkw):
    """Step the solver with the drive, record area-weighted coverages of the regions after every step, stop when
    every tracked region has passed stop_level (or at t_end). Returns the record (histories included)."""
    p = params(chem, seed=seed, **pkw)
    t0 = time.perf_counter()
    mc = voxel.VoxelMC(g, p)
    t_create = time.perf_counter() - t0
    f, regions = regions_fn(mc, g)
    w = f["weight"]
    track = {k: v for k, v in regions.items() if not isinstance(v, list)}
    lists = {k: v for k, v in regions.items() if isinstance(v, list)}
    t, n0 = drive_series(drive, dt, t_end)
    hist = {k: [] for k in track}
    hist_min_list = {k: [] for k in lists}
    hist_lists = {k: [] for k in lists}                                        # per member (e.g. per recess) histories
    hist_t = []
    snaps = {}
    area_below = []
    print(f"  {label}: {g.shape} {mc.n_faces} faces, {len(track)} regions, chem {chem}, drive {drive}, dt {dt}, M {M}, "
          f"create {t_create:.1f} s, {THREADS} threads", flush=True)
    for i in range(1, len(t)):
        nmid = float(np.interp(0.5 * (t[i] + t[i - 1]), t, n0))
        mc.step(max(0.0, nmid), dt, M)
        th = mc.face_coverage()
        hist_t.append(float(t[i]))
        for k, idx in track.items():
            hist[k].append(wmean(th, w, idx))
        for k, lst in lists.items():
            vals = np.array([wmean(th, w, idx) for idx in lst])
            hist_min_list[k].append(float(np.nanmin(vals)) if len(vals) else float("nan"))
            hist_lists[k].append([float(v) for v in vals])
        area_below.append(float(w[th < 0.9].sum() / w.sum()))
        if snapshot_every and i % snapshot_every == 0:
            snaps[float(t[i])] = th.copy()
        if verbose_every and i % verbose_every == 0:
            r = mc.rate()
            worst = min(v[-1] for v in hist.values())
            print(f"    t {t[i]:8.3f} s  worst region {worst:.4f}  area<0.9 {area_below[-1]:.4f}  {r['particles_per_s']:.0f} p/s  "
                  f"{r['segments_per_particle']:.0f} seg/p  wall {r['wall_s']:.0f} s", flush=True)
        done = all(v[-1] >= stop_level for v in hist.values()) and all(v[-1] >= stop_level for v in hist_min_list.values())
        if done and i >= min_steps:
            break
        if mc.wall_time + t_create > max_wall_s:                               # the 20-minute cap of the package
            print(f"    wall-clock cap {max_wall_s:.0f} s reached at t = {t[i]:.3f} s", flush=True)
            stopped_by_wall = True
            break
    th = mc.face_coverage()
    rate = mc.rate()
    stopped_by_wall = locals().get("stopped_by_wall", False)
    # the last snapshot is the final state
    snaps[float(hist_t[-1])] = th.copy()
    rec = dict(label=label, chem=chem, drive=drive, dt=dt, M=M, seed=seed, shape=list(g.shape), voxel_nm=g.voxel * 1e9,
               n_faces=mc.n_faces, mouth_area=mc.mouth_area, site_area=mc.site_area, threads=rate["threads"],
               wall_s=rate["wall_s"] + t_create, create_s=t_create, particles=rate["particles"], pps=rate["particles_per_s"],
               pps_thread=rate["particles_per_s_per_thread"], seg_per_particle=rate["segments_per_particle"],
               counts=mc.counts(), t_hist=hist_t, hist=hist, hist_min_list=hist_min_list, area_below_0p9=area_below,
               region_sizes={k: int(len(v)) for k, v in track.items()}, t_end_reached=hist_t[-1], stopped_by_wall=bool(stopped_by_wall),
               area_correction=bool(p.area_correction), normal_radius=int(p.normal_radius),
               meta={k: (v if not isinstance(v, (np.ndarray,)) else v.tolist()) for k, v in g.meta.items()})
    # time to 90 / 99 % of every region and of the slowest region ("everywhere")
    tt = np.array(hist_t)
    rec["t90"] = {k: voxel.time_to_level(tt, v, 0.9) for k, v in hist.items()}
    rec["t99"] = {k: voxel.time_to_level(tt, v, 0.99) for k, v in hist.items()}
    allmin = np.min(np.array([hist[k] for k in hist] + [hist_min_list[k] for k in hist_min_list]), axis=0)
    rec["hist_everywhere_min"] = [float(v) for v in allmin]
    rec["t90_everywhere"] = voxel.time_to_level(tt, allmin, 0.9)
    rec["t99_everywhere"] = voxel.time_to_level(tt, allmin, 0.99)
    rec["last_region_90"] = max(rec["t90"], key=lambda k: (rec["t90"][k] if np.isfinite(rec["t90"][k]) else 1e9))
    # per-member times of the list regions (per recess): where along the list the slowest member sits
    rec["t90_list"], rec["t99_list"], rec["slowest_member"] = {}, {}, {}
    for k, hl in hist_lists.items():
        H = np.array(hl)
        t90s = [voxel.time_to_level(tt, H[:, j], 0.9) for j in range(H.shape[1])]
        t99s = [voxel.time_to_level(tt, H[:, j], 0.99) for j in range(H.shape[1])]
        rec["t90_list"][k] = t90s; rec["t99_list"][k] = t99s
        fin = [v if np.isfinite(v) else -1.0 for v in t90s]
        j = int(np.argmax(fin))
        rec["slowest_member"][k] = dict(index=j, of=len(t90s), fraction=(j + 0.5) / len(t90s), t90=t90s[j], t99=t99s[j],
                                        last_t90=t90s[-1], last_t99=t99s[-1])
    rec["final_min_region"] = float(allmin[-1])
    rec["final_area_below_0p9"] = area_below[-1]
    rec["snapshots"] = snaps
    rec["faces"] = f
    rec["regions"] = track
    rec["solver"] = mc
    print(f"    done t {hist_t[-1]:.3f} s: t90/t99 everywhere {rec['t90_everywhere']:.3f}/{rec['t99_everywhere']:.3f} s, last region "
          f"{rec['last_region_90']}, {rate['wall_s']:.0f} s wall, {rate['particles_per_s_per_thread']:.0f} p/s/thread, "
          f"{rate['segments_per_particle']:.0f} seg/p", flush=True)
    return rec


def strip(rec):
    """JSON-serialisable copy (no faces, solver or snapshots; snapshots are stored separately where wanted)."""
    out = {k: v for k, v in rec.items() if k not in ("snapshots", "faces", "regions", "solver")}
    out["counts"] = {k: int(v) for k, v in rec["counts"].items()}
    out["meta"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in rec["meta"].items()}
    return out


# -------------------------------------------------------------------------------------------------------------------
# the cells
GAA_CELLS = {
    # IBM 2017 (Loubet et al., VLSI 2017): 3 sheets of 5 nm at 10 nm spacing, Lg 12 nm; W 40 nm (IBM, Semiconductor
    # Engineering figure caption); gap to the neighbouring stack 20 nm, trench height 40 nm above the top sheet and
    # 10 nm below the bottom sheet are this work's assumptions (stated in wp10).
    "gaa_ibm2017": dict(n_sheets=3, t_s=5e-9, g_s=10e-9, W=40e-9, Lg=12e-9, gap_x=20e-9, H_top=40e-9, g_bot=10e-9),
    # scaled: 4 sheets (Intel 18A has four ribbons), spacing 7 nm (the 7-8 nm minimum of Wong and Kakushima 2022),
    # W 25 nm (inside the 19-50 nm teardown range), Lg 12 nm, gap 15 nm
    "gaa_scaled4": dict(n_sheets=4, t_s=5e-9, g_s=7e-9, W=25e-9, Lg=12e-9, gap_x=15e-9, H_top=40e-9, g_bot=7e-9),
    # the IBM cell after 1 nm interfacial oxide and 2 nm HfO2 on every surface (no narrowing in the solver, so the
    # TiN step is run on the narrowed geometry explicitly): slots 4 nm, trench 6 nm, sheets 11 nm, gap 14 nm
    "gaa_ibm2017_after_hk": dict(n_sheets=3, t_s=11e-9, g_s=4e-9, W=46e-9, Lg=6e-9, gap_x=14e-9, H_top=37e-9, g_bot=4e-9),
}
GAA_VOXEL = {"gaa_ibm2017": 1.25e-9, "gaa_scaled4": 0.875e-9, "gaa_ibm2017_after_hk": 0.5e-9}   # 8 voxels per slot
CPP = 44e-9                                                   # contacted poly pitch (Loubet 2017: 44/48 nm)

NAND_CELLS = {
    # wp4 / wp9 cell with the channel-hole pillars added: 100 nm slit, 5 um stack, 29/51 nm word line/pitch
    # (-> 30/50 at 5 nm voxels), recesses 365 nm, holes 100 nm at 200 nm pitch, two staggered rows per side
    "nand_wp4_pillars": dict(w_s=100e-9, L_s=5e-6, g_r=29e-9, L_r=365e-9, p_r=51e-9, d_h=100e-9, pitch_h=200e-9, rows=2),
    # 2XX-layer cell: slit depth 12.7 um and word-line pitch 48 nm (TechInsights, YMTC 232L), slit 250 nm and
    # ~100 nm holes (Onto Innovation 2023), slit pitch > 2 um -> recesses 1 um, four staggered rows per side;
    # word line 30 / oxide 18 nm at 6 nm voxels (the split is an assumption: wp4's 29/22 ratio at the new pitch)
    "nand_2xx": dict(w_s=250e-9, L_s=12.7e-6, g_r=30e-9, L_r=1000e-9, p_r=48e-9, d_h=100e-9, pitch_h=200e-9, rows=4),
}
NAND_VOXEL = {"nand_wp4_pillars": 5e-9, "nand_2xx": 6e-9}


FIGURE_KEY = "gaa_ibm2017_TMA_constant_h1.25_M50000"
MAPS = FIG / "device_geometries_maps.json"      # face coverage snapshots for figures/gaa_coverage.png (part `maps`)
MAP_RUN = dict(name="gaa_ibm2017_after_hk", chem="TiCl4", drive="constant", M=50000)
NAND_NORMAL_RADIUS = 2      # cone kernel of 2 voxels: the 6-voxel recesses are not blended across (inventory 0.991 vs 0.953 at 5)


def build_gaa(name, voxel_size=None):
    return gaa_stack(h=voxel_size or GAA_VOXEL[name], **GAA_CELLS[name])


def build_nand(name, voxel_size=None):
    return nand_cell(h=voxel_size or NAND_VOXEL[name], **NAND_CELLS[name])


def gaa_wall_ratios(g, mc_site_area):
    """(cell ratio = site area / mouth area, wafer ratio with the flat ILD over the rest of the CPP)."""
    m = g.meta
    mouth = m["pitch_x"] * m["Lg_actual"]
    footprint = m["pitch_x"] * CPP
    wafer = (mc_site_area + m["pitch_x"] * (CPP - m["Lg_actual"])) / footprint
    return mc_site_area / mouth, wafer


def nand_wall_ratios(g, mc_site_area):
    """(cell ratio = site area / slit mouth area, wafer ratio for a wafer tiled at pitch w_s + 2 L_r)."""
    m = g.meta
    mouth = m["pitch_h_actual"] * m["w_actual"]
    footprint = m["pitch_h_actual"] * (m["w_actual"] + 2 * m["L_r_actual"])
    return mc_site_area / mouth, mc_site_area / footprint


# -------------------------------------------------------------------------------------------------------------------
def part_gaa(cache):
    """GAA cells: three chemistries, both drives, plus M doubling and a resolution change on the IBM cell."""
    runs = []
    for name in ("gaa_ibm2017", "gaa_scaled4"):
        for chem in ("TMA", "TDMAH", "TiCl4"):
            for drive in ("constant", "chamber"):
                runs.append((f"{name}_{chem}_{drive}_h{GAA_VOXEL[name]*1e9:.3g}_M50000", dict(name=name, chem=chem, drive=drive, M=50000)))
    for chem in ("TiCl4", "TDMAH"):
        for drive in ("constant", "chamber"):
            runs.append((f"gaa_ibm2017_after_hk_{chem}_{drive}_h0.5_M50000", dict(name="gaa_ibm2017_after_hk", chem=chem, drive=drive, M=50000)))
    runs.append(("gaa_ibm2017_TMA_constant_h1.25_M100000", dict(name="gaa_ibm2017", chem="TMA", drive="constant", M=100000)))
    runs.append(("gaa_ibm2017_TMA_constant_h0.625_M50000", dict(name="gaa_ibm2017", chem="TMA", drive="constant", M=50000, voxel_size=0.625e-9)))
    runs.append(("gaa_ibm2017_TMA_constant_h1.25_M50000_dt0.0025", dict(name="gaa_ibm2017", chem="TMA", drive="constant", M=50000, dt=0.0025)))
    for key, kw in runs:
        if key in cache:
            continue
        g = build_gaa(kw["name"], kw.get("voxel_size"))
        # flat Langmuir time 2.3 N0/(s n vbar/4) sets the scale; dt = that / 100, t_end generous
        c = CHEM[kw["chem"]]
        t_flat = 2.303 * c["N0"] / (c["s"] * N0_GAS * vbar_of(kw["chem"]) / 4)
        dt = kw.get("dt", max(t_flat / 100, 2.5e-4 if kw["drive"] == "constant" else 2.5e-3))
        t_end = 60 * t_flat if kw["drive"] == "constant" else 22.0
        # every face of the GAA cell is axis-aligned, so plain voxel faces (w = 1) are the exact area: no staircase
        rec = run_cell(g, kw["chem"], kw["drive"], dt, kw["M"], t_end, gaa_regions, label=key, area_correction=False,
                       snapshot_every=max(1, int(round(t_flat / 10 / dt))))
        cell, wafer = gaa_wall_ratios(g, rec["site_area"])
        rec["wall_ratio_cell"], rec["wall_ratio_wafer"] = cell, wafer
        rec["smooth_area"] = gaa_smooth_area(g.meta)
        rec["inventory_ratio"] = rec["site_area"] / rec["smooth_area"]
        rec["rule_pulse_s"] = rule_pulse(kw["chem"], wafer)
        rec["effusion_limit_s"] = effusion_limit(kw["chem"], cell)
        rec["t_flat_90_s"] = t_flat
        tt = np.array(rec["t_hist"])
        for lv in (0.9, 0.99):
            rec[f"exposure_{int(lv*100)}_Pa_s"] = P_PA * (rec["t90_everywhere"] if lv == 0.9 else rec["t99_everywhere"]) if kw["drive"] == "constant" else float("nan")
        # coverage at the rule's pulse: slowest region and top of the top sheet (step coverage = ratio)
        tr = rec["rule_pulse_s"]
        mn = np.interp(tr, tt, rec["hist_everywhere_min"], right=rec["hist_everywhere_min"][-1])
        top = np.interp(tr, tt, rec["hist"]["sheet0_top"], right=rec["hist"]["sheet0_top"][-1])
        rec["cov_min_at_rule"], rec["cov_top_at_rule"], rec["step_coverage_at_rule"] = float(mn), float(top), float(mn / top) if top > 0 else float("nan")
        out = strip(rec)
        cache[key] = out
        save_cache(cache)
        print(f"    wall-area ratio cell {cell:.2f}, wafer {wafer:.2f}; rule pulse {rec['rule_pulse_s']:.3f} s; step coverage at the rule "
              f"{rec['step_coverage_at_rule']:.3f}; inventory {rec['inventory_ratio']:.3f}", flush=True)


def part_maps(cache):
    """Face-coverage maps of the narrowed IBM cell under TiCl4 at 0.3 Pa (the run with the largest edge-to-centre
    gradient) at a third, two thirds and the whole of its t90, kept in a separate small file for the figure."""
    if MAPS.exists():
        return
    g = build_gaa(MAP_RUN["name"])
    c = CHEM[MAP_RUN["chem"]]
    t_flat = 2.303 * c["N0"] / (c["s"] * N0_GAS * vbar_of(MAP_RUN["chem"]) / 4)
    dt = max(t_flat / 100, 2.5e-4)
    rec = run_cell(g, MAP_RUN["chem"], MAP_RUN["drive"], dt, MAP_RUN["M"], 60 * t_flat, gaa_regions, label="maps", area_correction=False,
                   snapshot_every=1, verbose_every=0)
    t90 = rec["t90_everywhere"]
    times = sorted(rec["snapshots"])
    picks = [min(times, key=lambda q: abs(q - f * t90)) for f in (1 / 3, 2 / 3, 1.0)]
    f = rec["faces"]
    out = dict(run=MAP_RUN, voxel=g.voxel, shape=list(g.shape), t90=t90, meta={k: v for k, v in g.meta.items() if k in ("gap_k", "sheet_x", "trench_y")},
               x=[float(v) for v in f["x"]], y=[float(v) for v in f["y"]], z=[float(v) for v in f["z"]], direction=[int(v) for v in f["direction"]],
               snapshots={f"{q:.5f}": [round(float(v), 4) for v in rec["snapshots"][q]] for q in picks})
    MAPS.write_text(json.dumps(out))
    print(f"  maps at t = {', '.join(f'{q:.4f}' for q in picks)} s (t90 {t90:.4f} s) -> {MAPS}", flush=True)


def nand_runs(name, voxel_size, chems_drives, M, dt, t_end, extra=()):
    runs = []
    for chem, drive in chems_drives:
        te = t_end[chem] if isinstance(t_end, dict) else t_end
        if drive == "chamber":
            te = 22.0
        runs.append((f"{name}_{chem}_{drive}_h{voxel_size*1e9:.3g}_M{M}", dict(name=name, chem=chem, drive=drive, M=M, dt=dt, t_end=te, voxel_size=voxel_size)))
    runs += list(extra)
    return runs


def part_nand(cache, name, voxel_size, chems_drives, M, dt, t_end, extra=()):
    for key, kw in nand_runs(name, voxel_size, chems_drives, M, dt, t_end, extra):
        if key in cache:
            continue
        g = build_nand(kw["name"], kw["voxel_size"])
        rec = run_cell(g, kw["chem"], kw["drive"], kw["dt"], kw["M"], kw["t_end"], nand_regions, label=key, verbose_every=20,
                       normal_radius=NAND_NORMAL_RADIUS)
        cell, wafer = nand_wall_ratios(g, rec["site_area"])
        rec["wall_ratio_cell"], rec["wall_ratio_wafer"] = cell, wafer
        rec["smooth_area"] = nand_smooth_area(g.meta)
        rec["inventory_ratio"] = rec["site_area"] / rec["smooth_area"]
        rec["rule_pulse_s"] = rule_pulse(kw["chem"], wafer)
        rec["effusion_limit_s"] = effusion_limit(kw["chem"], cell)
        tt = np.array(rec["t_hist"])
        for lv, key90 in ((0.9, "t90_everywhere"), (0.99, "t99_everywhere")):
            rec[f"exposure_{int(lv*100)}_Pa_s"] = P_PA * rec[key90] if kw["drive"] == "constant" else float("nan")
        tr = rec["rule_pulse_s"]
        mn = np.interp(tr, tt, rec["hist_everywhere_min"], right=rec["hist_everywhere_min"][-1])
        top = np.interp(tr, tt, rec["hist"]["recess_top_deephalf"], right=rec["hist"]["recess_top_deephalf"][-1])
        rec["cov_min_at_rule"], rec["cov_top_at_rule"] = float(mn), float(top)
        rec["step_coverage_at_rule"] = float(mn / top) if top > 0 else float("nan")
        # slice profile at the end (coverage along the depth, for the figure)
        z, prof, wt = rec["solver"].slice_coverage()
        rec["final_profile"] = dict(z=[float(v) for v in z[wt > 0][::max(1, len(z) // 400)]], theta=[float(v) for v in prof[wt > 0][::max(1, len(z) // 400)]])
        cache[key] = strip(rec)
        save_cache(cache)
        print(f"    wall-area ratio cell {cell:.1f}, wafer {wafer:.1f}; rule {rec['rule_pulse_s']:.1f} s; effusion limit {rec['effusion_limit_s']:.2f} s; "
              f"inventory {rec['inventory_ratio']:.3f}", flush=True)


def part_nand5(cache):
    cd = [(c, d) for c in ("TMA", "TDMAH", "TiCl4") for d in ("constant", "chamber")]
    extra = [("nand_wp4_pillars_TMA_constant_h5_M200000", dict(name="nand_wp4_pillars", chem="TMA", drive="constant", M=200000, dt=0.1, t_end=30.0, voxel_size=5e-9)),
             ("nand_wp4_pillars_TMA_constant_h10_M100000", dict(name="nand_wp4_pillars", chem="TMA", drive="constant", M=100000, dt=0.1, t_end=30.0, voxel_size=10e-9))]
    part_nand(cache, "nand_wp4_pillars", 5e-9, cd, M=100000, dt=0.1, t_end=30.0, extra=extra)


def part_nand13(cache):
    cd = [("TMA", "constant"), ("TMA", "chamber"), ("TDMAH", "constant"), ("TiCl4", "constant"), ("TDMAH", "chamber"), ("TiCl4", "chamber")]
    extra = [("nand_2xx_TMA_constant_h6_M120000", dict(name="nand_2xx", chem="TMA", drive="constant", M=120000, dt=0.25, t_end=80.0, voxel_size=6e-9)),
             ("nand_2xx_TMA_constant_h8_M60000", dict(name="nand_2xx", chem="TMA", drive="constant", M=60000, dt=0.25, t_end=80.0, voxel_size=8e-9))]
    part_nand(cache, "nand_2xx", 6e-9, cd, M=60000, dt=0.25, t_end=dict(TMA=80.0, TDMAH=120.0, TiCl4=60.0), extra=extra)


def part_refs(cache):
    """Nested 1D recess model (aldconf.recess) for both NAND cells, hydraulic-diameter ('none') and parallel-plate
    ('plates') Knudsen diffusivities, constant 0.3 Pa and the chamber drive; TMA only (the comparison is about the
    transport, which is chemistry-independent up to N0/vbar scaling)."""
    refs = cache.setdefault("refs", {})
    td, nd = load_drive()
    for name, spec in NAND_CELLS.items():
        g = build_nand(name)
        m = g.meta
        for label, geom in (("nominal", RecessedSlit(spec["w_s"], spec["L_s"], spec["g_r"], spec["L_r"], spec["p_r"])),
                            ("rounded", RecessedSlit(m["w_actual"], m["L_actual"], m["g_r_actual"], m["L_r_actual"], m["p_r_actual"]))):
            for corr in ("none", "plates"):
                for chem in ("TMA", "TDMAH", "TiCl4"):
                    key = f"nested_{name}_{label}_{corr}_{chem}"
                    if key in refs:
                        continue
                    c = CHEM[chem]
                    pr = Precursor(sticking=c["s"], site_density=c["N0"], mass_amu=c["mass"])
                    t_dose = 60.0 if name == "nand_wp4_pillars" else 300.0
                    rB = simulate_recessed_constant(geom, pr, T, P_PA, t_dose, n_substeps=1200, n_samples=300, n_slit=200, n_recess=50, knudsen_correction=corr)
                    rA = simulate_recessed(geom, pr, T, td, nd, substeps_per_sample=5, n_slit=200, n_recess=50, knudsen_correction=corr)
                    refs[key] = dict(area_factor=geom.area_factor,
                                     B_slit_bottom_t90=time_to_coverage(rB.t, rB.hist_slit_bottom), B_recess_deep_t90=time_to_coverage(rB.t, rB.hist_recess_bottom_deep),
                                     B_recess_min_t90=time_to_coverage(rB.t, rB.hist_recess_bottom_min), B_recess_top_t90=time_to_coverage(rB.t, rB.hist_recess_bottom_top),
                                     B_recess_deep_t99=time_to_coverage(rB.t, rB.hist_recess_bottom_deep, 0.99),
                                     A_slit_bottom_t90=time_to_coverage(rA.t, rA.hist_slit_bottom), A_recess_deep_t90=time_to_coverage(rA.t, rA.hist_recess_bottom_deep),
                                     A_recess_min_t90=time_to_coverage(rA.t, rA.hist_recess_bottom_min), A_recess_top_t90=time_to_coverage(rA.t, rA.hist_recess_bottom_top),
                                     A_final_recess_deep=float(rA.hist_recess_bottom_deep[-1]),
                                     B_t=[float(v) for v in rB.t[::3]], B_recess_deep_hist=[float(v) for v in rB.hist_recess_bottom_deep[::3]],
                                     B_slit_bottom_hist=[float(v) for v in rB.hist_slit_bottom[::3]])
                    save_cache(cache)
                    print(f"{key}: const slit bottom {refs[key]['B_slit_bottom_t90']:.2f} s, deepest recess {refs[key]['B_recess_deep_t90']:.2f} s; "
                          f"chamber {refs[key]['A_slit_bottom_t90']:.2f} / {refs[key]['A_recess_deep_t90']:.2f} s", flush=True)


# -------------------------------------------------------------------------------------------------------------------
def write_csv(cache):
    rows = []
    for key, r in cache.items():
        if key == "refs":
            continue
        row = dict(run=key, cell=r["label"].split("_" + r["chem"])[0], chem=r["chem"], film=CHEM[r["chem"]]["film"], drive=r["drive"],
                   voxel_nm=r["voxel_nm"], dt_s=r["dt"], M=r["M"], seed=r["seed"], shape="x".join(str(v) for v in r["shape"]), n_faces=r["n_faces"],
                   mouth_area_m2=r["mouth_area"], site_area_m2=r["site_area"], smooth_area_m2=r.get("smooth_area"), inventory_ratio=r.get("inventory_ratio"),
                   wall_ratio_cell=r.get("wall_ratio_cell"), wall_ratio_wafer=r.get("wall_ratio_wafer"), rule_pulse_s=r.get("rule_pulse_s"),
                   effusion_limit_s=r.get("effusion_limit_s"), t90_everywhere_s=r["t90_everywhere"], t99_everywhere_s=r["t99_everywhere"],
                   exposure_90_Pa_s=r.get("exposure_90_Pa_s"), exposure_99_Pa_s=r.get("exposure_99_Pa_s"), last_region_to_90=r["last_region_90"],
                   cov_min_at_rule=r.get("cov_min_at_rule"), cov_top_at_rule=r.get("cov_top_at_rule"), step_coverage_at_rule=r.get("step_coverage_at_rule"),
                   final_min_region=r["final_min_region"], final_area_fraction_below_0p9=r["final_area_below_0p9"], t_end_reached_s=r["t_end_reached"],
                   particles=r["particles"], threads=r["threads"], wall_s=r["wall_s"], pps_thread=r["pps_thread"], seg_per_particle=r["seg_per_particle"],
                   escaped=r["counts"]["escaped"], stuck=r["counts"]["stuck"], left=r["counts"]["left"])
        sm = r.get("slowest_member", {}).get("recess_all_deephalf")
        if sm:
            row.update(slowest_recess_index=sm["index"], n_recesses=sm["of"], slowest_recess_depth_fraction=round(sm["fraction"], 3),
                       t90_slowest_recess=sm["t90"], t99_slowest_recess=sm["t99"], t90_deepest_recess=sm["last_t90"], t99_deepest_recess=sm["last_t99"])
        for k, v in r["t90"].items():
            row[f"t90_{k}"] = v
        for k, v in r["t99"].items():
            row[f"t99_{k}"] = v
        # thickness after N cycles at the rule's pulse: N x GPC x coverage (50 cycles)
        gpc = CHEM[r["chem"]]["gpc_nm"]
        row["film_50cyc_at_rule_min_nm"] = 50 * gpc * r.get("cov_min_at_rule", float("nan"))
        row["film_50cyc_at_rule_top_nm"] = 50 * gpc * r.get("cov_top_at_rule", float("nan"))
        rows.append(row)
    cols = []
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    rows.sort(key=lambda r: (r["cell"], r["chem"], r["drive"], r["run"]))
    with open(FIG / "device_geometries.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else (f"{r[c]:.6g}" if isinstance(r[c], float) else r[c])) for c in cols})
    # nested references as extra rows in a second block
    refs = cache.get("refs", {})
    if refs:
        with open(FIG / "device_geometries.csv", "a", newline="") as fh:
            fh.write("\n# nested 1D recess model references (aldconf.recess): key, area_factor, B = constant 0.3 Pa, A = chamber drive\n")
            w = csv.writer(fh)
            w.writerow(["key", "area_factor", "B_slit_bottom_t90", "B_recess_top_t90", "B_recess_deep_t90", "B_recess_min_t90", "B_recess_deep_t99",
                        "A_slit_bottom_t90", "A_recess_top_t90", "A_recess_deep_t90", "A_recess_min_t90", "A_final_recess_deep"])
            for k, v in refs.items():
                w.writerow([k] + [f"{v[c]:.4g}" for c in ("area_factor", "B_slit_bottom_t90", "B_recess_top_t90", "B_recess_deep_t90", "B_recess_min_t90", "B_recess_deep_t99",
                                                        "A_slit_bottom_t90", "A_recess_top_t90", "A_recess_deep_t90", "A_recess_min_t90", "A_final_recess_deep")])
    print(f"wrote {FIG / 'device_geometries.csv'} ({len(rows)} runs, {len(refs)} references)")


def plot_geometries():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    # GAA cells: cross-section across the sheet width (x-z at mid-Lg) and along the gate (y-z at mid-W)
    for ax, name in zip(axes[0], ("gaa_ibm2017", "gaa_scaled4", "gaa_ibm2017_after_hk")):
        g = build_gaa(name)
        m = g.meta
        j = (m["trench_y"][0] + m["trench_y"][1]) // 2
        sl = g.occupancy[:, j, :].T
        h = g.voxel * 1e9
        ax.imshow(sl, origin="upper", cmap="Greys", extent=[0, sl.shape[1] * h, sl.shape[0] * h, 0], interpolation="nearest", aspect="equal")
        ax.set_xlabel("x across the sheet width, nm (periodic)"); ax.set_ylabel("depth below the trench mouth z, nm")
        ax.set_title(f"{name}: {m['n_sheets']} sheets {m['t_s_actual']*1e9:.3g} nm, slots {m['g_s_actual']*1e9:.3g} nm, W {m['W_actual']*1e9:.3g} nm, "
                     f"Lg {m['Lg_actual']*1e9:.3g} nm\nvoxel {h:.3g} nm, {g.shape[0]}x{g.shape[1]}x{g.shape[2]} (slice at mid-Lg)", fontsize=8)
    # NAND cells: plan view through a recess (x-y) and the y-z section at the pillar row
    for ax, name in zip(axes[1][:2], ("nand_wp4_pillars", "nand_2xx")):
        g = build_nand(name)
        m = g.meta
        k = (m["recess_k"][0][0] + m["recess_k"][0][1]) // 2
        sl = g.occupancy[:, :, k]
        h = g.voxel * 1e9
        ax.imshow(sl.T, origin="lower", cmap="Greys", extent=[0, sl.shape[0] * h, 0, sl.shape[1] * h], interpolation="nearest", aspect="auto")
        ax.set_xlabel("x along the slit, nm (periodic; stretched)"); ax.set_ylabel("y across slit and recesses, nm")
        ax.set_title(f"{name}: plan view inside a recess\nslit {m['w_actual']*1e9:.0f} nm, recesses {m['L_r_actual']*1e9:.0f} nm, holes {m['d_h']*1e9:.0f} nm\n"
                     f"{m['n_recesses']} recesses of {m['g_r_actual']*1e9:.0f} nm at {m['p_r_actual']*1e9:.0f} nm pitch, {m['L_actual']*1e6:.1f} um deep", fontsize=8)
    ax = axes[1][2]
    g = build_nand("nand_2xx")
    m = g.meta
    h = g.voxel * 1e9
    i = int(round(m["pillar_xy"][0][0] / g.voxel)) % g.shape[0]
    kmax = m["recess_k"][5][1] + 3
    sl = g.occupancy[i, :, :kmax].T
    ax.imshow(sl, origin="upper", cmap="Greys", extent=[0, sl.shape[1] * h, sl.shape[0] * h, 0], interpolation="nearest", aspect="auto")
    ax.set_xlabel("y across slit and recesses, nm"); ax.set_ylabel("depth z, nm (first six word lines)")
    ax.set_title("nand_2xx: section through a pillar row at the top of the stack\n(white = gas: slit, recesses interrupted by the pillars)", fontsize=8)
    fig.suptitle("Device geometries as voxel occupancy arrays (black solid, white gas; mouth at z = 0)", fontsize=11)
    fig.savefig(FIG / "device_geometries.png", dpi=150)
    plt.close(fig)


def plot_gaa(cache):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    key = FIGURE_KEY
    if key not in cache:
        return
    r = cache[key]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 2.2]))
    # maps: coverage of the underside of the lowest sheet (direction -z, facing the bottom slot) at t90/3, 2 t90/3, t90
    if MAPS.exists():
        mp = json.loads(MAPS.read_text())
        h = mp["voxel"]; nx, ny = mp["shape"][0], mp["shape"][1]
        x = np.array(mp["x"]); y = np.array(mp["y"]); z = np.array(mp["z"]); d = np.array(mp["direction"])
        ka, kb = mp["meta"]["gap_k"][-1]
        sel = (d == 5) & (z > ka * h - 0.5 * h) & (z < kb * h + 0.5 * h)
        xi = np.round(x[sel] / h - 0.5).astype(int); yi = np.round(y[sel] / h - 0.5).astype(int)
        im = None
        for ax, q in zip(axes[0], sorted(mp["snapshots"], key=float)):
            th = np.array(mp["snapshots"][q])
            img = np.full((ny, nx), np.nan)
            img[yi, xi] = th[sel]
            im = ax.imshow(img, origin="lower", vmin=0, vmax=1, cmap="viridis", extent=[0, nx * h * 1e9, 0, ny * h * 1e9], aspect="equal")
            ax.set_title(f"narrowed IBM cell (4 nm slots), TiCl4 at 0.3 Pa: underside of the lowest sheet\nt = {float(q)*1e3:.0f} ms = {float(q)/mp['t90']:.2f} t90", fontsize=8)
            ax.set_xlabel("x along the sheet width, nm (period 60 nm)"); ax.set_ylabel("y along the gate, nm")
        if im is not None:
            fig.colorbar(im, ax=axes[0].tolist(), shrink=0.7, label="coverage", pad=0.01)
    else:
        for ax in axes[0]:
            ax.text(0.5, 0.5, "run --part maps for the coverage maps", ha="center", va="center", transform=ax.transAxes); ax.set_axis_off()
    # histories: regions for TMA constant; everywhere-min for three chemistries; chamber drive
    ax = axes[1][0]
    for k, v in r["hist"].items():
        if k.endswith("_mid") or k in ("sheet0_top", "floor_under_stack", "floor_between_stacks"):
            ax.plot(r["t_hist"], v, label=f"{k} (t90 {r['t90'][k]:.3f} s)")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlabel("dose time, s"); ax.set_ylabel("coverage"); ax.set_title("IBM-2017 cell, TMA at 0.3 Pa: regions", fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3); ax.set_xlim(0, max(r["t_hist"]))
    ax = axes[1][1]
    for chem, col in (("TMA", "tab:blue"), ("TDMAH", "tab:orange"), ("TiCl4", "tab:green")):
        for name, ls in (("gaa_ibm2017", "-"), ("gaa_scaled4", "--")):
            k = f"{name}_{chem}_constant_h{GAA_VOXEL[name]*1e9:.3g}_M50000"
            if k in cache:
                rr = cache[k]
                ax.plot(rr["t_hist"], rr["hist_everywhere_min"], color=col, ls=ls, label=f"{name} {chem}: t90 {rr['t90_everywhere']:.3f} s")
    for chem, col in (("TiCl4", "tab:red"), ("TDMAH", "tab:brown")):
        k = f"gaa_ibm2017_after_hk_{chem}_constant_h0.5_M50000"
        if k in cache:
            ax.plot(cache[k]["t_hist"], cache[k]["hist_everywhere_min"], color=col, label=f"after IL+HfO2 (4 nm slots), {chem}: t90 {cache[k]['t90_everywhere']:.3f} s")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlabel("dose time at 0.3 Pa, s"); ax.set_ylabel("slowest region coverage"); ax.set_xscale("log")
    ax.set_title("Slowest region, constant 0.3 Pa, three chemistries", fontsize=9); ax.legend(fontsize=6); ax.grid(alpha=0.3, which="both")
    ax = axes[1][2]
    td, nd = load_drive()
    ax2 = ax.twinx(); ax2.plot(td, nd / C0, color="0.7", lw=0.8); ax2.set_ylabel("centre-face density n / c0 (grey)"); ax2.set_ylim(0, 0.1)
    for chem, col in (("TMA", "tab:blue"), ("TDMAH", "tab:orange"), ("TiCl4", "tab:green")):
        k = f"gaa_ibm2017_{chem}_chamber_h1.25_M50000"
        if k in cache:
            rr = cache[k]
            ax.plot(rr["t_hist"], rr["hist_everywhere_min"], color=col, label=f"{chem}: t90 {rr['t90_everywhere']:.2f} s, t99 {rr['t99_everywhere']:.2f} s")
    k = "gaa_ibm2017_after_hk_TiCl4_chamber_h0.5_M50000"
    if k in cache:
        ax.plot(cache[k]["t_hist"], cache[k]["hist_everywhere_min"], color="tab:red", label=f"after IL+HfO2, TiCl4: t90 {cache[k]['t90_everywhere']:.2f} s")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlabel("time, s"); ax.set_ylabel("slowest region coverage"); ax.set_xlim(0, 6)
    ax.set_title("IBM-2017 cell under the chamber's centre-face history (AR-100 wafer, 1 % TMA)", fontsize=9); ax.legend(fontsize=7, loc="lower right"); ax.grid(alpha=0.3)
    fig.suptitle("GAA nanosheet stack: coverage maps and histories (lib/aldmc voxel Monte Carlo)", fontsize=11)
    fig.savefig(FIG / "gaa_coverage.png", dpi=150)
    plt.close(fig)


def plot_nand(cache):
    """figures/nand_cell_coverage.png: histories of the two NAND cells (slit floor, deepest recess, slowest recess) for
    the three chemistries at 0.3 Pa and under the chamber drive, the nested-model references, and the per-recess
    saturation profile along the depth."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    refs = cache.get("refs", {})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
    for row, (name, hkey, Mkey) in enumerate((("nand_wp4_pillars", "h5", "M100000"), ("nand_2xx", "h6", "M60000"))):
        for col, drive in enumerate(("constant", "chamber")):
            ax = axes[row][col]
            for chem, c in (("TMA", "tab:blue"), ("TDMAH", "tab:orange"), ("TiCl4", "tab:green")):
                k = f"{name}_{chem}_{drive}_{hkey}_{Mkey}"
                if k not in cache:
                    continue
                r = cache[k]
                ax.plot(r["t_hist"], r["hist"]["slit_floor"], color=c, ls=":", lw=1)
                ax.plot(r["t_hist"], r["hist"]["recess_deep_deephalf"], color=c, ls="--", lw=1)
                ax.plot(r["t_hist"], r["hist_min_list"]["recess_all_deephalf"], color=c, label=f"{chem}: slowest recess t90 {r['t90_everywhere']:.1f} s (deepest {r['t90']['recess_deep_deephalf']:.1f}, floor {r['t90']['slit_floor']:.1f})")
            if drive == "constant":
                for corr, ls, lab in (("none", "-.", "nested model (wp4 slit diffusivity)"), ("plates", (0, (1, 1)), "nested model (parallel-plate correction)")):
                    rk = f"nested_{name}_rounded_{corr}_TMA"
                    if rk in refs:
                        ax.plot(refs[rk]["B_t"], refs[rk]["B_recess_deep_hist"], color="k", ls=ls, lw=1, label=f"{lab}: deepest recess {refs[rk]['B_recess_deep_t90']:.1f} s")
            else:
                td, nd = load_drive()
                ax2 = ax.twinx(); ax2.plot(td, nd / C0, color="0.75", lw=0.8); ax2.set_ylim(0, 0.1); ax2.set_ylabel("centre-face n / c0 (grey)")
            ax.axhline(0.9, color="0.6", ls=":")
            if drive == "constant":
                xmax = [cache[f"{name}_{chem}_constant_{hkey}_{Mkey}"]["t_end_reached"] for chem in ("TMA", "TDMAH", "TiCl4") if f"{name}_{chem}_constant_{hkey}_{Mkey}" in cache]
                xmax += [1.4 * refs[f"nested_{name}_rounded_none_TMA"]["B_recess_deep_t90"]] if f"nested_{name}_rounded_none_TMA" in refs else []
                if xmax:
                    ax.set_xlim(0, max(xmax))
            ax.set_xlabel("time, s" if drive == "chamber" else "dose time at 0.3 Pa, s"); ax.set_ylabel("coverage")
            ax.set_title(f"{name}, {'constant 0.3 Pa' if drive == 'constant' else 'chamber centre-face history (AR-100 wafer)'}\n"
                         "dotted: slit floor, dashed: deepest recess, solid: slowest recess (deep halves)", fontsize=8)
            ax.legend(fontsize=6, loc="lower right"); ax.grid(alpha=0.3)
    # per-recess profiles
    ax = axes[0][2]
    prof = FIG / "device_geometries_recess_profile.csv"
    if prof.exists():
        rows = [l for l in prof.read_text().splitlines() if l and not l.startswith("#")]
        data = np.array([[float(v) for v in l.split(",")] for l in rows[1:]])
        ax.plot(data[:, 1], data[:, 2], "o-", ms=3, color="tab:blue", label="nand_wp4_pillars, TMA 0.3 Pa: t90 per recess (M 50 000)")
        ax.plot(data[:, 1], data[:, 3], "s-", ms=3, color="tab:cyan", label="t99 per recess")
    for name, hkey, Mkey, c in (("nand_2xx", "h6", "M60000", "tab:red"),):
        k = f"{name}_TMA_constant_{hkey}_{Mkey}"
        if k in cache and "t90_list" in cache[k]:
            r = cache[k]; g = build_nand(name); h = g.voxel
            z = [0.5 * (a + b) * h * 1e6 for a, b in g.meta["recess_k"]]
            ax.plot(z, r["t90_list"]["recess_all_deephalf"], "o-", ms=2, color=c, label=f"{name}, TMA 0.3 Pa: t90 per recess")
    ax.set_xlabel("recess depth below the slit mouth, um"); ax.set_ylabel("time to 90 / 99 % of the recess deep half, s")
    ax.set_title("Where the last point sits: saturation time along the stack", fontsize=9); ax.legend(fontsize=6); ax.grid(alpha=0.3)
    # final profile along the depth for the 2xx cell runs (slice coverage)
    ax = axes[1][2]
    for chem, c in (("TMA", "tab:blue"), ("TDMAH", "tab:orange"), ("TiCl4", "tab:green")):
        for drive, ls in (("constant", "-"), ("chamber", "--")):
            k = f"nand_2xx_{chem}_{drive}_h6_M60000"
            if k in cache and "final_profile" in cache[k]:
                fp = cache[k]["final_profile"]
                ax.plot(np.array(fp["z"]) * 1e6, fp["theta"], color=c, ls=ls, lw=1, label=f"{chem} {drive}: at t = {cache[k]['t_end_reached']:.1f} s")
    ax.set_xlabel("depth, um"); ax.set_ylabel("slice-mean coverage at the end of the run"); ax.set_title("nand_2xx: coverage along the depth at the end of each run", fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3)
    fig.suptitle("3D-NAND word-line cells with channel-hole pillars (lib/aldmc voxel Monte Carlo)", fontsize=11)
    fig.savefig(FIG / "nand_cell_coverage.png", dpi=150)
    plt.close(fig)


PARTS = dict(gaa=part_gaa, maps=part_maps, nand5=part_nand5, nand13=part_nand13, refs=part_refs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="", help="comma-separated: " + ",".join(PARTS))
    ap.add_argument("--figure", action="store_true")
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
    write_csv(cache)
    plot_geometries()
    plot_gaa(cache)
    plot_nand(cache)
    print("wrote figures/device_geometries.png, figures/gaa_coverage.png, figures/nand_cell_coverage.png")


if __name__ == "__main__":
    main()
