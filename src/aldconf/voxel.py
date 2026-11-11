"""ctypes binding of libaldmc (lib/aldmc): the 3D voxel test-particle Monte Carlo feature solver, with geometry
builders and a driver that mirrors aldconf.simulate_driven.

    from aldconf import voxel
    g = voxel.hole(d=100e-9, L=10e-6, voxel=100e-9 / 16)                 # occupancy array + voxel size
    mc = voxel.VoxelMC(g, voxel.Params(s=5e-3, N0=4e18, vbar=410.0, seed=1, threads=8))
    mc.step(n_mouth=3.79e19, dt=0.1, M=200_000)                           # one time step
    mc.bottom_coverage(); mc.slice_coverage(); mc.counts()
    res = voxel.simulate_driven_voxel(g, params, t_series, n0_series, M=200_000)   # like aldconf.simulate_driven

Physics, discretisation and the particle weighting are stated in lib/aldmc/include/aldmc.h and lib/aldmc/README.md:
molecular flow (straight flights, no gas-phase collisions), cosine entry over the open mouth, sticking s (1 - theta)
on the face's own coverage, loss r, diffuse (cosine) re-emission, escape through the mouth; faces carry N0 h^2 w sites
with the projected-area weight w of the estimated smooth normal; a step of M particles represents n vbar/4 A_mouth dt
molecules; the coverage update is the Langmuir relaxation of the Monte Carlo impingement. No film narrowing.

The library is looked up in this order: the file named by ALDMC_LIBRARY, the directory named by ALDMC_BUILD_DIR,
<repo>/lib/aldmc/build (also its Release/ subdirectory) and finally the system loader. `is_available()` says whether
it was found; `VoxelMC` raises AldMCNotBuilt otherwise with the build command.

Occupancy arrays are numpy uint8 of shape (nx, ny, nz), 1 = solid, 0 = gas; z increases into the feature and the
plane z = 0 is the open mouth. Lengths in metres.
"""
from __future__ import annotations

import ctypes
import os
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUILD_DIR = ROOT / "lib" / "aldmc" / "build"
BUILD_COMMAND = "cmake -S lib/aldmc -B lib/aldmc/build -DCMAKE_BUILD_TYPE=Release && cmake --build lib/aldmc/build"

UPDATES = {"exponential": 0, "deposit": 1}
EMISSIONS = {"smooth": 0, "facet": 1}
DIRECTIONS = {"+x": 0, "-x": 1, "+y": 2, "-y": 3, "+z": 4, "-z": 5}
ERRORS = {0: "ok", 1: "ALDMC_ERR_ARG", 2: "ALDMC_ERR_PARAM", 3: "ALDMC_ERR_STATE", 4: "ALDMC_ERR_ALLOC", 5: "ALDMC_ERR_INTERNAL"}


class AldMCError(RuntimeError):
    """A library call returned an error code; the message is aldmc_last_error()."""


class AldMCNotBuilt(AldMCError):
    """libaldmc was not found."""


class _Params(ctypes.Structure):
    _fields_ = [("s", ctypes.c_double), ("r", ctypes.c_double), ("N0", ctypes.c_double), ("vbar", ctypes.c_double),
                ("seed", ctypes.c_ulonglong), ("threads", ctypes.c_int), ("periodicX", ctypes.c_int),
                ("periodicY", ctypes.c_int), ("areaCorrection", ctypes.c_int), ("update", ctypes.c_int),
                ("normalRadius", ctypes.c_int), ("emission", ctypes.c_int), ("maxSegments", ctypes.c_longlong)]


def default_threads() -> int:
    """At most 8 worker threads by default (keeps the machine usable); ALDMC_THREADS overrides."""
    env = os.environ.get("ALDMC_THREADS")
    if env:
        return max(1, int(env))
    return max(1, min(8, os.cpu_count() or 1))


@dataclass
class Params:
    """The aldmc_params struct. Periodicity comes from the geometry (VoxelGeometry.periodic), not from here."""
    s: float = 5e-3
    r: float = 0.0
    N0: float = 4e18
    vbar: float = 410.0
    seed: int = 1
    threads: int = 0                 # 0 -> default_threads()
    area_correction: bool = True
    update: str = "exponential"
    normal_radius: int = 0           # 0 -> the library default (5 voxels, cone kernel)
    emission: str = "facet"          # re-emission about the "facet" normal (default) or the estimated "smooth" normal
    max_segments: int = 0            # 0 -> the library default (1e8)

    def to_struct(self, periodic: Tuple[bool, bool]) -> _Params:
        if self.update not in UPDATES:
            raise ValueError(f"update must be one of {sorted(UPDATES)}")
        if self.emission not in EMISSIONS:
            raise ValueError(f"emission must be one of {sorted(EMISSIONS)}")
        return _Params(s=float(self.s), r=float(self.r), N0=float(self.N0), vbar=float(self.vbar), seed=int(self.seed),
                       threads=int(self.threads) if self.threads > 0 else default_threads(),
                       periodicX=int(bool(periodic[0])), periodicY=int(bool(periodic[1])),
                       areaCorrection=int(bool(self.area_correction)), update=UPDATES[self.update],
                       normalRadius=int(self.normal_radius), emission=EMISSIONS[self.emission], maxSegments=int(self.max_segments))


# ---------------------------------------------------------------------------------------------------------------
# geometry builders
@dataclass
class VoxelGeometry:
    occupancy: np.ndarray                 # uint8 (nx, ny, nz), 1 solid
    voxel: float                          # m
    periodic: Tuple[bool, bool] = (False, False)
    meta: dict = field(default_factory=dict)

    @property
    def shape(self):
        return tuple(self.occupancy.shape)

    @property
    def depth(self) -> float:
        """Depth of the deepest gas voxel's far face (m)."""
        gas_k = np.where(self.occupancy.min(axis=(0, 1)) == 0)[0]
        return float((gas_k.max() + 1) * self.voxel) if len(gas_k) else 0.0

    def z_centres(self) -> np.ndarray:
        return (np.arange(self.occupancy.shape[2]) + 0.5) * self.voxel


def from_array(occupancy, voxel: float, periodic=(False, False), **meta) -> VoxelGeometry:
    occ = np.ascontiguousarray(np.asarray(occupancy) != 0).astype(np.uint8)
    if occ.ndim != 3:
        raise ValueError("occupancy must be a 3D array (nx, ny, nz)")
    return VoxelGeometry(occ, float(voxel), (bool(periodic[0]), bool(periodic[1])), dict(meta))


def _round_voxels(length: float, voxel: float, minimum: int = 1) -> int:
    return max(minimum, int(round(length / voxel)))


def hole(d: float, L: float, voxel: float, pad: int = 1, bottom_ratio: float = 1.0) -> VoxelGeometry:
    """Cylindrical hole of diameter d and depth L (frustum with bottom_ratio < 1), closed by a floor layer, with
    `pad` solid voxels around. In every slice the round(pi r^2) voxel centres nearest the axis are gas, so the digital
    cross-section (and the mouth area that sets the entering flux) equals the smooth one to within one voxel; a plain
    radius threshold gives digital discs 1-4 % off in area depending on the resolution."""
    nd = _round_voxels(d, voxel, 2)
    nL = _round_voxels(L, voxel, 1)
    n = nd + 2 * pad
    c = n / 2.0
    x = (np.arange(n) + 0.5 - c)[:, None]
    y = (np.arange(n) + 0.5 - c)[None, :]
    order = np.argsort((x ** 2 + y ** 2).ravel(), kind="stable")
    occ = np.ones((n, n, nL + 1), dtype=np.uint8)
    for k in range(nL):
        zc = (k + 0.5) / nL                                        # relative depth of the voxel centre
        radius = 0.5 * nd * (1.0 - (1.0 - bottom_ratio) * zc)      # in voxels
        count = int(round(np.pi * radius ** 2))
        layer = np.ones(n * n, dtype=np.uint8)
        layer[order[:count]] = 0
        occ[:, :, k] = layer.reshape(n, n)
    meta = dict(kind="hole" if bottom_ratio == 1.0 else "tapered_hole", d=d, L=L, bottom_ratio=bottom_ratio,
                voxels_across=nd, voxels_deep=nL, d_actual=nd * voxel, L_actual=nL * voxel)
    return VoxelGeometry(occ, float(voxel), (False, False), meta)


def tapered_hole(d: float, L: float, bottom_ratio: float, voxel: float, pad: int = 1) -> VoxelGeometry:
    """Hole whose diameter falls linearly from d at the mouth to bottom_ratio d at the bottom (aldconf.TaperedHole)."""
    if not (0.0 < bottom_ratio <= 1.0):
        raise ValueError("bottom_ratio must be in (0, 1]")
    return hole(d, L, voxel, pad=pad, bottom_ratio=bottom_ratio)


def slit(w: float, L: float, voxel: float, nx: int = 4, pad: int = 1) -> VoxelGeometry:
    """Slit of gap w (along y) and depth L, infinite along x (periodic over nx voxels), closed by a floor."""
    nw = _round_voxels(w, voxel, 1)
    nL = _round_voxels(L, voxel, 1)
    occ = np.ones((nx, nw + 2 * pad, nL + 1), dtype=np.uint8)
    occ[:, pad:pad + nw, :nL] = 0
    meta = dict(kind="slit", w=w, L=L, voxels_across=nw, voxels_deep=nL, w_actual=nw * voxel, L_actual=nL * voxel)
    return VoxelGeometry(occ, float(voxel), (True, False), meta)


def recessed_slit(w: float, L: float, g_r: float, L_r: float, p_r: float, voxel: float, nx: int = 4, pad: int = 1,
                  oxide_first: bool = True) -> VoxelGeometry:
    """The word-line structure of aldconf.recess.RecessedSlit: a slit of gap w and depth L whose two walls carry
    lateral recesses of height g_r (along z), lateral length L_r (along y) and vertical pitch p_r, periodic along
    the slit (x). The stack starts with a solid (oxide) layer of thickness p_r - g_r at the mouth when oxide_first.
    g_r and p_r are rounded to whole voxels; meta records the actual values."""
    nw = _round_voxels(w, voxel, 1)
    nL = _round_voxels(L, voxel, 1)
    ng = _round_voxels(g_r, voxel, 1)
    npitch = _round_voxels(p_r, voxel, ng + 1)
    nLr = _round_voxels(L_r, voxel, 1)
    ny = nw + 2 * nLr + 2 * pad
    occ = np.ones((nx, ny, nL + 1), dtype=np.uint8)
    y0 = pad + nLr
    occ[:, y0:y0 + nw, :nL] = 0                                # the slit
    ox = npitch - ng
    k = 0
    recesses = []
    while k < nL:
        start = k + (ox if oxide_first else 0)
        end = min(start + ng, nL)
        if start >= nL:
            break
        if end - start == ng:                                   # whole recesses only
            occ[:, pad:y0, start:end] = 0                       # left recess
            occ[:, y0 + nw:y0 + nw + nLr, start:end] = 0        # right recess
            recesses.append((start, end))
        k += npitch
    meta = dict(kind="recessed_slit", w=w, L=L, g_r=g_r, L_r=L_r, p_r=p_r, voxels_across=nw, voxels_deep=nL,
                w_actual=nw * voxel, L_actual=nL * voxel, g_r_actual=ng * voxel, p_r_actual=npitch * voxel,
                L_r_actual=nLr * voxel, n_recesses=len(recesses), recess_k=recesses, pad=pad, y_slit=(y0, y0 + nw),
                recess_end_y=(pad, y0 + nw + nLr))
    return VoxelGeometry(occ, float(voxel), (True, False), meta)


def open_box(nx: int, ny: int, nz: int, voxel: float) -> VoxelGeometry:
    """All gas, periodic laterally: no walls, every particle leaves (tests)."""
    return VoxelGeometry(np.zeros((nx, ny, nz), dtype=np.uint8), float(voxel), (True, True), dict(kind="open_box"))


def flat_surface(n: int, voxel: float, gap: int = 2) -> VoxelGeometry:
    """A flat floor `gap` voxels below the mouth, periodic laterally: the Langmuir reference surface (tests)."""
    occ = np.zeros((n, n, gap + 1), dtype=np.uint8)
    occ[:, :, gap] = 1
    return VoxelGeometry(occ, float(voxel), (True, True), dict(kind="flat"))


# ---------------------------------------------------------------------------------------------------------------
# library loading
def _candidates():
    names = {"darwin": ["libaldmc.dylib"], "win32": ["aldmc.dll", "libaldmc.dll"]}.get(sys.platform, ["libaldmc.so"])
    env_file = os.environ.get("ALDMC_LIBRARY")
    env_dir = os.environ.get("ALDMC_BUILD_DIR")
    if env_file:
        yield Path(env_file)
    dirs = []
    if env_dir:
        dirs += [Path(env_dir), Path(env_dir) / "Release"]
    if not env_file and not env_dir:
        dirs += [BUILD_DIR, BUILD_DIR / "Release"]
    for d in dirs:
        for n in names:
            yield d / n


_lib = None
_lib_path: Optional[Path] = None


def load_library(path: Optional[os.PathLike] = None) -> ctypes.CDLL:
    """Load libaldmc once and declare the prototypes. `path` overrides the search."""
    global _lib, _lib_path
    if _lib is not None and path is None:
        return _lib
    tried = []
    found = None
    if path is not None:
        found = Path(path)
    else:
        for c in _candidates():
            tried.append(str(c))
            if c.exists():
                found = c
                break
    if found is None and not (os.environ.get("ALDMC_LIBRARY") or os.environ.get("ALDMC_BUILD_DIR")):
        try:
            lib = ctypes.CDLL({"darwin": "libaldmc.dylib", "win32": "aldmc.dll"}.get(sys.platform, "libaldmc.so"))
            found = Path(lib._name)
        except OSError:
            raise AldMCNotBuilt("libaldmc not found (looked at: " + ", ".join(tried) + "); build it with: " + BUILD_COMMAND)
    elif found is None:
        raise AldMCNotBuilt("libaldmc not found (looked at: " + ", ".join(tried) + "); build it with: " + BUILD_COMMAND)
    else:
        lib = ctypes.CDLL(str(found))
    D, I, P, LL, ULL = ctypes.c_double, ctypes.c_int, ctypes.c_void_p, ctypes.c_longlong, ctypes.c_ulonglong
    PD, PI, PLL, PULL = ctypes.POINTER(D), ctypes.POINTER(I), ctypes.POINTER(LL), ctypes.POINTER(ULL)
    lib.aldmc_version.restype = ctypes.c_char_p; lib.aldmc_version.argtypes = []
    lib.aldmc_params_default.restype = None; lib.aldmc_params_default.argtypes = [ctypes.POINTER(_Params)]
    lib.aldmc_create.restype = P
    lib.aldmc_create.argtypes = [I, I, I, D, ctypes.POINTER(ctypes.c_ubyte), ctypes.POINTER(_Params)]
    lib.aldmc_destroy.restype = None; lib.aldmc_destroy.argtypes = [P]
    lib.aldmc_last_error.restype = ctypes.c_char_p; lib.aldmc_last_error.argtypes = [P]
    lib.aldmc_get_params.restype = I; lib.aldmc_get_params.argtypes = [P, ctypes.POINTER(_Params)]
    lib.aldmc_step.restype = I; lib.aldmc_step.argtypes = [P, D, D, LL]
    lib.aldmc_reset.restype = I; lib.aldmc_reset.argtypes = [P]
    lib.aldmc_n_faces.restype = I; lib.aldmc_n_faces.argtypes = [P]
    lib.aldmc_n_mouth_voxels.restype = I; lib.aldmc_n_mouth_voxels.argtypes = [P]
    for name in ("aldmc_mouth_area", "aldmc_site_area", "aldmc_time", "aldmc_bottom_coverage", "aldmc_face_coverage"):
        getattr(lib, name).restype = I; getattr(lib, name).argtypes = [P, PD]
    lib.aldmc_set_face_coverage.restype = I; lib.aldmc_set_face_coverage.argtypes = [P, PD]
    lib.aldmc_face_info.restype = I; lib.aldmc_face_info.argtypes = [P, PLL, PI, PD, PD, PD, PD]
    lib.aldmc_face_hits.restype = I; lib.aldmc_face_hits.argtypes = [P, PULL, PULL]
    lib.aldmc_slice_coverage.restype = I; lib.aldmc_slice_coverage.argtypes = [P, PD, PD]
    lib.aldmc_counts.restype = I; lib.aldmc_counts.argtypes = [P, PULL, PULL, PULL, PULL]
    lib.aldmc_last_counts.restype = I; lib.aldmc_last_counts.argtypes = [P, PULL, PULL, PULL, PULL]
    lib.aldmc_last_segments.restype = I; lib.aldmc_last_segments.argtypes = [P, PULL]
    _lib, _lib_path = lib, found
    return lib


def is_available() -> bool:
    try:
        load_library()
        return True
    except AldMCNotBuilt:
        return False


def library_path() -> Optional[Path]:
    load_library()
    return _lib_path


def version() -> str:
    return load_library().aldmc_version().decode()


def _ptr(a, ctype):
    return a.ctypes.data_as(ctypes.POINTER(ctype))


# ---------------------------------------------------------------------------------------------------------------
class VoxelMC:
    """A solver on one voxel geometry. Methods raise AldMCError on a non-zero return code."""

    def __init__(self, geometry: VoxelGeometry, params: Optional[Params] = None):
        self.lib = load_library()
        self.geometry = geometry
        self.params = params or Params()
        occ = np.ascontiguousarray(geometry.occupancy.transpose(2, 1, 0)).astype(np.uint8)   # (nz, ny, nx) -> (k ny + j) nx + i
        self._occ = occ
        nx, ny, nz = geometry.occupancy.shape
        self._nz = nz
        self._pstruct = self.params.to_struct(geometry.periodic)
        self._h = self.lib.aldmc_create(nx, ny, nz, float(geometry.voxel), _ptr(occ, ctypes.c_ubyte), ctypes.byref(self._pstruct))
        if not self._h:
            raise AldMCError(self.lib.aldmc_last_error(None).decode())
        self.n_faces = self.lib.aldmc_n_faces(self._h)
        self.wall_time = 0.0
        self.particles = 0
        self.segments = 0

    def __del__(self):
        self.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        h = getattr(self, "_h", None)
        if h:
            self.lib.aldmc_destroy(h)
            self._h = None

    def _check(self, code: int):
        if code != 0:
            raise AldMCError(f"{ERRORS.get(code, code)}: {self.lib.aldmc_last_error(self._h).decode()}")

    @property
    def threads(self) -> int:
        p = _Params(); self._check(self.lib.aldmc_get_params(self._h, ctypes.byref(p))); return p.threads

    @property
    def mouth_area(self) -> float:
        v = ctypes.c_double(); self._check(self.lib.aldmc_mouth_area(self._h, ctypes.byref(v))); return v.value

    @property
    def site_area(self) -> float:
        """Sum of the weighted face areas (m2); sites = N0 * site_area."""
        v = ctypes.c_double(); self._check(self.lib.aldmc_site_area(self._h, ctypes.byref(v))); return v.value

    @property
    def time(self) -> float:
        v = ctypes.c_double(); self._check(self.lib.aldmc_time(self._h, ctypes.byref(v))); return v.value

    def step(self, n_mouth: float, dt: float, M: int):
        """Launch M particles for mouth density n_mouth (1/m3) over dt (s); advances the coverage."""
        t0 = time.perf_counter()
        self._check(self.lib.aldmc_step(self._h, float(n_mouth), float(dt), int(M)))
        self.wall_time += time.perf_counter() - t0
        self.particles += int(M) if n_mouth > 0 and dt > 0 else 0
        seg = ctypes.c_ulonglong(); self.lib.aldmc_last_segments(self._h, ctypes.byref(seg)); self.segments += seg.value

    def reset(self):
        self._check(self.lib.aldmc_reset(self._h))
        self.wall_time = 0.0; self.particles = 0; self.segments = 0

    def face_coverage(self) -> np.ndarray:
        out = np.zeros(self.n_faces); self._check(self.lib.aldmc_face_coverage(self._h, _ptr(out, ctypes.c_double))); return out

    def set_face_coverage(self, theta):
        a = np.ascontiguousarray(np.asarray(theta, float))
        if a.shape != (self.n_faces,):
            raise ValueError("theta must have n_faces entries")
        self._check(self.lib.aldmc_set_face_coverage(self._h, _ptr(a, ctypes.c_double)))

    def face_hits(self):
        """Collisions and sticking events per face in the last step."""
        hits = np.zeros(self.n_faces, dtype=np.uint64); sticks = np.zeros(self.n_faces, dtype=np.uint64)
        self._check(self.lib.aldmc_face_hits(self._h, _ptr(hits, ctypes.c_ulonglong), _ptr(sticks, ctypes.c_ulonglong)))
        return hits, sticks

    def faces(self) -> dict:
        """Per-face voxel index, direction (0..5), weight and centre coordinates (m)."""
        n = self.n_faces
        vox = np.zeros(n, dtype=np.int64); d = np.zeros(n, dtype=np.int32); w = np.zeros(n)
        cx, cy, cz = np.zeros(n), np.zeros(n), np.zeros(n)
        self._check(self.lib.aldmc_face_info(self._h, _ptr(vox, ctypes.c_longlong), _ptr(d, ctypes.c_int), _ptr(w, ctypes.c_double),
                                             _ptr(cx, ctypes.c_double), _ptr(cy, ctypes.c_double), _ptr(cz, ctypes.c_double)))
        return dict(voxel=vox, direction=d, weight=w, x=cx, y=cy, z=cz)

    def slice_coverage(self):
        """(z_centres, weighted mean coverage per slice, weighted face area per slice in m2); slices without faces
        have area 0."""
        out = np.zeros(self._nz); wt = np.zeros(self._nz)
        self._check(self.lib.aldmc_slice_coverage(self._h, _ptr(out, ctypes.c_double), _ptr(wt, ctypes.c_double)))
        return self.geometry.z_centres(), out, wt

    def bottom_coverage(self) -> float:
        v = ctypes.c_double(); self._check(self.lib.aldmc_bottom_coverage(self._h, ctypes.byref(v))); return v.value

    def counts(self, last: bool = False) -> dict:
        e, s, l, x = (ctypes.c_ulonglong() for _ in range(4))
        fn = self.lib.aldmc_last_counts if last else self.lib.aldmc_counts
        self._check(fn(self._h, ctypes.byref(e), ctypes.byref(s), ctypes.byref(l), ctypes.byref(x)))
        return dict(escaped=e.value, stuck=s.value, lost=l.value, left=x.value)

    def rate(self) -> dict:
        """Wall time, particles and particles per second per thread so far."""
        th = self.threads
        pps = self.particles / self.wall_time if self.wall_time > 0 else float("nan")
        return dict(wall_s=self.wall_time, particles=self.particles, threads=th, particles_per_s=pps,
                    particles_per_s_per_thread=pps / th, segments_per_particle=self.segments / max(1, self.particles))


# ---------------------------------------------------------------------------------------------------------------
@dataclass
class VoxelDrivenResult:
    """Unpacks like aldconf.simulate_driven's (x, theta, hist_t, hist): x = slice centres (m), theta = slice
    coverage at the end, hist[:, i] = slice coverage at record_depths[i] after every sample."""
    x: np.ndarray
    theta: np.ndarray
    t: np.ndarray
    hist: np.ndarray
    bottom: np.ndarray            # floor coverage history
    slice_weight: np.ndarray      # weighted face area per slice (m2)
    snapshots: dict               # time -> slice coverage profile
    counts: dict
    rate: dict
    solver: Optional[VoxelMC] = field(default=None, repr=False)

    def __iter__(self):
        return iter((self.x, self.theta, self.t, self.hist))


def _resample(t_series, n0_series, dt):
    """Uniform grid of spacing dt over the series' span with the density interpolated (subsampling)."""
    t = np.asarray(t_series, float); n0 = np.asarray(n0_series, float)
    tt = np.arange(t[0], t[-1] + 1e-12, dt)
    if tt[-1] < t[-1] - 1e-9:
        tt = np.append(tt, t[-1])
    return tt, np.interp(tt, t, n0)


def simulate_driven_voxel(geometry: VoxelGeometry, params: Params, t_series, n0_series, M: int = 100_000,
                          substeps_per_sample: int = 1, record_depths=(0.0, 1.0), resample_dt: Optional[float] = None,
                          snapshot_times=(), stop_at_bottom: Optional[float] = None, t_stop: Optional[float] = None,
                          verbose: bool = False, solver: Optional[VoxelMC] = None) -> VoxelDrivenResult:
    """One half-cycle with a time-varying mouth density n0(t) [1/m3], the voxel twin of aldconf.simulate_driven.

    Each sample interval of the (optionally resampled) series is split into `substeps_per_sample` steps of M particles
    with the density at the sub-step mid-point. `record_depths` are relative depths (0 mouth, 1 bottom) whose slice
    coverage is recorded after every sample; the floor coverage history is `bottom`. `stop_at_bottom` ends the run
    once the floor coverage exceeds it (after the sample); `t_stop` ends it at that time."""
    if resample_dt is not None:
        t_series, n0_series = _resample(t_series, n0_series, resample_dt)
    t_series = np.asarray(t_series, float); n0_series = np.asarray(n0_series, float)
    mc = solver or VoxelMC(geometry, params)
    z, _, wt = mc.slice_coverage()
    with_faces = np.where(wt > 0)[0]
    kmax = int(with_faces.max()) if len(with_faces) else len(z) - 1
    idx = [min(kmax, int(round(d * kmax))) for d in record_depths]
    hist_t, hist, bottom = [], [], []
    snaps = {}; snap_left = sorted(snapshot_times)
    for i in range(1, len(t_series)):
        dt = (t_series[i] - t_series[i - 1]) / substeps_per_sample
        for k in range(substeps_per_sample):
            tt = t_series[i - 1] + (k + 0.5) * dt
            n0 = float(np.interp(tt, t_series, n0_series))
            mc.step(max(0.0, n0), dt, M)
        tnow = t_series[i]
        _, prof, _ = mc.slice_coverage()
        b = mc.bottom_coverage()
        hist_t.append(tnow); hist.append(prof[idx].copy()); bottom.append(b)
        while snap_left and tnow >= snap_left[0] - 1e-9:
            snaps[snap_left.pop(0)] = prof.copy()
        if verbose:
            r = mc.rate()
            print(f"  t = {tnow:7.3f} s  bottom {b:.4f}  mouth {prof[idx[0]]:.4f}  {r['particles_per_s']:.0f} p/s  "
                  f"{r['segments_per_particle']:.0f} seg/p  wall {r['wall_s']:.0f} s", flush=True)
        if stop_at_bottom is not None and b >= stop_at_bottom:
            break
        if t_stop is not None and tnow >= t_stop - 1e-9:
            break
    _, prof, wt = mc.slice_coverage()
    return VoxelDrivenResult(x=z, theta=prof, t=np.array(hist_t), hist=np.array(hist), bottom=np.array(bottom),
                             slice_weight=wt, snapshots=snaps, counts=mc.counts(), rate=mc.rate(), solver=mc)


def simulate_constant_voxel(geometry: VoxelGeometry, params: Params, n0: float, t_end: float, dt: float, M: int = 100_000,
                            **kw) -> VoxelDrivenResult:
    """Constant mouth density n0 for t_end seconds in steps of dt (two-point series, one step per sample)."""
    n = int(round(t_end / dt))
    t_series = np.linspace(0.0, n * dt, n + 1)
    return simulate_driven_voxel(geometry, params, t_series, np.full(n + 1, float(n0)), M=M, substeps_per_sample=1, **kw)


def time_to_level(t, series, level: float = 0.9) -> float:
    """First crossing of `level` (linear interpolation), nan if never; aldconf.recess.time_to_coverage's rule."""
    t = np.asarray(t, float); y = np.asarray(series, float)
    idx = np.where(y >= level)[0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]))


__all__ = ["Params", "VoxelGeometry", "VoxelMC", "VoxelDrivenResult", "AldMCError", "AldMCNotBuilt", "is_available",
           "load_library", "library_path", "version", "hole", "tapered_hole", "slit", "recessed_slit", "from_array",
           "open_box", "flat_surface", "simulate_driven_voxel", "simulate_constant_voxel", "time_to_level",
           "default_threads", "BUILD_COMMAND"]
