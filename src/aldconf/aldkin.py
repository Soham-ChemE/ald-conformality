"""ctypes binding of libaldkin (lib/aldkin): the multiscale ALD wall kinetics of the OpenFOAM boundary conditions as
a standalone library.

    from aldconf.aldkin import AldKin, Params, is_available
    k = AldKin(Params(s=5e-3, vbar=410, N0=4e18, c0=3.79e20, d=100e-9, AR=100, holeDensity=1.78e13, nCells=200))
    k.set_faces(radius=[0.00067])                 # one face at the wafer centre
    flux, keff = k.step(T_face=[0.02], dt=1e-3)    # advance one solver step with T = n/c0 on each face
    k.theta("bottom")                              # hole-bottom coverage per face

The library is looked up in this order: the file named by ALDKIN_LIBRARY, the directory named by ALDKIN_BUILD_DIR,
<repo>/lib/aldkin/build (also its Release/ subdirectory for multi-config generators) and finally the system loader.
`is_available()` says whether it was found; `AldKin` raises AldKinNotBuilt otherwise with the build command.

Quantities and units are those of include/aldkin.h: T_face is the gas density relative to c0 (set c0 = 1 to pass
absolute densities in 1/m3), dt in s, K_eff in m/s, the uptake flux in molecules/(m2 s). Scheme "euler" is the
boundary conditions' forward-Euler update (bit-identical to them when stepped once per solver step); "exponential"
is the exact relaxation used by aldconf.solver.
"""
from __future__ import annotations

import ctypes
import os
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUILD_DIR = ROOT / "lib" / "aldkin" / "build"
BUILD_COMMAND = "cmake -S lib/aldkin -B lib/aldkin/build -DCMAKE_BUILD_TYPE=Release && cmake --build lib/aldkin/build"

SCHEMES = {"euler": 0, "exponential": 1}
THETA = {"flat": 0, "mid": 1, "bottom": 2, "hole": 3, "blocked_flat": 4, "blocked_bottom": 5, "blocked_hole": 6}
PROFILE = {"AR": 0, "density": 1}
ERRORS = {0: "ok", 1: "ALDKIN_ERR_ARG", 2: "ALDKIN_ERR_PARAM", 3: "ALDKIN_ERR_STATE", 4: "ALDKIN_ERR_ALLOC", 5: "ALDKIN_ERR_INTERNAL"}


class AldKinError(RuntimeError):
    """A library call returned an error code; the message is aldkin_last_error()."""


class AldKinNotBuilt(AldKinError):
    """libaldkin was not found."""


class _Params(ctypes.Structure):
    _fields_ = [("s", ctypes.c_double), ("r", ctypes.c_double), ("kdes", ctypes.c_double), ("vbar", ctypes.c_double),
                ("N0", ctypes.c_double), ("c0", ctypes.c_double), ("d", ctypes.c_double), ("AR", ctypes.c_double),
                ("holeDensity", ctypes.c_double), ("nCells", ctypes.c_int), ("flatOnly", ctypes.c_int),
                ("scheme", ctypes.c_int),
                ("nu", ctypes.c_double), ("vbarB", ctypes.c_double), ("sB", ctypes.c_double), ("kB", ctypes.c_double),
                ("kEtch", ctypes.c_double), ("cB0", ctypes.c_double)]


@dataclass
class Params:
    """The aldkin_params struct. Defaults are the boundary conditions' (aldkin_params_default)."""
    s: float = 1.0
    r: float = 0.0
    kdes: float = 0.0
    vbar: float = 410.0
    N0: float = 4e18
    c0: float = 1.0
    d: float = 1e-7
    AR: float = 50.0
    holeDensity: float = 0.0
    nCells: int = 100
    flatOnly: bool = False
    scheme: str = "euler"
    # two-species option (library 1.1): a byproduct released at nu per adsorbed precursor, mean speed vbarB, that
    # re-adsorbs on free sites (sB), desorbs (kB, 1/s), etches covered sites (kEtch) and sits at cB0 [1/m3] at the
    # mouth; all 0 = the single-species path, bit for bit (aldconf.multispecies is the Python twin)
    nu: float = 0.0
    vbarB: float = 0.0
    sB: float = 0.0
    kB: float = 0.0
    kEtch: float = 0.0
    cB0: float = 0.0

    def to_struct(self) -> _Params:
        if self.scheme not in SCHEMES:
            raise ValueError(f"scheme must be one of {sorted(SCHEMES)}")
        d = asdict(self)
        d["flatOnly"] = int(bool(self.flatOnly))
        d["scheme"] = SCHEMES[self.scheme]
        d["nCells"] = int(self.nCells)
        return _Params(**d)

    @classmethod
    def from_struct(cls, p: _Params) -> "Params":
        inv = {v: k for k, v in SCHEMES.items()}
        return cls(s=p.s, r=p.r, kdes=p.kdes, vbar=p.vbar, N0=p.N0, c0=p.c0, d=p.d, AR=p.AR, holeDensity=p.holeDensity,
                   nCells=p.nCells, flatOnly=bool(p.flatOnly), scheme=inv[p.scheme],
                   nu=p.nu, vbarB=p.vbarB, sB=p.sB, kB=p.kB, kEtch=p.kEtch, cB0=p.cB0)


def _candidates():
    """Library files to try, in order. ALDKIN_LIBRARY (a file) or ALDKIN_BUILD_DIR (a directory) replace the default
    search in <repo>/lib/aldkin/build when set."""
    names = {"darwin": ["libaldkin.dylib"], "win32": ["aldkin.dll", "libaldkin.dll"]}.get(sys.platform, ["libaldkin.so"])
    env_file = os.environ.get("ALDKIN_LIBRARY")
    env_dir = os.environ.get("ALDKIN_BUILD_DIR")
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
    """Load libaldkin once and declare the prototypes. `path` overrides the search."""
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
    if found is None and not (os.environ.get("ALDKIN_LIBRARY") or os.environ.get("ALDKIN_BUILD_DIR")):
        # last resort: the system loader
        try:
            lib = ctypes.CDLL({"darwin": "libaldkin.dylib", "win32": "aldkin.dll"}.get(sys.platform, "libaldkin.so"))
            found = Path(lib._name)
        except OSError:
            raise AldKinNotBuilt("libaldkin not found (looked at: " + ", ".join(tried) + "); build it with: " + BUILD_COMMAND)
    elif found is None:
        raise AldKinNotBuilt("libaldkin not found (looked at: " + ", ".join(tried) + "); build it with: " + BUILD_COMMAND)
    else:
        lib = ctypes.CDLL(str(found))
    D, I, P = ctypes.c_double, ctypes.c_int, ctypes.c_void_p
    PD = ctypes.POINTER(ctypes.c_double)
    lib.aldkin_version.restype = ctypes.c_char_p
    lib.aldkin_version.argtypes = []
    lib.aldkin_params_default.restype = None
    lib.aldkin_params_default.argtypes = [ctypes.POINTER(_Params)]
    lib.aldkin_create.restype = P
    lib.aldkin_create.argtypes = [ctypes.POINTER(_Params)]
    lib.aldkin_clone.restype = P
    lib.aldkin_clone.argtypes = [P]
    lib.aldkin_destroy.restype = None
    lib.aldkin_destroy.argtypes = [P]
    lib.aldkin_last_error.restype = ctypes.c_char_p
    lib.aldkin_last_error.argtypes = [P]
    lib.aldkin_get_params.restype = I
    lib.aldkin_get_params.argtypes = [P, ctypes.POINTER(_Params)]
    lib.aldkin_set_faces.restype = I
    lib.aldkin_set_faces.argtypes = [P, I, PD]
    lib.aldkin_set_profile.restype = I
    lib.aldkin_set_profile.argtypes = [P, I, I, PD, PD]
    lib.aldkin_set_profiles.restype = I
    lib.aldkin_set_profiles.argtypes = [P, I, PD, PD, PD]
    lib.aldkin_n_faces.restype = I
    lib.aldkin_n_faces.argtypes = [P]
    lib.aldkin_n_cells.restype = I
    lib.aldkin_n_cells.argtypes = [P]
    lib.aldkin_face_pattern.restype = I
    lib.aldkin_face_pattern.argtypes = [P, PD, PD]
    lib.aldkin_step.restype = I
    lib.aldkin_step.argtypes = [P, PD, D, PD, PD]
    lib.aldkin_evaluate.restype = I
    lib.aldkin_evaluate.argtypes = [P, PD, PD, PD]
    lib.aldkin_theta.restype = I
    lib.aldkin_theta.argtypes = [P, I, PD]
    lib.aldkin_set_theta.restype = I
    lib.aldkin_set_theta.argtypes = [P, I, PD]
    lib.aldkin_keff.restype = I
    lib.aldkin_keff.argtypes = [P, PD]
    lib.aldkin_mixed_bc.restype = I
    lib.aldkin_mixed_bc.argtypes = [P, PD, PD, PD, PD]
    lib.aldkin_desorption_flux.restype = I
    lib.aldkin_desorption_flux.argtypes = [P, PD]
    lib.aldkin_hole_gas.restype = I
    lib.aldkin_hole_gas.argtypes = [P, I, PD]
    lib.aldkin_set_byproduct_mouth.restype = I
    lib.aldkin_set_byproduct_mouth.argtypes = [P, PD]
    lib.aldkin_byproduct_flux.restype = I
    lib.aldkin_byproduct_flux.argtypes = [P, PD]
    lib.aldkin_byproduct_balance.restype = I
    lib.aldkin_byproduct_balance.argtypes = [P, I, PD]
    lib.aldkin_hole_byproduct_gas.restype = I
    lib.aldkin_hole_byproduct_gas.argtypes = [P, I, PD]
    _lib, _lib_path = lib, found
    return lib


def is_available() -> bool:
    try:
        load_library()
        return True
    except AldKinNotBuilt:
        return False


def library_path() -> Optional[Path]:
    load_library()
    return _lib_path


def version() -> str:
    return load_library().aldkin_version().decode()


def default_params() -> Params:
    p = _Params()
    load_library().aldkin_params_default(ctypes.byref(p))
    return Params.from_struct(p)


def _arr(x, n: Optional[int] = None) -> np.ndarray:
    a = np.ascontiguousarray(np.asarray(x, dtype=np.float64).ravel())
    if n is not None and a.size != n:
        raise ValueError(f"expected {n} values, got {a.size}")
    return a


def _ptr(a: Optional[np.ndarray]):
    return None if a is None else a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


class AldKin:
    """One kinetics handle: a set of wall faces with their coverages."""

    def __init__(self, params: Params = Params(), _handle=None):
        self._lib = load_library()
        self._h = None
        if _handle is not None:
            self._h = _handle
            return
        self._h = self._lib.aldkin_create(ctypes.byref(params.to_struct()))
        if not self._h:
            raise AldKinError(self._lib.aldkin_last_error(None).decode())

    # -- lifetime ------------------------------------------------------------------------------------------------
    def close(self):
        if getattr(self, "_h", None):
            self._lib.aldkin_destroy(self._h)
            self._h = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def clone(self) -> "AldKin":
        h = self._lib.aldkin_clone(self._h)
        if not h:
            raise AldKinError(self._lib.aldkin_last_error(None).decode())
        return AldKin(_handle=h)

    def _check(self, rc: int, what: str):
        if rc != 0:
            raise AldKinError(f"{what}: {ERRORS.get(rc, rc)}: {self._lib.aldkin_last_error(self._h).decode()}")

    # -- setup ---------------------------------------------------------------------------------------------------
    @property
    def params(self) -> Params:
        p = _Params()
        self._check(self._lib.aldkin_get_params(self._h, ctypes.byref(p)), "aldkin_get_params")
        return Params.from_struct(p)

    @property
    def n_faces(self) -> int:
        return self._lib.aldkin_n_faces(self._h)

    @property
    def n_cells(self) -> int:
        return self._lib.aldkin_n_cells(self._h)

    def set_faces(self, radius: Optional[Sequence[float]] = None, n: Optional[int] = None):
        """Define the faces by their centre radii (m), or `n` faces at radius 0 when no profiles are used."""
        if radius is None:
            if n is None:
                raise ValueError("give radius or n")
            self._check(self._lib.aldkin_set_faces(self._h, int(n), None), "aldkin_set_faces")
        else:
            r = _arr(radius)
            self._check(self._lib.aldkin_set_faces(self._h, r.size, _ptr(r)), "aldkin_set_faces")

    def set_profile(self, which: str, r: Sequence[float], values: Sequence[float]):
        """Radial table for 'AR' or 'density' (independent radii); an empty table removes it."""
        rr, vv = _arr(r), _arr(values)
        if rr.size != vv.size:
            raise ValueError("r and values must have the same length")
        self._check(self._lib.aldkin_set_profile(self._h, PROFILE[which], rr.size, _ptr(rr) if rr.size else None,
                                                 _ptr(vv) if vv.size else None), "aldkin_set_profile")

    def set_profiles(self, r: Sequence[float], AR: Optional[Sequence[float]] = None, density: Optional[Sequence[float]] = None):
        """Both tables on one radius vector (None keeps the scalar for that quantity)."""
        rr = _arr(r)
        a = None if AR is None else _arr(AR, rr.size)
        dd = None if density is None else _arr(density, rr.size)
        self._check(self._lib.aldkin_set_profiles(self._h, rr.size, _ptr(rr), _ptr(a), _ptr(dd)), "aldkin_set_profiles")

    def face_pattern(self):
        """(AR, holeDensity) per face as used."""
        n = self.n_faces
        a, d = np.zeros(n), np.zeros(n)
        self._check(self._lib.aldkin_face_pattern(self._h, _ptr(a), _ptr(d)), "aldkin_face_pattern")
        return a, d

    # -- stepping ------------------------------------------------------------------------------------------------
    def step(self, T_face, dt: float):
        """Advance every face by dt given T_face = n/c0 per face; returns (uptake_flux, K_eff)."""
        n = self.n_faces
        T = _arr(T_face, n if n else None)          # with no faces the library reports the state error itself
        flux, keff = np.zeros(n), np.zeros(n)
        self._check(self._lib.aldkin_step(self._h, _ptr(T), float(dt), _ptr(flux), _ptr(keff)), "aldkin_step")
        return flux, keff

    def evaluate(self, T_face):
        """K_eff and flux for new face values without advancing; returns (uptake_flux, K_eff)."""
        n = self.n_faces
        T = _arr(T_face, n if n else None)
        flux, keff = np.zeros(n), np.zeros(n)
        self._check(self._lib.aldkin_evaluate(self._h, _ptr(T), _ptr(flux), _ptr(keff)), "aldkin_evaluate")
        return flux, keff

    # -- state ---------------------------------------------------------------------------------------------------
    def theta(self, which: str = "flat") -> np.ndarray:
        """'flat', 'mid', 'bottom', 'blocked_flat', 'blocked_bottom' (n_faces) or 'hole', 'blocked_hole' (n_faces x nCells)."""
        n, N = self.n_faces, self.n_cells
        full = which in ("hole", "blocked_hole")
        out = np.zeros(n * N if full else n)
        self._check(self._lib.aldkin_theta(self._h, THETA[which], _ptr(out)), "aldkin_theta")
        return out.reshape(n, N) if full else out

    def set_theta(self, which: str, values):
        n, N = self.n_faces, self.n_cells
        v = _arr(values, n * N if which in ("hole", "blocked_hole") else n)
        self._check(self._lib.aldkin_set_theta(self._h, THETA[which], _ptr(v)), "aldkin_set_theta")

    # -- two-species option ----------------------------------------------------------------------------------------
    def set_byproduct_mouth(self, nB_face=None):
        """Byproduct number density (1/m3, absolute) at every face and hole mouth; None restores the scalar cB0."""
        if nB_face is None:
            self._check(self._lib.aldkin_set_byproduct_mouth(self._h, None), "aldkin_set_byproduct_mouth")
        else:
            v = _arr(nB_face, self.n_faces)
            self._check(self._lib.aldkin_set_byproduct_mouth(self._h, _ptr(v)), "aldkin_set_byproduct_mouth")

    def byproduct_flux(self) -> np.ndarray:
        """Net byproduct flux into the chamber per unit wafer area from the last step/evaluate [molecules/(m2 s)]."""
        out = np.zeros(self.n_faces)
        self._check(self._lib.aldkin_byproduct_flux(self._h, _ptr(out)), "aldkin_byproduct_flux")
        return out

    def byproduct_balance(self, face: int) -> dict:
        """Hole bookkeeping rates per unit mouth area: produced, desorbed, escaped, adsorbed, etched."""
        out = np.zeros(5)
        self._check(self._lib.aldkin_byproduct_balance(self._h, int(face), _ptr(out)), "aldkin_byproduct_balance")
        return dict(zip(("produced", "desorbed", "escaped", "adsorbed", "etched"), out.tolist()))

    def hole_byproduct_gas(self, face: int) -> np.ndarray:
        """Absolute byproduct density (1/m3) along the hole of one face from the last step/evaluate."""
        out = np.zeros(self.n_cells)
        self._check(self._lib.aldkin_hole_byproduct_gas(self._h, int(face), _ptr(out)), "aldkin_hole_byproduct_gas")
        return out

    def keff(self) -> np.ndarray:
        out = np.zeros(self.n_faces)
        self._check(self._lib.aldkin_keff(self._h, _ptr(out)), "aldkin_keff")
        return out

    def mixed_bc(self, D_over_delta):
        """(refValue, refGrad, valueFraction) of the mixed boundary condition for D/delta per face."""
        n = self.n_faces
        Dd = _arr(D_over_delta, n)
        rv, rg, vf = np.zeros(n), np.zeros(n), np.zeros(n)
        self._check(self._lib.aldkin_mixed_bc(self._h, _ptr(Dd), _ptr(rv), _ptr(rg), _ptr(vf)), "aldkin_mixed_bc")
        return rv, rg, vf

    def desorption_flux(self) -> np.ndarray:
        out = np.zeros(self.n_faces)
        self._check(self._lib.aldkin_desorption_flux(self._h, _ptr(out)), "aldkin_desorption_flux")
        return out

    def hole_gas(self, face: int) -> np.ndarray:
        """Relative gas density n/n_mouth along the hole of one face from the last step/evaluate."""
        out = np.zeros(self.n_cells)
        self._check(self._lib.aldkin_hole_gas(self._h, int(face), _ptr(out)), "aldkin_hole_gas")
        return out


def drive(k: AldKin, t_series, T_series, substeps_per_sample: int = 1, record=("flat", "bottom")):
    """Step one or more faces through a sampled face-density history T(t) (same value on every face).

    Mirrors aldconf.solver.simulate_driven: each sample interval is split into `substeps_per_sample` steps and the
    drive is interpolated linearly at the mid-point of each sub-step. Returns (t_hist, {name: array(ntime, n_faces)})
    with the coverages after each sample interval."""
    t_series = np.asarray(t_series, float)
    T_series = np.asarray(T_series, float)
    n = k.n_faces
    hist_t, hist = [], {name: [] for name in record}
    for i in range(1, len(t_series)):
        dt = (t_series[i] - t_series[i - 1]) / substeps_per_sample
        for j in range(substeps_per_sample):
            tt = t_series[i - 1] + (j + 0.5) * dt
            Tf = float(np.interp(tt, t_series, T_series))
            k.step(np.full(n, Tf), dt)
        hist_t.append(t_series[i])
        for name in record:
            hist[name].append(k.theta(name).copy())
    return np.array(hist_t), {name: np.array(v) for name, v in hist.items()}
