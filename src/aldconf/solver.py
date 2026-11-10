"""Cycle-by-cycle feature-scale ALD model.

Per dose step, gas transport is quasi-steady (Ylilammi et al. 2018):

    d/dx ( D(x) dn/dx ) = (P/A)(x) * (vbar/4) * [ s (1 - theta(x)) + r ] * n(x),
    n(0) = n0,   dn/dx(L) = 0,

and the coverage advances with the gas profile frozen over a sub-step dt:

    theta <- 1 - (1 - theta) exp( - s (vbar/4) n dt / N0 ).

After each cycle the film thickness grows by gpc * theta at every x and the opening
narrows by walls * d(thickness), which updates D and P/A for the next cycle.
"""
from dataclasses import dataclass, field
from typing import Optional
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from .geometry import Slit, Hole
from .transport import mean_speed, knudsen_diffusivity, gas_diffusivity, bosanquet
from .kinetics import Precursor, number_density


@dataclass
class Process:
    geometry: object                 # Slit or Hole
    precursor: Precursor
    T: float = 573.15                # K
    p_partial: float = 1.0           # Pa, precursor partial pressure at the feature mouth
    p_total: Optional[float] = None  # Pa, total pressure for the Bosanquet correction; None = pure Knudsen
    t_dose: float = 0.1              # s
    cycles: int = 1
    n_cells: int = 200
    n_substeps: int = 200
    knudsen_factor: float = 1.0
    narrowing: bool = True

    @property
    def dose(self):
        """Exposure p * t in Pa s (1 Langmuir = 1.333e-4 Pa s)."""
        return self.p_partial * self.t_dose


@dataclass
class Result:
    x: np.ndarray
    theta: np.ndarray          # coverage at the end of the last dose
    thickness: np.ndarray      # cumulative film thickness (m)
    opening: np.ndarray        # local opening after growth (m)
    n_last: np.ndarray         # gas density at the end of the last dose
    process: Process = field(repr=False)

    @property
    def depth_over_gap(self):
        return self.x / self.process.geometry.gap


def _solve_gas(x, D_face, k, n0):
    """Steady diffusion with first-order sink k(x) [1/s]; Dirichlet at x=0, Neumann at x=L.

    Cell-centred finite volumes, D on faces, tridiagonal solve."""
    N = len(x)
    dx = x[1] - x[0]
    # Faces: D_face has N+1 entries; face 0 is the mouth, face N is the closed end.
    aW = D_face[:-1] / dx ** 2
    aE = D_face[1:] / dx ** 2
    aE[-1] = 0.0  # zero flux at closed end
    diag = aW + aE + k
    lower = -aW[1:]
    upper = -aE[:-1]
    rhs = np.zeros(N)
    # Dirichlet ghost at the mouth: n_ghost = 2 n0 - n_0 -> modifies first row
    diag[0] += aW[0]
    rhs[0] += 2.0 * aW[0] * n0
    ab = np.zeros((3, N))
    ab[0, 1:] = upper
    ab[1, :] = diag
    ab[2, :-1] = lower
    return solve_banded((1, 1), ab, rhs)


def simulate(proc: Process) -> Result:
    g = proc.geometry
    pr = proc.precursor
    L = g.length
    N = proc.n_cells
    x = (np.arange(N) + 0.5) * L / N
    opening = np.full(N, g.gap, dtype=float)
    thickness = np.zeros(N)
    theta = np.zeros(N)
    vbar = mean_speed(proc.T, pr.mass_amu)
    n0 = number_density(proc.p_partial, proc.T)
    dt = proc.t_dose / proc.n_substeps
    n = np.zeros(N)

    for _ in range(proc.cycles):
        theta[:] = 0.0  # each cycle starts on a surface the co-reactant has reset
        for _ in range(proc.n_substeps):
            d_h = g.hydraulic_diameter(opening)
            D = knudsen_diffusivity(d_h, proc.T, pr.mass_amu, proc.knudsen_factor)
            if proc.p_total is not None:
                D = bosanquet(D, gas_diffusivity(proc.T, proc.p_total, pr.mass_amu, pr.sigma))
            D_face = np.empty(N + 1)
            D_face[1:-1] = 2.0 * D[:-1] * D[1:] / (D[:-1] + D[1:])
            D_face[0] = D[0]
            D_face[-1] = D[-1]
            k = g.surface_to_volume(opening) * (vbar / 4.0) * (pr.sticking * (1.0 - theta) + pr.recombination)
            n = _solve_gas(x, D_face, k, n0)
            theta = 1.0 - (1.0 - theta) * np.exp(-pr.sticking * (vbar / 4.0) * n * dt / pr.site_density)
        growth = pr.gpc * theta
        thickness += growth
        if proc.narrowing:
            opening = np.maximum(opening - g.walls() * growth, 0.05 * g.gap)
    return Result(x=x, theta=theta, thickness=thickness, opening=opening, n_last=n, process=proc)


def step_coverage(res: Result, depth: Optional[float] = None, entrance_fraction: float = 0.02):
    """Thickness at `depth` (default: feature bottom) over thickness near the mouth."""
    ref = res.thickness[: max(1, int(entrance_fraction * len(res.x)))].mean()
    if depth is None:
        val = res.thickness[-1]
    else:
        val = np.interp(depth, res.x, res.thickness)
    return float(val / ref) if ref > 0 else 0.0


def penetration_depth(res: Result, level: float = 0.5):
    """Depth where the profile falls to `level` of its entrance value (PD50 by default)."""
    ref = res.thickness[: max(1, int(0.02 * len(res.x)))].mean()
    y = res.thickness / ref if ref > 0 else res.thickness
    below = np.where(y < level)[0]
    if len(below) == 0:
        return float(res.x[-1])
    i = below[0]
    if i == 0:
        return 0.0
    # linear interpolation between cells i-1 and i
    x0, x1, y0, y1 = res.x[i - 1], res.x[i], y[i - 1], y[i]
    return float(x0 + (level - y0) * (x1 - x0) / (y1 - y0))


def dose_for_coverage(proc: Process, target: float = 0.95, t_lo: float = 1e-4, t_hi: float = 1e3):
    """Dose time (s) at which step coverage at the bottom reaches `target`. Returns nan if unreachable."""
    def f(log_t):
        p = Process(**{**proc.__dict__, "t_dose": float(np.exp(log_t))})
        return step_coverage(simulate(p)) - target
    lo, hi = np.log(t_lo), np.log(t_hi)
    if f(hi) < 0:
        return float("nan")
    if f(lo) > 0:
        return t_lo
    return float(np.exp(brentq(f, lo, hi, xtol=1e-3)))
