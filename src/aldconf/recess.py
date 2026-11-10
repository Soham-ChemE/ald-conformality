"""Nested feature model: a vertical slit whose walls carry lateral recesses (3D-NAND word-line cavities fed from
a gate-line slit).

Geometry (per unit width into the page): a slit of gap w_s and depth L_s; each wall carries recesses of gap g_r
(the word-line height), lateral length L_r and vertical pitch p_r, so a fraction f_r = g_r / p_r of the slit wall
is recess mouth and 1 - f_r is flat wall. Every recess is a 1D slit sub-problem (Knudsen diffusivity 2 g_r vbar/3,
Langmuir sink on its two faces, closed end) fed by the slit gas density at its mouth. The slit wall's net uptake
per unit wall area per unit gas density is

    K_wall(x) = (1 - f_r) (vbar/4) [s (1 - theta_s(x)) + r] + f_r j_r(x),

with j_r(x) the recess mouth flux per unit mouth area per unit mouth density from the quasi-steady sub-problem:
the same construction as the chamber's patternedSurface boundary condition one level down. Gas is quasi-steady at
both levels (Ylilammi et al. 2018); coverages are advanced per substep. The recess end wall (area g_r per unit
width) is neglected, as are film narrowing and the recess-to-recess interaction through the channel holes.

Limits: f_r = 0 is the plain slit; L_r -> 0 (recess interior at the mouth density) is a slit with its wall uptake
multiplied by the area factor A_f = (p_r - g_r + 2 L_r) / p_r, i.e. `AreaEquivalentSlit`.
"""
from dataclasses import dataclass
import numpy as np

from .geometry import Slit
from .transport import mean_speed, knudsen_diffusivity
from .kinetics import Precursor
from .solver import _solve_gas


@dataclass
class RecessedSlit:
    width: float            # slit gap w_s (m)
    depth: float            # slit depth L_s (m)
    recess_gap: float       # g_r (m); 0 = no recesses
    recess_length: float    # L_r (m)
    recess_pitch: float     # p_r (m)

    def __post_init__(self):
        if self.recess_gap < 0 or self.recess_gap > self.recess_pitch:
            raise ValueError("recess_gap must lie in [0, recess_pitch]")

    @property
    def recessed_fraction(self):
        return self.recess_gap / self.recess_pitch

    @property
    def area_factor(self):
        """Wall area per unit area of slit wall: flat part plus the two faces of each recess (end wall excluded)."""
        return (self.recess_pitch - self.recess_gap + 2.0 * self.recess_length) / self.recess_pitch

    def wall_area_per_width(self):
        """Total wall area per unit width into the page (both slit walls with their recesses)."""
        return 2.0 * self.depth * self.area_factor

    def wall_area_ratio(self, slit_pitch):
        """A_wall / A_wafer for a wafer tiled with these slits at lateral pitch `slit_pitch` (patterned area only)."""
        return self.wall_area_per_width() / slit_pitch

    def plain_slit(self):
        return Slit(self.width, self.depth)


@dataclass
class AreaEquivalentSlit(Slit):
    """Slit whose wall sink per unit volume is multiplied by `area_factor` (same total wall area as a RecessedSlit
    whose recesses are so short that their interior sits at the mouth density)."""
    area_factor: float = 1.0

    def surface_to_volume(self, opening):
        return self.area_factor * 2.0 / np.asarray(opening)


@dataclass
class RecessResult:
    x: np.ndarray                 # slit depth coordinate (cell centres)
    y: np.ndarray                 # recess depth coordinate (cell centres)
    theta_slit: np.ndarray        # flat-wall coverage along the slit at the end
    theta_recess: np.ndarray      # (n_slit, n_recess) recess coverage at the end
    t: np.ndarray                 # history times
    hist_slit_bottom: np.ndarray  # flat-wall coverage at the slit bottom vs time
    hist_slit_mouth: np.ndarray
    hist_recess_bottom_top: np.ndarray     # recess-bottom coverage for the recess at the slit mouth
    hist_recess_bottom_mid: np.ndarray     # ... at mid depth
    hist_recess_bottom_deep: np.ndarray    # ... at the slit bottom
    hist_recess_bottom_min: np.ndarray     # minimum over all recesses
    snapshots: dict               # time -> (theta_slit, theta_recess)


def _thomas_batch(lower, diag, upper, rhs):
    """Solve a batch of tridiagonal systems; arrays (B, N) with lower[:, 0] and upper[:, -1] unused."""
    B, N = diag.shape
    cp = np.zeros((B, N)); dp = np.zeros((B, N))
    cp[:, 0] = upper[:, 0] / diag[:, 0]; dp[:, 0] = rhs[:, 0] / diag[:, 0]
    for i in range(1, N):
        m = diag[:, i] - lower[:, i] * cp[:, i - 1]
        cp[:, i] = upper[:, i] / m if i < N - 1 else 0.0
        dp[:, i] = (rhs[:, i] - lower[:, i] * dp[:, i - 1]) / m
    out = np.zeros((B, N)); out[:, -1] = dp[:, -1]
    for i in range(N - 2, -1, -1):
        out[:, i] = dp[:, i] - cp[:, i] * out[:, i + 1]
    return out


def _solve_recesses(D_r, dy, k_r):
    """Quasi-steady relative density in every recess (mouth density 1, closed end) and the mouth flux per unit
    mouth density [m/s]. k_r: (n_slit, n_recess) first-order sink."""
    B, N = k_r.shape
    a = D_r / dy ** 2
    diag = 2.0 * a + k_r
    diag[:, -1] = a + k_r[:, -1]          # closed end: only the inner face conducts
    diag[:, 0] += a                       # Dirichlet ghost at the mouth
    lower = np.full((B, N), -a); upper = np.full((B, N), -a)
    rhs = np.zeros((B, N)); rhs[:, 0] = 2.0 * a
    nrel = _thomas_batch(lower, diag, upper, rhs)
    j = D_r * (1.0 - nrel[:, 0]) / (0.5 * dy)
    return nrel, j


def simulate_recessed(geom: RecessedSlit, pr: Precursor, T: float, t_series, n0_series, substeps_per_sample=5,
                      n_slit=200, n_recess=50, snapshot_times=(), knudsen_factor=1.0) -> RecessResult:
    """One half-cycle with a time-varying gas density n0(t) [1/m3] at the slit mouth.

    `t_series`/`n0_series` as in aldconf.simulate_driven (a constant density is the two-point series
    [0, t_dose], [n0, n0] with substeps_per_sample = number of substeps)."""
    t_series = np.asarray(t_series, float); n0_series = np.asarray(n0_series, float)
    Ns, Nr = n_slit, n_recess
    x = (np.arange(Ns) + 0.5) * geom.depth / Ns
    vbar = mean_speed(T, pr.mass_amu)
    s, r, N0 = pr.sticking, pr.recombination, pr.site_density
    slit = geom.plain_slit()
    D_s = knudsen_diffusivity(slit.hydraulic_diameter(geom.width), T, pr.mass_amu, knudsen_factor)
    D_face = np.full(Ns + 1, D_s)
    sv = 2.0 / geom.width
    f_r = geom.recessed_fraction
    have_recess = geom.recess_gap > 0 and geom.recess_length > 0
    if have_recess:
        y = (np.arange(Nr) + 0.5) * geom.recess_length / Nr
        dy = geom.recess_length / Nr
        D_r = knudsen_diffusivity(2.0 * geom.recess_gap, T, pr.mass_amu, knudsen_factor)
        sv_r = 2.0 / geom.recess_gap
    else:
        y = np.zeros(0); Nr = 0
    theta_s = np.zeros(Ns); theta_r = np.zeros((Ns, Nr))
    hist_t = []; hb = []; hm = []; rt = []; rm = []; rd = []; rmin = []
    snaps = {}; snap_left = sorted(snapshot_times)
    for i in range(1, len(t_series)):
        dt = (t_series[i] - t_series[i - 1]) / substeps_per_sample
        for k in range(substeps_per_sample):
            tt = t_series[i - 1] + (k + 0.5) * dt
            n0 = float(np.interp(tt, t_series, n0_series))
            if have_recess:
                k_r = sv_r * (vbar / 4.0) * (s * (1.0 - theta_r) + r)
                nrel, j_r = _solve_recesses(D_r, dy, k_r)
                K_wall = (1.0 - f_r) * (vbar / 4.0) * (s * (1.0 - theta_s) + r) + f_r * j_r
            else:
                K_wall = (vbar / 4.0) * (s * (1.0 - theta_s) + r)
            if n0 > 0:
                n_s = _solve_gas(x, D_face, sv * K_wall, n0)
            else:
                n_s = np.zeros(Ns)
            theta_s = 1.0 - (1.0 - theta_s) * np.exp(-s * (vbar / 4.0) * n_s * dt / N0)
            if have_recess:
                theta_r = 1.0 - (1.0 - theta_r) * np.exp(-s * (vbar / 4.0) * (n_s[:, None] * nrel) * dt / N0)
        tnow = t_series[i]
        hist_t.append(tnow); hb.append(theta_s[-1]); hm.append(theta_s[0])
        if have_recess:
            rt.append(theta_r[0, -1]); rm.append(theta_r[Ns // 2, -1]); rd.append(theta_r[-1, -1]); rmin.append(theta_r[:, -1].min())
        else:
            rt.append(np.nan); rm.append(np.nan); rd.append(np.nan); rmin.append(np.nan)
        while snap_left and tnow >= snap_left[0] - 1e-9:
            snaps[snap_left.pop(0)] = (theta_s.copy(), theta_r.copy())
    return RecessResult(x=x, y=y, theta_slit=theta_s, theta_recess=theta_r, t=np.array(hist_t),
                        hist_slit_bottom=np.array(hb), hist_slit_mouth=np.array(hm),
                        hist_recess_bottom_top=np.array(rt), hist_recess_bottom_mid=np.array(rm),
                        hist_recess_bottom_deep=np.array(rd), hist_recess_bottom_min=np.array(rmin), snapshots=snaps)


def simulate_recessed_constant(geom, pr, T, p_partial, t_dose, n_substeps=400, **kw):
    """Constant mouth density p_partial [Pa] for t_dose [s]; history sampled every t_dose / n_samples."""
    from .kinetics import number_density
    n0 = number_density(p_partial, T)
    n_samples = kw.pop("n_samples", 200)
    sub = max(1, n_substeps // n_samples)
    t_series = np.linspace(0.0, t_dose, n_samples + 1)
    return simulate_recessed(geom, pr, T, t_series, np.full(n_samples + 1, n0), substeps_per_sample=sub, **kw)


def time_to_coverage(t, theta, level=0.9):
    """First time the coverage history crosses `level` (linear interpolation); nan if never."""
    t = np.asarray(t); theta = np.asarray(theta)
    idx = np.where(theta >= level)[0]
    if len(idx) == 0:
        return float("nan")
    i = idx[0]
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - theta[i - 1]) * (t[i] - t[i - 1]) / (theta[i] - theta[i - 1]))
