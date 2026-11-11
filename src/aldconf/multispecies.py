"""Two-species feature model: the precursor P and a reaction byproduct B that can re-adsorb.

Physics. At the chamber pressures of this work (100 to 1,000 Pa) the mean free path is 10 to 120 um, so inside a
100 nm feature the Knudsen number is 100 to 1,000: molecules collide with the walls, not with each other. A byproduct
therefore cannot block the precursor in the gas phase; each species diffuses independently with its own Knudsen
diffusivity D_K = d_h vbar/3 (vbar ~ 1/sqrt(M)). What a byproduct can do is surface-side: re-adsorb on a free site and
block it until it desorbs, or react with the film (etching). For TMA/H2O the byproduct is methane, which does not
re-adsorb at these temperatures, so the single-species model stands; for chloride chemistries (TiCl4/NH3, HfCl4/H2O)
HCl re-adsorption and etching are documented and this model quantifies what they do inside a feature.

Model (one half-cycle of P at the feature mouth). Coverages: theta_P, sites taken by the precursor (irreversible),
theta_B, sites blocked by the byproduct (reversible); f = 1 - theta_P - theta_B is the free fraction. With a_P =
s_P (vbar_P/4) n_P / N0, a_B = s_B (vbar_B/4) n_B / N0 and e = k_etch (vbar_B/4) n_B / N0 (all 1/s, local):

    d theta_P/dt = a_P f - e theta_P
    d theta_B/dt = a_B f - k_B theta_B

Gas, quasi-steady for both species (Ylilammi et al. 2018 for the precursor; the byproduct relaxes on the same
L^2/D_K ~ 10 us scale), closed end, Dirichlet mouth densities n_P0(t) and n_B0 = mouth_fraction n_P0:

    (1/A) d/dx (A D_P dn_P/dx) = (P/A) (vbar_P/4) [ s_P f + r ] n_P
    (1/A) d/dx (A D_B dn_B/dx) = (P/A) (vbar_B/4) [ s_B f + k_etch theta_P ] n_B
                                 - (P/A) [ nu s_P (vbar_P/4) f n_P + k_B N0 theta_B ]

nu is the number of B released per adsorbed P (3 CH4 per TMA over a full Al2O3 cycle, Puurunen 2005; 4 HCl per TiCl4
over a TiN cycle, Ritala et al. 1995, of which about 0.5 leave in the TiCl4 half-step, Tiznado and Zaera 2006).
The coverage pair is advanced per substep with the gas frozen by the exact solution of the linear 2 x 2 system
(matrix exponential in closed form), so with s_B = k_etch = 0 the precursor update is the single-species
theta <- 1 - (1 - theta) exp(-a_P dt) to round-off and the gas equation for P is unchanged: the single-species solver
is reproduced. k_B -> infinity empties the blocked sites instantly and reproduces it too.

Byproduct bookkeeping (gas side, per substep, with the frozen coverages the gas was solved with): produced (nu x P
uptake) + desorbed = escaped through the mouth + re-adsorbed + consumed by etching. The quasi-steady finite-volume
solve conserves this identity to round-off; `TwoSpeciesResult.balance` holds the time integrals.

Not modelled: the etch product (TiCl_x) is not tracked as a third species; B produced on the flat wafer surface
outside the feature is swept away by the chamber flow (that is the mouth_fraction = 0 default; a nonzero value stands
for byproduct accumulating in the chamber); gas-phase collisions between species (Knudsen number >> 1).
"""
from dataclasses import dataclass, field
from typing import Optional
import numpy as np

from .solver import Process, _solve_gas, _initial_opening, _area_terms
from .transport import mean_speed, knudsen_diffusivity, gas_diffusivity, bosanquet
from .solver import _geometry_factor
from .kinetics import number_density
from .recess import time_to_coverage


@dataclass
class Byproduct:
    """The second species.

    name: label
    mass_amu: molecular mass (sets vbar_B and D_K,B)
    nu: molecules of B released per adsorbed P (stoichiometry of the half-reaction)
    sticking: s_B, probability that a B collision with a free site blocks it
    k_des: k_B, first-order desorption rate of a blocked site (1/s); 0 = irreversible blocking
    k_etch: probability that a B collision with a precursor-covered site removes it (0 = no etching)
    mouth_fraction: byproduct density at the feature mouth as a fraction of the precursor mouth density
        (0: the chamber sweeps the byproduct away; > 0: byproduct accumulated in the chamber)
    sigma: hard-sphere diameter for the transition-regime (Bosanquet) correction when Process.p_total is set
    """
    name: str = "CH4"
    mass_amu: float = 16.04
    nu: float = 3.0
    sticking: float = 0.0
    k_des: float = 0.0
    k_etch: float = 0.0
    mouth_fraction: float = 0.0
    sigma: float = 3.8e-10

    @property
    def active(self):
        """Whether B acts on the surface at all; if not, the precursor sees the single-species model."""
        return self.sticking > 0.0 or self.k_etch > 0.0


# Standard byproducts. Masses: CH4 16.04, HCl 36.46 g/mol.
METHANE = Byproduct("CH4", 16.04, nu=3.0)
HCL_TICL4_HALFSTEP = Byproduct("HCl", 36.46, nu=0.5)    # Cl:Ti 4 -> 3.5 after the TiCl4 step (Tiznado and Zaera 2006)
HCL_TICL4_CYCLE = Byproduct("HCl", 36.46, nu=4.0)       # whole-cycle stoichiometry 6 TiCl4 + 8 NH3 -> 6 TiN + 24 HCl + N2


@dataclass
class TwoSpeciesResult:
    x: np.ndarray
    theta_P: np.ndarray            # precursor-derived coverage at the end
    theta_B: np.ndarray            # blocked-site coverage at the end
    n_P: np.ndarray                # precursor gas density at the end (1/m3)
    n_B: np.ndarray                # byproduct gas density at the end (1/m3)
    t: np.ndarray                  # history times
    hist_P: np.ndarray             # (ntime, ndepth) theta_P at the recorded relative depths
    hist_B: np.ndarray             # (ntime, ndepth) theta_B
    balance: dict                  # byproduct bookkeeping, molecules per unit mouth area (see module docstring)
    thickness: np.ndarray          # cumulative film (m), gpc * theta_P per cycle
    process: Process = field(repr=False)
    byproduct: Byproduct = field(repr=False)

    @property
    def theta_free(self):
        return 1.0 - self.theta_P - self.theta_B

    def time_to_bottom_coverage(self, level=0.9):
        """First time the bottom precursor coverage reaches `level` (nan if never)."""
        return time_to_coverage(self.t, self.hist_P[:, -1], level)

    def time_to_mouth_coverage(self, level=0.9):
        return time_to_coverage(self.t, self.hist_P[:, 0], level)


# ---------------------------------------------------------------------------------------------------------------
# exact advance of the linear coverage pair at frozen gas
# ---------------------------------------------------------------------------------------------------------------
def _e1(z):
    """(exp(z) - 1)/z, 1 at z = 0."""
    z = np.asarray(z, dtype=float)
    small = np.abs(z) < 1e-12
    zs = np.where(small, 1.0, z)
    return np.where(small, 1.0, np.expm1(zs) / zs)


def _e1_prime(z):
    """d/dz [(exp(z) - 1)/z] = (exp(z) (z - 1) + 1)/z^2, 1/2 at z = 0; series below |z| = 1e-3."""
    z = np.asarray(z, dtype=float)
    small = np.abs(z) < 1e-3
    zs = np.where(small, 1.0, z)
    direct = (np.exp(zs) * (zs - 1.0) + 1.0) / (zs * zs)
    series = 0.5 + z / 3.0 + z * z / 8.0 + z ** 3 / 30.0
    return np.where(small, series, direct)


def advance_coverages(theta_P, theta_B, a_P, a_B, e, k_B, dt):
    """Exact solution over dt of  d/dt [thP, thB] = M [thP, thB] + [a_P, a_B]  with
    M = [[-(a_P + e), -a_P], [-a_B, -(a_B + k_B)]]  (all rates frozen), element-wise over arrays.

    theta(dt) = theta0 + phi (M theta0 + b),  phi = int_0^dt exp(M s) ds = psi(l1) I + (M - l1 I) psi[l1, l2],
    psi(l) = (exp(l dt) - 1)/l and psi[l1, l2] its divided difference (Cayley-Hamilton for a 2 x 2 matrix; the
    eigenvalues l1 >= l2 are real because the off-diagonal product a_P a_B >= 0). Singular M (k_B = e = 0, or a cell
    without gas) is handled: psi(0) = dt and the divided difference goes to its derivative."""
    thP = np.asarray(theta_P, dtype=float); thB = np.asarray(theta_B, dtype=float)
    a_P = np.asarray(a_P, dtype=float); a_B = np.asarray(a_B, dtype=float); e = np.asarray(e, dtype=float)
    m11 = -(a_P + e); m12 = -a_P; m21 = -a_B; m22 = -(a_B + k_B)
    tr = m11 + m22
    det = a_P * k_B + e * a_B + e * k_B                 # m11 m22 - m12 m21 without cancellation (all terms >= 0)
    disc = np.sqrt(np.maximum(tr * tr - 4.0 * det, 0.0))
    l2 = 0.5 * (tr - disc)                              # the eigenvalue of larger magnitude, accurate
    l1 = np.where(l2 != 0.0, det / np.where(l2 != 0.0, l2, 1.0), 0.0)   # the smaller one from det/l2, no cancellation
    psi1 = dt * _e1(l1 * dt); psi2 = dt * _e1(l2 * dt)
    dl = (l1 - l2) * dt
    near = np.abs(dl) < 1e-6
    dls = np.where(near, 1.0, l1 - l2)
    dd_direct = (psi1 - psi2) / dls
    dd_deriv = dt * dt * _e1_prime(0.5 * (l1 + l2) * dt)
    dd = np.where(near, dd_deriv, dd_direct)
    r1 = m11 * thP + m12 * thB + a_P
    r2 = m21 * thP + m22 * thB + a_B
    phi11 = psi1 + (m11 - l1) * dd
    phi12 = m12 * dd
    phi21 = m21 * dd
    phi22 = psi1 + (m22 - l1) * dd
    thP_new = thP + phi11 * r1 + phi12 * r2
    thB_new = thB + phi21 * r1 + phi22 * r2
    # round-off guard only (the exact solution stays in the simplex)
    thP_new = np.clip(thP_new, 0.0, 1.0)
    thB_new = np.clip(thB_new, 0.0, 1.0 - thP_new)
    return thP_new, thB_new


# ---------------------------------------------------------------------------------------------------------------
# the half-cycle
# ---------------------------------------------------------------------------------------------------------------
def _diffusivities(g, opening, proc, mass, sigma):
    N = len(opening)
    D = knudsen_diffusivity(g.hydraulic_diameter(opening), proc.T, mass, proc.knudsen_factor * _geometry_factor(g))
    if proc.p_total is not None:
        D = bosanquet(D, gas_diffusivity(proc.T, proc.p_total, mass, sigma))
    D_face = np.empty(N + 1)
    D_face[1:-1] = 2.0 * D[:-1] * D[1:] / (D[:-1] + D[1:])
    D_face[0] = D[0]; D_face[-1] = D[-1]
    return D_face


def _half_cycle(proc: Process, bp: Byproduct, t_series, n0_series, substeps_per_sample, record_depths, opening,
                theta_P, theta_B, balance):
    """Advance one half-cycle on the given opening from the given coverages (modified in place); appends the
    byproduct bookkeeping to `balance`. Returns (theta_P, theta_B, n_P, n_B, hist_t, hist_P, hist_B)."""
    g = proc.geometry; pr = proc.precursor
    N = proc.n_cells
    L = g.length
    x = (np.arange(N) + 0.5) * L / N
    dx = x[1] - x[0]
    vbar_P = mean_speed(proc.T, pr.mass_amu)
    vbar_B = mean_speed(proc.T, bp.mass_amu)
    s, r, N0 = pr.sticking, pr.recombination, pr.site_density
    sB, kB, ke, nu = bp.sticking, bp.k_des, bp.k_etch, bp.nu
    t_series = np.asarray(t_series, float); n0_series = np.asarray(n0_series, float)
    idx = [min(N - 1, int(round(d * (N - 1)))) for d in record_depths]
    D_face_P = _diffusivities(g, opening, proc, pr.mass_amu, pr.sigma)
    D_face_B = _diffusivities(g, opening, proc, bp.mass_amu, bp.sigma)
    sv = g.surface_to_volume(opening)
    A_face, A_cell = _area_terms(g, opening)
    A_vol = g.cross_section(opening) if A_cell is None else A_cell          # cell cross-section for the bookkeeping
    A_mouth = A_vol[0] if A_face is None else A_face[0]
    n_P = np.zeros(N); n_B = np.zeros(N)
    hist_t = []; hP = []; hB = []
    for i in range(1, len(t_series)):
        dt = (t_series[i] - t_series[i - 1]) / substeps_per_sample
        for k in range(substeps_per_sample):
            tt = t_series[i - 1] + (k + 0.5) * dt
            n0 = float(np.interp(tt, t_series, n0_series))
            nB0 = bp.mouth_fraction * n0
            free = 1.0 - theta_P - theta_B
            kP = sv * (vbar_P / 4.0) * (s * free + r)
            n_P = _solve_gas(x, D_face_P, kP, n0, A_face, A_cell) if n0 > 0 else np.zeros(N)
            uptake = s * (vbar_P / 4.0) * n_P * free                  # P adsorption flux per unit wall area
            if bp.active or nu > 0.0 or nB0 > 0.0:
                src = sv * (nu * uptake + kB * N0 * theta_B)           # molecules/(m3 s)
                kBsink = sv * (vbar_B / 4.0) * (sB * free + ke * theta_P)
                if nB0 > 0.0 or np.any(src > 0.0):
                    n_B = _solve_gas(x, D_face_B, kBsink, nB0, A_face, A_cell, source=src)
                else:
                    n_B = np.zeros(N)
            else:
                src = np.zeros(N); n_B = np.zeros(N)
            a_P = s * (vbar_P / 4.0) * n_P / N0
            a_B = sB * (vbar_B / 4.0) * n_B / N0
            e = ke * (vbar_B / 4.0) * n_B / N0
            # bookkeeping with the frozen coverages the gas was solved with (per unit mouth area, this substep)
            if nu > 0.0 or bp.active or nB0 > 0.0:
                vol = A_vol * dx
                balance["produced"] += float(np.sum(sv * nu * uptake * vol)) * dt / A_mouth
                balance["desorbed"] += float(np.sum(sv * kB * N0 * theta_B * vol)) * dt / A_mouth
                balance["adsorbed"] += float(np.sum(sv * sB * (vbar_B / 4.0) * free * n_B * vol)) * dt / A_mouth
                balance["etched"] += float(np.sum(sv * ke * (vbar_B / 4.0) * theta_P * n_B * vol)) * dt / A_mouth
                balance["escaped"] += float(D_face_B[0] * 2.0 * (n_B[0] - nB0) / dx) * dt
                balance["P_uptake"] += float(np.sum(sv * uptake * vol)) * dt / A_mouth
            if bp.active:
                theta_P, theta_B = advance_coverages(theta_P, theta_B, a_P, a_B, e, kB, dt)
            else:
                theta_P = 1.0 - (1.0 - theta_P) * np.exp(-s * (vbar_P / 4.0) * n_P * dt / N0)
        hist_t.append(t_series[i]); hP.append(theta_P[idx].copy()); hB.append(theta_B[idx].copy())
    return theta_P, theta_B, n_P, n_B, np.array(hist_t), np.array(hP), np.array(hB)


def _new_balance():
    return {"produced": 0.0, "desorbed": 0.0, "adsorbed": 0.0, "etched": 0.0, "escaped": 0.0, "P_uptake": 0.0}


def balance_closure(balance):
    """Relative closure of (produced + desorbed) - (escaped + adsorbed + etched); 0 when nothing was produced."""
    inn = balance["produced"] + balance["desorbed"]
    out = balance["escaped"] + balance["adsorbed"] + balance["etched"]
    return 0.0 if inn == 0.0 else (inn - out) / inn


def simulate_two_species(proc: Process, bp: Byproduct = Byproduct(), record_depths=(0.0, 1.0),
                         n_samples: Optional[int] = None) -> TwoSpeciesResult:
    """Constant mouth density proc.p_partial for proc.t_dose, proc.cycles cycles (as aldconf.simulate: each cycle
    starts on a reset surface, the film grows by gpc * theta_P and the opening narrows if proc.narrowing).
    The coverage history is recorded every t_dose / n_samples (default: every substep)."""
    g = proc.geometry
    N = proc.n_cells
    x = (np.arange(N) + 0.5) * g.length / N
    opening = _initial_opening(g, x)
    thickness = np.zeros(N)
    n0 = number_density(proc.p_partial, proc.T)
    n_samples = proc.n_substeps if n_samples is None else int(n_samples)
    sub = max(1, proc.n_substeps // n_samples)
    t_series = np.linspace(0.0, proc.t_dose, n_samples + 1)
    n0_series = np.full(n_samples + 1, n0)
    balance = _new_balance()
    for _ in range(proc.cycles):
        theta_P = np.zeros(N); theta_B = np.zeros(N)
        theta_P, theta_B, n_P, n_B, ht, hP, hB = _half_cycle(proc, bp, t_series, n0_series, sub, record_depths,
                                                             opening, theta_P, theta_B, balance)
        growth = proc.precursor.gpc * theta_P
        thickness += growth
        if proc.narrowing:
            opening = np.maximum(opening - g.walls() * growth, 0.05 * _initial_opening(g, x))
    balance["closure"] = balance_closure(balance)
    return TwoSpeciesResult(x=x, theta_P=theta_P, theta_B=theta_B, n_P=n_P, n_B=n_B, t=ht, hist_P=hP, hist_B=hB,
                            balance=balance, thickness=thickness, process=proc, byproduct=bp)


def simulate_two_species_driven(proc: Process, bp: Byproduct, t_series, n0_series, substeps_per_sample: int = 5,
                                record_depths=(0.0, 1.0)) -> TwoSpeciesResult:
    """One half-cycle with a time-varying precursor density n0(t) [1/m3] at the mouth (the reactor's wafer-surface
    density history), as aldconf.simulate_driven; the byproduct mouth density follows it as mouth_fraction * n0(t)."""
    g = proc.geometry
    N = proc.n_cells
    x = (np.arange(N) + 0.5) * g.length / N
    opening = _initial_opening(g, x)
    balance = _new_balance()
    theta_P = np.zeros(N); theta_B = np.zeros(N)
    theta_P, theta_B, n_P, n_B, ht, hP, hB = _half_cycle(proc, bp, t_series, n0_series, substeps_per_sample,
                                                         record_depths, opening, theta_P, theta_B, balance)
    balance["closure"] = balance_closure(balance)
    return TwoSpeciesResult(x=x, theta_P=theta_P, theta_B=theta_B, n_P=n_P, n_B=n_B, t=ht, hist_P=hP, hist_B=hB,
                            balance=balance, thickness=proc.precursor.gpc * theta_P, process=proc, byproduct=bp)


def knudsen_number(p_total, T, d, sigma=3.8e-10):
    """Mean free path over the feature opening (hard-sphere carrier of diameter sigma; N2 3.8 A)."""
    from .transport import mean_free_path
    return mean_free_path(T, p_total, sigma) / d
