"""Gas transport in the free-molecular and transition regimes. SI units."""
import numpy as np

K_B = 1.380649e-23
AMU = 1.66053906660e-27


def mean_speed(T, mass_amu):
    """Mean thermal speed sqrt(8 k T / (pi m))."""
    m = mass_amu * AMU
    return np.sqrt(8.0 * K_B * T / (np.pi * m))


def knudsen_diffusivity(d_h, T, mass_amu, factor=1.0):
    """Knudsen diffusivity D_K = d_h * vbar / 3 for a long channel of hydraulic diameter d_h.

    Exact for a cylinder (Knudsen's long-tube transmission 4 d / 3 L). For a wide slit (d_h = 2 w) the coefficient is
    not exact: molecules that leave a plate at grazing angles travel arbitrarily far before the next wall collision,
    so the free-molecular diffusivity of parallel plates grows logarithmically with the channel length (Clausing;
    Berman 1965; Steckelmacher 1986), see `slit_knudsen_factor`. `factor` multiplies D_K and is 1 by default; the
    geometry's own factor (Slit(knudsen_correction="plates")) is applied by the solver on top of it.
    """
    return factor * d_h * mean_speed(T, mass_amu) / 3.0


def slit_knudsen_factor(gap, length):
    """Length-dependent correction of the slit Knudsen diffusivity: D_slit / (2 w vbar / 3) = (3/8) ln(L / w).

    For free-molecular flow between parallel plates of gap w and length L >> w the transmission probability is
    W = (w / L) [ln(L / w) + c] (Clausing's slit result; tabulated by Berman, J. Appl. Phys. 36, 3356 (1965); reviewed by
    Steckelmacher, Rep. Prog. Phys. 49, 1083 (1986)); a test-particle Monte Carlo of the infinitely wide slit made
    for this repository (paper/wp/wp9b_slit_transport.md) gives W L / w - ln(L / w) = +0.10 for L / w = 10 to 200 and
    W = 0.6844 at L / w = 1, Clausing's value. With Fick's law over the channel, D = W L vbar / 4, so
    D_slit = (w vbar / 4) [ln(L / w) + c] against the hydraulic-diameter value 2 w vbar / 3, a ratio of
    (3/8) [ln(L / w) + c]. The constant is dropped here: with c = 0 the one-dimensional slit model matches Machball's
    ballistic Trench within 5 % at aspect ratios 20 to 100 (t90 of the bottom coverage at 0.3 Pa, s = 5e-3), whereas
    the full c = 0.10 puts it 6 % fast at AR 50, the continuum-versus-ballistic offset also seen for holes. The
    hydraulic-diameter value is recovered at L / w = e^(8/3) = 14.4; below that the factor is less than 1 and the
    slit is reaction-limited anyway. Returns 1 for L <= w.
    """
    ratio = np.asarray(length, dtype=float) / np.asarray(gap, dtype=float)
    return np.where(ratio > 1.0, 0.375 * np.log(np.maximum(ratio, 1.0)), 1.0)


def mean_free_path(T, p_total, sigma):
    """Hard-sphere mean free path k T / (sqrt(2) pi sigma^2 p)."""
    return K_B * T / (np.sqrt(2.0) * np.pi * sigma ** 2 * p_total)


def gas_diffusivity(T, p_total, mass_amu, sigma):
    """Self/carrier diffusivity D = lambda * vbar / 3 (hard-sphere, kinetic theory)."""
    return mean_free_path(T, p_total, sigma) * mean_speed(T, mass_amu) / 3.0


def bosanquet(D_K, D_gas):
    """Transition-regime interpolation 1/D = 1/D_K + 1/D_gas."""
    return 1.0 / (1.0 / D_K + 1.0 / D_gas)
