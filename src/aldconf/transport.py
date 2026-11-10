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

    Exact for a cylinder. For a wide slit the exact coefficient carries a logarithmic
    geometry correction (Ylilammi et al. 2018); `factor` absorbs it and is 1 by default.
    """
    return factor * d_h * mean_speed(T, mass_amu) / 3.0


def mean_free_path(T, p_total, sigma):
    """Hard-sphere mean free path k T / (sqrt(2) pi sigma^2 p)."""
    return K_B * T / (np.sqrt(2.0) * np.pi * sigma ** 2 * p_total)


def gas_diffusivity(T, p_total, mass_amu, sigma):
    """Self/carrier diffusivity D = lambda * vbar / 3 (hard-sphere, kinetic theory)."""
    return mean_free_path(T, p_total, sigma) * mean_speed(T, mass_amu) / 3.0


def bosanquet(D_K, D_gas):
    """Transition-regime interpolation 1/D = 1/D_K + 1/D_gas."""
    return 1.0 / (1.0 / D_K + 1.0 / D_gas)
