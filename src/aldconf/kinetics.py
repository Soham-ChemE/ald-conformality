"""Surface kinetics: irreversible Langmuir adsorption, optional radical recombination."""
from dataclasses import dataclass
import numpy as np
from .transport import mean_speed, K_B


@dataclass
class Precursor:
    """A reactant pulse.

    name: label
    mass_amu: molecular mass
    sticking: sticking coefficient on an open site, s (irreversible Langmuir)
    recombination: wall loss probability per collision regardless of coverage, r (0 for thermal ALD;
        ~1e-5 to 1e-1 for plasma radicals, Arts et al. 2019/2021)
    site_density: adsorption sites per m^2 that one saturated half-cycle fills, N0
    gpc: film thickness added per fully saturated cycle (m)
    sigma: hard-sphere diameter for the transition-regime correction (m)
    """
    name: str = "TMA"
    mass_amu: float = 72.09
    sticking: float = 1e-3
    recombination: float = 0.0
    site_density: float = 4.0e18
    gpc: float = 0.10e-9
    sigma: float = 6.0e-10


def impingement_flux(n, T, mass_amu):
    """Wall impingement flux n vbar / 4 (m^-2 s^-1)."""
    return n * mean_speed(T, mass_amu) / 4.0


def number_density(p_partial, T):
    return p_partial / (K_B * T)
