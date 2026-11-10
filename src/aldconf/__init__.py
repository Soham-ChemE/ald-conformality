"""Feature-scale ALD / PEALD conformality model."""
from .geometry import Slit, Hole
from .transport import mean_speed, knudsen_diffusivity, gas_diffusivity, bosanquet
from .kinetics import Precursor, impingement_flux
from .solver import Process, simulate, simulate_driven, step_coverage, penetration_depth, dose_for_coverage
from .fit import fit_sticking
from .surrogate import load_surrogate, rule_pulse, code_factors, decode_factors

__all__ = [
    "Slit", "Hole", "mean_speed", "knudsen_diffusivity", "gas_diffusivity", "bosanquet",
    "Precursor", "impingement_flux", "Process", "simulate", "simulate_driven", "step_coverage",
    "penetration_depth", "dose_for_coverage", "fit_sticking",
    "load_surrogate", "rule_pulse", "code_factors", "decode_factors",
]
