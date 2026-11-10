"""Extract the sticking coefficient from a measured thickness-vs-depth profile."""
import numpy as np
from scipy.optimize import minimize_scalar
from .solver import Process, simulate
from .kinetics import Precursor


def fit_sticking(x_meas, thick_meas, proc: Process, s_bounds=(1e-6, 1.0)):
    """Least-squares fit of log10(s) so the model profile matches `thick_meas(x_meas)`.

    The measured profile is normalised by its entrance value; so is the model, which removes
    the growth-per-cycle scale and leaves the shape, which is what carries s (Arts et al. 2019).
    Returns (s_fit, rms_residual)."""
    x_meas = np.asarray(x_meas, float)
    y_meas = np.asarray(thick_meas, float)
    y_meas = y_meas / np.mean(y_meas[: max(1, len(y_meas) // 50)])

    def cost(log_s):
        pr = Precursor(**{**proc.precursor.__dict__, "sticking": 10.0 ** log_s})
        p = Process(**{**proc.__dict__, "precursor": pr})
        res = simulate(p)
        y_mod = res.thickness / np.mean(res.thickness[: max(1, len(res.x) // 50)])
        y_i = np.interp(x_meas, res.x, y_mod)
        return float(np.mean((y_i - y_meas) ** 2))

    r = minimize_scalar(cost, bounds=(np.log10(s_bounds[0]), np.log10(s_bounds[1])), method="bounded",
                        options={"xatol": 1e-3})
    return 10.0 ** r.x, float(np.sqrt(r.fun))
