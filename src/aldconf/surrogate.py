"""Supply-limit pulse rule and the numpy-only Gaussian-process surrogate of the patterned-wafer DOE.

The surrogate is read from ``figures/doe_gp_python.json`` (written by ``analysis/doe_gp.py``), which holds, for each
response, an ordinary-least-squares linear mean on the coded factors plus a Gaussian process on its residuals with an
ARD squared-exponential kernel and a white-noise term (the scikit-learn equivalent of MATLAB fitrgp with
'BasisFunction','linear' and 'KernelFunction','ardsquaredexponential'). Prediction here needs numpy only, so the same
numbers can be reproduced in a browser from the JSON.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

KB = 1.380649e-23            # J/K
P_STP, T_STP = 101325.0, 273.15
SCCM_TO_M3S = 1.6667e-8      # m3/s per sccm at STP
R_SHOWERHEAD = 0.15          # m, 300 mm wafer-sized showerhead

DEFAULT_JSON = Path(__file__).resolve().parents[2] / "figures" / "doe_gp_python.json"

FACTOR_LO = dict(x=0.005, Q=150.0, p=100.0, g=10.0)
FACTOR_HI = dict(x=0.03, Q=600.0, p=1000.0, g=25.0)
CODING = {
    "x": "x = 0.005*(0.03/0.005)^((c_x+1)/2)   (TMA fraction, log-coded)",
    "Q": "Q_sccm = 150*4^((c_Q+1)/2)           (carrier flow, log-coded)",
    "p": "p_Pa = 100*10^((c_p+1)/2)            (pressure, log-coded)",
    "g": "g_mm = 10 + 7.5*(c_g+1)              (gap, linear)",
    "inverse": "c = 2*log(v/lo)/log(hi/lo) - 1 for x, Q, p; c_g = (g_mm - 10)/7.5 - 1",
}


def code_factors(x, Q_sccm, p_Pa, g_mm):
    """Real factors -> coded [-1, 1] factors (log coding for x, Q, p; linear for g). Extrapolation allowed."""
    cx = 2 * np.log(np.asarray(x, float) / FACTOR_LO["x"]) / np.log(FACTOR_HI["x"] / FACTOR_LO["x"]) - 1
    cq = 2 * np.log(np.asarray(Q_sccm, float) / FACTOR_LO["Q"]) / np.log(FACTOR_HI["Q"] / FACTOR_LO["Q"]) - 1
    cp = 2 * np.log(np.asarray(p_Pa, float) / FACTOR_LO["p"]) / np.log(FACTOR_HI["p"] / FACTOR_LO["p"]) - 1
    cg = (np.asarray(g_mm, float) - FACTOR_LO["g"]) / 7.5 - 1
    return np.stack(np.broadcast_arrays(cx, cq, cp, cg), axis=-1)


def decode_factors(c):
    """Coded factors (..., 4) -> dict of real factors x, Q_sccm, p_Pa, g_mm."""
    c = np.asarray(c, float)
    return dict(
        x=FACTOR_LO["x"] * (FACTOR_HI["x"] / FACTOR_LO["x"]) ** ((c[..., 0] + 1) / 2),
        Q_sccm=FACTOR_LO["Q"] * 4.0 ** ((c[..., 1] + 1) / 2),
        p_Pa=FACTOR_LO["p"] * 10.0 ** ((c[..., 2] + 1) / 2),
        g_mm=FACTOR_LO["g"] + 7.5 * (c[..., 3] + 1),
    )


# ----------------------------------------------------------------------------------------------------------------------
# Supply-limit rule
# ----------------------------------------------------------------------------------------------------------------------
def showerhead_velocity(Q_sccm, p_Pa, T_K, radius_m=R_SHOWERHEAD):
    """Mean gas velocity (m/s) below a circular showerhead of the given radius at p and T."""
    vol = Q_sccm * SCCM_TO_M3S * (P_STP / p_Pa) * (T_K / T_STP)
    return vol / (np.pi * radius_m ** 2)


def rule_pulse(x, Q_sccm, p_Pa, T_K, AR, d_m, hole_density, N0=4e18, factor=1.2):
    """Supply-limit pulse rule: t_pulse = factor*(1 + A_wall/A_wafer)*N0/(U*c0).

    x: precursor mole fraction at the showerhead; Q_sccm: carrier flow; p_Pa, T_K: chamber pressure and temperature;
    AR, d_m: hole aspect ratio and diameter; hole_density: holes per m2 of wafer. A_wall/A_wafer = density*pi*d*AR*d.
    Returns dict(t_supply_s, t_pulse_s, wall_ratio, U_m_s, c0) with c0 in molecules/m3.
    """
    wall_ratio = hole_density * np.pi * d_m * AR * d_m
    U = showerhead_velocity(Q_sccm, p_Pa, T_K)
    c0 = x * p_Pa / (KB * T_K)
    t_supply = (1 + wall_ratio) * N0 / (U * c0)
    return dict(t_supply_s=float(t_supply), t_pulse_s=float(factor * t_supply), wall_ratio=float(wall_ratio),
                U_m_s=float(U), c0=float(c0))


# ----------------------------------------------------------------------------------------------------------------------
# GP surrogate from JSON
# ----------------------------------------------------------------------------------------------------------------------
def _ard_rbf(A, B, const, ls):
    """k(a, b) = const * exp(-0.5 * sum(((a-b)/ls)^2)) for A (n,4), B (m,4)."""
    d = (A[:, None, :] - B[None, :, :]) / ls
    return const * np.exp(-0.5 * np.sum(d * d, axis=-1))


class _Response:
    def __init__(self, spec):
        self.name = spec["name"]
        self.transform = spec.get("transform", "identity")
        self.lin = np.asarray(spec["linear"]["coef"], float)               # intercept + 4
        self.const = float(spec["kernel"]["constant_value"])
        self.ls = np.asarray(spec["kernel"]["length_scale"], float)
        self.noise = float(spec["kernel"]["noise_level"])
        self.X = np.asarray(spec["X_train"], float)
        self.alpha = np.asarray(spec["alpha"], float)
        self.K_inv = np.asarray(spec["K_inv"], float) if "K_inv" in spec else None
        self.y_mean = float(spec["y_normalisation"]["mean"])
        self.y_std = float(spec["y_normalisation"]["std"])

    def predict(self, C, return_sd=False):
        C = np.atleast_2d(np.asarray(C, float))
        mean_lin = self.lin[0] + C @ self.lin[1:]
        Ks = _ard_rbf(C, self.X, self.const, self.ls)
        mu = mean_lin + self.y_std * (Ks @ self.alpha) + self.y_mean
        if not return_sd:
            return mu
        if self.K_inv is None:
            raise ValueError("JSON has no K_inv; standard deviation unavailable")
        var = self.const + self.noise - np.einsum("ij,jk,ik->i", Ks, self.K_inv, Ks)
        sd = self.y_std * np.sqrt(np.maximum(var, 0.0))
        return mu, sd


class Surrogate:
    """Three GP surrogates (log10 t90_worst, edge_lag_pct, utilization) evaluated from the JSON with numpy only."""

    def __init__(self, spec):
        self.spec = spec
        self.responses = {r["name"]: _Response(r) for r in spec["responses"]}
        self.loo_r2 = {r["name"]: r.get("loo_r2") for r in spec["responses"]}

    def predict_coded(self, C):
        """Coded factors (n, 4) -> dict of arrays t90_s, t90_sd_factor, edge_lag_pct, utilization."""
        C = np.atleast_2d(np.asarray(C, float))
        lt, sd = self.responses["log10_t90_worst"].predict(C, return_sd=True)
        return dict(t90_s=10 ** lt, t90_sd_factor=10 ** sd,
                    edge_lag_pct=self.responses["edge_lag_pct"].predict(C),
                    utilization=self.responses["utilization"].predict(C))

    def predict(self, x, Q, p, g):
        """Real factors (TMA fraction, sccm, Pa, mm), scalars -> dict of floats; arrays -> dict of arrays."""
        C = code_factors(x, Q, p, g).reshape(-1, 4)
        out = self.predict_coded(C)
        if np.ndim(x) == 0 and np.ndim(Q) == 0 and np.ndim(p) == 0 and np.ndim(g) == 0:
            return {k: float(v[0]) for k, v in out.items()}
        return out

    __call__ = predict

    def grid(self, n=21):
        """Full-factorial coded grid n^4 in [-1, 1]^4, in the same order as MATLAB ndgrid(:) (first factor fastest)."""
        gl = np.linspace(-1, 1, n)
        A, B, C_, D = np.meshgrid(gl, gl, gl, gl, indexing="ij")
        return np.column_stack([A.ravel(order="F"), B.ravel(order="F"), C_.ravel(order="F"), D.ravel(order="F")])

    def pareto_front(self, max_lag=5.0, n=21):
        """Pareto front over the n^4 grid with predicted edge lag <= max_lag: minimum precursor per saturating pulse
        (1/utilization, relative to the stoichiometric minimum) for each pulse time. Returns (front, all_feasible) as
        structured-array-like dicts of columns t90_s, util, precursor_over_minimum, x, Q_sccm, p_Pa, g_mm, tma_sccm_s."""
        G = self.grid(n)
        pr = self.predict_coded(G)
        feas = pr["edge_lag_pct"] <= max_lag
        real = decode_factors(G[feas])
        t90 = pr["t90_s"][feas]
        util = pr["utilization"][feas]
        cols = dict(t90_s=t90, util=util, precursor_over_minimum=1.0 / util, x=real["x"], Q_sccm=real["Q_sccm"],
                    p_Pa=real["p_Pa"], g_mm=real["g_mm"], tma_sccm_s=real["x"] * real["Q_sccm"] * t90)
        order = np.argsort(t90, kind="stable")
        best = np.inf
        keep = []
        for i in order:
            if cols["precursor_over_minimum"][i] < best - 1e-12:
                best = cols["precursor_over_minimum"][i]
                keep.append(i)
        keep = np.asarray(keep, int)
        front = {k: v[keep] for k, v in cols.items()}
        return front, cols


def load_surrogate(path=DEFAULT_JSON):
    """Load figures/doe_gp_python.json and return a Surrogate whose .predict(x, Q, p, g) gives
    dict(t90_s, t90_sd_factor, edge_lag_pct, utilization)."""
    with open(path) as f:
        spec = json.load(f)
    return Surrogate(spec)
