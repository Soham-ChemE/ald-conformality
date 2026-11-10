"""scikit-learn port of analysis/doe_gp.m: Gaussian-process surrogates of the 25-run patterned-wafer DOE.

For each response (log10 t90_worst, edge_lag_pct, utilization) an ordinary-least-squares linear model on the four
coded factors is fitted first and a GaussianProcessRegressor with ConstantKernel*RBF(ARD, 4 length scales) +
WhiteKernel is fitted to its residuals. This two-stage fit is the equivalent of MATLAB fitrgp's
'BasisFunction','linear' with 'KernelFunction','ardsquaredexponential' (explicit linear mean plus a zero-mean GP on
the remainder; MATLAB estimates the basis coefficients jointly, here they are the OLS values, which is what the GP
marginal likelihood converges to when the kernel explains little of the linear trend).

Outputs: figures/doe_gp_python_report.txt (Python and MATLAB numbers side by side), figures/doe_gp_python.json
(everything a browser needs to evaluate the three surrogates), figures/doe_pareto.csv and figures/doe_pareto.png.

Run:  PYTHONPATH=src .venv/bin/python analysis/doe_gp.py
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from aldconf.surrogate import CODING, Surrogate, decode_factors, rule_pulse  # noqa: E402

RESPONSES = ("log10_t90_worst", "edge_lag_pct", "utilization")
N_RESTARTS, SEED = 8, 0
# Length-scale bounds on the coded factors. The face-centred design has three levels one coded unit apart, so
# variation on scales shorter than ~0.3 is not identifiable from it; without the lower bound the edge-lag GP (a noisy
# response) takes a length scale of 0.04 in c_p, i.e. it memorises the points, and its leave-one-out R^2 drops. MATLAB's
# fitrgp starts from length scale = std(X) and its optimizer stays in the smooth basin; the bound reproduces that.
LS_BOUNDS = (0.3, 30.0)
LS_BOUNDS_LOOSE = (1e-2, 1e2)
warnings.filterwarnings("ignore", category=ConvergenceWarning)


# ----------------------------------------------------------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------------------------------------------------------
def fit_response(X, y, n_restarts=N_RESTARTS, seed=SEED, ls_bounds=LS_BOUNDS):
    """OLS linear mean + GP(ARD RBF + white noise) on the residuals. Returns (lin_coef[5], gp)."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    A = np.column_stack([np.ones(len(X)), X])
    lin, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ lin
    kernel = (ConstantKernel(1.0, (1e-4, 1e3)) * RBF(length_scale=[1.0] * 4, length_scale_bounds=ls_bounds)
              + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-10, 1e1)))
    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=n_restarts, random_state=seed)
    gp.fit(X, resid)
    return lin, gp


def predict(lin, gp, C, return_std=False):
    C = np.atleast_2d(np.asarray(C, float))
    mean_lin = lin[0] + C @ lin[1:]
    if return_std:
        mu, sd = gp.predict(C, return_std=True)
        return mean_lin + mu, sd
    return mean_lin + gp.predict(C)


def loo(X, y, ls_bounds=LS_BOUNDS):
    """True leave-one-out: 25 refits, each with its own OLS mean and GP hyperparameters."""
    n = len(y)
    yhat = np.zeros(n)
    for i in range(n):
        idx = np.ones(n, bool)
        idx[i] = False
        lin, gp = fit_response(X[idx], y[idx], ls_bounds=ls_bounds)
        yhat[i] = predict(lin, gp, X[i:i + 1])[0]
    r2 = 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    rmse = np.sqrt(np.mean((y - yhat) ** 2))
    return r2, rmse


def fit_all(T):
    X = T[["c_x", "c_Q", "c_p", "c_g"]].to_numpy(float)
    ys = {"log10_t90_worst": np.log10(T.t90_worst.to_numpy(float)),
          "edge_lag_pct": T.edge_lag_pct.to_numpy(float),
          "utilization": T.utilization.to_numpy(float)}
    return X, ys, {k: fit_response(X, v) for k, v in ys.items()}


# ----------------------------------------------------------------------------------------------------------------------
# JSON export
# ----------------------------------------------------------------------------------------------------------------------
def response_spec(name, X, lin, gp, loo_r2, loo_rmse):
    k = gp.kernel_
    L = gp.L_
    K_inv = np.linalg.solve(L.T, np.linalg.solve(L, np.eye(len(X))))
    return {
        "name": name,
        "transform": "log10" if name.startswith("log10") else "identity",
        "linear": {"basis": "1, c_x, c_Q, c_p, c_g", "coef": lin.tolist()},
        "kernel": {"type": "ConstantKernel*RBF(ARD) + WhiteKernel on the OLS residuals",
                   "constant_value": float(k.k1.k1.constant_value),
                   "length_scale": np.asarray(k.k1.k2.length_scale, float).tolist(),
                   "noise_level": float(k.k2.noise_level)},
        "X_train": X.tolist(),
        "alpha": np.asarray(gp.alpha_, float).ravel().tolist(),
        "K_inv": K_inv.tolist(),
        "y_normalisation": {"mean": float(gp._y_train_mean), "std": float(gp._y_train_std),
                            "note": "GP fitted to residuals r = y - linear; sklearn normalize_y: r_n = (r - mean)/std"},
        "loo_r2": float(loo_r2), "loo_rmse": float(loo_rmse),
        "prediction": "y = linear(c) + mean + std * sum_j alpha_j * const*exp(-0.5*sum_k ((c_k - X_jk)/ls_k)^2); "
                      "var = std^2 * (const + noise - k*^T K_inv k*)",
    }


def verify_json(spec, fits, rng):
    """Re-implement the prediction from the JSON alone (aldconf.surrogate, numpy only) and compare with sklearn."""
    sur = Surrogate(json.loads(json.dumps(spec)))
    C = rng.uniform(-1, 1, size=(10, 4))
    worst = 0.0
    for name, (lin, gp) in fits.items():
        mu_sk, sd_sk = predict(lin, gp, C, return_std=True)
        mu_js, sd_js = sur.responses[name].predict(C, return_sd=True)
        worst = max(worst, np.abs(mu_sk - mu_js).max(), np.abs(sd_sk - sd_js).max())
    return worst


# ----------------------------------------------------------------------------------------------------------------------
# MATLAB reference numbers (parsed from figures/doe_gp_report.txt so nothing is retyped)
# ----------------------------------------------------------------------------------------------------------------------
def parse_matlab_report(path):
    out = {"loo": {}, "points": {}, "mc": None}
    if not path.exists():
        return out
    txt = path.read_text()
    for m in re.finditer(r"^(\w+): LOO R\^2 = ([\d.]+), LOO RMSE = ([\d.]+)", txt, re.M):
        out["loo"][m.group(1)] = (float(m.group(2)), float(m.group(3)))
    pat = (r"^(Centre point|Fastest with edge lag <= 5 %|Same precursor supply[^:]*): x=([\d.]+) %, Q=(\d+) sccm, "
           r"p=(\d+) Pa, g=([\d.]+) mm -> t90 ([\d.]+) s(?: \(sd factor ([\d.]+)\))?, edge lag (-?[\d.]+) %, "
           r"utilization (\d+) %")
    keys = {"Centre point": "centre", "Fastest with edge lag <= 5 %": "fastest"}
    for m in re.finditer(pat, txt, re.M):
        key = keys.get(m.group(1), "same_supply")
        out["points"][key] = dict(x_pct=float(m.group(2)), Q=float(m.group(3)), p=float(m.group(4)),
                                  g=float(m.group(5)), t90=float(m.group(6)),
                                  sd=float(m.group(7)) if m.group(7) else None, lag=float(m.group(8)),
                                  util=float(m.group(9)))
    m = re.search(r"t_pulse = ([\d.]+) s, 90 % band ([\d.]+) to ([\d.]+) s", txt)
    if m:
        out["mc"] = tuple(float(g) for g in m.groups())
    return out


# ----------------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------------
def fmt_point(label, real, t90, sd, lag, util, matlab=None):
    s = (f"{label}: x={100 * real['x']:.2f} %, Q={real['Q_sccm']:.0f} sccm, p={real['p_Pa']:.0f} Pa, "
         f"g={real['g_mm']:.1f} mm -> t90 {t90:.2f} s")
    if sd is not None:
        s += f" (sd factor {sd:.2f})"
    s += f", edge lag {lag:.1f} %, utilization {100 * util:.0f} %"
    if matlab:
        s += (f"\n    MATLAB: x={matlab['x_pct']:.2f} %, Q={matlab['Q']:.0f} sccm, p={matlab['p']:.0f} Pa, "
              f"g={matlab['g']:.1f} mm -> t90 {matlab['t90']:.2f} s, edge lag {matlab['lag']:.1f} %, "
              f"utilization {matlab['util']:.0f} %   (t90 difference {100 * (t90 / matlab['t90'] - 1):+.1f} %)")
    return s


def main():
    import sklearn
    T = pd.read_csv(ROOT / "figures/doe_results.csv")
    X, ys, fits = fit_all(T)
    ml = parse_matlab_report(ROOT / "figures/doe_gp_report.txt")
    lines = [f"Python {sys.version.split()[0]}, scikit-learn {sklearn.__version__}, numpy {np.__version__}: "
             f"GP surrogates on {len(T)} DOE runs (OLS linear mean + GP on residuals, ARD RBF + white noise, "
             f"n_restarts_optimizer={N_RESTARTS}, random_state={SEED})",
             "Leave-one-out (25 refits per response), Python vs MATLAB fitrgp (figures/doe_gp_report.txt):",
             f"  {'response':18s} {'R2 py':>8s} {'R2 MATLAB':>10s} {'RMSE py':>9s} {'RMSE MATLAB':>12s}"]
    loo_stats = {}
    for name in RESPONSES:
        r2, rmse = loo(X, ys[name])
        loo_stats[name] = (r2, rmse)
        mr2, mrmse = ml["loo"].get(name, (np.nan, np.nan))
        lines.append(f"  {name:18s} {r2:8.3f} {mr2:10.3f} {rmse:9.4f} {mrmse:12.4f}")
    r2_loose, rmse_loose = loo(X, ys["edge_lag_pct"], ls_bounds=LS_BOUNDS_LOOSE)
    lines.append(f"  Length-scale bounds on the coded factors are {LS_BOUNDS} (three-level design, levels 1.0 apart). "
                 f"The t90 and utilization LOO statistics do not depend on them (t90's c_p length scale sits at the "
                 f"lower bound, 0.18 when unbounded, with the same LOO R^2 and RMSE to the digits shown); edge lag does: "
                 f"with bounds {LS_BOUNDS_LOOSE} its c_p length scale collapses to 0.04 and the LOO R^2 falls to "
                 f"{r2_loose:.3f} (RMSE {rmse_loose:.4f}).")
    lines.append("")
    lines.append("Fitted hyperparameters (coded factors c_x, c_Q, c_p, c_g):")
    for name in RESPONSES:
        lin, gp = fits[name]
        k = gp.kernel_
        lines.append(f"  {name}: linear = {np.array2string(lin, precision=4)}, const = {k.k1.k1.constant_value:.4g}, "
                     f"length scales = {np.array2string(np.asarray(k.k1.k2.length_scale), precision=3)}, "
                     f"noise = {k.k2.noise_level:.3g}, LML = {gp.log_marginal_likelihood_value_:.2f}")

    # grid search, 21^4 over the coded cube
    sur_tmp = Surrogate({"responses": [response_spec(n, X, *fits[n], *loo_stats[n]) for n in RESPONSES]})
    G = sur_tmp.grid(21)
    lin_t, gp_t = fits["log10_t90_worst"]
    t_pred, t_sd = predict(lin_t, gp_t, G, return_std=True)
    lag_pred = predict(*fits["edge_lag_pct"], G)
    util_pred = predict(*fits["utilization"], G)
    feas = lag_pred <= 5
    i0 = int(np.argmin(np.sum(np.abs(G), axis=1)))
    ib = int(np.argmin(np.where(feas, t_pred, np.inf)))
    xq = decode_factors(G)["x"] * decode_factors(G)["Q_sccm"]
    same = np.abs(np.log(xq / xq[i0])) < 0.05
    cand = np.flatnonzero(same & feas)
    ib2 = int(cand[np.argmin(t_pred[cand])])
    real = decode_factors(G)
    pick = lambda i: {k: float(v[i]) for k, v in real.items()}
    lines.append("")
    lines.append(fmt_point("Centre point", pick(i0), 10 ** t_pred[i0], 10 ** t_sd[i0], lag_pred[i0], util_pred[i0],
                           ml["points"].get("centre")))
    lines.append(fmt_point("Fastest with edge lag <= 5 %", pick(ib), 10 ** t_pred[ib], 10 ** t_sd[ib], lag_pred[ib],
                           util_pred[ib], ml["points"].get("fastest")))
    lines.append(fmt_point("Same precursor supply as centre (x*Q within 5 % in log), edge lag <= 5 %, fastest",
                           pick(ib2), 10 ** t_pred[ib2], None, lag_pred[ib2], util_pred[ib2],
                           ml["points"].get("same_supply")))

    # Monte Carlo on N0 through the supply rule at the centre point
    rng = np.random.default_rng(1)
    N0 = 4e18 * 1.25 ** rng.standard_normal(10000)
    c = pick(i0)
    base = rule_pulse(c["x"], c["Q_sccm"], c["p_Pa"], 573.15, 100, 100e-9, 1.78e13, N0=1.0)
    t_rule = base["t_pulse_s"] * N0
    med, lo, hi = np.median(t_rule), np.percentile(t_rule, 5), np.percentile(t_rule, 95)
    lines.append("")
    lines.append(f"Monte Carlo (N0 lognormal x/÷1.25, 10000 draws, seed 1) on the rule at the centre point "
                 f"(x={100 * c['x']:.3f} %, Q={c['Q_sccm']:.0f} sccm, p={c['p_Pa']:.0f} Pa, T=573.15 K, wall ratio "
                 f"{base['wall_ratio']:.1f}, U={base['U_m_s']:.4f} m/s): t_pulse = {med:.1f} s, 90 % band {lo:.1f} to {hi:.1f} s")
    if ml["mc"]:
        lines.append(f"    MATLAB (rounded constants x=1.22 %, p=300 Pa, U=0.0501 m/s): t_pulse = {ml['mc'][0]:.1f} s, "
                     f"90 % band {ml['mc'][1]:.1f} to {ml['mc'][2]:.1f} s")

    # JSON export and verification
    spec = {
        "description": "GP surrogates of the 25-run patterned-wafer DOE (figures/doe_results.csv), fitted by "
                       "analysis/doe_gp.py with scikit-learn. Each response = OLS linear mean on the coded factors + "
                       "GP(ARD RBF + white noise) on the residuals. Evaluate with aldconf.surrogate or in JavaScript.",
        "software": {"python": sys.version.split()[0], "scikit-learn": sklearn.__version__, "numpy": np.__version__},
        "factors": ["c_x", "c_Q", "c_p", "c_g"],
        "coding": CODING,
        "factor_ranges": {"x": [0.005, 0.03], "Q_sccm": [150, 600], "p_Pa": [100, 1000], "g_mm": [10, 25]},
        "responses": [response_spec(n, X, *fits[n], *loo_stats[n]) for n in RESPONSES],
        "t90_from_log10": "t90_s = 10^log10_t90_worst; sd factor = 10^sd",
    }
    worst = verify_json(spec, fits, np.random.default_rng(0))
    assert worst < 1e-6, worst
    (ROOT / "figures/doe_gp_python.json").write_text(json.dumps(spec, indent=1))
    lines.append("")
    lines.append(f"JSON check: numpy re-implementation from figures/doe_gp_python.json vs sklearn predict at 10 random "
                 f"points, worst |difference| in mean and sd = {worst:.1e}")

    # Pareto front: pulse time vs precursor per saturating pulse (1/utilization), edge lag <= 5 %
    sur = Surrogate(json.loads((ROOT / "figures/doe_gp_python.json").read_text()))
    front, allf = sur.pareto_front(max_lag=5.0, n=21)
    pd.DataFrame(front).to_csv(ROOT / "figures/doe_pareto.csv", index=False, float_format="%.6g")
    t0, u0 = 10 ** t_pred[i0], util_pred[i0]
    t2, u2 = 10 ** t_pred[ib2], util_pred[ib2]
    tb, ub = 10 ** t_pred[ib], util_pred[ib]
    # best precursor use the front offers at or below the centre's pulse time
    at_t0 = front["precursor_over_minimum"][front["t90_s"] <= t0].min()
    j_end = int(np.argmin(front["precursor_over_minimum"]))
    j_fast = 0
    lines.append("")
    lines.append(f"Pareto front ({len(front['t90_s'])} points from {int(feas.sum())} feasible of {len(G)} grid points, "
                 f"figures/doe_pareto.csv, figures/doe_pareto.png):")
    lines.append(f"  At the centre point's pulse time ({t0:.1f} s) the centre uses {1 / u0:.3f}x the stoichiometric "
                 f"minimum of TMA per saturating pulse ({100 * u0:.0f} % utilization); the front reaches "
                 f"{at_t0:.3f}x ({100 / at_t0:.0f} % utilization) at the same or shorter pulse, a saving of "
                 f"{100 * (1 - at_t0 * u0):.0f} % of the precursor per wafer pass at no cost in time.")
    lines.append(f"  The front runs from {front['t90_s'][j_fast]:.2f} s at {front['precursor_over_minimum'][j_fast]:.3f}x "
                 f"(x={100 * front['x'][j_fast]:.2f} %, Q={front['Q_sccm'][j_fast]:.0f} sccm, p={front['p_Pa'][j_fast]:.0f} Pa, "
                 f"g={front['g_mm'][j_fast]:.1f} mm) to {front['t90_s'][j_end]:.2f} s at {front['precursor_over_minimum'][j_end]:.3f}x "
                 f"(x={100 * front['x'][j_end]:.2f} %, Q={front['Q_sccm'][j_end]:.0f} sccm, p={front['p_Pa'][j_end]:.0f} Pa, "
                 f"g={front['g_mm'][j_end]:.1f} mm): the last {100 * (1 - front['precursor_over_minimum'][j_end] / at_t0):.0f} % "
                 f"of precursor saving beyond the centre-time optimum costs {front['t90_s'][j_end] / t0:.1f}x the pulse time.")
    lines.append(f"  Same-supply optimum sits at {t2:.2f} s, {1 / u2:.3f}x; fastest feasible at {tb:.2f} s, {1 / ub:.3f}x.")
    lines.append("  TMA delivered per saturating pulse (x*Q*t90, sccm*s) is in the csv column tma_sccm_s.")

    make_figure(T, front, allf, (t0, u0), (t2, u2), (tb, ub))

    report = "\n".join(lines) + "\n"
    (ROOT / "figures/doe_gp_python_report.txt").write_text(report)
    print(report)


def make_figure(T, front, allf, centre, same, fastest):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    blue, orange, aqua, ink, ink2, grid, cloud = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#d9d8d4", "#ecebe8"
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ticks = [2, 3, 5, 10, 20, 30, 50]
    for a in ax:
        a.set_xscale("log")
        a.xaxis.set_major_locator(FixedLocator(ticks))
        a.xaxis.set_minor_locator(NullLocator())
        a.xaxis.set_major_formatter(ScalarFormatter())
        a.grid(True, which="major", color=grid, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(colors=ink2)
        a.set_xlabel("saturating pulse, worst point on the wafer t90 (s)", color=ink)
    # panel 1: pulse time vs precursor use
    a = ax[0]
    a.scatter(allf["t90_s"], allf["precursor_over_minimum"], s=2, color=cloud, label="feasible grid, edge lag <= 5 %",
              zorder=1, rasterized=True)
    a.scatter(T.t90_worst, 1 / T.utilization, s=34, facecolor="white", edgecolor=blue, lw=1.6, label="25 DOE runs",
              zorder=3)
    a.plot(front["t90_s"], front["precursor_over_minimum"], color=orange, lw=2, label="Pareto front (GP fit)", zorder=4)
    a.scatter(*zip((centre[0], 1 / centre[1]), (same[0], 1 / same[1])), marker="D", s=60, color=aqua, edgecolor="white",
              lw=1, zorder=5)
    a.annotate("centre point", (centre[0], 1 / centre[1]), xytext=(8, 10), textcoords="offset points", color=ink, fontsize=9)
    a.annotate("same supply, fastest", (same[0], 1 / same[1]), xytext=(8, -14), textcoords="offset points", color=ink,
               fontsize=9)
    a.set_ylabel("TMA per saturating pulse / stoichiometric minimum (= 1 / utilization)", color=ink, fontsize=9)
    a.set_title("Pulse time vs precursor use: the front is the cheapest recipe for each pulse time", fontsize=9,
                color=ink, loc="left")
    a.legend(fontsize=8, frameon=False, loc="upper right")
    # panel 2: utilization vs pulse time
    a = ax[1]
    a.scatter(allf["t90_s"], 100 * allf["util"], s=2, color=cloud, zorder=1, rasterized=True)
    a.scatter(T.t90_worst, 100 * T.utilization, s=34, facecolor="white", edgecolor=blue, lw=1.6, zorder=3,
              label="25 DOE runs")
    a.plot(front["t90_s"], 100 * front["util"], color=orange, lw=2, zorder=4, label="Pareto front (GP fit)")
    for (t, u), lab, dy in ((centre, "centre point", -16), (same, "same supply, fastest", 10),
                            (fastest, "fastest, edge lag <= 5 %", 10)):
        a.scatter(t, 100 * u, marker="D", s=60, color=aqua, edgecolor="white", lw=1, zorder=5)
        a.annotate(lab, (t, 100 * u), xytext=(8, dy), textcoords="offset points", color=ink, fontsize=9)
    a.set_ylabel("precursor utilization up to saturation (%)", color=ink)
    a.set_ylim(40, 100)
    a.set_title("Utilization vs pulse time, with the baseline and the same-supply optimum", fontsize=9, color=ink,
                loc="left")
    a.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/doe_pareto.png", dpi=150, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
