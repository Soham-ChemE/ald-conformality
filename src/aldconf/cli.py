"""Command line: supply-limit rule, GP surrogate prediction, Pareto front.

    python -m aldconf rule --x 0.01 --Q 300 --p 300 --T 573 --AR 100 --d 100e-9 --density 1.78e13
    python -m aldconf surrogate --x 0.012 --Q 300 --p 316 --g 17.5
    python -m aldconf pareto --max-lag 5
"""
from __future__ import annotations

import argparse
import sys

from .surrogate import DEFAULT_JSON, load_surrogate, rule_pulse


def build_parser():
    p = argparse.ArgumentParser(prog="python -m aldconf", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("rule", help="supply-limit pulse rule t = factor*(1 + A_wall/A_wafer)*N0/(U c0)")
    r.add_argument("--x", type=float, required=True, help="precursor mole fraction at the showerhead")
    r.add_argument("--Q", type=float, required=True, help="carrier flow (sccm)")
    r.add_argument("--p", type=float, required=True, help="pressure (Pa)")
    r.add_argument("--T", type=float, default=573.15, help="temperature (K)")
    r.add_argument("--AR", type=float, required=True, help="hole aspect ratio")
    r.add_argument("--d", type=float, required=True, help="hole diameter (m)")
    r.add_argument("--density", type=float, required=True, help="holes per m2 of wafer")
    r.add_argument("--N0", type=float, default=4e18, help="site density (/m2)")
    r.add_argument("--factor", type=float, default=1.2, help="rule constant")

    s = sub.add_parser("surrogate", help="GP surrogate prediction from figures/doe_gp_python.json")
    s.add_argument("--x", type=float, required=True, help="TMA fraction (0.005-0.03 fitted range)")
    s.add_argument("--Q", type=float, required=True, help="carrier flow (sccm, 150-600)")
    s.add_argument("--p", type=float, required=True, help="pressure (Pa, 100-1000)")
    s.add_argument("--g", type=float, required=True, help="showerhead gap (mm, 10-25)")
    s.add_argument("--json", default=str(DEFAULT_JSON), help="surrogate JSON path")

    f = sub.add_parser("pareto", help="Pareto front: minimum precursor per saturating pulse for each pulse time")
    f.add_argument("--max-lag", type=float, default=5.0, help="edge-lag constraint (%%)")
    f.add_argument("--n", type=int, default=21, help="grid points per factor")
    f.add_argument("--json", default=str(DEFAULT_JSON), help="surrogate JSON path")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    out = sys.stdout
    if a.cmd == "rule":
        r = rule_pulse(a.x, a.Q, a.p, a.T, a.AR, a.d, a.density, N0=a.N0, factor=a.factor)
        out.write(f"A_wall/A_wafer = {r['wall_ratio']:.2f}\n")
        out.write(f"U = {r['U_m_s']:.4f} m/s below the 300 mm showerhead at {a.p:g} Pa, {a.T:g} K\n")
        out.write(f"c0 = {r['c0']:.4g} molecules/m3 (x = {a.x:g})\n")
        out.write(f"t_supply = (1 + A_wall/A_wafer) N0 / (U c0) = {r['t_supply_s']:.3f} s\n")
        out.write(f"t_pulse = {a.factor:g} x t_supply = {r['t_pulse_s']:.3f} s\n")
    elif a.cmd == "surrogate":
        sur = load_surrogate(a.json)
        pr = sur.predict(a.x, a.Q, a.p, a.g)
        lo, hi = pr["t90_s"] / pr["t90_sd_factor"], pr["t90_s"] * pr["t90_sd_factor"]
        out.write(f"x = {a.x:g}, Q = {a.Q:g} sccm, p = {a.p:g} Pa, gap = {a.g:g} mm\n")
        out.write(f"t90 worst point = {pr['t90_s']:.3f} s (GP sd factor {pr['t90_sd_factor']:.3f}: "
                  f"{lo:.3f} to {hi:.3f} s)\n")
        out.write(f"edge lag = {pr['edge_lag_pct']:.2f} %\n")
        out.write(f"utilization = {100 * pr['utilization']:.1f} %  "
                  f"(TMA per saturating pulse = {1 / pr['utilization']:.3f} x stoichiometric minimum)\n")
        r2 = sur.loo_r2
        out.write("fit LOO R2: " + ", ".join(f"{k} {v:.3f}" for k, v in r2.items() if v is not None) + "\n")
    elif a.cmd == "pareto":
        sur = load_surrogate(a.json)
        front, _ = sur.pareto_front(max_lag=a.max_lag, n=a.n)
        cols = ["t90_s", "util", "precursor_over_minimum", "x", "Q_sccm", "p_Pa", "g_mm", "tma_sccm_s"]
        heads = ["t90_s", "util", "precursor/min", "x", "Q_sccm", "p_Pa", "g_mm", "tma_sccm_s"]
        out.write(f"Pareto front, edge lag <= {a.max_lag:g} %, {a.n}^4 grid, {len(front['t90_s'])} points "
                  f"(precursor/min = TMA per saturating pulse over the stoichiometric minimum = 1/util)\n")
        out.write("".join(f"{c:>14s}" for c in heads) + "\n")
        for i in range(len(front["t90_s"])):
            out.write("".join(f"{front[c][i]:14.4g}" for c in cols) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
