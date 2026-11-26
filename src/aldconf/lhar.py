"""Sticking coefficient from the front slope of a lateral-channel (PillarHall-type) saturation profile, with the effects that
the plain slope method s = 13.9 slope^2 leaves out: narrowing of the gap by the film over the cycles, the measurement spot,
and the slit transport (hydraulic diffusion, or a Monte Carlo of the real 500 nm cavity without gas or with N2 at 100, 300
or 1000 Pa and 300 C). The slopes come from a table of the complete forward model (lhar_slope_table.csv, built by
docs/paper/build_slope_table.py): depth 150 to 300 gaps, filled fraction 0 to 0.4, spot 0 or 11 gaps, s 0.0025 to 0.06.

    >>> from aldconf.lhar import sticking_from_slope
    >>> r = sticking_from_slope(slope=0.019, z50=208, fill=0.19, transport="gas300")
    >>> r["s"], r["s_lo"], r["s_hi"], r["sensitivity"]

The result says how much the slope still depends on s (0.5 for the ideal slope method) and whether the slope pins s down at
all; at 40 % filling it does not, for any transport model.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import numpy as np

TRANSPORTS = ("hydraulic", "free", "gas100", "gas300", "gas1000")


def narrowing_limit(fill: float) -> float:
    """Upper bound of slope x half-thickness depth for a film that fills the fraction `fill` = 2 t / H of the gap."""
    if not 0 < fill <= 1: raise ValueError("fill must be in (0, 1]")
    return (1 - fill / 2) / fill


def relative_sensitivity(fill: float, z50: float, slope_unfilled: float) -> float:
    """How much of the slope's sensitivity to s survives channel narrowing, relative to an unfilled channel (1 = all of it).
    One number decides: u = fill x z50 / (1 - fill/2) x slope_unfilled, the width of the narrowing smear over the width of
    the intrinsic front; the complete model follows 1 / (1 + (u / 0.92)^3) within 0.07 (rms, 567 table points, three
    transport families). 80 % of the sensitivity is kept for u < 0.58, half at u = 0.92."""
    if fill <= 0: return 1.0
    u = fill * z50 / (1 - fill / 2) * slope_unfilled; return float(1.0 / (1.0 + (u / 0.92) ** 3))


def max_fill_for(z50: float, slope_unfilled: float, keep: float = 0.8) -> float:
    """Largest filled fraction 2 t / H that keeps the fraction `keep` of the slope's sensitivity to s."""
    if not 0 < keep < 1: raise ValueError("keep must be in (0, 1)")
    u = 0.92 * (1 / keep - 1) ** (1 / 3); a = u / (z50 * slope_unfilled); return float(a / (1 + a / 2))


@lru_cache(maxsize=1)
def _table():
    import csv
    rows = list(csv.DictReader(open(Path(__file__).with_name("lhar_slope_table.csv"))))
    T = {}
    for r in rows:
        T.setdefault((r["transport"], float(r["spot_gaps"])), {}).setdefault((float(r["z50"]), float(r["fill"])), []).append((float(r["s"]), float(r["slope"]) if r["slope"] else np.nan))
    return T


def slope_curve(z50: float, fill: float, transport: str = "gas300", spot_gaps: float = 11.0):
    """(s, slope) arrays of the complete model at this depth and filled fraction (bilinear in depth and fill, spot 0 or 11 gaps)."""
    if transport not in TRANSPORTS: raise ValueError(f"transport must be one of {TRANSPORTS}")
    if spot_gaps not in (0.0, 11.0): raise ValueError("spot_gaps must be 0 or 11 (5.5 um on a 500 nm gap)")
    tab = _table()[(transport, float(spot_gaps))]; zs = sorted({k[0] for k in tab}); fs = sorted({k[1] for k in tab})
    if not (zs[0] <= z50 <= zs[-1] and fs[0] <= fill <= fs[-1]): raise ValueError(f"outside the table: depth {zs[0]:g}-{zs[-1]:g} gaps, fill {fs[0]:g}-{fs[-1]:g}")
    i = min(np.searchsorted(zs, z50, side="right") - 1, len(zs) - 2); j = min(np.searchsorted(fs, fill, side="right") - 1, len(fs) - 2)
    u = (z50 - zs[i]) / (zs[i + 1] - zs[i]); v = (fill - fs[j]) / (fs[j + 1] - fs[j]); s = np.array([a for a, _ in tab[(zs[i], fs[j])]]); out = np.zeros(len(s))
    for (a, b, w) in ((i, j, (1 - u) * (1 - v)), (i + 1, j, u * (1 - v)), (i, j + 1, (1 - u) * v), (i + 1, j + 1, u * v)):
        if w > 1e-12: out = out + w * np.log(np.array([c for _, c in tab[(zs[a], fs[b])]]))
    ok = np.isfinite(out)
    if ok.sum() < 3: raise ValueError("this depth is beyond the simulated exposures for this transport model")
    return s[ok], np.exp(out[ok])


def sticking_from_slope(slope: float, z50: float, fill: float, transport: str = "gas300", spot_gaps: float = 11.0, rel_slope_error: float = 0.13) -> dict:
    """Invert the complete model. `slope` is d(thickness/plateau)/d(x/H) at half thickness, `z50` the half-thickness depth in
    gaps, `fill` = 2 t / H. `rel_slope_error` is the one-sigma relative error of the slope (0.13 is the scatter of repeated
    runs in Yim 2020). Returns s, its interval (None where the slope no longer limits s on that side within the simulated
    range), the sensitivity d ln(slope) / d ln(s) at the solution, the model slope at the largest simulated s (the ceiling),
    and the plain slope-method value for comparison."""
    s, sl = slope_curve(z50, fill, transport, spot_gaps); ls, lsl = np.log(s), np.log(sl); mono = np.maximum.accumulate(lsl)
    def inv(x):
        if x >= mono[-1] or x <= mono[0]: return None
        return float(np.exp(np.interp(x, mono, ls)))
    x = np.log(slope); s0 = inv(x); lo = inv(x - rel_slope_error); hi = inv(x + rel_slope_error)
    sens = None
    if s0 is not None:
        k = int(np.clip(np.searchsorted(s, s0), 1, len(s) - 1)); sens = float((lsl[k] - lsl[k - 1]) / (ls[k] - ls[k - 1]))
    frac = float(slope / sl[-1]); note = []
    if s0 is None: note.append("measured slope is above the model ceiling: no s reproduces it" if x >= mono[-1] else "measured slope is below the model at the smallest simulated s")
    if s0 is not None and hi is None: note.append("no upper bound on s within the simulated range")
    if s0 is not None and lo is None: note.append("no lower bound on s within the simulated range")
    return dict(s=s0, s_lo=lo, s_hi=hi, sensitivity=sens, fraction_of_ceiling=frac, slope_at_largest_s=float(sl[-1]), s_range=(float(s[0]), float(s[-1])),
                slope_method_s=13.9 * slope ** 2, narrowing_limit_fraction=(slope * z50 / narrowing_limit(fill)) if fill > 0 else None, transport=transport, note="; ".join(note))
