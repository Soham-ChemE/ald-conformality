"""Parsers for the patternedSurface BC output (thetaFlat, thetaBottom, thetaMid, Keff, value, and the per-face
pattern ARField, holeDensityField) per wafer face per time."""
import re
from pathlib import Path
import numpy as np

NR = 150

def _block(txt, name="wafer"):
    m = re.search(r"\n\s+%s\s*\{(.*?)\n\s+\}" % name, txt, re.S)
    if m: return m.group(1)
    m = re.search(r"\n\s+%s\s*\{([^}]*)\}" % name, txt, re.S)
    return m.group(1) if m else ""

def _field(blk, key, n=NR, default=0.0):
    m = re.search(r"\b%s\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)" % key, blk, re.S)
    if m: return np.fromstring(m.group(2), sep="\n")
    m = re.search(r"\b%s\s+uniform\s+([0-9.eE+-]+)" % key, blk)
    return np.full(n, float(m.group(1))) if m else np.full(n, default)

def _scalar(blk, key, default=np.nan):
    m = re.search(r"\b%s\s+([0-9.eE+-]+);" % key, blk)
    return float(m.group(1)) if m else default

def load_patterned(case):
    """Returns t, fields, params, pulse. fields: thetaFlat, thetaBottom, thetaMid, Keff, value (ntime x 150), plus
    ARField and holeDensityField (ntime x 150; filled from the scalar AR / holeDensity for cases written before the
    per-face pattern existed). params: the scalar entries s, vbar, N0, c0, d, AR, holeDensity, nCells."""
    case = Path(case)
    times = sorted((float(d.name), d) for d in case.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times])
    keys = ("thetaFlat", "thetaBottom", "thetaMid", "Keff", "value")
    out = {k: np.zeros((len(t), NR)) for k in keys + ("ARField", "holeDensityField")}
    params = {}
    for i, (_, d) in enumerate(times):
        blk = _block((d / "T").read_text())
        if not params and "patternedSurface" in blk:
            params = {k: _scalar(blk, k) for k in ("s", "vbar", "N0", "c0", "d", "AR", "holeDensity", "nCells")}
        for k in keys: out[k][i] = _field(blk, k)
        out["ARField"][i] = _field(blk, "ARField", default=params.get("AR", np.nan))
        out["holeDensityField"][i] = _field(blk, "holeDensityField", default=params.get("holeDensity", np.nan))
    pulse = float(re.search(r"\(\(0 1\) \(([0-9.]+) 1\)", (case / "0/T").read_text()).group(1))
    return t, out, params, pulse
