"""Sticking coefficient per published LHAR sample per method, from the authors' own tabulated front slope and position
(Yim et al. 2020, Table 1; data/lhar_profiles/yim2020_table1.csv).
Methods (all: s = A x slope^2, slope = |d theta / d(x/H)| at half coverage):
  arts        A = 13.9          (Arts et al. 2019; continuum, D_K = 2 H vbar / 3)
  gonsalves   A = 11.1          (Gonsalves et al. 2024, Kn >> 1 branch; Ylilammi model, same D_K)
  corr_front  A = 13.9 f(x50)   (this work: D_K x (3/8) ln(L_c/H), cut-off L_c = front position)
  corr_pillar A = 13.9 f(L_p)   (cut-off = pillar mean free path 1/(n d): 170 / 520 / 2080 um for layouts v1a / v1b / v2a)
  ballistic   A_ball(s, x50)    (Machball Trench calibration, figures/ballistic_slope_calibration.csv, inverted per sample)
Interval per sample: the published slope has two significant figures; the half-unit rounding of the last digit is propagated.
Writes figures/lhar_slope_table.csv and prints group statistics.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src"))
def f(ar): return 0.375 * np.log(ar)
PILLAR_D = 4e-6; PITCH = {"v1a": 28e-6, "v1b": 49e-6, "v2a": 98e-6}      # Yim 2020 ESI Table S2 (hexagonal arrays)
def pillar_mfp(layout): return (np.sqrt(3) / 2 * PITCH[layout] ** 2) / PILLAR_D   # 1/(n d): mean flight between pillar hits
L_PILLAR = pillar_mfp("v1b")
cal = pd.read_csv(ROOT / "figures/ballistic_slope_calibration.csv"); b = cal[cal.model == "ballistic"]; h = cal[cal.model == "1d_hydraulic"]
print(f"pillar mean free path {L_PILLAR*1e6:.0f} um; continuum check: A(1D hydraulic, this code) = {h.A.median():.2f} (range {h.A.min():.2f}-{h.A.max():.2f})")
S = np.array(sorted(b.s.unique()))
ZR = {s0: (b[b.s == s0].z50.min(), b[b.s == s0].z50.max()) for s0 in S}
def invert_ballistic(slope, z50):
    """s such that the ballistic slope at z50 equals the measured one, interpolating in log s between the calibration
    curves that cover z50. Returns (s, flag); flag is True if z50 or the slope lies outside the calibrated range."""
    use = [s0 for s0 in S if ZR[s0][0] <= z50 <= ZR[s0][1]]
    if len(use) < 2: return float("nan"), True
    sl = []
    for s0 in use:
        d = b[b.s == s0].sort_values("z50"); sl.append(float(np.interp(z50, d.z50, d.slope)))
    sl = np.array(sl); use = np.array(use)
    if slope < sl[0] or slope > sl[-1]:
        i = 0 if slope < sl[0] else len(use) - 2
        e = np.log(use[i + 1] / use[i]) / np.log(sl[i + 1] / sl[i]); return float(use[i] * (slope / sl[i]) ** e), True
    return float(np.exp(np.interp(np.log(slope), np.log(sl), np.log(use)))), False
t = pd.read_csv(ROOT / "data/lhar_profiles/yim2020_table1.csv"); t["source"] = "Yim2020"; t["species"] = "TMA"; t["slope_rel_err"] = np.nan
# Arts 2019 Fig. 2, digitised here (data/lhar_profiles/arts2019_fig2_*.csv): front fitted between theta 0.3 and 0.7
def arts_front(name):
    a = pd.read_csv(ROOT / f"data/lhar_profiles/arts2019_fig2_{name}.csv"); sel = (a.normalized_thickness > 0.3) & (a.normalized_thickness < 0.7)
    k = np.polyfit(a.x_over_H[sel], a.normalized_thickness[sel], 1); return (0.5 - k[1]) / k[0], -k[0]
extra = []
for name, sp, relerr, smp in (("tma_limited", "TMA", 0.10, 101), ("h2o_limited", "H2O", 0.12, 102)):   # slope uncertainty incl. plateau normalisation (TMA) and fit-window sensitivity (H2O)
    z50, sl = arts_front(name)
    extra.append(dict(series="Arts", sample=smp, pulse_purge_s="dose-" + name, H_nm=500, pillar_layout="v1b", cycles=400, magnification=0, x50_tilde=z50, slope_dtheta_dxtilde=sl, source="Arts2019", species=sp, slope_rel_err=relerr))
t = pd.concat([t, pd.DataFrame(extra)], ignore_index=True); rows = []
for _, r in t.iterrows():
    H = r.H_nm * 1e-9; sl = r.slope_dtheta_dxtilde; z = r.x50_tilde
    dig = 10 ** (np.floor(np.log10(sl)) - 1); lo, hi = sl - dig / 2, sl + dig / 2
    if np.isfinite(r.slope_rel_err): lo, hi = sl * (1 - r.slope_rel_err), sl * (1 + r.slope_rel_err)
    out = dict(source=r.source, species=r.species, sample=int(r["sample"]), series=r.series, H_nm=r.H_nm, cycles=r.cycles, magnification=r.magnification, pulse=r.pulse_purge_s.split("-")[0], layout=r.pillar_layout, x50_tilde=z, slope=sl)
    for name, A in (("arts", 13.9), ("gonsalves", 11.1), ("corr_front", 13.9 * f(z)), ("corr_pillar", 13.9 * f(pillar_mfp(r.pillar_layout) / H))):
        out[f"s_{name}"] = A * sl ** 2; out[f"s_{name}_lo"] = A * lo ** 2; out[f"s_{name}_hi"] = A * hi ** 2
    sb, ex = invert_ballistic(sl, z); out["s_ballistic"] = sb; out["s_ballistic_lo"] = invert_ballistic(lo, z)[0]; out["s_ballistic_hi"] = invert_ballistic(hi, z)[0]; out["ballistic_extrapolated"] = ex
    out["factor_front"] = f(z); out["factor_pillar"] = f(pillar_mfp(r.pillar_layout) / H); out["pillar_mfp_um"] = pillar_mfp(r.pillar_layout) * 1e6
    rows.append(out)
T = pd.DataFrame(rows); T.to_csv(ROOT / "figures/lhar_slope_table.csv", index=False)
pd.set_option("display.width", 250); print(T[["source", "sample", "series", "H_nm", "cycles", "magnification", "x50_tilde", "slope", "s_arts", "s_gonsalves", "s_corr_front", "s_corr_pillar", "s_ballistic", "ballistic_extrapolated"]].to_string(float_format=lambda v: f"{v:.3g}"))
# reference set: 500 nm gap, 50x objective (5-6 um spot), 500 cycles, 0.1 s pulse: the like-for-like group
T["s_ballistic_relation"] = np.nan
ref = T[(T.source == "Yim2020") & (T.H_nm == 500) & (T.magnification == 50) & (T.cycles == 500) & (T.pulse == "0.1")]
print(f"\nreference group (500 nm, 50x, 500 cycles, 0.1 s TMA): samples {sorted(ref['sample'].tolist())}")
for m in ("arts", "gonsalves", "corr_front", "corr_pillar", "ballistic"):
    v = ref[f"s_{m}"]; g = np.exp(np.log(v).mean()); sd = np.exp(np.log(v).std(ddof=1))
    print(f"  {m:12s} geometric mean {g:.2e}  x/÷ {sd:.2f}  range {v.min():.2e}-{v.max():.2e}  ratio to SFG 4e-3: {g/4e-3:.2f}")
