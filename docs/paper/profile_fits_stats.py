"""Every summary statistic of the profile-by-profile fits (figures/profile_fits*.csv): two fit windows, two film thicknesses
(the thicker S_IIb of the main analysis and the thinner S_I), three transport kernels, with and without narrowing, and five
sets of fronts, so that no combination is reported selectively.
Sets: all 16 fronts; 13 measurements (the two channels of an EDS scan averaged in ln s); the 3 EDS scans (channels averaged);
the 10 reflectometry fronts; the 9 reflectometry fronts at a 500 nm gap.
Slopes of ln s against the filled fraction, three ways:
  slope             ordinary least squares on the best values (fits without an upper limit enter at their best value);
  slope_uncensored  the same without the fits whose 16-84 % interval reaches an end of the simulated range;
  slope_censored    maximum likelihood with normal scatter in which a fit without an upper limit is a right-censored
                    observation at its 16th percentile and a fit without a lower limit a left-censored one at its 84th
                    (asymptotic standard error from a finite-difference Hessian of the log-likelihood, without a small-sample
                    correction, so it is about 7 % smaller than the least-squares convention at n = 16).
loo_min / loo_max: range of the least-squares slope when the fronts of one sample at a time are left out (all 16 only);
slope_without_1000_cycle_sample: the slope without the three fronts of the 1000-cycle sample, the only ones above f = 0.22.
Also: counts of fits at an end of the range and without an upper / lower limit, standard deviation of ln s, geometric mean,
mean rms, and (rows with narrowing) the number of fronts whose rms is lower with narrowing than without.
Writes figures/profile_fits_stats.csv.   .venv/bin/python docs/paper/profile_fits_stats.py"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
import warnings; warnings.filterwarnings("ignore", category=RuntimeWarning)
ROOT = Path(__file__).resolve().parents[2]; F = ROOT / "figures"; rows = []
def ols(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.ptp(x) == 0: return np.nan, np.nan
    if len(x) == 3: return float(np.polyfit(x, y, 1)[0]), np.nan
    c, cov = np.polyfit(x, y, 1, cov=True); return float(c[0]), float(np.sqrt(cov[0, 0]))
def tobit(x, y, lo, up, right, left):
    """y exact where not censored; right-censored at lo (true value above), left-censored at up (true value below)."""
    x, y = np.asarray(x, float), np.asarray(y, float); ex = ~(right | left)
    if ex.sum() < 3 or not (right | left).any(): return ols(x, y)
    def nll(p):
        a, b, ls = p; sg = np.exp(ls); m = a + b * x; v = -norm.logpdf(y[ex], m[ex], sg).sum()
        if right.any(): v -= norm.logsf(lo[right], m[right], sg).sum()
        if left.any(): v -= norm.logcdf(up[left], m[left], sg).sum()
        return v
    b0 = np.polyfit(x, y, 1); r = minimize(nll, [b0[1], b0[0], np.log(max(np.std(y - np.polyval(b0, x)), 0.1))], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxiter=20000))
    p0 = r.x; h = 1e-4; Hm = np.zeros((3, 3))                      # finite-difference Hessian of the negative log-likelihood at the optimum
    for i in range(3):
        for j in range(3):
            ei, ej = np.eye(3)[i] * h, np.eye(3)[j] * h; Hm[i, j] = (nll(p0 + ei + ej) - nll(p0 + ei - ej) - nll(p0 - ei + ej) + nll(p0 - ei - ej)) / (4 * h * h)
    try: se = float(np.sqrt(np.linalg.inv(Hm)[1, 1]))
    except Exception: se = np.nan
    return float(p0[1]), se
for film, fs in (("S_IIb", ""), ("S_I", "_filmSI")):
    for win, suf in ((0.3, "_cut0.3"), (0.1, "")):
        f = F / f"profile_fits{suf}{fs}.csv"
        if not f.exists(): continue
        D = pd.read_csv(f); D["lns"] = np.log(D.s_fit); D["scan"] = D.id.str.replace("_left|_right", "", regex=True)
        for fam in ("hydraulic", "gas", "ballistic"):
            a = D[(D.transport == fam) & D.narrowing].set_index("id"); b = D[(D.transport == fam) & ~D.narrowing].set_index("id"); better = (a.rms < b.rms)
            for nar, g in ((False, b), (True, a)):
                g = g.reset_index(); avg = g.groupby("scan").agg(fill=("fill", "mean"), lns=("lns", "mean"), method=("method", "first"), H_nm=("H_nm", "first"), cu=("no_upper_limit", "max"), cl=("no_lower_limit", "max")).reset_index()
                R_ = g.method == "reflectometry"; R5 = R_ & (g.H_nm == 500)
                sets = {"all 16 fronts": g, "13 measurements (EDS channels averaged)": None, "3 EDS scans (channels averaged)": None, "10 reflectometry fronts": g[R_], "9 reflectometry fronts, 500 nm gap": g[R5]}
                for name, sub in sets.items():
                    if sub is None:
                        q = avg if name.startswith("13") else avg[avg.method == "EDS"]; x, y, cens = q.fill.values, q.lns.values, (q.cu | q.cl).values; cs = (np.nan, np.nan); loo = (np.nan, np.nan); m11 = (np.nan, np.nan)
                    else:
                        x, y, cens = sub.fill.values, sub.lns.values, (sub.no_upper_limit | sub.no_lower_limit).values; cs = tobit(x, y, np.log(sub.s_p16.values), np.log(sub.s_p84.values), sub.no_upper_limit.values, sub.no_lower_limit.values)
                        l = [ols(sub[sub["sample"] != smp].fill, sub[sub["sample"] != smp].lns)[0] for smp in sub["sample"].unique()] if name.startswith("all") else [np.nan]; loo = (np.nanmin(l), np.nanmax(l))
                        if name.startswith("all"): m11 = ols(sub[sub["sample"] != "Yim11"].fill, sub[sub["sample"] != "Yim11"].lns)
                        else: m11 = (np.nan, np.nan)
                    sl, se = ols(x, y); slu, seu = ols(x[~cens], y[~cens])
                    rows.append(dict(film=film, window=win, transport={"ballistic": "gas-free", "gas": "gas 300 Pa"}.get(fam, fam), narrowing=nar, fronts=name, n=len(x), at_range_end=int(sub.at_grid_edge.sum()) if sub is not None else np.nan, no_upper_limit=int(sub.no_upper_limit.sum()) if sub is not None else int(q.cu.sum()), no_lower_limit=int(sub.no_lower_limit.sum()) if sub is not None else int(q.cl.sum()),
                                     slope=round(sl, 2), slope_se=round(se, 2), slope_uncensored=round(slu, 2), slope_uncensored_se=round(seu, 2), n_uncensored=int((~cens).sum()), slope_censored=round(cs[0], 2), slope_censored_se=round(cs[1], 2), loo_min=round(loo[0], 2), loo_max=round(loo[1], 2), slope_without_1000_cycle_sample=round(m11[0], 2), slope_without_1000_cycle_sample_se=round(m11[1], 2),
                                     sd_ln_s=round(float(np.std(y, ddof=1)), 2), s_geomean=float(f"{np.exp(np.mean(y)):.3g}"), mean_rms=round(float(sub.rms.mean()), 4) if sub is not None else np.nan, n_rms_lower_with_narrowing=int(better[sub.id].sum()) if (sub is not None and nar) else np.nan))
R = pd.DataFrame(rows); R.to_csv(F / "profile_fits_stats.csv", index=False); pd.set_option("display.width", 300); pd.set_option("display.max_rows", 300)
print(R[R.fronts.isin(["all 16 fronts", "10 reflectometry fronts", "9 reflectometry fronts, 500 nm gap"])].drop(columns=["at_range_end", "s_geomean"]).to_string(index=False))
