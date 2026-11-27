"""Every summary statistic of the profile-by-profile fits (figures/profile_fits*.csv), for both fit windows, the three
transport kernels, with and without narrowing, and five sets of fronts, so that no combination is reported selectively.
Sets: all 16 fronts; 13 independent measurements (the two channels of an EDS scan averaged in ln s); the 3 EDS scans (channels
averaged); the 10 reflectometry fronts; the 9 reflectometry fronts at a 500 nm gap.
Columns: number of fits at an end of the simulated range, number whose 16-84 % interval has no upper / no lower limit,
slope of ln s against the filled fraction with its standard error (ordinary least squares; `slope_uncensored` leaves out fits
whose interval reaches an end of the range), standard deviation of ln s, geometric mean, mean rms, and for the rows with
narrowing the number of fronts whose rms is lower with narrowing than without.
Writes figures/profile_fits_stats.csv.   .venv/bin/python docs/paper/profile_fits_stats.py"""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; F = ROOT / "figures"; rows = []
def ols(x, y):
    if len(x) < 3 or np.ptp(x) == 0: return np.nan, np.nan
    c, cov = np.polyfit(x, y, 1, cov=True) if len(x) > 3 else (np.polyfit(x, y, 1), None); se = float(np.sqrt(cov[0, 0])) if cov is not None else np.nan; return float(c[0]), se
for win, suf in ((0.3, "_cut0.3"), (0.1, "")):
    D = pd.read_csv(F / f"profile_fits{suf}.csv"); D["lns"] = np.log(D.s_fit); D["scan"] = D.id.str.replace("_left|_right", "", regex=True)
    for fam in ("hydraulic", "gas", "ballistic"):
        a = D[(D.transport == fam) & D.narrowing].set_index("id"); b = D[(D.transport == fam) & ~D.narrowing].set_index("id"); better = (a.rms < b.rms)
        for nar, g in ((False, b), (True, a)):
            g = g.reset_index(); avg = g.groupby("scan").agg(fill=("fill", "mean"), lns=("lns", "mean"), method=("method", "first"), H_nm=("H_nm", "first"), cens=("no_upper_limit", "max"), cens2=("no_lower_limit", "max")).reset_index()
            sets = {"all 16 fronts": (g.fill, g.lns, g.no_upper_limit | g.no_lower_limit, g), "13 measurements (EDS channels averaged)": (avg.fill, avg.lns, avg.cens | avg.cens2, None), "3 EDS scans (channels averaged)": (avg[avg.method == "EDS"].fill, avg[avg.method == "EDS"].lns, None, None),
                    "10 reflectometry fronts": (g[g.method == "reflectometry"].fill, g[g.method == "reflectometry"].lns, (g.no_upper_limit | g.no_lower_limit)[g.method == "reflectometry"], g[g.method == "reflectometry"]), "9 reflectometry fronts, 500 nm gap": (g[(g.method == "reflectometry") & (g.H_nm == 500)].fill, g[(g.method == "reflectometry") & (g.H_nm == 500)].lns, None, g[(g.method == "reflectometry") & (g.H_nm == 500)])}
            for name, (x, y, cens, sub) in sets.items():
                sl, se = ols(x.values, y.values); slu, seu = ols(x.values[~cens.values], y.values[~cens.values]) if cens is not None else (np.nan, np.nan)
                rows.append(dict(window=win, transport={"ballistic": "gas-free", "gas": "gas 300 Pa"}.get(fam, fam), narrowing=nar, fronts=name, n=len(x), at_range_end=int(sub.at_grid_edge.sum()) if sub is not None else np.nan, no_upper_limit=int(sub.no_upper_limit.sum()) if sub is not None else np.nan, no_lower_limit=int(sub.no_lower_limit.sum()) if sub is not None else np.nan,
                                 slope=round(sl, 2), slope_se=round(se, 2), slope_uncensored=round(slu, 2), slope_uncensored_se=round(seu, 2), n_uncensored=int((~cens.values).sum()) if cens is not None else np.nan, sd_ln_s=round(float(np.std(y, ddof=1)), 2), s_geomean=float(f"{np.exp(np.mean(y)):.3g}"),
                                 mean_rms=round(float(sub.rms.mean()), 4) if sub is not None else np.nan, n_rms_lower_with_narrowing=int(better[sub.id].sum()) if (sub is not None and nar) else np.nan))
R = pd.DataFrame(rows); R.to_csv(F / "profile_fits_stats.csv", index=False); pd.set_option("display.width", 280)
print(R[R.fronts.isin(["all 16 fronts", "13 measurements (EDS channels averaged)", "10 reflectometry fronts", "3 EDS scans (channels averaged)"])].to_string(index=False))
