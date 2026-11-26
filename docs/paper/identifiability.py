"""Phase 2: what a lateral-channel saturation profile can say about the sticking coefficient.

For a profile measured at 2 um steps on a 500 nm gap (one point per 4 gaps, 0 to 400 gaps) with independent Gaussian noise
sigma on thickness/plateau, and the complete forward model (transport kernel + channel narrowing + 5.5 um spot):
  1. information  the Cramer-Rao bound on ln(s) from the full profile, with the exposure and the plateau height as nuisance
                  parameters, against the filled fraction, for each transport model and s; and the same from the slope alone
                  (slope error / sensitivity);
  2. ambiguity    a profile generated with one transport model is fitted with another (s, exposure and plateau free): the rms
                  misfit says whether the data could tell the models apart, the fitted s is the bias if they cannot;
  3. rule         the loss of sensitivity as a function of one number, u = f x z50 x slope_0 (width of the narrowing smear over
                  the width of the intrinsic front), with a fitted closed form;
  4. designs      the bound on ln(s) for candidate experiments (cycle number at fixed gap), for single profiles and pairs.
Writes figures/identifiability_information.csv, identifiability_ambiguity.csv, identifiability_summary.json and
docs/paper/figures/fig13_identifiability_map.png.
  PYTHONPATH=src .venv/bin/python docs/paper/identifiability.py
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar, curve_fit
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
FIG = ROOT / "figures"; X = np.arange(2.0, 400, 4.0); SPOT = 11.0; Z50 = 215.0; SIGMA = 0.03
KERN = {"hydraulic": ("hyd", 1.0), "free": ("walls90_pillars49", 1.0), "gas300": ("gas300", 0.92)}
FILLS = (0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)

def prof(run, e, fill, n): return mt.cumulative(run, e, fill, X, SPOT, n)

def expo(run, fill, n, z50=Z50): return mt.slope_at_depth(run, z50, fill, SPOT, n)

rows = []
for name, (kn, n) in KERN.items():
    K = mt.kernel(kn); ss = list(K)
    for fill in FILLS:
        for k in range(1, len(ss) - 1):
            s = ss[k]; sl, e = expo(K[s], fill, n)
            if not np.isfinite(sl): continue
            y = prof(K[s], e, fill, n)
            dls = (prof(K[ss[k + 1]], e, fill, n) - prof(K[ss[k - 1]], e, fill, n)) / np.log(ss[k + 1] / ss[k - 1])       # d y / d ln s at fixed exposure
            dle = (prof(K[s], e * 1.03, fill, n) - prof(K[s], e / 1.03, fill, n)) / (2 * np.log(1.03))                    # d y / d ln exposure
            J = np.c_[dls, dle, y]; F = J.T @ J / SIGMA ** 2; cov = np.linalg.inv(F + 1e-12 * np.eye(3))
            sl_lo, sl_hi = expo(K[ss[k - 1]], fill, n)[0], expo(K[ss[k + 1]], fill, n)[0]
            sens = np.log(sl_hi / sl_lo) / np.log(ss[k + 1] / ss[k - 1]) if np.isfinite(sl_lo) and np.isfinite(sl_hi) else np.nan
            rows.append(dict(transport=name, fill=fill, s=s, exposure=round(float(e), 3), slope=round(float(sl), 5), sensitivity=round(float(sens), 4),
                             sigma_ln_s_profile=round(float(np.sqrt(cov[0, 0])), 4), sigma_ln_s_profile_exposure_known=round(float(np.sqrt(np.linalg.inv(F[np.ix_([0, 2], [0, 2])])[0, 0])), 4),
                             sigma_ln_s_slope=round(float(0.13 / sens), 3) if sens > 1e-3 else np.inf))
    print(name, "information done", flush=True)
I = pd.DataFrame(rows); I.to_csv(FIG / "identifiability_information.csv", index=False)

# ---- transport ambiguity: truth in one model, fit with another ----
amb = []
def best_fit(Kfit, nfit, y, fill):
    out = []
    for s, run in Kfit.items():
        def c(le):
            m = prof(run, 10 ** le, fill, nfit); a = float(m @ y / (m @ m)); return float(np.mean((a * m - y) ** 2))
        q = minimize_scalar(c, bounds=(np.log10(run[0][3]), np.log10(run[0][-2])), method="bounded", options={"xatol": 1e-4}); out.append((s, q.fun))
    out = np.array(out); j = int(np.argmin(out[:, 1])); jj = min(max(j, 1), len(out) - 2); p = np.polyfit(np.log(out[jj - 1:jj + 2, 0]), out[jj - 1:jj + 2, 1], 2)
    sb = float(np.exp(np.clip(-p[1] / (2 * p[0]), np.log(out[jj - 1, 0]), np.log(out[jj + 1, 0])))) if p[0] > 0 else float(out[j, 0])
    return sb, float(np.sqrt(out[j, 1])), bool(j in (0, len(out) - 1))
for truth, s_true in (("gas300", 0.0158), ("gas300", 0.008), ("free", 0.025), ("hydraulic", 0.008)):
    Kt, nt = mt.kernel(KERN[truth][0]), KERN[truth][1]
    for fill in (0.0, 0.05, 0.1, 0.2, 0.4):
        sl, e = expo(Kt[s_true], fill, nt)
        if not np.isfinite(sl): continue
        y = prof(Kt[s_true], e, fill, nt)
        for fitm in KERN:
            if fitm == truth: continue
            sb, rms, edge = best_fit(mt.kernel(KERN[fitm][0]), KERN[fitm][1], y, fill)
            amb.append(dict(truth=truth, s_true=s_true, fill=fill, fitted_with=fitm, s_fitted=float(f"{sb:.3g}"), bias_factor=round(sb / s_true, 2), rms_misfit=round(rms, 4), at_grid_edge=edge))
    print(truth, s_true, "ambiguity done", flush=True)
A = pd.DataFrame(amb); A.to_csv(FIG / "identifiability_ambiguity.csv", index=False)

# ---- one-number rule: sensitivity against u = f x z1 x slope_0, z1 = z50 / (1 - f/2) ----
T = pd.read_csv(FIG / "slope_table.csv").dropna(subset=["slope"]); T = T[T.spot_gaps == 11.0]; pts = []
for (tr, z50), g in T.groupby(["transport", "z50"]):
    base = g[g.fill == 0.0].set_index("s").slope
    for fill, h in g.groupby("fill"):
        h = h.sort_values("s"); s = h.s.values; sl = h.slope.values
        for k in range(1, len(s) - 1):
            if s[k] not in base.index or s[k - 1] not in base.index or s[k + 1] not in base.index: continue
            sens = np.log(sl[k + 1] / sl[k - 1]) / np.log(s[k + 1] / s[k - 1]); sens0 = np.log(base[s[k + 1]] / base[s[k - 1]]) / np.log(s[k + 1] / s[k - 1])
            pts.append(dict(transport=tr, z50=z50, fill=fill, s=s[k], u=fill * z50 / (1 - fill / 2) * base[s[k]], ratio=sens / sens0, sens=sens, sens0=sens0))
P = pd.DataFrame(pts); Pf = P[(P.fill > 0) & (P.sens0 > 0.2)]
g_ = lambda u, a, p: 1.0 / (1.0 + (u / a) ** p); (a_, p_), _ = curve_fit(g_, Pf.u.values, Pf.ratio.values, p0=(1.5, 2.0)); res = Pf.ratio.values - g_(Pf.u.values, a_, p_)
u_half = float(a_); u_80 = float(a_ * (0.25) ** (1 / p_))
summ = dict(noise_sigma=SIGMA, z50=Z50, rule=dict(form="sensitivity / sensitivity_unfilled = 1 / (1 + (u / a)^p), u = f x z50 / (1 - f/2) x slope_unfilled", a=round(float(a_), 3), p=round(float(p_), 3), rms_residual=round(float(np.sqrt(np.mean(res ** 2))), 3), n_points=int(len(Pf)),
                                                   u_for_80_percent=round(u_80, 3), u_for_50_percent=round(u_half, 3), mean_unfilled_sensitivity={k: round(float(v), 3) for k, v in P[P.fill == 0].groupby("transport").sens.mean().items()}))
# ---- candidate experiments: bound on ln s from one profile and from pairs (independent noise, own exposure and plateau each) ----
des = {}
for name in ("gas300", "hydraulic"):
    g = I[(I.transport == name) & (np.isclose(I.s, 0.0158 if name == "gas300" else 0.008))].set_index("fill")
    des[name] = {f"{int(round(f * 2500))} cycles ({int(round(f * 100))} % filled)": dict(sigma_ln_s_profile=float(g.loc[f, "sigma_ln_s_profile"]), factor_1sigma=round(float(np.exp(g.loc[f, "sigma_ln_s_profile"])), 2), sigma_ln_s_slope=float(g.loc[f, "sigma_ln_s_slope"])) for f in g.index if f > 0}
summ["designs_500nm_gap"] = des
summ["ambiguity_examples"] = A[(A.truth == "gas300") & (A.s_true == 0.0158)].to_dict("records")
(FIG / "identifiability_summary.json").write_text(json.dumps(summ, indent=1))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5, "savefig.dpi": 200})
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5)); col = {"hydraulic": "k", "free": "C0", "gas300": "C3"}
for tr, g in P[P.fill > 0].groupby("transport"): ax[0].plot(g.u, g.ratio, ".", ms=3, color=col.get(tr, "C2"), alpha=.6, label=tr if tr in col else None)
uu = np.geomspace(0.05, 12, 200); ax[0].plot(uu, g_(uu, a_, p_), "k-", lw=1.2, label=f"1 / (1 + (u/{a_:.2f})$^{{{p_:.2f}}}$)"); ax[0].set_xscale("log"); ax[0].set_xlabel("u = f · z$_{50}$/(1 − f/2) · slope without narrowing"); ax[0].set_ylabel("sensitivity relative to an unfilled channel"); ax[0].set_title("(a) one number decides"); ax[0].legend(frameon=False)
for tr, c in col.items():
    g = I[(I.transport == tr) & np.isclose(I.s, {"hydraulic": 0.008, "free": 0.025, "gas300": 0.0158}[tr])].sort_values("fill")
    ax[1].plot(g.fill * 100, np.exp(g.sigma_ln_s_profile), "-o", ms=3, color=c, label=f"{tr}: full profile"); ax[1].plot(g.fill * 100, np.exp(np.minimum(g.sigma_ln_s_slope, 5)), ":s", ms=3, color=c, label=f"{tr}: slope alone")
ax[1].set_yscale("log"); ax[1].set_xlabel("gap filled by the film (%)"); ax[1].set_ylabel("1σ uncertainty factor on s"); ax[1].set_title(f"(b) best possible precision, noise {SIGMA}"); ax[1].legend(frameon=False, fontsize=6.5)
g = A[(A.truth == "gas300") & (A.s_true == 0.0158)]
for fm, c in (("hydraulic", "k"), ("free", "C0")):
    h = g[g.fitted_with == fm]; ax[2].plot(h.fill * 100, h.rms_misfit, "-o", ms=3, color=c, label=f"fitted with {fm} (s off by ×{h.bias_factor.iloc[0]:.2g} to ×{h.bias_factor.iloc[-1]:.2g})")
ax[2].axhline(SIGMA, color="grey", ls=":", lw=1); ax[2].text(1, SIGMA * 1.05, "noise", fontsize=7, color="grey"); ax[2].set_xlabel("gap filled by the film (%)"); ax[2].set_ylabel("rms misfit of the wrong transport model"); ax[2].set_title("(c) can the profile tell the transport models apart?"); ax[2].legend(frameon=False, fontsize=6.5)
fig.tight_layout(); fig.savefig(ROOT / "docs/paper/figures/fig13_identifiability_map.png"); plt.close(fig)
print(json.dumps(summ["rule"], indent=1)); print(json.dumps(des, indent=1)); print(A.to_string(index=False))
print(I[np.isclose(I.s, 0.0158) & (I.transport == "gas300")].to_string(index=False))
