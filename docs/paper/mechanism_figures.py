"""Figures and summary numbers of the mechanism study (draft_v3), from figures/mechanism_*.csv/json and the raw Monte Carlo
frames in figures/pillar_mc/ (untracked; figures/pillar_mc_manifest.csv lists them). Writes
docs/paper/figures/fig9_transport_candidates.png, fig10_narrowing_law.png, fig11_heldout.png, fig12_identifiability.png,
figures/mechanism_transport.json (front diagnostics of every transport variant) and figures/mechanism_narrowing_check.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/mechanism_figures.py
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).parent))
import mechanism_test as mt
FIG = ROOT / "figures"; OUT = ROOT / "docs/paper/figures"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "legend.fontsize": 7.5, "figure.dpi": 150, "savefig.dpi": 200})
S0, Z0 = 0.0158, 215.0

# ---- 1. transport candidates at one sticking value and one front depth ----
VARIANTS = [("walls90_pillars49", "real cavity, cosine walls, no gas", "C0", "-"), ("spec20", "20 % specular", "C1", ":"), ("spec50", "50 % specular", "C1", "--"),
            ("retro20", "20 % back-scattered", "C2", ":"), ("retro50", "50 % back-scattered", "C2", "--"), ("cos2", "cos$^2$ lobe", "C4", ":"), ("cos4", "cos$^4$ lobe", "C4", "--"),
            ("cos8", "cos$^8$ lobe", "C4", "-."), ("gas100", "N$_2$ 100 Pa", "C3", ":"), ("gas300", "N$_2$ 300 Pa", "C3", "-"), ("gas1000", "N$_2$ 1000 Pa", "C3", "--"), ("hyd", "hydraulic diffusion model", "k", "-")]
diag = {}; curves = {}
for name, label, c, ls in VARIANTS:
    K = mt.kernel(name)
    if S0 not in K: continue
    run = K[S0]; E, TH, z = run
    z50 = np.array([mt.front(run, e, 0)[0] for e in E]); ok = np.isfinite(z50)
    if ok.sum() < 5 or np.nanmax(z50) < Z0: continue
    e = float(np.interp(Z0, z50[ok], E[ok])); th = mt.at(E, TH, e); sl = mt.front(run, e, 0)[1]
    diag[name] = dict(label=label.replace("$", ""), exposure_Pa_s=round(e, 2), A=round(S0 / sl ** 2, 1), slope=round(float(sl), 4), tail_plus_60=round(float(np.interp(Z0 + 60, z, th)), 4), tail_plus_100=round(float(np.interp(Z0 + 100, z, th)), 4))
    curves[name] = (z, th)
(FIG / "mechanism_transport.json").write_text(json.dumps(dict(s=S0, z50=Z0, variants=diag), indent=1))
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5), gridspec_kw=dict(width_ratios=[1, 1, 0.8]))
for name, label, c, ls in VARIANTS:
    if name not in curves: continue
    z, th = curves[name]; j = 0 if name in ("walls90_pillars49", "spec20", "spec50", "retro20", "retro50", "cos2", "cos4", "cos8", "hyd") else None
    if name in ("walls90_pillars49", "hyd") or j == 0: ax[0].plot(z, th, color=c, ls=ls, lw=1.2, label=label)
    if name in ("walls90_pillars49", "hyd", "gas100", "gas300", "gas1000"): ax[1].plot(z, th, color=c, ls=ls, lw=1.2, label=label)
for a, t in ((ax[0], "(a) wall-scattering law"), (ax[1], "(b) carrier-gas collisions")):
    a.set_xlim(120, 460); a.set_ylim(0, 1.02); a.set_xlabel("distance / gap"); a.set_ylabel("coverage"); a.set_title(t); a.legend(frameon=False, loc="upper right", fontsize=6.5)
names = [n for n, *_ in VARIANTS if n in diag]; y = np.arange(len(names))
ax[2].barh(y, [diag[n]["tail_plus_60"] for n in names], color=[v[2] for v in VARIANTS if v[0] in diag]); ax[2].set_yticks(y); ax[2].set_yticklabels([v[1] for v in VARIANTS if v[0] in diag], fontsize=7); ax[2].invert_yaxis()
ax[2].set_xlabel("coverage 60 gaps ahead of the front"); ax[2].set_title("(c) coverage ahead of the front")
fig.tight_layout(); fig.savefig(OUT / "fig9_transport_candidates.png"); plt.close(fig)

# ---- 2. the narrowing law against the full multi-cycle continuum model, and against the measured series ----
from aldconf import Slit, Precursor, Process, simulate
def prof(H, cyc, s, narrow):
    pr = Precursor(sticking=s, site_density=4e18, mass_amu=72.09, gpc=1e-10)
    r = simulate(Process(Slit(H, 250e-6 * H / 500e-9, knudsen_correction="none"), pr, T=573.15, p_partial=25, t_dose=0.1, cycles=cyc, n_cells=400, n_substeps=60, knudsen_factor=1.0, narrowing=narrow))
    return r.x / H, r.thickness / r.thickness[:10].mean()
rows = []; fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
for (H, cyc), c in zip(((500e-9, 250), (500e-9, 500), (500e-9, 1000), (100e-9, 500), (2000e-9, 500)), ("C0", "C1", "C2", "C3", "C4")):
    for s in (0.01, 0.03):
        z, y1 = prof(H, 1, s, False); z, yN = prof(H, cyc, s, True); fill = 2 * cyc * 1e-10 / H
        ys = np.mean([np.interp(z / f, z, y1) for f in np.linspace(max(1 - fill, 0.02), 1, 41)], 0); ys /= ys[:10].mean()
        def fs(y):
            i = int(np.argmax(y < 0.5)); z50 = z[i - 1] + (y[i - 1] - .5) / (y[i - 1] - y[i]) * (z[i] - z[i - 1]); cr = lambda l: np.interp(l, y[::-1], z[::-1]); return z50, 0.2 / (cr(.4) - cr(.6))
        a, b = fs(yN), fs(ys); rows.append(dict(H_nm=H * 1e9, cycles=cyc, s=s, fill=fill, z50_full=a[0], slope_full=a[1], z50_law=b[0], slope_law=b[1], max_abs_diff=float(np.abs(yN - ys).max())))
        if s == 0.03: ax[0].plot(z, yN, color=c, lw=1.2, label=f"{H*1e9:.0f} nm, {cyc} cycles"); ax[0].plot(z, ys, color=c, lw=1.2, ls=":")
pd.DataFrame(rows).to_csv(FIG / "mechanism_narrowing_check.csv", index=False)
ax[0].set_xlim(0, 220); ax[0].set_xlabel("distance / initial gap"); ax[0].set_ylabel("normalised thickness"); ax[0].set_title("(a) full multi-cycle model (lines) and the law (dotted)"); ax[0].legend(frameon=False)
H = mt.heldout(); H["fill"] = np.minimum(2 * H.film_nm / H.H_nm, 1.0); H["z1"] = H.z50 / (1 - H.fill / 2); H["nn"] = H.measured_slope * H.fill * H.z1
ff = np.linspace(0.03, 1, 100); ax[1].plot(ff, 1 / (ff * (1 - 0)) * (1 - ff / 2), "k-", lw=1, label="step-front value: slope x depth = (1 - f/2)/f")
for mag, mk, lab in ((50, "o", "reflectometry, 5-6 um spot"), (10, "s", "25 um spot"), (5, "D", "50 um spot"), (0, "^", "EDS (Haimi 2021)"), (100, "v", "Arts 2019")):
    g = H[H.magnification == mag]; ax[1].plot(g.fill, g.measured_slope * g.z50, mk, ms=5, label=lab)
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("fraction of the gap filled, f = 2 t / H"); ax[1].set_ylabel("measured slope x half-thickness depth"); ax[1].set_title("(b) measured fronts against the step-front value"); ax[1].legend(frameon=False, loc="upper right", fontsize=6.5)
d = H[(H.H_nm == 500) & (H.magnification.isin([50, 0]))].sort_values("cycles"); cur = json.loads((FIG / "mechanism_slope_curves.json").read_text()); cmp_ = pd.read_csv(FIG / "mechanism_comparison.csv")
for model, c in (("gas", "C3"), ("gas+narrowing+spot", "C0")):
    g = cmp_[(cmp_.model == model) & cmp_["sample"].isin(d["sample"])].sort_values("cycles"); ax[2].plot(g.cycles, g.predicted_slope, "x" if "narrowing" in model else "+", color=c, ms=7, label=("with" if "narrowing" in model else "without") + " narrowing (predicted)")
ax[2].plot(d.cycles, d.measured_slope, "ko", ms=4, label="measured"); ax[2].set_xlabel("ALD cycles (gap 500 nm)"); ax[2].set_ylabel("front slope at half thickness"); ax[2].set_title("(c) cycle series, gas model at 300 Pa"); ax[2].legend(frameon=False); ax[2].set_ylim(0, 0.03)
fig.tight_layout(); fig.savefig(OUT / "fig10_narrowing_law.png"); plt.close(fig)

# ---- 3. held-out slopes ----
summ = json.loads((FIG / "mechanism_summary.json").read_text()); fig, ax = plt.subplots(1, 2, figsize=(8.6, 4.0), sharex=True, sharey=True)
for a, models, t in ((ax[0], ("hydraulic", "ballistic", "gas"), "(a) transport model alone"), (ax[1], ("hydraulic+narrowing+spot", "ballistic+narrowing+spot", "gas+narrowing+spot"), "(b) with channel narrowing and spot size")):
    for m, c, mk in zip(models, ("k", "C0", "C3"), ("o", "s", "^")):
        g = cmp_[(cmp_.model == m) & cmp_.predicted_slope.notna()]; q = summ["models"][m]
        a.plot(g.measured_slope, g.predicted_slope, mk, color=c, ms=4.5, mfc="none", label=f"{ {'ballistic': 'gas-free', 'gas': 'gas, 300 Pa'}.get(m.split('+')[0], m.split('+')[0]) }: s = {q['s']:.2g}, rms ln error {q['rms_ln']} (n = {q['n']})")
    a.plot([0.003, 0.07], [0.003, 0.07], "k-", lw=0.7); a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("measured front slope"); a.set_title(t); a.legend(frameon=False, loc="upper left")
    a.set_xticks([0.005, 0.01, 0.02, 0.05]); a.set_xticklabels(["0.005", "0.01", "0.02", "0.05"]); a.set_yticks([0.005, 0.01, 0.02, 0.05]); a.set_yticklabels(["0.005", "0.01", "0.02", "0.05"]); a.minorticks_off()
ax[0].set_ylabel("predicted front slope (s frozen from the calibration profiles)"); fig.tight_layout(); fig.savefig(OUT / "fig11_heldout.png"); plt.close(fig)

# ---- 4. identifiability: slope against s for representative samples, complete model ----
fig, ax = plt.subplots(1, 2, figsize=(8.8, 3.6)); per = cur["gas+narrowing+spot"]; Hh = H.set_index("sample"); SIG = summ["replicate_sigma_ln_slope_times_depth"]
for smp, lab, c in (("Yim9", "2000 nm, 500 cycles (5 % filled)", "C4"), ("Yim10", "500 nm, 250 cycles (10 %)", "C0"), ("Yim4", "500 nm, 500 cycles (19 %)", "C1"), ("Yim11", "500 nm, 1000 cycles (40 %)", "C2"), ("Arts", "Arts 2019, deep front (18 %)", "C5")):
    if smp not in per: continue
    ax[0].plot(per[smp]["s"], per[smp]["slope"], "-", color=c, lw=1.3, label=lab); m = Hh.loc[smp, "measured_slope"]; ax[0].axhspan(m * np.exp(-SIG), m * np.exp(SIG), color=c, alpha=0.10, lw=0)
ss = np.geomspace(0.003, 0.06, 50); ax[0].plot(ss, np.sqrt(ss / 13.9), "k:", lw=1, label="slope method: $\\sqrt{s/13.9}$")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("sticking coefficient s"); ax[0].set_ylabel("front slope (complete gas model, 300 Pa)"); ax[0].set_title("(a) what the slope still says about s"); ax[0].legend(frameon=False, fontsize=6.5)
ax[0].set_xticks([0.003, 0.01, 0.03, 0.06]); ax[0].set_xticklabels(["0.003", "0.01", "0.03", "0.06"]); ax[0].set_yticks([0.005, 0.01, 0.02, 0.05]); ax[0].set_yticklabels(["0.005", "0.01", "0.02", "0.05"]); ax[0].minorticks_off()
for m, c, mk in zip(mt.COMPLETE, ("k", "C0", "C3"), ("o", "s", "^")):
    idf = pd.DataFrame(summ["ceiling"][m]).T; ax[1].plot(idf.fill.astype(float), idf.measured_over_ceiling.astype(float), mk, ms=5, color=c, mfc="none", label={"ballistic": "gas-free", "gas": "gas, 300 Pa"}.get(m.split("+")[0], m.split("+")[0]))
idf = pd.DataFrame(summ["ceiling"]["gas+narrowing+spot"]).T
for smp, r in idf.iterrows():
    if smp in ("Yim9", "Yim10", "Yim11", "Yim7", "Arts", "HaimiA", "Yim6"): ax[1].annotate(smp.replace("Yim", "#"), (float(r.fill), float(r.measured_over_ceiling)), textcoords="offset points", xytext=(4, 3), fontsize=7)
ax[1].axhline(1.0, color="k", ls=":", lw=1); ax[1].set_xscale("log"); ax[1].set_xlabel("fraction of the gap filled"); ax[1].set_ylabel("measured slope / largest model slope"); ax[1].set_title("(b) measured slope over the model slope at the largest s"); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT / "fig12_identifiability.png"); plt.close(fig)
print(json.dumps(diag, indent=1)); print(pd.DataFrame(rows).round(4).to_string())
