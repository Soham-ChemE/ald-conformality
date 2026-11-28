"""Two robustness checks of the profile-by-profile fits of profile_fits.py (16 TMA/water fronts, s and exposure free).

  gap       the 500 nm chips refitted with the measured gap in place of the nominal one: 460 nm, the fitted mean air gap of
            a nominal 500 nm PillarHall chip (Ikonen et al., Phys. Status Solidi A 223, e202501022, 2026; smallest 454 nm,
            largest 480 nm). Depths and the measurement spot are rescaled to the new gap and the filling becomes 2 t/460 nm.
            The 2000 nm sample keeps its nominal gap. Hydraulic and carrier-gas models, both fit windows.
  reading   systematic reading errors of the digitised curves: the plateau level taken 3 % higher or lower, and the channel
            entrance shifted by 5 gap heights (2.5 um) either way. Hydraulic model, window to 0.3 of the plateau.
For each variant: the slope of ln s against the filled fraction over the 16 fronts, without and with channel narrowing.
  .venv/bin/python docs/paper/robustness_fits.py gap        -> figures/robustness_gap_fits.csv, figures/robustness_gap.csv
  .venv/bin/python docs/paper/robustness_fits.py reading    -> figures/robustness_reading_fits.csv, figures/robustness_reading.csv
  .venv/bin/python docs/paper/robustness_fits.py censored   -> figures/robustness_censored.csv (from the two fit files)
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(Path(__file__).parent))
import profile_fits as pf
import mechanism_test as mt
F = ROOT / "figures"; mode = sys.argv[1]
if mode == "censored":       # censored (Tobit) slopes of the fits with narrowing, a fit without an upper limit entering as a lower bound at its 16th percentile
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()): import profile_fits_stats as st
    out = []
    for f in ("gap", "reading"):
        D = pd.read_csv(F / f"robustness_{f}_fits.csv")
        for (v, win, fam, nar), g in D.groupby(["variant", "window", "transport", "narrowing"], sort=False):
            if not nar: continue
            o = st.ols(g.fill.values, np.log(g.s_fit.values)); c = st.tobit(g.fill.values, np.log(g.s_fit.values), np.log(g.s_p16.values), np.log(g.s_p84.values), g.no_upper_limit.values.astype(bool), np.zeros(len(g), bool))
            out.append(dict(check=f, variant=v, window=win, transport=fam, count=len(g), no_upper_limit=int(g.no_upper_limit.sum()), dlns_df=round(o[0], 2), se=round(o[1], 2), dlns_df_censored=round(c[0], 2), se_censored=round(c[1], 2)))
    C = pd.DataFrame(out); C.to_csv(F / "robustness_censored.csv", index=False); pd.set_option("display.width", 200); print(C.to_string(index=False)); sys.exit(0)
H = mt.heldout().set_index("sample"); profiles = []
for f in sorted((ROOT / "data/blind_revealed_v2").glob("yim2020_s*.csv")):
    k = int(f.stem.split("_s")[1]); d = pd.read_csv(f); name = f"Yim{k}"
    if name in H.index: profiles.append(dict(id=f"opt_s{k:02d}", sample=name, method="reflectometry", z=d.x_over_H.values, y=d.thickness_over_plateau.values, spot_um=mt.SPOT_UM[int(H.loc[name, "magnification"])]))
for c in json.loads((ROOT / "docs/blindtest/v2/posthoc/eds_cycle_series_curves.json").read_text()):
    profiles.append(dict(id=f"eds_{c['cycles']}_{c['channel']}", sample={250: "Yim10", 500: "Yim2", 1000: "Yim11"}[c["cycles"]], method="EDS", z=np.array(c["z"]), y=np.array(c["y"]), spot_um=1.0))
def kernel(fam, r, fill, z50):
    if r.H_nm != 500: return mt.kernel_for(fam, r, fill)
    if fam == "hydraulic": return "hyd_deep"
    return "gas300_deep" if z50 / (1 - fill / 2) > 250 else mt.kernel_for(fam, r, fill)
def run(variants, fams, windows):
    rows = []
    for win in windows:
        pf.CUT = win
        for P in profiles:
            r = H.loc[P["sample"]]
            for vname, gap, yscale, shift in variants:
                g = gap if r.H_nm == 500 else r.H_nm; k = r.H_nm / g; fill = min(2 * r.film_nm / g, 1.0); z = P["z"] * k + shift; y = P["y"] / yscale; spot = P["spot_um"] * 1000 / g
                for fam in fams:
                    K = pf.kern(kernel(fam, r, fill, r.z50 * k))
                    for narrow in (False, True):
                        q = pf.fit(z, y, fill, spot, K, mt.GAP_EXPONENT[fam], narrow); rows.append(dict(variant=vname, window=win, id=P["id"], sample=P["sample"], method=P["method"], fill=round(fill, 3), transport=fam, narrowing=narrow, s_fit=q["s_fit"], s_p16=q["s_p16"], s_p84=q["s_p84"], rms=q["rms"], no_upper_limit=q["no_upper_limit"]))
            print(win, P["id"], "done", flush=True)
    D = pd.DataFrame(rows); out = []
    for (v, win, fam, nar), g in D.groupby(["variant", "window", "transport", "narrowing"], sort=False):
        for name, q in (("all 16 fronts", g), ("10 reflectometry fronts", g[g.method == "reflectometry"])):
            c, cov = np.polyfit(q.fill, np.log(q.s_fit), 1, cov=True)
            out.append(dict(variant=v, window=win, transport=fam, narrowing=nar, fronts=name, count=len(q), dlns_df=round(float(c[0]), 2), se=round(float(np.sqrt(cov[0, 0])), 2), sd_ln_s=round(float(np.log(q.s_fit).std(ddof=1)), 3), s_geomean=float(f"{np.exp(np.log(q.s_fit).mean()):.3g}"), mean_rms=round(float(q.rms.mean()), 4), no_upper_limit=int(q.no_upper_limit.sum())))
    return D, pd.DataFrame(out)
if mode == "gap":
    D, S = run([("nominal gap 500 nm", 500.0, 1.0, 0.0), ("measured gap 460 nm", 460.0, 1.0, 0.0)], ("hydraulic", "gas"), (0.3, 0.1))
else:
    D, S = run([("as read", 500.0, 1.0, 0.0), ("plateau 3 % higher", 500.0, 1.03, 0.0), ("plateau 3 % lower", 500.0, 0.97, 0.0), ("entrance 5 gap heights deeper", 500.0, 1.0, -5.0), ("entrance 5 gap heights shallower", 500.0, 1.0, 5.0)], ("hydraulic",), (0.3,))
D.to_csv(F / f"robustness_{mode}_fits.csv", index=False); S.to_csv(F / f"robustness_{mode}.csv", index=False); pd.set_option("display.width", 220); print(S[S.fronts == "all 16 fronts"].to_string(index=False))
