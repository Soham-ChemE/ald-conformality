"""Build docs/assets/data-mech.js (window.ALD_MECH) for the hero's slit steps and the page's "Inside the data" section.
Sources: the single-exposure coverage frames of the pillar-cavity Monte Carlo (figures/pillar_mc/, untracked; listed with hashes in
figures/pillar_mc_manifest.csv), the 1D hydraulic kernel (docs/paper/mechanism_test.py), the two measured profiles used for the
calibration, and the summary files of the mechanism test.
  PYTHONPATH=src .venv/bin/python docs/build_mech_assets.py
"""
import sys, json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "docs/paper"))
import mechanism_test as mt
S = (0.004, 0.0063, 0.01, 0.0158, 0.025, 0.04)                 # sticking values simulated for every kernel below
Z = np.arange(0, 404, 4.0)                                     # depth in gap heights
E = np.round(6.0 * (np.arange(1, 13) / 12.0) ** 2, 4)          # exposures in Pa s (every kernel's ladder reaches 6)
KERNELS = {"free": ("walls90_pillars49", "real cavity, no carrier gas (free molecular)", 0, 1.0), "gas100": ("gas100", "real cavity, N2 at 100 Pa", 100, 0.92),
           "gas300": ("gas300", "real cavity, N2 at 300 Pa", 300, 0.92), "gas1000": ("gas1000", "real cavity, N2 at 1000 Pa", 1000, 0.92),
           "hyd": ("hyd", "1D hydraulic-diameter diffusion model", None, 1.0)}
out = {"z": Z.tolist(), "E": E.tolist(), "s": list(S), "kernels": {}, "gap_nm": 500, "gpc_nm": 0.1,
       "source": "figures/pillar_mc/<case>_s<s>.csv (figures/pillar_mc_manifest.csv); hydraulic: docs/paper/mechanism_test.py"}
for key, (case, label, p, n) in KERNELS.items():
    K = mt.kernel(case); th = []
    for s in S:
        Ek, TH, z = K[s]; assert Ek[-1] >= E[-1] - 1e-9, (case, s, Ek[-1])
        th.append([[round(float(v), 3) for v in np.interp(Z, z, mt.at(Ek, TH, e))] for e in E])
    out["kernels"][key] = dict(label=label, pressure_Pa=p, gap_exponent=n, case=case, theta=th)
x, y, fill = mt.measured_profile("zenodo"); out["measured"] = {"zenodo": dict(label="Zenodo 3974438, mean of 8 line scans (500 cycles)", x=[round(float(v), 1) for v in x], y=[round(float(v), 4) for v in y], fill=round(float(fill), 3), file="data/zenodo_3974438")}
x, y, fill = mt.measured_profile("ylilammi_al2o3"); out["measured"]["ylilammi"] = dict(label="Ylilammi et al. 2018, Fig. 6 (500 cycles)", x=[round(float(v), 1) for v in x], y=[round(float(v), 4) for v in y], fill=round(float(fill), 3), file="data/lhar_profiles/ylilammi2018_fig6_al2o3.csv")
summ = json.loads((ROOT / "figures/mechanism_summary.json").read_text()); tr = json.loads((ROOT / "figures/mechanism_transport.json").read_text())
out["summary"] = dict(rms_common={m: v["rms_ln_common"] for m, v in summ["models"].items()}, s_calibrated={m: v["s"] for m, v in summ["models"].items()},
                      replicate_sigma=summ["replicate_sigma_ln_slope_times_depth"], transport={k: dict(tail=v["tail_plus_60"], exposure=v["exposure_Pa_s"]) for k, v in tr["variants"].items()},
                      files=["figures/mechanism_summary.json", "figures/mechanism_transport.json"])
js = "window.ALD_MECH=" + json.dumps(out, separators=(",", ":")) + ";\n"
(ROOT / "docs/assets/data-mech.js").write_text(js); print("wrote docs/assets/data-mech.js", len(js), "bytes")
