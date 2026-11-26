"""POST HOC: does the penetration scale more weakly with the gap at higher carrier pressure (which would raise the narrowing
limit)? Pairs of runs at 400 and 500 nm gaps, same exposure, s = 0.02, N2 at 1400 and 3000 Pa (lib/pillarmc, frames in
figures/pillar_mc/scale_p*_h*_s0.02.csv). Writes docs/blindtest/posthoc/gap_scaling_pressure.csv."""
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[3]; rows = []
def z50s(f, H):
    d = np.loadtxt(f, delimiter=",", comments="#"); E, TH = d[:, 0], d[:, 1:]; z = (np.arange(TH.shape[1]) + 0.5) * H; out = []
    for th in TH:
        i = int(np.argmax(th < 0.5)); out.append(np.nan if i == 0 else z[i - 1] + (th[i - 1] - 0.5) / (th[i - 1] - th[i]) * H)
    return E, np.array(out)
for p in (1400, 3000):
    E1, a = z50s(ROOT / f"figures/pillar_mc/scale_p{p}_h0.5_s0.02.csv", 0.5); E2, b = z50s(ROOT / f"figures/pillar_mc/scale_p{p}_h0.4_s0.02.csv", 0.4)
    for e in (4, 6, 8, 10, 12):
        u, v = np.interp(e, E1, a), np.interp(e, E2, b); rows.append(dict(pressure_Pa=p, exposure_Pa_s=e, x50_um_gap500=round(float(u), 1), x50_um_gap400=round(float(v), 1), exponent=round(float(np.log(v / u) / np.log(0.8)), 3)))
d = pd.DataFrame(rows); d.to_csv(ROOT / "docs/blindtest/posthoc/gap_scaling_pressure.csv", index=False); print(d.groupby("pressure_Pa").exponent.agg(["mean", "std"]).round(3))
