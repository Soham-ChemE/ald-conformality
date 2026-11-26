"""Multi-cycle Monte Carlo of the PillarHall cavity with the gap narrowing as the film grows (Phase 1, check of the narrowing law
against particle transport). Cycles are followed in groups: each group runs one exposure in lib/pillarmc on bare (reset) walls
with the film of all earlier cycles on both plates (tfile=), and adds group_size x gpc x coverage(x) to the film. The per-cycle
exposure and the sticking coefficient are fixed, as in the experiments. Pillars keep their radius; the film is capped just
below half the gap.
  PYTHONPATH=src .venv/bin/python docs/paper/multicycle_mc.py run <case>       (hours; raw output in figures/pillar_mc_multicycle/, untracked)
  PYTHONPATH=src .venv/bin/python docs/paper/multicycle_mc.py compare          (figures/multicycle_mc_check.csv and .json)
The profile difference is taken where the law is defined by the simulated cavity (x <= depth x (1 - fill)^n); beyond it the law
would need coverage from outside the 400-gap cavity, and without gas the closed end collects a coverage rise that is not a front.
"""
import sys, subprocess, json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; EXE = ROOT / "lib/pillarmc/pillarmc"; OUT = ROOT / "figures/pillar_mc_multicycle"; OUT.mkdir(exist_ok=True)
GPC_UM = 1.0e-4                                                                 # 0.1 nm per cycle on each plate
# case: gap um, cavity depth in gaps, pillar spacing um, width um, extra options, s, per-cycle exposure Pa s, cycles, group size, checkpoints
CASES = {"gas300": (0.5, 400, 49, 90, ["p=300"], 0.0158, 4.6, 1000, 25, (250, 500, 1000)),
         "free": (0.5, 400, 49, 90, [], 0.025, 3.9, 1000, 25, (250, 500, 1000)),
         "h100_gas300": (0.1, 400, 49, 90, ["p=300"], 0.0158, 4.6, 500, 10, (250, 500)),
         "gas300_s004": (0.5, 400, 49, 90, ["p=300"], 0.004, 4.6, 500, 25, (250, 500))}
NSTEP, M = 50, 100000
def run(case):
    H, nb, a, W, extra, s, E, cycles, gs, cps = CASES[case]; t = np.zeros(nb); hist = []; d = OUT / case; d.mkdir(exist_ok=True)
    for g in range(cycles // gs):
        f = d / f"group{g:03d}.csv"; tf = d / f"film{g:03d}.txt"; np.savetxt(tf, t, fmt="%.8f")
        if not f.exists():
            subprocess.check_call([str(EXE), f"{s:g}", f"{E:g}", str(NSTEP), str(M), str(H), str(nb * H), str(a), "2", str(W), str(g + 1), "8", str(f), f"tfile={tf}"] + extra, stderr=subprocess.DEVNULL)
        th = np.loadtxt(f, delimiter=",", comments="#")[-1, 1:]; hist.append(th); t = np.minimum(t + gs * GPC_UM * th, 0.495 * H)
        if (g + 1) * gs in cps: np.savetxt(d / f"film_after_{(g + 1) * gs}.txt", t, fmt="%.8f")
        print(case, "cycles", (g + 1) * gs, "film at mouth nm", round(t[2] * 1000, 2), flush=True)
    np.save(d / "theta_groups.npy", np.array(hist))
def fs(z, y):
    y = y / y[1:5].mean(); i = int(np.argmax(y < 0.5)); x50 = z[i - 1] + (y[i - 1] - 0.5) / (y[i - 1] - y[i]); c = lambda l: (lambda k: z[k - 1] + (y[k - 1] - l) / (y[k - 1] - y[k]))(int(np.argmax(y < l))); return x50, 0.2 / (c(0.4) - c(0.6)), y
def compare():
    rows = []; out = {}
    for case, (H, nb, a, W, extra, s, E, cycles, gs, cps) in CASES.items():
        d = OUT / case
        if not (d / "group000.csv").exists(): continue
        z = np.arange(nb) + 0.5; th1 = np.loadtxt(d / "group000.csv", delimiter=",", comments="#")[-1, 1:]; x1, s1, _ = fs(z, th1); out[case] = dict(s=s, exposure=E, single=dict(x50=round(x1, 1), slope=round(s1, 5)), checkpoints={})
        for cp in cps:
            ff = d / f"film_after_{cp}.txt"
            if not ff.exists(): continue
            t = np.loadtxt(ff); xm, sm, ym = fs(z, t); fill = min(2 * cp * GPC_UM / H, 0.99)
            for n in (1.0, 0.92):
                law = np.mean([np.interp(z / f ** n, z, th1, right=0.0) for f in np.linspace(max(1 - fill, 0.02), 1, 41)], 0); xl, sl, yl = fs(z, law)
                rows.append(dict(case=case, s=s, cycles=cp, fill=round(fill, 3), gap_exponent=n, x50_mc=round(xm, 1), x50_law=round(xl, 1), x50_err_pct=round(100 * (xl / xm - 1), 2), slope_mc=round(sm, 5), slope_law=round(sl, 5), slope_err_pct=round(100 * (sl / sm - 1), 1), max_abs_diff=round(float(np.abs(yl - ym)[z <= nb * max(1 - fill, 0.02) ** n].max()), 4), slope_single=round(s1, 5)))
            out[case]["checkpoints"][str(cp)] = dict(z=z[::4].tolist(), mc=[round(float(v), 4) for v in ym[::4]])
    dd = pd.DataFrame(rows); dd.to_csv(ROOT / "figures/multicycle_mc_check.csv", index=False); (ROOT / "figures/multicycle_mc_check.json").write_text(json.dumps(out)); print(dd.to_string(index=False))
if __name__ == "__main__":
    run(sys.argv[2]) if sys.argv[1] == "run" else compare()
