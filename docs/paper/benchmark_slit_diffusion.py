"""Benchmark of the Monte Carlo's carrier-gas model against kinetic theory for a slit: Siewert and Valougeorgis, Eur. J. Mech.
B/Fluids 23, 645 (2004), McCormack model, rigid spheres, diffuse plates, diffusion problem (Tables 7-9), heavy species flow
rate U_C,2 at light-species mole fraction C = 0.9 for Ne-Ar and He-Xe, widths 2a = 0.01 to 100 viscosity mean free paths.
Conversion (derived from their Eqs. 1-22, 27a, 39a; checked against the binary diffusion coefficient at 2a = 100):
    D_slit = (H / 2) v2 U_C,2 / C,   v2 = (2 k T / m_heavy)^(1/2),   H = 2a l0,   l0 n d1^2 = 0.2270 (Ne-Ar), 0.1403 (He-Xe) at C = 0.9
The Monte Carlo (`pillarmc slitdiffusion`) follows a dilute heavy tracer in the pure light gas, so the comparison carries the
difference between 10 % heavy species and the tracer limit (their mean free path changes by about 10 % between the two).
Also listed: the hydraulic Knudsen diffusivity in series with the gas diffusivity (Bosanquet), the transport of the published
slope method. Writes figures/benchmark_slit_diffusion.csv.
  PYTHONPATH=src .venv/bin/python docs/paper/benchmark_slit_diffusion.py
"""
import subprocess
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]; EXE = ROOT / "lib/pillarmc/pillarmc"; KB, T, H = 1.380649e-23, 300.0, 0.5e-6
PAIRS = {"Ne-Ar": dict(m_light=20.183, m_heavy=39.948, d1=2.58, ratio=1.406, l0nd2=0.2270, U={0.01: 2.59225, 0.1: 1.50406, 0.5: 0.829410, 1: 0.577961, 2: 0.366857, 5: 0.174230, 10: 0.0922516, 100: 0.00969010}),
         "He-Xe": dict(m_light=4.0026, m_heavy=131.30, d1=2.18, ratio=2.226, l0nd2=0.1403, U={0.01: 2.82660, 0.1: 1.75516, 0.5: 1.09273, 1: 0.833313, 2: 0.592313, 5: 0.324434, 10: 0.183512, 100: 0.0204040})}
rows = []
for name, g in PAIRS.items():
    d12 = g["d1"] * (1 + g["ratio"]) / 2
    for a2, U in g["U"].items():
        n = a2 * g["l0nd2"] / (H * (g["d1"] * 1e-10) ** 2); p = n * KB * T; nmol = 6000 if a2 <= 1 else 3000
        out = subprocess.check_output([str(EXE), "slitdiffusion", f"{p:.6g}", f"{H * 1e6:g}", str(nmol), str(T), str(g["m_heavy"]), str(g["m_light"]), f"{d12:.5f}"], text=True).split()
        v = {out[i]: float(out[i + 1]) for i in range(0, len(out), 2)}; v2 = np.sqrt(2 * KB * T / (g["m_heavy"] * 1.66053906660e-27)); unit = H * v2 / 2
        rows.append(dict(pair=name, width_2a=a2, pressure_Pa=round(p, 2), kinetic_theory=round(U / 0.9, 4), monte_carlo=round(v["D_slit"] / unit, 4), mc_rel_err=v["rel_err"], hydraulic_bosanquet=round(v["D_bosanquet_hydraulic"] / unit, 4),
                         mc_over_theory=round(v["D_slit"] / unit / (U / 0.9), 3), hydraulic_over_theory=round(v["D_bosanquet_hydraulic"] / unit / (U / 0.9), 3)))
        print(rows[-1], flush=True)
d = pd.DataFrame(rows)
# Composition effect in the continuum limit. With 10 % heavy species the two diffusing fluxes obey Graham's law at diffuse walls
# (no net momentum to the wall, Kramers and Kistemaker 1943), which adds a counterflow that the tracer limit does not have:
# heavy flux = Fick flux x [1 - x2 (sqrt m2 - sqrt m1) / (x1 sqrt m1 + x2 sqrt m2)]. It vanishes in the free-molecular limit.
for name, g in PAIRS.items():
    r1, r2 = np.sqrt(g["m_light"]), np.sqrt(g["m_heavy"]); d.loc[d.pair == name, "tracer_over_mixture_continuum"] = round(1 / (1 - 0.1 * (r2 - r1) / (0.9 * r1 + 0.1 * r2)), 3)
d.to_csv(ROOT / "figures/benchmark_slit_diffusion.csv", index=False)

# The case of the experiments: TMA tracer in N2 at 300 C in a 500 nm slit (program defaults), against the hydraulic formula.
rows = []
for p in (30, 100, 300, 1000, 3000, 30000):
    out = subprocess.check_output([str(EXE), "slitdiffusion", str(p), "0.5", "4000"], text=True).split(); v = {out[i]: float(out[i + 1]) for i in range(0, len(out), 2)}
    rows.append(dict(pressure_Pa=p, gap_over_mean_free_path=round(v["H_over_lambda"], 4), D_monte_carlo_m2_s=v["D_slit"], mc_rel_err=v["rel_err"], D_hydraulic_bosanquet_m2_s=v["D_bosanquet_hydraulic"], mc_over_hydraulic=round(v["D_slit"] / v["D_bosanquet_hydraulic"], 3)))
pd.DataFrame(rows).to_csv(ROOT / "figures/slit_diffusion_tma_n2.csv", index=False); print(pd.DataFrame(rows).to_string(index=False))
