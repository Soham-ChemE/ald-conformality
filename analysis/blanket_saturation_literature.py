"""Blanket-wafer TMA saturation of this chamber model against published Al2O3 ALD saturation data (work package 3, part A).

Two exposures are distinguished throughout, because the literature quotes both without always saying which:
  feed exposure     p_TMA(feed) x pulse time, what a recipe sets;
  surface exposure  the time integral of the TMA partial pressure at the surface, what the Langmuir kinetics see.
For our blanket run (reactor/showerhead_pulse_langmuir_p0.5: 300 Pa, 573 K, 300 sccm, 1 % TMA, s = 5e-3, N0 = 4e18)
both are computed here from the written coverage and face values; for each literature source the exposure the
experiment applied is tabulated as stated, and the planar 90 %-coverage exposure that the source's own sticking
coefficient or reaction cross-section implies is computed as
    E90 = ln(10) N0 / (s F1),   F1 = 1/sqrt(2 pi m kT)  (impingement rate per Pa),
or E90 = ln(10)/(sigma F1) when a cross-section is given (no site density needed). Exposure scales as 1/s, so the
sticking coefficient this model would need to reproduce a source is s_needed = s_ours x E90_ours / E_source.
Output: figures/blanket_saturation_literature.csv and .png.  Run: .venv/bin/python analysis/blanket_saturation_literature.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "reactor"))
from postprocess_langmuir import load

KB, AMU = 1.380649e-23, 1.66054e-27
M_TMA = 72.09 * AMU
def F1(T): return 1.0 / np.sqrt(2 * np.pi * M_TMA * KB * T)          # molecules / m2 / s / Pa
def E90_s(s, N0, T): return np.log(10) * N0 / (s * F1(T))            # Pa s
def E90_sigma(sigma_cm2, T): return np.log(10) / (sigma_cm2 * 1e-4 * F1(T))
L_PER_PAS = 1.0 / 1.33322e-4                                          # 1 L = 1e-6 Torr s = 1.333e-4 Pa s

# ---- our run -------------------------------------------------------------------------------------------------
case = ROOT / "reactor/showerhead_pulse_langmuir_p0.5"
t, theta, face, s_ours, vbar, c0, N0_ours, pulse = load(case)
T_OURS = 573.15
p_feed = c0 * KB * T_OURS
k90 = np.argmax(theta[:, 0] >= 0.9); t90 = t[k90]
E_feed = p_feed * t90
E_surf = np.trapezoid(face[:k90 + 1, 0] * p_feed, t[:k90 + 1])
kw = np.argmax(theta[:, 148] >= 0.9); t90w = t[kw]
E_feed_w = p_feed * t90w; E_surf_w = np.trapezoid(face[:kw + 1, 148] * p_feed, t[:kw + 1])
E90_kin_ours = E90_s(s_ours, N0_ours, T_OURS)

rows = [
    dict(source="this work, chamber model, blanket wafer (feed exposure, centre)", doi="this repository",
         T_C=300, reactor="axisymmetric showerhead, 300 Pa, 300 sccm N2, 1 % TMA", kind="feed p x t90",
         exposure_Pa_s=E_feed, s=s_ours, N0=N0_ours, note=f"p_feed {p_feed:.2f} Pa, t90 {t90:.3f} s (worst on wafer {E_feed_w:.2f} Pa s at {t90w:.3f} s); run reactor/showerhead_pulse_langmuir_p0.5"),
    dict(source="this work, chamber model, blanket wafer (surface exposure, centre)", doi="this repository",
         T_C=300, reactor="same run", kind="surface int p_s dt to theta = 0.9",
         exposure_Pa_s=E_surf, s=s_ours, N0=N0_ours, note=f"mean p_s/p_feed {E_surf/E_feed:.3f}; worst on wafer {E_surf_w:.3f} Pa s; equals the kinetic E90 {E90_kin_ours:.3f} Pa s"),
    dict(source="Ylilammi, Ylivaara, Puurunen 2018, applied (PillarHall, fitted entrance pressure x pulse)", doi="10.1063/1.5028178",
         T_C=300, reactor="Picosun R-150, 300 Pa, 150 sccm N2 per line, 0.1 s pulse", kind="applied: p_A0 x t_pulse",
         exposure_Pa_s=147.0 * 0.1, s=5.72e-3, N0=5e18, note="Table I: pA0 = 147 Pa (fitted), gamma = 0.00572, K = 219 /Pa; flat region of the channel saturated (constant thickness)"),
    dict(source="Ylilammi 2018, implied by fitted sticking coefficient", doi="10.1063/1.5028178",
         T_C=300, reactor="same", kind="kinetic E90 from s", exposure_Pa_s=E90_s(5.72e-3, 5e18, 573.15), s=5.72e-3, N0=5e18,
         note="site density q = 5e18 /m2 as in their worked example (not listed in Table I)"),
    dict(source="Arts et al. 2019, applied TMA dose, TMA-limited LHAR profile", doi="10.1116/1.5093620",
         T_C=220, reactor="Oxford OpAL, set-point 275 C (substrate ~220 C)", kind="applied: reactor p x t (stated)",
         exposure_Pa_s=190e-3 * 133.322, s=None, N0=None, note="~190 mTorr s TMA dose (text); planar surface saturated, penetration ~220 um"),
    dict(source="Arts et al. 2019, applied TMA dose, H2O-limited profile", doi="10.1116/1.5093620",
         T_C=220, reactor="same", kind="applied: reactor p x t (stated)", exposure_Pa_s=1200e-3 * 133.322, s=None, N0=None,
         note="~1200 mTorr s TMA dose (text)"),
    dict(source="Arts et al. 2019, implied by s0 = 2e-3 (upper)", doi="10.1116/1.5093620",
         T_C=220, reactor="same", kind="kinetic E90 from s", exposure_Pa_s=E90_s(2e-3, 4e18, 493.15), s=2e-3, N0=4e18,
         note="Table I: s0,TMA = (0.5-2)e-3 from the profile-front slope; N0 = 4e18 assumed here"),
    dict(source="Arts et al. 2019, implied by s0 = 5e-4 (lower)", doi="10.1116/1.5093620",
         T_C=220, reactor="same", kind="kinetic E90 from s", exposure_Pa_s=E90_s(5e-4, 4e18, 493.15), s=5e-4, N0=4e18, note="as above"),
    dict(source="Vandalon and Kessels 2019 (SFG), implied by reaction cross-section, evaluated at 300 C", doi="10.1021/acs.langmuir.9b01600",
         T_C=300, reactor="SFG vacuum chamber, 100 C measurement; 2017 JVST A 35, 05C313 (10.1116/1.4993597): TMA cross-section independent of temperature 100-300 C",
         kind="kinetic E90 from sigma", exposure_Pa_s=E90_sigma(6.5e-18, 573.15), s=None, N0=None,
         note="sigma_TMA = 6.5 +/- 0.6 e-18 cm2 (steady growth, Table 1); no site density needed; Arts Table I quotes the same work as s0 = (4 +/- 1)e-3 at 300 C"),
    dict(source="Vandalon and Kessels 2019, same cross-section evaluated at 100 C", doi="10.1021/acs.langmuir.9b01600",
         T_C=100, reactor="same", kind="kinetic E90 from sigma", exposure_Pa_s=E90_sigma(6.5e-18, 373.15), s=None, N0=None, note="F1 scales as T^-1/2"),
    dict(source="Gakis et al. 2018/2019, impingement exposure of the wafer during a 25 ms TMA valve opening (read from Fig. 12, 2018)", doi="10.1016/j.cherd.2018.02.031",
         T_C=300, reactor="Ultratech Fiji F200, 9.6 Pa Ar, 180 sccm, 200 mm wafer", kind="measured-model: surface impingement int p dt",
         exposure_Pa_s=0.0124 * np.sqrt(2 * np.pi * 0.07209 * 8.314 * 573.15), s=3.44e-3, N0=4.1e18,
         note="colour scale 0.0118-0.0130 mol/m2 (reading +/- 0.0002), x sqrt(2 pi M R T) = 46.45 Pa s per mol/m2 -> 0.55-0.60 Pa s; this exposure gave GPC 0.975-0.998 A/cycle (Fig. 10, 2019, digitised +/- 0.003) against 1.00 at 60 ms: >= 97 % saturated"),
    dict(source="Gakis et al. 2019, implied by fitted s0 = 0.004 x 86 % reaction probability", doi="10.1016/j.ces.2018.09.037",
         T_C=300, reactor="same", kind="kinetic E90 from s", exposure_Pa_s=E90_s(3.44e-3, 4.1e18, 573.15), s=3.44e-3, N0=4.1e18,
         note="N0 = 4.1e18 from their GPC 1.0 A/cycle at 3500 kg/m3 (one TMA per Al)"),
    dict(source="PillarHall record used in this repository (Zenodo 3974438), fitted mouth exposure", doi="10.5281/zenodo.3974438",
         T_C=300, reactor="Picosun, 500 cycles x 0.1 s", kind="applied: fitted p x t (this repo)", exposure_Pa_s=7.5, s=8.5e-3, N0=4e18,
         note="coupling/fit_pillarhall.py: 75 Pa x 0.1 s, s = 8.5e-3 (x/1.5)"),
    dict(source="PillarHall fit, implied by the fitted s", doi="10.5281/zenodo.3974438",
         T_C=300, reactor="same", kind="kinetic E90 from s", exposure_Pa_s=E90_s(8.5e-3, 4e18, 573.15), s=8.5e-3, N0=4e18, note=""),
    dict(source="Soto and Tysoe 1991 (UHV, room temperature) as summarised by Puurunen 2005", doi="10.1063/1.1940727",
         T_C=25, reactor="UHV, flat substrate", kind="kinetic E90 from s", exposure_Pa_s=E90_s(0.01, 4e18, 298.15), s=0.01, N0=4e18,
         note="reaction probability 0.01 for TMA on OH (Puurunen 2005 Section IV, citing JVST A 9, 2686); N0 = 4e18 assumed"),
    dict(source="Dillon et al. 1995 (porous alumina, FTIR) as summarised by Puurunen 2005", doi="10.1063/1.1940727",
         T_C=227, reactor="UHV FTIR cell, high-surface-area alumina membrane, 0.01 Torr", kind="kinetic E90 from first-order rate constant (transport-limited substrate)",
         exposure_Pa_s=np.log(10) / 0.0076 * 1.333, s=None, N0=None,
         note="ka = 0.0076 /s at 0.01 Torr and 227 C (Puurunen 2005, citing Surf. Sci. 322, 230): E90 = ln10/ka x p = 403 Pa s; porous sample, so this is diffusion into the membrane, not planar kinetics; effective s would be 1.8e-6"),
]
df = pd.DataFrame(rows)
df["exposure_L"] = df.exposure_Pa_s * L_PER_PAS
# s this model would need to reproduce each kinetic value (exposure ~ 1/s); applied exposures are upper bounds and get none
df["s_needed_for_this_model"] = np.where(df.kind.str.startswith("kinetic"), s_ours * E_surf / df.exposure_Pa_s, np.nan)
df.to_csv(ROOT / "figures/blanket_saturation_literature.csv", index=False)
print(f"ours: s {s_ours:g}, N0 {N0_ours:.3g}, p_feed {p_feed:.3f} Pa, t90 {t90:.3f} s -> feed exposure {E_feed:.3f} Pa s = {E_feed*L_PER_PAS:.0f} L; "
      f"surface exposure {E_surf:.3f} Pa s = {E_surf*L_PER_PAS:.0f} L; kinetic E90 {E90_kin_ours:.3f} Pa s")
print(df[["source", "T_C", "kind", "exposure_Pa_s", "exposure_L", "s_needed_for_this_model"]].to_string())

# ---- figure: saturating exposure against temperature ----------------------------------------------------------
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
fig, ax = plt.subplots(figsize=(8.2, 5.2))
app = df[df.kind.str.startswith("applied")]
kin = df[df.kind.str.startswith("kinetic") & ~df.source.str.startswith("Dillon") & ~df.source.str.startswith("PillarHall") & ~df.source.str.startswith("this work")]
ax.scatter(app.T_C, app.exposure_Pa_s, marker="s", s=70, facecolor="none", edgecolor=BLUE, lw=1.8, label="exposure the experiment applied (feed or reactor p x t)")
ax.scatter(kin.T_C, kin.exposure_Pa_s, marker="o", s=60, color=ORANGE, label="planar 90 % exposure implied by the reported s or sigma")
d = df[df.source.str.startswith("Dillon")]
ax.scatter(d.T_C, d.exposure_Pa_s, marker="v", s=70, color=YELLOW, label="Dillon 1995 (porous membrane, transport-limited)")
g = df[df.kind.str.startswith("measured-model")]
ax.scatter(g.T_C, g.exposure_Pa_s, marker="D", s=60, color=AQUA, label="Gakis 2018/2019: impingement exposure that saturated to >= 97 %")
ph = df[df.source.str.startswith("PillarHall")]
ax.scatter(ph.T_C, ph.exposure_Pa_s, marker="x", s=60, color="#555555", label="PillarHall fit in this repository (applied / implied)")
ax.scatter([300], [E_feed], marker="*", s=220, color=BLUE, zorder=5, label=f"this chamber model: feed exposure {E_feed:.2f} Pa s ({E_feed*L_PER_PAS:.0f} L)")
ax.scatter([300], [E_surf], marker="*", s=220, color=ORANGE, zorder=5, label=f"this chamber model: surface exposure {E_surf:.3f} Pa s ({E_surf*L_PER_PAS:.0f} L)")
Ts = np.linspace(20, 320, 50)
for s_, ls in ((1e-3, ":"), (5e-3, "--"), (1e-2, "-.")):
    ax.plot(Ts, E90_s(s_, 4e18, Ts + 273.15), ls, color="#888888", lw=1)
    ax.text(130, E90_s(s_, 4e18, 403.15) * 1.12, f"Langmuir E90, s = {s_:g}", fontsize=7, color="#555555", va="bottom")
# labels
for _, r in df.iterrows():
    if r.source.startswith("this work"): continue
    lab = r.source.split(",")[0].split(" (")[0].replace(" et al.", "").replace(" and Kessels", "").replace("Ylilammi 2018", "Ylilammi")
    if lab.startswith("PillarHall record"): lab = "PillarHall (this repo), applied"
    if lab.startswith("Gakis et al. 2018/2019") or lab.startswith("Gakis 2018"): lab = "Gakis 2018/2019, 25 ms"
    dx, dy = (5, 3) if r.T_C < 290 else (7, -9 if "Gakis 2018" in lab else 3)
    ax.annotate(lab, (r.T_C, r.exposure_Pa_s), xytext=(dx, dy), textcoords="offset points", fontsize=6.5, color="#444444")
ax.set_yscale("log"); ax.set_xlim(0, 360); ax.set_ylim(0.02, 1000)
ax.set_xlabel("temperature (C)"); ax.set_ylabel("TMA exposure for saturation (Pa s)")
ax2 = ax.secondary_yaxis("right", functions=(lambda e: e * L_PER_PAS, lambda l: l / L_PER_PAS)); ax2.set_ylabel("Langmuir")
ax.set_title("Blanket Al2O3 ALD: saturating TMA exposure, this model vs published data (grey: Langmuir E90 at N0 = 4e18 /m2)", fontsize=9)
ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=6.5, loc="upper left")
fig.tight_layout(); fig.savefig(ROOT / "figures/blanket_saturation_literature.png", dpi=150)
print("saved", ROOT / "figures/blanket_saturation_literature.png")
