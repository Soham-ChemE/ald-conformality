"""Second chemistry in the chamber: HfO2 from TDMAH (Hf[N(CH3)2]4, M = 354.7 g/mol) on the 3D-NAND-like AR-100
wafer. Builds reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah as a copy of
reactor/showerhead_pulse_patterned_AR100_long (1 % precursor, 300 Pa, 573 K, 300 sccm, 20 s pulse, 22 s run,
frozen carrier flow in 0/U and 0/phi) with the precursor-specific entries changed:
  vbar  sqrt(8 R T / (pi M)) for M = 354.7 g/mol at 573 K,
  s     7.0e-3   initial sticking probability of the Hf amide precursor on the hydroxylated surface
                 (Reiter et al. 2024, Solid-State Electron. 211, 108816, Table 3, calibrated to Lam 3D-NAND
                 HfO2 step-coverage data; no direct TDMAH measurement exists, see paper/wp/wp4_geometry_chemistry.md),
  N0    2.8e18   Hf atoms per m2 per saturated cycle = GPC x rho / M(HfO2) x N_A with GPC 0.10 nm (Blaschke et
                 al. 2020: 1.08 A at 200 C, 0.93 A at 250 C; Gieraltowska et al. 2023: 1.2 A at 180-240 C) and
                 the bulk monoclinic density 9.68 g/cm3 (Delabie et al. 2005),
  D     hard-sphere Chapman-Enskog binary diffusivity in N2, scaled from the chamber's TMA value 8.4e-3 m2/s by
        the ratio of the two Chapman-Enskog estimates (same N2 diameter, so only the precursor's mass and
        diameter enter); TDMAH diameter from its liquid density by the rule used by Ylilammi et al. 2018 for TMA,
  f     Danckwerts inlet value fraction U / (U + D / delta), delta = 0.25 mm (first cell height).
c0 (3.79e20 /m3 = 1 % of 300 Pa at 573 K), d, AR, holeDensity, nCells unchanged.
  .venv/bin/python reactor/chemistry/make_tdmah_case.py            # build and print the parameters
  .venv/bin/python reactor/chemistry/make_tdmah_case.py --run      # build, run (one scalarTransportFoam call) and postprocess
  .venv/bin/python reactor/chemistry/make_tdmah_case.py --post     # postprocess an existing run
  .venv/bin/python reactor/chemistry/make_tdmah_case.py --run --d-sensitivity
      # also builds and runs ..._tdmah_Dtma: the same case with the TMA diffusivity (8.4e-3) and inlet fraction, to
      # separate the effect of the lower gas diffusivity from that of s, N0 and vbar
"""
import re, shutil, subprocess, sys, time, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "reactor/showerhead_pulse_patterned_AR100_long"
CASE = ROOT / "reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah"
CASE_DTMA = ROOT / "reactor/chemistry/showerhead_pulse_patterned_AR100_tdmah_Dtma"   # --d-sensitivity variant
LAUNCHER = Path.home() / "Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam"
sys.path.insert(0, str(ROOT / "coupling")); sys.path.insert(0, str(ROOT / "src"))

KB, R, NA, AMU = 1.380649e-23, 8.314462618, 6.02214076e23, 1.66053906660e-27
T, P, U, DELTA = 573.15, 300.0, 0.0501, 0.25e-3
C0 = 0.01 * P / (KB * T)                      # 3.79e20 /m3
WALL = 1.78e13 * np.pi * 100e-9 * 100 * 100e-9  # 55.9
D_TMA_CHAMBER = 8.4e-3                        # m2/s, constant/transportProperties of the TMA cases
# precursor data: M (g/mol), hard-sphere diameter (m), s, N0 (1/m2)
M_N2, SIGMA_N2 = 28.0134, 3.798e-10           # Poling, Prausnitz, O'Connell 2001, Appendix B (Lennard-Jones sigma)
M_TMA, SIGMA_TMA = 72.09, 5.91e-10            # Ylilammi et al. 2018, J. Appl. Phys. 123, 205301 (d_TMA = 591 pm)
M_TDMAH, RHO_TDMAH_LIQ = 354.70, 1.098e3      # kg/m3 at 25 C (vendor property compilation, ChemicalBook CAS 19782-68-4; no DOI)
S_TDMAH, GPC_TDMAH, RHO_HFO2, M_HFO2 = 7.0e-3, 0.10e-9, 9.68e3, 210.49


def vbar(M, T_K=T):
    return np.sqrt(8.0 * R * T_K / (np.pi * M * 1e-3))


def hard_sphere_diameter(M, rho):
    """d = (6 M / (pi rho N_A))^(1/3): the molecule's share of the condensed-phase volume as a sphere."""
    return (6.0 * (M * 1e-3) / (np.pi * rho * NA)) ** (1.0 / 3.0)


def chapman_enskog_hs(M1, M2, sigma1, sigma2):
    """Hard-sphere binary diffusivity D12 = 3/16 sqrt(2 pi k^3 T^3 / mu) / (p pi sigma12^2) (collision integral 1)."""
    mu = (M1 * M2 / (M1 + M2)) * AMU
    sigma12 = 0.5 * (sigma1 + sigma2)
    return 3.0 / 16.0 * np.sqrt(2.0 * np.pi * KB ** 3 * T ** 3 / mu) / (P * np.pi * sigma12 ** 2)


def parameters():
    sigma_tdmah = hard_sphere_diameter(M_TDMAH, RHO_TDMAH_LIQ)
    D_tma_ce = chapman_enskog_hs(M_TMA, M_N2, SIGMA_TMA, SIGMA_N2)
    D_tdmah_ce = chapman_enskog_hs(M_TDMAH, M_N2, sigma_tdmah, SIGMA_N2)
    D_tdmah = D_TMA_CHAMBER * D_tdmah_ce / D_tma_ce
    N0 = GPC_TDMAH * RHO_HFO2 / (M_HFO2 * 1e-3) * NA
    f = U / (U + D_tdmah / DELTA)
    return dict(vbar_tma=vbar(M_TMA), vbar_tdmah=vbar(M_TDMAH), sigma_tdmah=sigma_tdmah, D_tma_ce=D_tma_ce,
                D_tdmah_ce=D_tdmah_ce, D_tdmah=D_tdmah, N0=N0, f=f, s=S_TDMAH, c0=C0,
                rule_s=1.2 * (1 + WALL) * N0 / (U * C0), rule_tma_s=1.2 * (1 + WALL) * 4e18 / (U * C0),
                supply_time_s=N0 / (U * C0))


def build(par, case=None):
    case = case or CASE
    if case.exists(): shutil.rmtree(case)
    case.mkdir(parents=True)
    for sub in ("0", "constant", "system"): shutil.copytree(SRC / sub, case / sub)
    t = case / "0/T"; s = t.read_text()
    bc = (f"type patternedSurface; s {par['s']:g}; r 0; vbar {par['vbar_tdmah']:.0f}; N0 {par['N0']:.2e}; c0 3.79e20; "
          "d 100e-9; AR 100; holeDensity 1.78e13; nCells 200; value uniform 0;")
    s, n1 = re.subn(r"    wafer      \{ type patternedSurface;[^\n]*\n", "    wafer      { " + bc + " }\n", s)
    s, n2 = re.subn(r"uniformValueFraction constant [0-9.e-]+;", f"uniformValueFraction constant {par['f']:.6g};", s)
    s, n3 = re.subn(r"// T is the TMA mole fraction[^\n]*", "// T is the TDMAH mole fraction normalised to its showerhead value during the pulse (dimensionless carrier).", s)
    s, n4 = re.subn(r"    // exactly U\*c0\*A during the pulse\. As a mixed BC: value fraction f = U/\(U \+ D/delta\) = [^\n]*",
                    f"    // exactly U*c0*A during the pulse. As a mixed BC: value fraction f = U/(U + D/delta) = 0.0501/(0.0501 + {par['D_tdmah'] / DELTA:.3g}) = {par['f']:.4g}.", s)
    s, n5 = re.subn(r"    // TMA: vbar 410 m/s at 573 K; N0 4e18 sites/m2 \(about one TMA-saturated Al2O3 monolayer\); c0 = TMA number density\n    // in the feed at T = 1, taken as 3 Pa partial pressure at 573 K -> 3.79e20 1/m3.\n",
                    f"    // TDMAH (M 354.7): vbar {par['vbar_tdmah']:.0f} m/s at 573 K; s {par['s']:g} and N0 {par['N0']:.2e} Hf/m2 from the literature (see\n    // reactor/chemistry/make_tdmah_case.py); c0 = precursor number density in the feed at T = 1, 3 Pa at 573 K -> 3.79e20 1/m3.\n", s)
    assert n1 == n2 == n3 == n4 == n5 == 1, (n1, n2, n3, n4, n5); t.write_text(s)
    tp = case / "constant/transportProperties"; d = tp.read_text()
    d, m1 = re.subn(r"DT\s+[0-9.e-]+;", f"DT              {par['D_tdmah']:.5g};", d)
    d, m2 = re.subn(r"// TMA in N2 at 300 Pa, 573 K\. Chapman-Enskog estimate scaled from ~8e-6 m2/s at 1 atm, 300 K:\n// D = 8e-6 \* \(101325/300\) \* \(573/300\)\^1\.75 = 8\.4e-3 m2/s\. Treated as a parameter; see paper/assumptions\.md\.",
                    f"// TDMAH in N2 at 300 Pa, 573 K: the TMA value 8.4e-3 m2/s scaled by the ratio of hard-sphere Chapman-Enskog\n// estimates (TDMAH diameter {par['sigma_tdmah'] * 1e12:.0f} pm from its liquid density, TMA 591 pm): {par['D_tdmah']:.3g} m2/s. See reactor/chemistry/make_tdmah_case.py.", d)
    assert m1 == m2 == 1, (m1, m2); tp.write_text(d)
    return case


def run(case=None):
    case = case or CASE
    # one command through the launcher; the bundled bash segfaults when a shell is spawned inside the environment
    t0 = time.time()
    rc = subprocess.run([str(LAUNCHER), "-c", f"cd '{case}' && scalarTransportFoam > log.scalarTransportFoam 2>&1"]).returncode
    dt = time.time() - t0
    (case / "runtime.txt").write_text(f"scalarTransportFoam exit {rc}, wall time {dt:.0f} s\n")
    return rc, dt


def t90_of(t, series, level=0.9):
    i = np.argmax(series >= level); return float(t[i]) if series.max() >= level else float("nan")


# rule-only chemistries (no chamber run): s, GPC, density and the derived N0 with sources (DOIs in paper/wp/wp4_geometry_chemistry.md)
CHEMISTRIES = [
    # name, precursor, M (g/mol), s, GPC (m), film density (kg/m3), film formula mass per metal atom (g/mol), T (C), source note
    ("Al2O3", "TMA", 72.09, 5e-3, 0.10e-9, None, None, 300, "chamber runs: N0 4e18 (one TMA-saturated monolayer), s 5e-3 literature-range"),
    ("HfO2", "TDMAH", 354.70, 7.0e-3, 0.10e-9, 9.68e3, 210.49, 300, "s Reiter 2024 (Hf amide, calibrated); GPC Blaschke 2020 / Gieraltowska 2023; rho bulk Delabie 2005"),
    ("HfO2 at 250 C", "TDMAH", 354.70, 7.0e-3, 0.093e-9, 9.68e3, 210.49, 250, "GPC 0.93 A at 250 C (Blaschke 2020); U c0 at fixed sccm and fraction is temperature-independent"),
    ("TiN", "TiCl4", 189.68, 0.10, 0.030e-9, 4.23e3, 61.87, 325, "s Ylilammi 2018 (TiCl4/H2O, 110 C); GPC 0.30 A and rho 4.23 g/cm3 (XRR) Kinnunen 2017 at 325 C"),
    ("W", "WF6", 297.83, 0.40, 0.25e-9, 19.25e3, 183.84, 300, "s Elam 2001 (WF6 on Si-terminated W, 573 K); GPC 2.5 A Klaus 2000; rho from bcc a = 316.5 pm"),
]


def rule_table(par):
    """Rule-only predictions t = 1.2 (1 + 55.9) N0 / (U c0) at 1 % precursor in the baseline chamber, per chemistry."""
    import csv
    rows = []
    for name, prec, M, s_, gpc, rho, Mf, Tc, note in CHEMISTRIES:
        N0 = 4e18 if rho is None else gpc * rho / (Mf * 1e-3) * NA
        vb = vbar(M, Tc + 273.15)
        rows.append(dict(film=name, precursor=prec, M_g_mol=M, T_C=Tc, s=s_, GPC_nm=gpc * 1e9, N0_per_m2=N0, vbar_m_s=vb,
                         uptake_coeff_s_vbar_4_m_s=s_ * vb / 4.0, supply_time_N0_over_Uc0_s=N0 / (U * C0),
                         rule_AR100_wafer_s=1.2 * (1 + WALL) * N0 / (U * C0), rule_blanket_s=1.2 * N0 / (U * C0), note=note))
    with open(ROOT / "figures/chemistry_rule_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(f"{r['film']:14s} {r['precursor']:6s} s {r['s']:<7g} GPC {r['GPC_nm']:.3f} nm  N0 {r['N0_per_m2']:.2e}  vbar {r['vbar_m_s']:.0f} m/s  s vbar/4 {r['uptake_coeff_s_vbar_4_m_s']:.2f} m/s  "
              f"supply time {r['supply_time_N0_over_Uc0_s']:.3f} s  rule AR-100 wafer {r['rule_AR100_wafer_s']:.2f} s  (blanket supply x1.2: {r['rule_blanket_s']:.3f} s)")
    return rows


def postprocess(par):
    from patterned_lib import load_patterned
    out = {}
    cases = [("TMA", SRC), ("TDMAH", CASE)]
    if CASE_DTMA.exists() and any(d.name == "22" for d in CASE_DTMA.iterdir()): cases.append(("TDMAH_Dtma", CASE_DTMA))
    for name, case in cases:
        t, f, p, pulse = load_patterned(case)
        tb, tf = f["thetaBottom"], f["thetaFlat"]
        out[name] = dict(params=p, pulse_s=pulse, bottom_t90_centre=t90_of(t, tb[:, 0]), flat_t90_centre=t90_of(t, tf[:, 0]),
                         bottom_t90_worst=float(np.nanmax([t90_of(t, tb[:, j]) for j in range(149)])),
                         bottom_t90_edge=t90_of(t, tb[:, 148]), flat_t90_edge=t90_of(t, tf[:, 140]),
                         surface_density_rel_at_5s=float(f["value"][np.argmin(abs(t - 5.0)), 0]),
                         bottom_theta_end_of_pulse=float(tb[np.argmin(abs(t - pulse)), 0]))
    for k in out: out[k]["rule_s"] = par["rule_tma_s"] if k == "TMA" else par["rule_s"]
    for k in out: out[k]["ratio_to_rule"] = out[k]["bottom_t90_centre"] / out[k]["rule_s"]
    names = list(out)
    out["parameters"] = {k: float(v) for k, v in par.items()}
    (CASE / "results.json").write_text(json.dumps(out, indent=1, default=float))
    for k in names:
        o = out[k]
        print(f"{k}: bottom t90 centre {o['bottom_t90_centre']:.2f} s (worst {o['bottom_t90_worst']:.2f}, edge {o['bottom_t90_edge']:.2f}), flat t90 {o['flat_t90_centre']:.2f} s, "
              f"rule {o['rule_s']:.2f} s, ratio {o['ratio_to_rule']:.3f}; surface density at 5 s {o['surface_density_rel_at_5s']:.4f} c0")
    return out


if __name__ == "__main__":
    par = parameters()
    print("parameters:", {k: f"{v:.4g}" for k, v in par.items()}, flush=True)
    sens = "--d-sensitivity" in sys.argv
    par_d = {**par, "D_tdmah": D_TMA_CHAMBER, "f": U / (U + D_TMA_CHAMBER / DELTA)}
    if "--post" not in sys.argv:
        build(par); print("built", CASE.name, flush=True)
        if sens: build(par_d, CASE_DTMA); print("built", CASE_DTMA.name, "(TMA diffusivity)", flush=True)
    if "--run" in sys.argv:
        rc, dt = run(); print(f"run: exit {rc}, {dt / 60:.1f} min", flush=True)
        if sens: rc, dt = run(CASE_DTMA); print(f"run {CASE_DTMA.name}: exit {rc}, {dt / 60:.1f} min", flush=True)
    if "--run" in sys.argv or "--post" in sys.argv:
        postprocess(par)
    rule_table(par)
