"""Two-species feature kinetics with byproduct re-adsorption (aldconf.multispecies) on the AR-100 channel hole.

Physics first: at 100 to 1,000 Pa the mean free path is 10 to 120 um, so in a 100 nm hole the Knudsen number is 100
to 1,000 and a byproduct cannot block the precursor in the gas phase; each species diffuses independently with its own
Knudsen diffusivity (D_K ~ 1/sqrt(M)). A byproduct acts only at the surface: re-adsorption on free sites (blocking
until it desorbs) or etching. Methane (TMA/H2O) does not re-adsorb at these temperatures, so case (a) is the
single-species result; HCl (TiCl4/NH3) does, and cases (b)-(d) quantify what re-adsorption, etching and chamber
accumulation of HCl do to the saturating pulse of a 100 nm x 10 um hole.

Cases (every run also appears in figures/multispecies.csv; the figure is figures/multispecies.png):
  (a) TMA with CH4 (nu 3, s_B 0): identical to the single-species solver; and a hypothetical bound s_B 1e-3, k_B 1e3 1/s.
  (b) TiCl4 with HCl: nu 0.5 (HCl released in the TiCl4 half-step, Cl:Ti 4 -> 3.5, Tiznado and Zaera 2006) and nu 4
      (whole-cycle stoichiometry, Ritala et al. 1995, as a bound); s_B 1e-3 to 0.4 (Elam et al. 2000 reactive sticking
      1e-3 on alpha-Al2O3(0001) at 298 K to the Knudsen-cell initial uptake 0.38 on alumina, Lostier et al. 2026; no
      value exists on TiN at 600 K); k_B 0.1 to 1e3 1/s (the 300-650 K desorption span of Elam et al. 2000 read with
      a 1e13 1/s prefactor at 573 K); k_etch 0 and bounds 1e-4, 1e-3 per collision (HCl does not etch TiO2, Leem et al.
      2014; Ti removal as TiCl3/TiCl4 by surface Cl in the TiN kinetic Monte Carlo of Kim et al. 2022, no rate).
  (c) time to 90 % bottom coverage and the coverage profile at the rule's pulse, at the chamber's centre-face density
      history (reactor/showerhead_pulse_patterned_AR100_long, the TMA run; for TiCl4 this history is a stand-in, since
      no TiCl4 chamber run exists) and at a constant 0.3 Pa mouth density.
  (d) byproduct mouth density 0, 1, 3, 10 % of the precursor mouth density (chamber accumulation).
Also: the substep convergence of the TiCl4 case (the byproduct rates set the substep, not the precursor's).
Light: 1D models only; a few minutes on one core.
"""
import sys, csv, time
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "coupling"))
from aldconf import Hole, Precursor, Process, simulate, simulate_driven
from aldconf.multispecies import Byproduct, simulate_two_species, simulate_two_species_driven, knudsen_number
from aldconf.transport import mean_speed, mean_free_path, knudsen_diffusivity
from aldconf.recess import time_to_coverage
from patterned_lib import load_patterned

CASE = ROOT / "reactor/showerhead_pulse_patterned_AR100_long"
D, AR, T = 100e-9, 100, 573.15
HOLE = Hole(D, AR * D)
TMA = Precursor("TMA", 72.09, sticking=5e-3, site_density=4e18, sigma=5.91e-10)
TICL4 = Precursor("TiCl4", 189.68, sticking=0.10, site_density=1.24e18, sigma=7.04e-10)
RULE = {"TMA": 14.38, "TiCl4": 4.44}            # rule pulses on the AR-100 wafer at 1 % (figures/chemistry_rule_table.csv)
P_CONST = 0.3                                     # Pa, constant-mouth benchmark (as the Machball check)
DT = 1e-3                                         # substep for every two-species run (see the convergence block)
DEPTHS = np.linspace(0, 1, 41)


def t90(res):
    return res.time_to_bottom_coverage(0.9)


def const_run(pr, bp, t_dose, dt=DT, depths=DEPTHS):
    n = int(round(t_dose / dt))
    proc = Process(HOLE, pr, T=T, p_partial=P_CONST, t_dose=t_dose, n_cells=200, n_substeps=n, narrowing=False)
    return simulate_two_species(proc, bp, record_depths=depths, n_samples=min(n, 2000))


def driven_run(pr, bp, t, n0, dt=DT, depths=DEPTHS):
    sub = max(1, int(round((t[1] - t[0]) / dt)))
    proc = Process(HOLE, pr, T=T, p_partial=1.0, n_cells=200)
    return simulate_two_species_driven(proc, bp, t, n0, substeps_per_sample=sub, record_depths=depths)


def profile_at(res, tt):
    i = int(np.argmin(np.abs(res.t - tt)))
    return res.hist_P[i], res.hist_B[i], float(res.t[i])


def main():
    t0 = time.time()
    rows = []; profiles = {}
    # --- the physics numbers ---------------------------------------------------------------------------------------
    print("Knudsen numbers of a 100 nm feature (N2 carrier, sigma 3.8 A, 573 K):")
    for p in (100.0, 300.0, 1000.0):
        print(f"  {p:6.0f} Pa: mean free path {mean_free_path(T, p, 3.8e-10) * 1e6:6.1f} um, Kn = {knudsen_number(p, T, D):7.0f}")
    for a, b in ((TMA, Byproduct('CH4', 16.04)), (TICL4, Byproduct('HCl', 36.46))):
        print(f"  D_K({b.name}) / D_K({a.name}) = {mean_speed(T, b.mass_amu) / mean_speed(T, a.mass_amu):.2f} (vbar {mean_speed(T, b.mass_amu):.0f} vs {mean_speed(T, a.mass_amu):.0f} m/s)")

    # --- chamber drive ------------------------------------------------------------------------------------------------
    t, f, params, pulse = load_patterned(CASE)
    n0_series = f["value"][:, 0] * params["c0"]
    chamber_t90 = time_to_coverage(t, f["thetaBottom"][:, 0])
    print(f"chamber (two-way) bottom t90 at the centre: {chamber_t90:.2f} s; drive: {len(t)} samples, pulse {pulse:g} s")

    def add(case, pr, bp, drive, res, t90_ref=None, dt=DT, note=""):
        tb = t90(res); tm = res.time_to_mouth_coverage(0.9)
        thP_rule, thB_rule, _ = profile_at(res, RULE[pr.name])
        rows.append(dict(case=case, precursor=pr.name, byproduct=bp.name, nu=bp.nu, s_B=bp.sticking, k_B=bp.k_des,
                         k_etch=bp.k_etch, mouth_fraction=bp.mouth_fraction, drive=drive, dt_s=dt,
                         t90_bottom_s=tb, t90_mouth_s=tm,
                         delta_t90_percent=(np.nan if t90_ref is None or not np.isfinite(tb) else 100 * (tb / t90_ref - 1)),
                         theta_bottom_at_rule=thP_rule[-1], theta_blocked_max_at_rule=thB_rule.max(),
                         theta_blocked_bottom_end=res.theta_B[-1], theta_P_bottom_end=res.theta_P[-1],
                         escaped_fraction=(res.balance["escaped"] / (res.balance["produced"] + res.balance["desorbed"])
                                           if res.balance["produced"] > 0 else np.nan),
                         balance_closure=res.balance["closure"], note=note))
        return tb

    # --- (a) TMA / methane ----------------------------------------------------------------------------------------------
    print("\n(a) TMA with methane, AR-100 hole")
    for drive in ("constant 0.3 Pa", "chamber centre face"):
        if drive.startswith("constant"):
            proc1 = Process(HOLE, TMA, T=T, p_partial=P_CONST, t_dose=20.0, n_cells=200, n_substeps=int(20.0 / DT), narrowing=False)
            # single-species reference history on the same grid through simulate_driven (constant series)
            x, th1, ht1, h1 = simulate_driven(proc1, np.linspace(0, 20.0, 2001), np.full(2001, P_CONST / (1.380649e-23 * T)), substeps_per_sample=int(round(0.01 / DT)), record_depths=DEPTHS)
            run = lambda bp: const_run(TMA, bp, 20.0)
        else:
            proc1 = Process(HOLE, TMA, T=T, p_partial=1.0, n_cells=200)
            x, th1, ht1, h1 = simulate_driven(proc1, t, n0_series, substeps_per_sample=int(round((t[1] - t[0]) / DT)), record_depths=DEPTHS)
            run = lambda bp: driven_run(TMA, bp, t, n0_series)
        t_single = time_to_coverage(ht1, h1[:, -1])
        res0 = run(Byproduct("CH4", 16.04, nu=3.0, sticking=0.0))
        diff = np.abs(res0.theta_P - th1).max()
        tb0 = add("a", TMA, res0.byproduct, drive, res0, t_single, note=f"max |theta - single-species| {diff:.1e}")
        resb = run(Byproduct("CH4", 16.04, nu=3.0, sticking=1e-3, k_des=1e3))
        tbb = add("a", TMA, resb.byproduct, drive, resb, t_single, note="hypothetical bound, no literature re-adsorption of CH4")
        print(f"  {drive}: single-species t90 {t_single:.3f} s; two-species s_B = 0: {tb0:.3f} s (max |dtheta| {diff:.1e}); bound s_B 1e-3, k_B 1e3: {tbb:.3f} s ({100 * (tbb / t_single - 1):+.2f} %), escaped {resb.balance['escaped'] / resb.balance['produced']:.4f}")
        profiles[("a", drive)] = (res0, resb)

    # --- (b)+(c) TiCl4 / HCl ----------------------------------------------------------------------------------------------
    print("\n(b) TiCl4 with HCl, AR-100 hole, constant 0.3 Pa mouth")
    T_DOSE = 30.0
    base_c = const_run(TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.0), T_DOSE)
    t_base_c = add("b", TICL4, base_c.byproduct, "constant 0.3 Pa", base_c, note="no re-adsorption (single-species)")
    print(f"  no HCl action: bottom t90 {t_base_c:.3f} s")
    grid_s = (1e-3, 1e-2, 0.1, 0.4); grid_k = (0.1, 1.0, 10.0, 100.0, 1000.0)
    tb_grid = {}
    for nu in (0.5, 4.0):
        for sB in grid_s:
            line = []
            for kB in grid_k:
                res = const_run(TICL4, Byproduct("HCl", 36.46, nu=nu, sticking=sB, k_des=kB), T_DOSE)
                tb = add("b", TICL4, res.byproduct, "constant 0.3 Pa", res, t_base_c)
                tb_grid[(nu, sB, kB)] = tb; line.append(tb)
            print(f"  nu {nu:3.1f}, s_B {sB:5.0e}: t90 = " + ", ".join(f"{v:6.2f}" for v in line) + f" s for k_B = {grid_k} 1/s")
    print("  with etching (k_etch per HCl collision on a Ti site):")
    for nu in (0.5, 4.0):
        for sB in (1e-2, 0.1):
            for kB in (1.0, 100.0):
                for ke in (1e-4, 1e-3):
                    res = const_run(TICL4, Byproduct("HCl", 36.46, nu=nu, sticking=sB, k_des=kB, k_etch=ke), T_DOSE)
                    tb = add("b-etch", TICL4, res.byproduct, "constant 0.3 Pa", res, t_base_c, note="etch bound")
                    print(f"  nu {nu:3.1f}, s_B {sB:5.0e}, k_B {kB:6.1f}, k_etch {ke:5.0e}: t90 {tb:7.2f} s (blocking only {tb_grid[(nu, sB, kB)]:.2f} s); bottom theta_P at 30 s {res.theta_P[-1]:.3f}, blocked {res.theta_B[-1]:.3f}")
    # substep convergence
    print("  substep convergence (nu 0.5, s_B 0.1, k_B 1):")
    for dt in (4e-3, 2e-3, 1e-3, 5e-4):
        res = const_run(TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=1.0), T_DOSE, dt=dt)
        tb = add("b-dt", TICL4, res.byproduct, "constant 0.3 Pa", res, t_base_c, dt=dt, note="substep convergence")
        print(f"    dt {dt * 1e3:4.1f} ms: t90 {tb:.4f} s, closure {res.balance['closure']:.1e}")

    print("\n(c) TiCl4 with HCl at the chamber's centre-face density history (TMA run as a stand-in)")
    base_d = driven_run(TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.0), t, n0_series)
    t_base_d = add("c", TICL4, base_d.byproduct, "chamber centre face", base_d, note="no re-adsorption (single-species); TMA chamber history as stand-in")
    print(f"  no HCl action: bottom t90 {t_base_d:.3f} s; bottom coverage at the TiCl4 rule pulse ({RULE['TiCl4']} s): {profile_at(base_d, RULE['TiCl4'])[0][-1]:.3f}")
    driven = {}
    for nu in (0.5, 4.0):
        for sB in (1e-2, 0.1, 0.4):
            for kB in (1.0, 100.0):
                res = driven_run(TICL4, Byproduct("HCl", 36.46, nu=nu, sticking=sB, k_des=kB), t, n0_series)
                tb = add("c", TICL4, res.byproduct, "chamber centre face", res, t_base_d, note="TMA chamber history as stand-in")
                driven[(nu, sB, kB, 0.0)] = res
                print(f"  nu {nu:3.1f}, s_B {sB:5.0e}, k_B {kB:6.1f}: t90 {tb:7.3f} s ({100 * (tb / t_base_d - 1):+.1f} %), theta_P bottom at the rule pulse {profile_at(res, RULE['TiCl4'])[0][-1]:.3f}, at 14.4 s {profile_at(res, 14.38)[0][-1]:.4f}")
    for ke in (1e-4, 1e-3):
        res = driven_run(TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=1.0, k_etch=ke), t, n0_series)
        tb = add("c-etch", TICL4, res.byproduct, "chamber centre face", res, t_base_d, note="etch bound; TMA chamber history as stand-in")
        driven[(0.5, 0.1, 1.0, ke)] = res
        print(f"  etch k_etch {ke:5.0e} (nu 0.5, s_B 0.1, k_B 1): t90 {tb:7.3f} s, bottom theta_P at 20 s {res.theta_P[-1]:.3f}")

    # --- (d) byproduct at the mouth -----------------------------------------------------------------------------------------
    print("\n(d) byproduct mouth density (chamber accumulation), fraction of the precursor mouth density")
    sweep = {}
    for label, pr, bp0, runner, ref in (
            ("TiCl4 const, k_B 1", TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=1.0), lambda bp: const_run(TICL4, bp, T_DOSE), t_base_c),
            ("TiCl4 const, k_B 100", TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=100.0), lambda bp: const_run(TICL4, bp, T_DOSE), t_base_c),
            ("TiCl4 driven, k_B 1", TICL4, Byproduct("HCl", 36.46, nu=0.5, sticking=0.1, k_des=1.0), lambda bp: driven_run(TICL4, bp, t, n0_series), t_base_d),
            ("TMA const, bound", TMA, Byproduct("CH4", 16.04, nu=3.0, sticking=1e-3, k_des=1e3), lambda bp: const_run(TMA, bp, 20.0), None)):
        line = []
        for frac in (0.0, 0.01, 0.03, 0.10):
            bp = Byproduct(bp0.name, bp0.mass_amu, nu=bp0.nu, sticking=bp0.sticking, k_des=bp0.k_des, mouth_fraction=frac)
            res = runner(bp)
            drive = "chamber centre face" if "driven" in label else "constant 0.3 Pa"
            tb = add("d", pr, bp, drive, res, ref if ref is not None else None, note=label)
            line.append(tb); sweep[(label, frac)] = tb
        print(f"  {label}: t90 = " + ", ".join(f"{v:7.3f}" for v in line) + " s at 0, 1, 3, 10 %")

    # --- outputs -----------------------------------------------------------------------------------------------------------
    with open(ROOT / "figures/multispecies.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\nwrote figures/multispecies.csv ({len(rows)} runs) in {time.time() - t0:.0f} s")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    # (1) TMA histories
    ax = axes[0, 0]
    for (case, drive), (r0, rb) in profiles.items():
        ls = "-" if drive.startswith("constant") else "--"
        ax.plot(r0.t, r0.hist_P[:, -1], ls, color="k", label=f"TMA, CH4 s_B = 0 ({drive})")
        ax.plot(rb.t, rb.hist_P[:, -1], ls, color="tab:red", lw=1, label=f"TMA, CH4 s_B 1e-3, k_B 1e3 (bound)")
    ax.axhline(0.9, color="0.6", ls=":"); ax.set_xlim(0, 22); ax.set_xlabel("time (s)"); ax.set_ylabel("bottom coverage theta_P")
    ax.set_title("(a) TMA/H2O: methane does not re-adsorb; the bound changes nothing", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # (2) TiCl4 t90 vs k_B
    ax = axes[0, 1]
    colors = {1e-3: "tab:blue", 1e-2: "tab:green", 0.1: "tab:orange", 0.4: "tab:red"}
    for nu, mk in ((0.5, "o-"), (4.0, "s--")):
        for sB in grid_s:
            ax.plot(grid_k, [tb_grid[(nu, sB, k)] for k in grid_k], mk, color=colors[sB], label=f"nu {nu:g}, s_B {sB:g}")
    ax.axhline(t_base_c, color="k", lw=1, label="no HCl action")
    ax.set_xscale("log"); ax.set_xlabel("HCl desorption rate k_B (1/s)"); ax.set_ylabel("time to 90 % bottom coverage (s)")
    ax.set_title("(b) TiCl4/NH3 hole at 0.3 Pa: HCl re-adsorption range", fontsize=9); ax.legend(fontsize=6, ncol=2); ax.grid(alpha=0.3)
    # (3) profiles at the rule pulse, driven
    ax = axes[1, 0]
    pP, pB, tt = profile_at(base_d, RULE["TiCl4"])
    ax.plot(DEPTHS * AR, pP, "k-", label=f"no HCl action, theta_P at {tt:.2f} s")
    for key, c in (((0.5, 0.1, 1.0, 0.0), "tab:orange"), ((0.5, 0.4, 1.0, 0.0), "tab:red"), ((4.0, 0.1, 1.0, 0.0), "tab:purple"), ((0.5, 0.1, 1.0, 1e-3), "tab:brown")):
        res = driven[key]; pP, pB, tt = profile_at(res, RULE["TiCl4"])
        lab = f"nu {key[0]:g}, s_B {key[1]:g}, k_B {key[2]:g}" + (f", k_etch {key[3]:g}" if key[3] else "")
        ax.plot(DEPTHS * AR, pP, "-", color=c, label=lab + " (theta_P)")
        ax.plot(DEPTHS * AR, pB, ":", color=c, label=lab + " (blocked)")
    ax.set_xlabel("depth / diameter"); ax.set_ylabel("coverage"); ax.set_title(f"(c) TiCl4 driven by the chamber centre-face history: profiles at the TiCl4 rule pulse ({RULE['TiCl4']} s)", fontsize=8)
    ax.legend(fontsize=6); ax.grid(alpha=0.3)
    # (4) mouth fraction sweep
    ax = axes[1, 1]
    for label, c in (("TiCl4 const, k_B 1", "tab:orange"), ("TiCl4 const, k_B 100", "tab:green"), ("TiCl4 driven, k_B 1", "tab:blue"), ("TMA const, bound", "tab:red")):
        fr = [0.0, 0.01, 0.03, 0.10]
        ax.plot([100 * q for q in fr], [sweep[(label, q)] for q in fr], "o-", color=c, label=label)
    ax.set_xlabel("byproduct mouth density, % of the precursor mouth density"); ax.set_ylabel("time to 90 % bottom coverage (s)")
    ax.set_title("(d) byproduct accumulated in the chamber (nu 0.5, s_B 0.1 for TiCl4; TMA bound s_B 1e-3, k_B 1e3)", fontsize=8); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(ROOT / "figures/multispecies.png", dpi=150); print("saved figures/multispecies.png")


if __name__ == "__main__":
    main()
