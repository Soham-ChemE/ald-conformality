# Work package 12: two-species feature kinetics with byproduct re-adsorption

At the chamber pressures of this work (100 to 1,000 Pa of N2 at 573 K) the mean free path is 12 to 123 um, so
inside a 100 nm feature the Knudsen number is 120 to 1,200: molecules collide with the walls, not with each other. A
reaction byproduct therefore cannot block the precursor in the gas phase, and each species diffuses independently with
its own Knudsen diffusivity D_K = d vbar/3, proportional to 1/sqrt(M) (methane 2.1 x faster than TMA, HCl 2.3 x faster
than TiCl4). What a byproduct can do is surface-side: re-adsorb on a free site and block it until it desorbs, or react
with the film (etching, documented for HCl on some oxides and nitrides). For TMA/H2O the byproduct is methane, which
does not re-adsorb at these temperatures, so the single-species model of the earlier work packages stands and this
package shows it; for chloride chemistries (TiCl4/NH3 giving HCl, HfCl4/H2O) HCl re-adsorption is documented
(Tiznado and Zaera 2006; Snyder et al. 2006; Leem et al. 2014) and the two-species model matters there.

Repository state: no commits made; nothing outside the file list below was edited. README.md, CLAIMS.md,
paper/assumptions.md, paper/draft.md, main.tex and docs/ are untouched (their text is proposed at the end). lib/aldmc,
src/aldconf/voxel.py and tests/test_voxel_mc.py were not touched (another package works there).

## What was built

| File | Purpose |
|---|---|
| `src/aldconf/multispecies.py` (new) | `Byproduct` (name, mass, nu, sticking s_B, k_des, k_etch, mouth_fraction), `simulate_two_species` (constant mouth density, cycles and narrowing as `simulate`), `simulate_two_species_driven` (a mouth density history, as `simulate_driven`), `TwoSpeciesResult` (theta_P, theta_B, both gas profiles, histories at chosen depths, the byproduct bookkeeping, `time_to_bottom_coverage`), `advance_coverages` (closed-form exact advance of the linear coverage pair), `knudsen_number`; the named byproducts `METHANE`, `HCL_TICL4_HALFSTEP` (nu 0.5), `HCL_TICL4_CYCLE` (nu 4) |
| `src/aldconf/solver.py` | `_solve_gas` gained an optional `source` argument (volumetric right-hand side); `None` leaves the equation and its arithmetic unchanged (the 95 earlier tests pass unchanged) |
| `src/aldconf/__init__.py` | exports `Byproduct`, `TwoSpeciesResult`, `simulate_two_species`, `simulate_two_species_driven` |
| `tests/test_multispecies.py` (new) | 8 tests: s_B = 0 (any nu) and nu = 0 reproduce `simulate` to 1e-9 (measured 0 and 1e-16); k_B = 1e12 reproduces it to 1e-9; the driven variant equals `simulate_driven` to 1e-12; byproduct mass balance closes to 1e-6 (measured 4e-16 to 7e-16) with and without etching and with a mouth byproduct density; the surface inventory agrees with the gas-side bookkeeping to first order in the substep; the exact 2 x 2 advance against a fine explicit integration (1e-5) and its degenerate cases; HCl re-adsorption lengthens the TiCl4 pulse and etching lowers the bottom coverage further; a mouth byproduct density blocks sites with nu = 0; the Knudsen numbers |
| `coupling/multispecies.py` (new) | the cases (a)-(d) below plus the substep convergence; writes `figures/multispecies.csv` (96 runs) and `figures/multispecies.png`; 12 minutes on one core |
| `figures/multispecies.csv`, `figures/multispecies.png` (new) | one row per run: case, precursor, byproduct, nu, s_B, k_B, k_etch, mouth fraction, drive, dt, bottom and mouth t90, change in per cent, bottom coverage at the rule pulse, blocked coverage, escaped fraction, balance closure, note |
| `lib/aldkin/include/aldkin.h`, `src/aldkin.cpp` | version 1.1.0: `aldkin_params` gained `nu`, `vbarB`, `sB`, `kB`, `kEtch`, `cB0` (appended; all 0 by default = the single-species path, bit for bit); selectors `ALDKIN_THETA_BLOCKED_FLAT/_BOTTOM/_HOLE`; `aldkin_set_byproduct_mouth` (per-face byproduct density for a host that transports it), `aldkin_byproduct_flux` (net byproduct flux into the chamber per unit wafer area), `aldkin_byproduct_balance` (the hole's produced/desorbed/escaped/re-adsorbed/etched rates), `aldkin_hole_byproduct_gas`. The two-species pass (`updateTwo`) runs only when `sB > 0` or `kEtch > 0`; the single-species pass is the previous code with the previous operand order (the optional `thB` pointer of `solveHole` is null there) |
| `lib/aldkin/tests/test_aldkin.cpp`, `CMakeLists.txt` | 5 new ctest cases (16 in all): bit-identity of the inactive fields for both models and both schemes; the k_B -> infinity limit (1e14: coverages within 1e-9, blocked below 1e-9); the mass balance after every step to 1e-9 with and without etching and with a mouth density, the per-face flux as the open-area blend, the per-face mouth override; re-adsorption slows the bottom and etching more, Euler and exponential within 0.01; the error paths |
| `lib/aldkin/README.md` | the new fields, functions and tests; the dictionary keywords left as the next step |
| `src/aldconf/aldkin.py`, `tests/test_aldkin.py` | the binding exposes the six fields (`Params`), `theta("blocked_flat" / "blocked_bottom" / "blocked_hole")`, `set_byproduct_mouth`, `byproduct_flux`, `byproduct_balance`, `hole_byproduct_gas`; 3 new tests (12 in all): round trip and inactive bit-identity, the library against `aldconf.multispecies` on the same hole (coverages within 1e-9 over 2,000 steps, byproduct gas within 1e-9 relative, bookkeeping closed), the flat-face mouth density |

The OpenFOAM boundary conditions still call `aldkin_params_default` and set the single-species fields, so they run
the single-species path unchanged (regression below). Reading `nu`, `vbarB`, `sB`, `kB`, `kEtch`, `cB0` from the
dictionary and feeding a transported HCl field to `aldkin_set_byproduct_mouth` is the stated next step; it was not
a five-line change that keeps the written fields unchanged, so it was left out.

## Model

One half-cycle of the precursor P at the feature mouth. Coverages per cell: theta_P (sites taken by P, irreversible),
theta_B (sites blocked by the byproduct B, reversible); free fraction f = 1 - theta_P - theta_B. With the local rates
a_P = s_P (vbar_P/4) n_P/N0, a_B = s_B (vbar_B/4) n_B/N0 and e = k_etch (vbar_B/4) n_B/N0:

    d theta_P/dt = a_P f - e theta_P
    d theta_B/dt = a_B f - k_B theta_B

Gas, quasi-steady for both species (the byproduct relaxes on the same L^2/D_K ~ 10 us scale as the precursor), closed
end, Dirichlet mouth densities n_P0(t) and n_B0 = mouth_fraction x n_P0:

    (1/A) d/dx (A D_P dn_P/dx) = (P/A) (vbar_P/4) [s_P f + r] n_P
    (1/A) d/dx (A D_B dn_B/dx) = (P/A) (vbar_B/4) [s_B f + k_etch theta_P] n_B - (P/A) [nu s_P (vbar_P/4) f n_P + k_B N0 theta_B]

Two tridiagonal solves per substep (the second with the source on the right-hand side). The coverage pair is advanced
with the gas frozen by the exact solution of its linear 2 x 2 system, theta(dt) = theta0 + phi (M theta0 + b) with
phi = integral of exp(M s), written in closed form from the two real eigenvalues (the smaller one taken as det/l2 so
that k_B -> infinity causes no cancellation). With s_B = k_etch = 0 the precursor update is exactly the single-species
theta <- 1 - (1 - theta) exp(-a_P dt) and the precursor gas equation is unchanged. Bookkeeping, gas side, per substep
with the coverages the gas was solved with: produced (nu x P uptake) + desorbed = escaped through the mouth +
re-adsorbed + consumed by etching; the finite-volume solve conserves it to round-off (closure 1e-16 to 1e-14 in every
run). The surface inventory N0 theta_B agrees with the gas-side integral to first order in the substep; the byproduct
rates (a_B up to 50 1/s for HCl at s_B 0.05) set the substep, not the precursor's 0.3 1/s, so every two-species run
uses 1 ms substeps (convergence table below). Not modelled: the etch product as a third species; byproduct produced
on the flat wafer surface outside the feature (swept away by the flow unless mouth_fraction > 0); gas-phase
collisions between species (Knudsen number >> 1).

Stoichiometry. nu = 3 methane per TMA over an Al2O3 cycle (ligand exchange, part in each half-step; Puurunen 2005);
nu = 4 HCl per TiCl4 over a TiN cycle (6 TiCl4 + 8 NH3 -> 6 TiN + 24 HCl + N2, Ritala et al. 1995), of which about 0.5
leave during the TiCl4 half-step (Cl:Ti falls from 4 to about 3.5 after the TiCl4 exposure, Tiznado and Zaera 2006)
and the rest during the NH3 step. The TiCl4 half-cycle therefore runs at nu = 0.5 (anchored) with nu = 4 as the bound
in which the whole cycle's HCl is released in the hole during the TiCl4 pulse.

## Parameters and their standing

| Quantity | Value(s) run | Standing | Source |
|---|---|---|---|
| Feature, temperature, drive | 100 nm x 10 um hole (AR 100), 573 K; constant 0.3 Pa mouth (the Machball-check benchmark) and the chamber's centre-face density history of `reactor/showerhead_pulse_patterned_AR100_long` (1 % TMA, 300 Pa, 300 sccm, 20 s pulse) | as in every earlier package | wp4, wp8 |
| TMA | s 5e-3, N0 4e18, M 72.09; CH4 M 16.04, nu 3 | chamber-run values | Puurunen 2005 (nu) |
| TiCl4 | s 0.10, N0 1.24e18, M 189.68; HCl M 36.46 | s is for TiCl4/H2O on TiO2 at 110 C (no TiCl4/NH3 value published); N0 from 0.30 A per cycle and 4.23 g/cm3 | Ylilammi et al. 2018; Kinnunen et al. 2017 (as wp4) |
| nu (TiCl4 half-step) | 0.5 anchored; 4 bound | Cl:Ti 4 -> 3.5 after the TiCl4 step; whole-cycle 4 | Tiznado and Zaera 2006; Ritala et al. 1995 |
| s_B (HCl) | 1e-3, 1e-2, 0.1, 0.4 | range: no value on TiN or TiO2 at 573 K exists; 1e-3 is the reactive sticking of HCl on alpha-Al2O3(0001) at 298 K at low coverage, 0.38 +/- 0.05 the Knudsen-cell initial uptake coefficient on alumina at 296 K | Elam et al. 2000; Lostier et al. 2026 (Table 2) |
| k_B (HCl desorption) | 0.1, 1, 10, 100, 1,000 1/s | range: HCl desorbs from alpha-Al2O3(0001) over 300-650 K in TPD, which with a 1e13 1/s prefactor is E_des 0.8-1.8 eV (Redhead, 2 K/s) and k_B = 1e-3 to 1e6 1/s at 573 K; the run range covers the middle; k_B <= 0.1 is the strongest-bound-site end | Elam et al. 2000 |
| k_etch (HCl on Ti sites) | 0; bounds 1e-4, 1e-3 per collision | no rate exists: HCl does not etch TiO2 in TiCl4/H2O ALD (it reduces the adsorption sites instead); the TiN kinetic Monte Carlo lists Ti removal as TiCl3/TiCl4 by surface Cl without a rate usable here | Leem et al. 2014; Kim et al. 2022 |
| CH4 re-adsorption | s_B 1e-3, k_B 1e3 1/s | hypothetical bound; no literature re-adsorption of methane on Al2O3 at 573 K | - |
| mouth fraction | 0, 1, 3, 10 % | stands for byproduct accumulated in the chamber; a parameter | - |

The literature-anchored items are the stoichiometries, the Knudsen numbers and diffusivity ratios, and the TMA/TiCl4
parameters already in use; every HCl kinetic number is a range or a bound and is labelled so in the tables.

## Results

Knudsen numbers of the 100 nm hole (N2 carrier, hard-sphere 3.8 A, 573 K): mean free path 123 / 41 / 12 um and
Kn = 1,233 / 411 / 123 at 100 / 300 / 1,000 Pa. D_K(CH4)/D_K(TMA) = 2.12 (vbar 870 vs 410 m/s);
D_K(HCl)/D_K(TiCl4) = 2.28 (577 vs 253 m/s).

### (a) TMA/H2O with methane (literature-anchored: methane does not re-adsorb)

| drive | single-species bottom t90 | two-species, CH4 nu 3, s_B 0 | bound: s_B 1e-3, k_B 1e3 1/s | methane escaping the hole (bound) |
|---|---|---|---|---|
| constant 0.3 Pa mouth | 15.600 s | 15.600 s, max coverage difference 0 | 15.601 s (+0.01 %) | 99.95 % |
| chamber centre-face history | 14.005 s | 14.005 s, max coverage difference 0 | 14.006 s (+0.01 %) | 99.95 % |

The single-species result is reproduced exactly (the byproduct gas is solved and 1e-9 of the sites are transiently
blocked at the bound, which changes nothing). The 15.600 s here against 15.616 s in the library test is the 10 ms
history grid of this script against the 1 ms grid there.

### (b) TiCl4/NH3 with HCl at a constant 0.3 Pa mouth density (ranges)

Bottom t90 without HCl action (single-species): 7.754 s. With re-adsorption (nu, s_B, k_B; no etching), bottom t90 in s:

| nu | s_B | k_B 0.1 | k_B 1 | k_B 10 | k_B 100 | k_B 1,000 1/s |
|---|---|---|---|---|---|---|
| 0.5 (anchored) | 1e-3 | 23.51 (+203 %) | 8.59 (+11 %) | 7.75 | 7.75 | 7.75 |
| 0.5 | 1e-2 | 24.34 (+214 %) | 8.65 (+12 %) | 7.75 | 7.75 | 7.75 |
| 0.5 | 0.1 | 24.51 (+216 %) | 10.06 (+30 %) | 7.93 (+2 %) | 7.75 | 7.75 |
| 0.5 | 0.4 | 27.67 (+257 %) | 12.00 (+55 %) | 8.51 (+10 %) | 7.79 (+0.5 %) | 7.75 |
| 4 (bound) | 1e-3 | 24.73 | 8.41 (+8 %) | 7.75 | 7.75 | 7.75 |
| 4 | 1e-2 | 27.32 | 8.97 (+16 %) | 7.78 | 7.75 | 7.75 |
| 4 | 0.1 | not within 30 s | 13.36 (+72 %) | 8.42 (+9 %) | 7.77 | 7.75 |
| 4 | 0.4 | not within 30 s | 20.39 (+163 %) | 10.16 (+31 %) | 7.97 (+3 %) | 7.75 |

Reading: k_B is the controlling parameter. At k_B >= 100 1/s (E_des below about 1.25 eV at 573 K) HCl is a spectator
whatever its sticking; at k_B = 1 1/s the anchored half-step stoichiometry costs 11 to 55 % of the pulse depending on
s_B; at k_B = 0.1 1/s the blocked sites outlive the pulse and the hole fills three times slower, almost
independently of s_B (the sites are blocked faster than they are freed in any case). The sticking coefficient matters
only in the middle of the range.

With etching (bound; k_etch per HCl collision on a Ti site), bottom t90 and the state at 30 s:

| nu | s_B | k_B 1/s | k_etch 1e-4 | k_etch 1e-3 | blocking only (above) |
|---|---|---|---|---|---|
| 0.5 | 1e-2 | 1 | 10.44 s | 13.98 s | 8.65 s |
| 0.5 | 1e-2 | 100 | 9.32 s | 12.46 s | 7.75 s |
| 0.5 | 0.1 | 1 | 11.44 s | 14.17 s | 10.06 s |
| 0.5 | 0.1 | 100 | 9.74 s | 13.11 s | 7.75 s |
| 4 | 1e-2 | 1 | never: bottom theta_P 0.045, blocked 0.94 at 30 s | never: 0.001 / 0.98 | 8.97 s |
| 4 | 1e-2 | 100 | never: 0.000 / 0.37 | never: 0.000 / 0.37 | 7.75 s |
| 4 | 0.1 | 1 | never: 0.161 / 0.84 | never: 0.011 / 0.99 | 13.36 s |
| 4 | 0.1 | 100 | never: 0.000 / 0.85 | never: 0.000 / 0.85 | 7.77 s |

Etching acts even when blocking does not (k_B 100): every HCl collision on a coated site has a chance to undo it, and
the HCl flux in the hole is nu times the precursor uptake. At nu = 0.5 an etch probability of 1e-4 costs 20 to 35 %
and 1e-3 costs 70 to 80 % of the pulse; at nu = 4 the bottom never reaches 90 %: the hole settles into an
etch-limited steady state in which the bottom cells are mostly blocked or bare. That bound is the regime a chlorine
chemistry must avoid, and it is where a measured HCl etch rate on TiN at 573 K would be decisive.

Substep convergence (nu 0.5, s_B 0.1, k_B 1): bottom t90 10.111 / 10.074 / 10.059 / 10.052 s at 4 / 2 / 1 / 0.5 ms;
the 1 ms runs above are converged to 0.1 %. Byproduct balance closure 1e-16 to 1e-14 in every run (column
`balance_closure`).

### (c) The saturating pulse and the profile at the rule's pulse, chamber-driven and constant

Chamber-driven TiCl4 uses the TMA run's centre-face density history as a stand-in (no TiCl4 chamber run exists; the
TiCl4 wafer, with 3.2 x fewer sites, would saturate sooner and its surface density would rise sooner), so only the
changes are meaningful there, not the absolute times. Bottom t90 without HCl action: 9.622 s driven, 7.754 s constant.

| nu | s_B | k_B 1/s | driven bottom t90 | change | bottom theta_P at the TiCl4 rule pulse 4.44 s | at the TMA rule pulse 14.4 s | constant-mouth t90 | change |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 1e-2 | 1 | 10.72 s | +11 % | 0.000 | 0.997 | 8.65 s | +12 % |
| 0.5 | 1e-2 | 100 | 9.62 s | 0 | 0.000 | 1.000 | 7.75 s | 0 |
| 0.5 | 0.1 | 1 | 11.76 s | +22 % | 0.000 | 0.991 | 10.06 s | +30 % |
| 0.5 | 0.1 | 100 | 9.62 s | 0 | 0.000 | 1.000 | 7.75 s | 0 |
| 0.5 | 0.4 | 1 | 13.11 s | +36 % | 0.000 | 0.959 | 12.00 s | +55 % |
| 0.5 | 0.4 | 100 | 9.67 s | +0.5 % | 0.000 | 1.000 | 7.79 s | +0.5 % |
| 4 | 1e-2 | 1 | 10.89 s | +13 % | 0.012 | 0.997 | 8.97 s | +16 % |
| 4 | 1e-2 | 100 | 9.62 s | 0 | 0.000 | 1.000 | 7.78 s | +0.4 % |
| 4 | 0.1 | 1 | 14.06 s | +46 % | 0.058 | 0.922 | 13.36 s | +72 % |
| 4 | 0.1 | 100 | 9.65 s | +0.3 % | 0.000 | 1.000 | 7.77 s | +0.3 % |
| 4 | 0.4 | 1 | 16.72 s | +74 % | 0.118 | 0.698 | 20.39 s | +163 % |
| 4 | 0.4 | 100 | 9.84 s | +2.3 % | 0.000 | 1.000 | 7.97 s | +3 % |
| 0.5, k_etch 1e-4 | 0.1 | 1 | 12.51 s | +30 % | - | - | 11.44 s | +48 % |
| 0.5, k_etch 1e-3 | 0.1 | 1 | 13.74 s | +43 % | - | - | 14.17 s | +83 % |

The relative cost of re-adsorption is 20 to 40 % smaller under the chamber history than at a constant mouth density
(the rising surface density at the end of the chamber pulse refills the released sites faster), the same direction as
the supply-limited taper result of wp4. The profiles at the rule pulse (figure, panel c) show the mechanism: with
blocking the precursor front is broader and reaches deeper at a given time (the blocked near-mouth sites let TiCl4
pass), but the sites behind the front must be released and refilled, which is where the extra time goes; the blocked
coverage peaks just ahead of the precursor front. The TiCl4 rule pulse of 4.44 s (wp4) is the TiCl4 wafer's own
supply time at 1 %; under the TMA history the hole is not yet saturated at 4.44 s even without HCl, which is the
stand-in's limitation, not a result.

### (d) Byproduct at the mouth (chamber accumulation; a parameter)

Bottom t90 for a byproduct mouth density of 0 / 1 / 3 / 10 % of the precursor mouth density:

| case | 0 % | 1 % | 3 % | 10 % |
|---|---|---|---|---|
| TiCl4 constant 0.3 Pa, nu 0.5, s_B 0.1, k_B 1 | 10.059 s | 10.111 (+0.5 %) | 10.215 (+1.6 %) | 10.568 (+5.1 %) |
| TiCl4 constant, nu 0.5, s_B 0.1, k_B 100 | 7.749 s | 7.749 | 7.749 | 7.749 |
| TiCl4 chamber-driven, nu 0.5, s_B 0.1, k_B 1 | 11.759 s | 11.808 (+0.4 %) | 11.905 (+1.2 %) | 12.234 (+4.0 %) |
| TMA constant, CH4 bound s_B 1e-3, k_B 1e3 | 15.601 s | 15.601 | 15.601 | 15.601 |

HCl accumulated in the chamber at 10 % of the precursor density costs 4 to 5 % of the pulse when k_B is 1 1/s and
nothing when it is 100 1/s: inside the feature the locally produced HCl (up to 0.5 to 4 molecules per adsorbed TiCl4,
at a density that rivals the precursor's near the front) outweighs what the chamber brings in. Purge design for HCl
is a chamber-side matter (corrosion, particles) rather than a feature-side one, within this model.

### Summary of the pulse-time change by case

| case | anchoring | pulse-time change |
|---|---|---|
| (a) TMA/CH4, s_B 0 | literature (methane does not re-adsorb) | 0 (identical to the single-species solver) |
| (a) TMA/CH4 bound s_B 1e-3, k_B 1e3 | hypothetical bound | +0.01 % |
| (b) TiCl4/HCl, nu 0.5, k_B >= 100 1/s | nu anchored; k_B range | 0 to +0.5 % |
| (b) TiCl4/HCl, nu 0.5, k_B 1 1/s | nu anchored; s_B and k_B ranges | +11 % (s_B 1e-3) to +55 % (s_B 0.4) |
| (b) TiCl4/HCl, nu 0.5, k_B 0.1 1/s | strongest-bound-site end of the range | +203 to +257 % |
| (b) TiCl4/HCl, nu 4 (whole-cycle HCl in the half-step) | bound | up to +163 % at k_B 1; no saturation within 30 s at k_B 0.1 and s_B >= 0.1 |
| (b) etching k_etch 1e-4 / 1e-3, nu 0.5 | bound | +20 to +35 % / +70 to +80 %; at nu 4 the bottom never saturates |
| (c) same under the chamber history | stand-in history | +11 to +74 % (blocking), +30 / +43 % (etch bounds); 20 to 40 % smaller relative cost than at constant density |
| (d) HCl at the mouth 1 / 3 / 10 % | parameter | +0.5 / +1.6 / +5 % at k_B 1; 0 at k_B 100 |

## Library

ctest (`lib/aldkin/build`): 16 tests, all passed (11 before this package + 5 two-species cases); the inactive-field
bit-identity test compares coverages, K_eff and flux with `==` over 400 steps for the patterned and flat-only models
under both schemes. pytest `tests/test_aldkin.py`: 12 passed (9 + 3), among them the library against
`aldconf.multispecies` on the same TiCl4/HCl hole: precursor and blocked coverages within 1e-9 over 2,000 steps of
1 ms, byproduct gas within 1e-9 relative, bookkeeping closed to 1e-9. During the build a modelling inconsistency was
caught by that twin test: the library's precursor sink first used 1 - theta instead of the free fraction
1 - theta - theta_B, which stalled the front behind the Python result; the single-species path was not affected (the
pointer is null there) and the fix is tested.

Regression: `./reproduce.sh --regress` through the launcher at the end of this package (cmake + ctest for lib/aldkin,
16 passed, and lib/aldmc, 10 passed; wmake of the adapters against the 1.1 library; steady flow, 2000 iterations; the
2 s regression pulse; `tests/test_regression_reference.py`): "3 passed in 2.79s". lib/aldmc was built in a scratch
directory (`ALDMC_BUILD_DIR`) so as not to disturb the other package's build directory; reproduce.sh still rewrites
`lib/aldmc/log.cmake`.

## Figures

- `figures/multispecies.png`: (a) TMA bottom-coverage histories, single-species against the two-species model with
  methane at s_B = 0 and at the hypothetical bound, constant mouth (solid) and chamber-driven (dashed): the curves lie
  on top of each other. (b) TiCl4/NH3 hole at 0.3 Pa: time to 90 % bottom coverage against the HCl desorption rate
  k_B for s_B 1e-3 to 0.4 and nu 0.5 (circles) and 4 (squares); the line is the single-species 7.75 s. (c) Coverage
  along the hole at the TiCl4 rule pulse (4.44 s) under the chamber's centre-face history: precursor coverage (solid)
  and blocked coverage (dotted) for no HCl action, nu 0.5 at s_B 0.1 and 0.4, nu 4 at s_B 0.1 and the etch bound.
  (d) Time to 90 % against the byproduct mouth density (0 to 10 % of the precursor's) for the TiCl4 cases at k_B 1 and
  100 1/s, driven and constant, and the TMA bound.
- `figures/multispecies.csv`: the 96 runs.

## Proposed paper text

### 3.13 Byproducts: what a second species does and does not do

At the pressures of this work the mean free path is 12 to 123 um, so inside a 100 nm feature the Knudsen number is
120 to 1,200 and molecules collide with the walls, not with each other: a reaction byproduct cannot block the
precursor in the gas phase, and each species diffuses independently with its own Knudsen diffusivity (methane 2.1 x
faster than TMA, HCl 2.3 x faster than TiCl4). What a byproduct can do is surface-side: re-adsorb on a free site and
block it until it desorbs, or etch the film. We extended the feature solver to two species (aldconf.multispecies and
the library option): the byproduct is released at nu molecules per adsorbed precursor, solved quasi-steadily with a
source on the right-hand side, re-adsorbs on free sites with sticking s_B, desorbs at k_B and removes a coated site
with probability k_etch per collision; the precursor sees the free fraction 1 - theta_P - theta_B. With s_B = k_etch = 0,
or k_B -> infinity, the single-species solver is reproduced to round-off, and the byproduct balance (produced +
desorbed = escaped + re-adsorbed + etched) closes to 1e-14 in every run. For TMA/H2O the byproduct is methane
(3 per TMA over a cycle [Puurunen 2005]), which does not re-adsorb at these temperatures; the two-species result is
the single-species one, and even a hypothetical re-adsorption at s_B = 1e-3 with 1 ms residence changes the AR-100
pulse by 0.01 % while 99.95 % of the methane leaves the hole (Figure, multispecies.png, panel a). Every result of
this paper therefore stands as computed.

For TiCl4/NH3 the byproduct is HCl, 4 per TiCl4 over a cycle [Ritala et al. 1995], about 0.5 of them in the TiCl4
half-step [Tiznado and Zaera 2006], and its re-adsorption is documented on TiN, silica and TiO2 [Tiznado and Zaera
2006; Snyder et al. 2006; Leem et al. 2014]; no sticking coefficient or desorption rate exists for HCl on TiN at
573 K, so we ran the range that the HCl/alumina surface science spans (s_B 1e-3 to 0.4 [Elam et al. 2000; Lostier et
al. 2026]; k_B 0.1 to 1,000 1/s, the 300-650 K TPD span of [Elam et al. 2000] read at 573 K). The desorption rate
decides (panel b, table): at k_B >= 100 1/s HCl is a spectator whatever its sticking; at k_B = 1 1/s the anchored
half-step stoichiometry lengthens the 100 nm x 10 um hole's saturating pulse by 11 to 55 % (7.75 s to 8.6-12.0 s at a
constant 0.3 Pa mouth) and by 11 to 36 % under the chamber's density history; at k_B = 0.1 1/s the blocked sites
outlive the pulse and the hole fills three times slower. Releasing the whole cycle's four HCl in the half-step (a
bound) doubles these costs, and an etch probability of 1e-3 per HCl collision (a bound; HCl does not etch TiO2 [Leem
et al. 2014], and the TiN kinetic Monte Carlo of [Kim et al. 2022] lists Ti removal by surface Cl without a rate)
costs 70 to 80 % at nu 0.5 and prevents saturation at nu 4: the hole settles into an etch-limited steady state. HCl
accumulated in the chamber at 10 % of the precursor density adds 4 to 5 % at k_B = 1 1/s and nothing at 100 1/s,
because the HCl produced inside the feature outweighs what the chamber brings in (panel d). The profiles at the rule
pulse (panel c) show the mechanism: blocked sites near the mouth let the precursor pass, so the front is broader and
deeper at a given time, but the sites behind it must be released and refilled, and that is where the time goes.

What this means for a chloride recipe is a measurement, not a number: the supply rule (Section 3.4) is untouched
when k_B >= 100 1/s, carries a 10 to 50 % chemistry margin at k_B ~ 1 1/s, and fails at k_B <= 0.1 1/s or with
etching, so a TiN or HfO2 tool model needs the HCl desorption rate on the growing film at the deposition
temperature (one TPD of HCl from ALD TiN would fix it) before the wall-area rule can be applied with the same
confidence as for Al2O3. The library carries the option behind defaults that reproduce the single-species boundary
conditions bit for bit; a transported HCl field in the chamber is the next step.

Discussion bullet:
- Byproducts: at Knudsen numbers of 100 to 1,000 a byproduct cannot block the precursor in the gas, only on the
  surface. For TMA/H2O (methane) the single-species model is exact to the purposes of this paper; for chlorides the
  HCl desorption rate at the deposition temperature decides between no effect (k_B >= 100 1/s), a 10 to 50 % longer
  pulse (k_B ~ 1 1/s) and failure of the supply rule (k_B <= 0.1 1/s or etching), and that rate is not in the
  literature for TiN or HfO2 at 573 K; the two-species solver and the library option (aldkin 1.1) are ready for it.

## README text (feature-solver row of "What is in the box")

Python: Knudsen diffusion + Langmuir uptake in slits and holes, constant or time-driven mouth density, fitting, scaling
tests; a two-species variant (`aldconf.multispecies`) adds a reaction byproduct that diffuses on its own, re-adsorbs
on free sites and may etch, and reduces to the single-species solver to round-off when it does neither. At the
Knudsen numbers of a 100 nm feature (100 to 1,000) methane from TMA changes nothing, while HCl from TiCl4 lengthens
the pulse by 0 to 55 % across the literature range of its desorption rate (`coupling/multispecies.py`).

## CLAIMS rows (CLAIMS.md format)

```
| Knudsen number of a 100 nm feature at the chamber pressures | mean free path 123 / 41 / 12 um, Kn 1,233 / 411 / 123 at 100 / 300 / 1,000 Pa (N2, 3.8 A, 573 K); D_K(CH4)/D_K(TMA) 2.12, D_K(HCl)/D_K(TiCl4) 2.28 | figures/multispecies.csv (script header output), src/aldconf/multispecies.py | PYTHONPATH=src .venv/bin/python coupling/multispecies.py |
| Two-species solver limits | s_B = 0 (any nu) and nu = 0 reproduce aldconf.simulate to 1e-9 (measured 0); k_B = 1e12 to 1e-9; byproduct balance produced + desorbed = escaped + re-adsorbed + etched closes to 1e-6 (measured 1e-16 to 1e-14) | tests/test_multispecies.py | pytest |
| TMA/H2O with methane (nu 3) on the AR-100 hole | bottom t90 15.600 s constant 0.3 Pa and 14.005 s chamber-driven, identical to the single-species solver; hypothetical re-adsorption s_B 1e-3, k_B 1e3 1/s: +0.01 %, 99.95 % of the methane escapes | figures/multispecies.csv (case a), figures/multispecies.png | PYTHONPATH=src .venv/bin/python coupling/multispecies.py; literature: Puurunen 2005, doi 10.1063/1.1940727 (stoichiometry) |
| TiCl4/NH3 with HCl (nu 0.5 in the TiCl4 half-step) on the AR-100 hole, constant 0.3 Pa | bottom t90 7.75 s without HCl action; 7.75-7.79 s at k_B >= 100 1/s; 8.59 / 8.65 / 10.06 / 12.00 s at k_B 1 1/s for s_B 1e-3 / 1e-2 / 0.1 / 0.4; 23.5-27.7 s at k_B 0.1 1/s | figures/multispecies.csv (case b), figures/multispecies.png | same; literature: Tiznado and Zaera 2006, doi 10.1021/jp062019f (nu 0.5); Ritala et al. 1995, doi 10.1149/1.2050083 (4 HCl per cycle); Elam et al. 2000, doi 10.1016/S0039-6028(99)01247-9 and Lostier et al. 2026, doi 10.48550/arXiv.2605.14826 (s_B and k_B ranges, HCl on alumina) |
| Same with the whole cycle's HCl released in the half-step (nu 4, bound) | bottom t90 8.41-20.39 s at k_B 1 1/s; no saturation within 30 s at k_B 0.1 1/s and s_B >= 0.1 | figures/multispecies.csv (case b) | same |
| HCl etching bound (k_etch 1e-4 / 1e-3 per collision on a Ti site) | nu 0.5: t90 9.3-11.4 s / 12.5-14.2 s against 7.75-10.06 s blocking only; nu 4: the bottom never reaches 90 % (etch-limited steady state, bottom theta_P 0.00-0.16 at 30 s) | figures/multispecies.csv (case b-etch) | same; literature: Leem et al. 2014, doi 10.5012/bkcs.2014.35.4.1195 (no etching of TiO2); Kim et al. 2022, doi 10.1016/j.commatsci.2022.111620 (Ti removal by surface Cl, no rate) |
| TiCl4/HCl under the chamber's centre-face history (TMA run as stand-in) | bottom t90 9.62 s without HCl action; +11 / +22 / +36 % (nu 0.5, k_B 1, s_B 1e-2 / 0.1 / 0.4), 0 to +0.5 % at k_B 100; +46 / +74 % at nu 4, s_B 0.1 / 0.4, k_B 1; etch bounds +30 / +43 %; relative cost 20-40 % below the constant-mouth case | figures/multispecies.csv (case c), figures/multispecies.png | same |
| HCl accumulated in the chamber (mouth density 1 / 3 / 10 % of the precursor's) | +0.5 / +1.6 / +5.1 % constant, +0.4 / +1.2 / +4.0 % driven at k_B 1 1/s; 0 at k_B 100 1/s; TMA bound unchanged | figures/multispecies.csv (case d) | same |
| Two-species substep convergence | bottom t90 10.111 / 10.074 / 10.059 / 10.052 s at 4 / 2 / 1 / 0.5 ms (nu 0.5, s_B 0.1, k_B 1) | figures/multispecies.csv (case b-dt) | same |
| Library two-species option (aldkin 1.1) | inactive fields bit-identical to the single-species handle (coverages, K_eff, flux ==, both models, both schemes); k_B 1e14 within 1e-9; balance to 1e-9 after every step; against aldconf.multispecies within 1e-9 over 2,000 steps | lib/aldkin/tests/test_aldkin.cpp, tests/test_aldkin.py | ctest --test-dir lib/aldkin/build; pytest tests/test_aldkin.py |
```

## Run log for paper/assumptions.md

```
## Night 4, part 3: two-species kinetics (src/aldconf/multispecies.py, coupling/multispecies.py, lib/aldkin 1.1; figures/multispecies.png, .csv)
Physics first: at 100-1,000 Pa the mean free path is 123-12 um (N2, 3.8 A, 573 K), Kn = 1,233-123 in a 100 nm hole,
so a byproduct cannot block the precursor in the gas phase; each species diffuses with its own Knudsen diffusivity
(CH4 2.12 x TMA, HCl 2.28 x TiCl4) and a byproduct acts only at the surface (re-adsorption, etching). The feature
solver gained a two-species variant: precursor P (irreversible, sticking on the free fraction 1 - theta_P - theta_B),
byproduct B released at nu per adsorbed P, quasi-steady gas for both (two tridiagonal solves per substep, the second
with the wall source nu x uptake + k_B N0 theta_B on the right-hand side; _solve_gas gained an optional source
argument, None leaves the earlier arithmetic unchanged), re-adsorption s_B, desorption k_B, etch probability k_etch
per collision on a coated site, byproduct mouth density as a fraction of the precursor's; the coverage pair advanced
per substep by the exact closed-form solution of its linear 2 x 2 system at frozen gas (eigenvalue of smaller
magnitude from det/l2, no cancellation at k_B -> infinity). Limits tested: s_B = 0 or nu = 0 reproduces simulate
(difference 0), k_B = 1e12 to 1e-9, the driven variant equals simulate_driven to 1e-12; the byproduct balance
produced + desorbed = escaped + re-adsorbed + etched closes to 1e-16 to 1e-14 (gas side, frozen coverages); the
surface inventory agrees with it to first order in the substep, and the HCl rates (a_B up to 50 1/s) set the substep
at 1 ms (t90 10.111 / 10.074 / 10.059 / 10.052 s at 4 / 2 / 1 / 0.5 ms). Stoichiometries: 3 CH4 per TMA per cycle
(Puurunen 2005); 4 HCl per TiCl4 per cycle (Ritala 1995), 0.5 in the TiCl4 half-step (Cl:Ti 4 -> 3.5, Tiznado and
Zaera 2006); nu 4 run as the bound. HCl kinetics: no value on TiN or TiO2 at 573 K; s_B 1e-3 (reactive sticking on
alpha-Al2O3(0001) at 298 K, Elam 2000) to 0.4 (Knudsen-cell initial uptake 0.38 on alumina at 296 K, Lostier 2026);
k_B 0.1-1,000 1/s from the 300-650 K TPD span of Elam 2000 (E_des 0.8-1.8 eV, k 1e-3-1e6 1/s at 573 K); k_etch 0
with bounds 1e-4, 1e-3 (HCl does not etch TiO2, Leem 2014; Ti removal by surface Cl in the TiN kMC of Kim 2022, no
rate). Results, AR-100 hole: (a) TMA/CH4 at s_B 0 identical to the single-species solver (15.600 s constant 0.3 Pa,
14.005 s chamber-driven), the hypothetical bound s_B 1e-3, k_B 1e3 +0.01 % with 99.95 % of the methane escaping;
(b) TiCl4/HCl constant 0.3 Pa: 7.754 s without HCl action; k_B >= 100: 0 to +0.5 %; k_B 1: 8.59 / 8.65 / 10.06 /
12.00 s for s_B 1e-3 / 1e-2 / 0.1 / 0.4 (nu 0.5), 8.41-20.39 s (nu 4); k_B 0.1: 23.5-27.7 s (nu 0.5), no saturation
in 30 s at nu 4 and s_B >= 0.1; etch bounds 1e-4 / 1e-3: +20-35 % / +70-80 % at nu 0.5, no saturation at nu 4
(etch-limited steady state); (c) under the chamber's centre-face history (the TMA run as a stand-in; only changes
are meaningful): 9.622 s without HCl action, +11 / +22 / +36 % (nu 0.5, k_B 1, s_B 1e-2 / 0.1 / 0.4), +46 / +74 %
(nu 4, s_B 0.1 / 0.4), 0 to +2.3 % at k_B 100, etch bounds +30 / +43 %; the relative cost is 20-40 % below the
constant-mouth one (the rising late-pulse density refills released sites); profiles at the rule pulse show a
broader, deeper front with blocking and the blocked coverage peaking ahead of it; (d) HCl at the mouth at 1 / 3 /
10 % of the precursor density: +0.5 / +1.6 / +5.1 % constant and +0.4 / +1.2 / +4.0 % driven at k_B 1, nothing at
k_B 100, nothing for the TMA bound. Library: aldkin 1.1.0 with nu, vbarB, sB, kB, kEtch, cB0 appended to
aldkin_params (defaults 0 = the single-species path), blocked-coverage selectors, aldkin_set_byproduct_mouth,
aldkin_byproduct_flux, aldkin_byproduct_balance, aldkin_hole_byproduct_gas; exponential scheme uses the same
closed-form 2 x 2 advance, Euler one forward step of both; 5 new ctest cases (16 passed): inactive fields ==
bit-identical for both models and schemes, k_B 1e14 within 1e-9, balance to 1e-9 per step, blocking and etching
ordering with Euler/exponential within 0.01, error paths; the binding exposes the fields and 3 new tests compare the
library with aldconf.multispecies within 1e-9 over 2,000 steps (a free-fraction inconsistency in the library's
precursor sink was caught by this test and fixed; the single-species path was not affected). pytest: 106 passed
(95 + 8 + 3; tests/test_voxel_mc.py excluded while the other package builds). ./reproduce.sh --regress through the
launcher with the 1.1 library: ctest 16 passed (aldkin) and 10 passed (aldmc), wmake, steady flow, regression pulse,
"3 passed" (the boundary conditions run the single-species path unchanged). The boundary conditions do not read
the new fields yet (next step: dictionary keywords and a transported HCl field). Caveats: every HCl kinetic number
is a range or a bound, labelled in the tables; the TiCl4 sticking is the TiCl4/H2O value; the chamber-driven TiCl4
runs use the TMA run's density history; the etch product is not tracked; byproduct produced on the flat wafer
surface is not fed back to the feature unless mouth_fraction > 0.
```

## Sources

- Puurunen, Journal of Applied Physics 97, 121301 (2005), doi 10.1063/1.1940727. TMA/H2O surface chemistry: ligand
  exchange releases methane in both half-steps, three CH4 per Al over a cycle.
- Ritala, Leskelä, Rauhala, Haussalo, Journal of The Electrochemical Society 142, 2731 (1995), doi 10.1149/1.2050083.
  TiN from TiCl4 and NH3: 6 TiCl4 + 8 NH3 -> 6 TiN + 24 HCl + N2 (4 HCl per TiCl4); chlorine content falls with
  deposition temperature.
- Tiznado, Zaera, Journal of Physical Chemistry B 110, 13491 (2006), doi 10.1021/jp062019f. XPS of TiCl4/NH3 TiN ALD:
  Cl:Ti about 3.5 after the TiCl4 exposure (about 0.5 HCl released per Ti in the TiCl4 half-step), most of the
  remaining Cl removed by ammonia; residual Cl attributed to re-adsorption of the HCl byproduct.
- Snyder, McCool, DiCarlo, Tripp, DeSisto, Thin Solid Films 514, 97 (2006), doi 10.1016/j.tsf.2006.03.013. Infrared
  study of TiCl4/NH3 TiN ALD on silica: re-adsorption of the HCl byproduct during the initial cycles.
- Leem, Park, Li, Zhou, Jin, Shin, Min, Bulletin of the Korean Chemical Society 35, 1195 (2014),
  doi 10.5012/bkcs.2014.35.4.1195. Role of HCl in TiCl4/H2O TiO2 ALD: HCl does not etch the grown TiO2 but reduces
  the adsorption sites for the next step.
- Kim, An, Oh, Jung, Kim, Nam, Han, Computational Materials Science 213, 111620 (2022),
  doi 10.1016/j.commatsci.2022.111620. Atomistic kinetic Monte Carlo of TiCl4/NH3 TiN ALD with DFT reaction energies:
  surface Cl generates HCl; Ti removal as TiCl3/TiCl4 by surface Cl (abstract and indexed text; the full text was not
  accessible from here, so no barrier from it is used).
- Elam, Nelson, Tolbert, George, Surface Science 450, 64 (2000), doi 10.1016/S0039-6028(99)01247-9. HCl on
  alpha-Al2O3(0001): reactive sticking coefficient about 1e-3 at 298 K at the lowest coverages; desorption over
  300-650 K from a broad distribution of sites.
- Lostier, Segev, Kislev, Schwartz Roitman, Locoge, Romanias, arXiv:2605.14826 (2026), doi 10.48550/arXiv.2605.14826.
  Knudsen-cell uptake of HCl on mineral surrogates at 296 K, Table 2: initial uptake coefficient (geometric area)
  0.38 +/- 0.05 on Al2O3, 0.43 on CaCO3, 0.052 on quartz, 0.014 on fumed silica; saturation 1.9e14 molecules/cm2 on
  Al2O3, 33 % reversibly adsorbed.
- Ylilammi, Ylivaara, Puurunen, Journal of Applied Physics 123, 205301 (2018), doi 10.1063/1.5028178, and Kinnunen
  et al., Nuclear Instruments and Methods B 406, 152 (2017), doi 10.1016/j.nimb.2016.12.032: the TiCl4 sticking (0.10,
  TiCl4/H2O on TiO2 at 110 C) and the TiN site density (0.30 A per cycle, 4.23 g/cm3), as in wp4.
