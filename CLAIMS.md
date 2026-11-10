# Claims audit

Every quantitative statement in the README and the paper traces to a file in this repository. Nothing here is an
estimate from memory. Rows marked "run" are reproduced by the named script; "literature" rows cite the source.

| Claim | Value | Source in repo | How to reproduce |
|---|---|---|---|
| Blanket wafer saturating TMA pulse, baseline chamber | 0.35 s centre, 0.38 s edge | figures/reactor_langmuir_coverage.png | reactor/Allrun, langmuir cases |
| Patterned AR-100 wafer (56x wall area) saturating pulse at 1 % TMA | 14.05 s centre, 14.1 s worst | figures/pulse_time_map.csv (run 7) | reactor/sweep_pattern_map.py |
| Design rule constant | 1.05 to 1.27 x supply limit, 12 patterned cases | figures/pulse_time_map.csv | reactor/sweep_pattern_map.py |
| Blanket wafers vs supply limit | 1.5 to 4.3 x | figures/pulse_time_map.csv | same |
| Hole-model vs Machball | 1.20/1.29 s (AR 20), 4.45/4.46 (AR 50), 16.2/15.6 (AR 100) | tests/test_machball_crosscheck.py | pytest |
| Embedded BC vs standalone solver | within 1 % at all times | paper/assumptions.md part 5 | reactor/showerhead_pulse_patterned_check |
| PillarHall fit | s_TMA = 8.5e-3 (x/÷1.5), exposure 7.5 Pa s, rms 0.038 | figures/pillarhall_fit.png | coupling/fit_pillarhall.py |
| Argonne fluidized bed saturation vs supply limit | tau_sat = 1.12 (0.92-1.31) | figures/validation_supply_rule.png | coupling/validate_supply_rule.py; literature: Lu et al. JVST A 40 (2022), OSTI 1902316 |
| Argonne rotary drum utilization vs loading | 24.7 % (2 g) to 70.7 % (75 g) | same | literature: Coile et al. JVST A 38, 052403 (2020), OSTI 1660751, Table 2 |
| H2O half-cycle | blanket 1.74 s; patterned AR-100 13.4 s | figures/coupling_patterned_h2o_vs_tma.png | reactor/make_h2o_cases.py |
| PEALD recombination ceiling | AR 100 unreachable at r >= 3e-3 (5 s) | figures/peald_ceiling.png | coupling/peald_ceiling.py |
| Temporal vs spatial | 29.9 vs 4 s at 56x area (k = 30) | figures/temporal_vs_spatial.png | coupling/temporal_vs_spatial.py (assumed zone flux) |

Known limitations and open items are listed in paper/assumptions.md and the paper's Discussion. Sticking (5e-3) and
site density (4e18 /m2) in the chamber runs are literature-range inputs; the PillarHall fit gives 8.5e-3 and Arts et
al. 2019 give 0.5-2e-3. The patterned-wafer conclusions are supply-limited and insensitive to this; the blanket-wafer
numbers are not.
