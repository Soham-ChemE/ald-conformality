# Invention note: pattern-aware, exhaust-confirmed pulse control for ALD (draft, 2026-11-10)

Status: a concept with simulation evidence from this work. Not checked by a patent attorney. Not filed.

## The idea in one paragraph
An ALD tool controller that (1) predicts, from the wafer's patterned-area ratio A = A_wall/A_wafer and the precursor
delivery rate, the pulse time at which the last feature on the wafer saturates, using the supply-limit rule
t = 1.2 (1 + A) N0/(U c0) or the full chamber-feature twin; (2) measures the precursor breakthrough at the exhaust
during the pulse (QMS, RGA, FTIR or optical absorption downstream of the wafer); (3) uses the first cycle's
breakthrough time to infer the effective A of the wafer actually in the chamber and sets the remaining hundreds of
cycles, both half-cycles, from it; (4) terminates each pulse at the confirmed breakthrough plus a model-derived
margin that covers the lag between exhaust breakthrough and feature-bottom saturation; (5) flags drift (ampoule
depletion, delivery-line faults, wrong wafer type) when the breakthrough time departs from the prediction.

## Evidence from the twin (coupling/exhaust_breakthrough.py, figures/exhaust_breakthrough.png)
- The exhaust precursor fraction stays at ~5 % of the feed while the patterned wafer is absorbing, then rises
  sharply when the wafer saturates. The 50 % breakthrough time tracks the time the hole bottoms reach 90 % coverage
  across all 12 patterned cases (AR 20-100, 0.5-10 % TMA): lag 0.2 to 0.7 s, breakthrough 2 to 20 % late, never
  early. So the exhaust is a conservative end-point for the deepest features.
- The early-pulse exhaust level does NOT encode the loading (5.5 % of feed for every patterned wafer): while the
  wafer is supply-limited the exhaust only sees the bypass fraction. Loading must come from design data or from the
  breakthrough time itself.
- Cycle-1 calibration works: inferring A from the 1 % TMA breakthrough gives A_est within 6 % of the true value for
  AR 20/50/100 and predicts the saturating pulse at 0.5, 3 and 10 % TMA within 0 to 9 %, and the H2O half-cycle on
  the same wafer within 6 % (14.2 predicted vs 13.4 s). One calibration cycle out of several hundred.
- The breakthrough curve has structure: a slow ramp while the flat regions saturate (1.95 s on the AR-100 wafer)
  followed by the steep rise when the holes finish; a two-stage model of the curve separates flat-area and
  feature-area saturation.

## Prior art found (2026-11-10, web and Google Patents via search; not exhaustive)
- Downstream QMS sensing of ALD half-cycles and "precursor dose end-point control" (Wei Lei, UMD thesis 2006,
  Rubloff group): QMS on tungsten ALD, end-point control on HfO2 in a furnace; no feature or pattern model.
- In-situ QMS/pressure breakthrough to confirm saturation on high-surface-area powders (Argonne fluidized-bed ALD;
  fixed-bed mesoporous silica with QMS breakthrough curves). Monitoring, not control; powders, not patterned wafers.
- Lam TW202113142A (2019): in-situ control of film properties by feedback on ampoule temperature, pressure and flow;
  extended dose time in selected cycles; no exhaust saturation signal, no pattern input.
- Navabi and Takoudis, arXiv 2602.18565 (2026): physics-informed Bayesian active learning tunes pulse times across
  runs with a Langmuir model in the GP kernel; no reactor transport, no pattern or loading input, no in-pulse signal.
- 3D-NAND gap-fill patents (US 11289374, US 12394670): state that large internal area demands high dose per cycle and
  that time and partial pressure trade off; no control method.
- Christofides group (UCLA, 2019-2024): multiscale CFD + kinetic Monte Carlo for ALD chamber design; no feature
  conformality, no exhaust end-point control.
- A Chinese filing on "ALD process management based on digital twinning" with a physics-informed neural network
  (assignee and claims not read). Must be read before any filing.
- Generic virtual-metrology run-to-run control patents (e.g. US 8437870B2) cover wafer-level recipe adjustment from
  predicted metrology, not in-pulse end-pointing from exhaust signals with a pattern model.

## What looks new (to be confirmed by a proper search of CPC C23C16/455 and C23C16/52)
The combination: pattern-aware physics prediction of the saturating pulse for the deepest feature + exhaust
breakthrough as the in-pulse confirmation + a model-derived lag/margin mapping exhaust breakthrough to feature-bottom
saturation + cycle-1 loading calibration transferred to the other half-cycle and all later cycles + drift detection
against the prediction.

## Next steps
1. Google Patents full-text search in C23C16/455 and C23C16/52 for: "end point" "atomic layer deposition" exhaust
   OR "residual gas" OR "mass spectrometer"; "loading" OR "pattern density" "dose"; read the Chinese digital-twin
   filing's claims.
2. Simulation evidence to strengthen: (a) a two-zone showerhead version where the edge zone is dosed independently
   using the predicted edge lag; (b) drift scenario: ampoule depletion (c0 falls 20 %) detected from the
   breakthrough shift within one cycle; (c) wrong-wafer scenario.
3. IP route: Penn Center for Innovation handles student inventions; Penn's policy on student-owned IP depends on
   use of Penn resources and funding (this work used a personal computer and open-source tools, but Soham is an
   enrolled student). Ask PCI before publishing the control section; a preprint is prior art against a later filing.
   A US provisional costs little and buys 12 months.
4. Keep the physics paper separate from the control claim, or publish the control section only after a provisional.
