# Design notes for the live page

The engineering rationale behind `docs/index.html`, `docs/assets/scene.js`, `docs/assets/body.js` and `docs/build_assets.py`.
The page uses the same design system as Soham Kavathekar's other project pages (`site.css`, `site.js`, `three-kit.js`,
`viz.js`, copied unchanged from the sibling repositories; their DESIGN.md files record the research behind the stack: a
pinned hero with one scrubbed GSAP timeline, Lenis synced to ScrollTrigger, three.js 0.186 with bloom and depth of field on
the hero only, one shared forward renderer scissor-clipped to every body figure, per-theme palettes and rendering constants,
render-on-demand, static PNG fallbacks for reduced motion and software renderers, Lighthouse-checked). What is specific to
this repository is below.

## 1. What the page has to say

One number carries the project: a blanket wafer saturates in 0.345 s and the patterned wafer tuned from it needs 14.05 s.
Everything else is the argument for why that number is a one-line rule, how it was checked, which knobs move it and which do
not. The page follows the paper's order (model, rule, validation, inside the hole, DOE, water/plasma/spatial, robustness) with
two additions the paper only reports in text: a calculator that runs the rule and the Gaussian-process surrogate in the
browser, and the non-uniform layout runs.

## 2. The hero: a 16-second pulse you scrub

The stage is the actual case `reactor/showerhead_pulse_patterned_AR100_long`: a 300 mm wafer in a 400 mm chamber, 15 mm gap
(drawn 3× so it reads), showerhead above, outlet ring at r = 200 mm. Three things are data, not decoration:

- **Gas.** 2,600 instanced particles fill the gap at fixed positions; every frame each one takes the TMA mole fraction at its
  (r, z) from the exported field (40 radial × 10 vertical bins, 29 time frames from 0.05 to 22 s, linearly interpolated in
  time) and is coloured and scaled by it. At 5 s the gas over the wafer sits at 9 to 16 % of the showerhead value: the holes are
  draining it, and the particles show it.
- **Wafer top.** A polar mesh whose vertex colour is the flat-surface coverage at that radius and time (221 times × 30 radii).
- **The cut-open hole.** A half cylinder beside the chamber (100 nm × 10 µm drawn at 1:200,000) coloured along its height by
  the coverage-against-depth history of one hole at the wafer centre, computed by the standalone feature solver driven by the
  chamber's face concentration (`aldconf.simulate_driven`; it reproduces the embedded boundary condition's 14.05 s). A ring
  marks the 50 % front.

Scroll maps to time (a short hold at 0, then 16 s over the track); the readout shows time, flat and hole-bottom coverage, and
monolayers delivered to the patterned area, which reaches 1.2 when the bottoms saturate: the rule, read off the scene. The
camera path: overview, low across the wafer, toward the hole, settle with both in frame. Under reduced motion the scene sits at
16 s.

## 3. Body scenes (`assets/body.js`), all from `assets/data.js`

| Placeholder | Scene | Data |
|---|---|---|
| `stack` | the chamber exploded (stage, wafer, floor, outlet, wall, showerhead) and a magnified patch of the hole array | geometry constants |
| `map` | saturating pulse against TMA fraction, blanket and AR 20/50/100, each with its supply-limit ribbon | `figures/pulse_time_map.csv` |
| `wiwnu` | hole-bottom coverage across the wafer for a pulse stopped at 8, 10, 12, 14 s | hero coverage history |
| `layout` | saturation time against radius for the uniform wafer, the dense-centre layout and the AR ramp, with the patterns | `figures/pattern_layout.csv` |
| `rule` | precursor utilization against total reactive area: chamber runs, fluidized bed, rotary drum | `coupling/validate_supply_rule.py` numbers |
| `pillar` | the PillarHall saturation profile, the fitted model and two alternative sticking coefficients | Zenodo 3974438 + `aldconf` at the fitted parameters |
| `hole` | the cut-open hole and its coverage-against-depth ribbon at the sticky step's time | hero hole history |
| `doe` | the JMP response surface over TMA fraction and carrier flow with the 25 runs as spheres | `figures/doe_results.csv`, `doe_rsm_report.txt` |
| `coo` | precursor delivered per saturating pulse against pulse time with the surrogate's Pareto front | `figures/doe_pareto.csv` |
| `water` | water coverage and outlet signal for four desorption constants | `reactor/showerhead_h2o_blanket_kdes*` |
| `peald` | largest aspect ratio reachable in 5 s against wall recombination probability | `aldconf`, as in `coupling/peald_ceiling.py` |
| `spatial` | temporal against spatial cycle time from the rule | arithmetic in the browser from the rule's constants |
| `robust` | the saturating pulse under every removed idealization | edge-ring, temperature, hole-size cases; water case |

## 4. The calculator

`#try` evaluates the rule from the pattern (d, AR, hole density), the showerhead supply (TMA fraction, carrier flow, pressure,
temperature) and N₀, with U the showerhead face velocity at chamber conditions. The surrogate is the scikit-learn Gaussian
process from `analysis/doe_gp.py`, exported to `figures/doe_gp_python.json` (linear mean, ARD RBF kernel, training inputs,
alpha, K⁻¹, normalisation) and evaluated in plain JavaScript; the same arithmetic as `python -m aldconf surrogate`. It warns
when the inputs leave the design box or the pattern differs from the one the DOE used.

## 5. Build and QA

`docs/build_assets.py --data` writes `assets/data.js` (about 300 KB) from the OpenFOAM time directories (not committed), the
committed `figures/*.csv` and the `aldconf` solver; `--snap`, `--hero-png`, `--qa` and `--lighthouse` drive
`docs/tools/page_qa.py`, a Playwright-for-Python port of the sibling repositories' `page-qa.mjs` (headless Chromium with
software WebGL and `window.__forceGL`, so frame timing is not representative; the Lighthouse run has WebGL disabled and
measures the page's static path, as on the other pages: 100 / 100 / 100 / 100 on the desktop preset, `docs/qa/README.md`).
One change to the shared runtime, in this repository's copy of `site.js`: the software-renderer probe now asks for a WebGL
context with `failIfMajorPerformanceCaveat`, which a software renderer refuses at once, instead of creating one and reading
its renderer string; creating a SwiftShader context cost 2.4 s of main thread before first paint in the audit environment.
Pinned CDN builds: three.js 0.186.1, GSAP 3.15.0, Lenis 1.3.27.
