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

## 6. Night three: eleven more scenes and the layout-to-recipe tool

Section 11 ("The reviewer's questions, answered with runs") restages the seven follow-up packages as scenes drawn by the
shared renderer from `window.ALD.n3` (`assets/body2.js`): `sticking` (saturating pulse against the sticking coefficient on
log axes, supply limits as thin ribbons, the published range as a translucent slab), `stepcov` (step coverage and
thickness after 50 cycles against pulse time), `cycle` (the four segments of the Al₂O₃ cycle as stacked boxes for three
recipes), `wph` (wafers per hour against wall-area ratio with the chamber run as a point), `exposure` (this model among
published saturating exposures, with Langmuir lines), `gakis` (coverage across the 200 mm wafer of the Fiji F200
approximation against the digitised measured points with error bars), `wafer3d` (the 3D chamber's 2,048-face wafer as a
triangulated height map, colour and height the deviation from the wedge, the pump port as an arc), `azimuth` (saturation
time against azimuth at three radii), `taper` (four cut-open frustums coloured by coverage at 10 s), `recess` (a word-line
stack: slit walls and 800 instanced recess segments coloured by coverage at 9 s) and `chemistry` (the rule per chemistry
with the TDMAH run). Each card links the audit figure (`assets/fig_*.png`, copied and quantised from `figures/` by
`build_assets.py --figs`). Section 10 is the layout-to-recipe tool (`assets/layout-tool.js`, markup from
`partials/layout-tool.html`): zones in, per-zone local rule out, the lateral-feed margin calibrated on the two non-uniform
runs and labelled as such, an SVG wafer map, CSV in and out; it has its own Playwright checks (`partials/README.md`).
The data file grew to about 380 KB; the 3D-chamber block is every third wafer face with its Delaunay triangulation.

## 7. The device cells as live scenes, and the GAA preset

The two device-geometry results of wp10 (`paper/wp/wp10_device_geometries.md`) were static figure cards in section 11; they
are now scenes drawn from the solver's own faces. `build_assets.py` (`device_data`, cached per-face snapshots in
`figures/device_geometries_scene.json`, rerun by `coupling/device_geometries.py`'s functions with the paper's parameters and
seeds: the GAA run reproduces 0.558 s, the NAND run 11.30 s) exports every exposed face as a tile: position in half-voxel
units and the face direction packed in three or four bytes, coverage quantised to one byte per snapshot, all base64 of
packed typed arrays (212 KB for the three sets; `data.js` is now about 585 KB). `gaacell` holds the IBM 2017 cell under TMA
(8,400 faces at 1.25 nm) and the narrowed cell after 3 nm of oxide under TiCl₄ (its 36,864 faces of 0.5 nm merged 2 × 2 to
1 nm tiles, from `figures/device_geometries_maps.json`), each at a third, two thirds and all of its time to 90 % everywhere;
`nandcell` is the 5 µm word-line cell with pillars cut through its row-0 pillar centre line (the faces of a two-voxel slab,
46,440 of 1.19 million, merged to 6,205 tiles: plates in groups of four along the recess, slit walls per oxide band, end
walls per recess) at 3, 6, 9 s and the end of the run, plus the area-weighted coverage of the pillar wall per recess and
side. The scenes follow the GPU solver's view: each GAA face is an unlit instanced plane with its normal into the gas and
front-face culling, so the camera sees the inside of the cavity from any angle; slots are drawn three times their height
and the trench above the stack at half height (stated in the caption) so the slot interiors are visible, and the cell is
drawn twice along its periodic axis. The NAND slab is drawn as boxes lying along the slit (mouth at the left, 4× wider
across the slit than deep, each face 100 nm long along the slit), with cylinder segments for the pillars and a line at the
slowest recess (86 of 100 at 4.3 µm). A `.vizctl` row under each card switches the snapshot (and, for the GAA card, the
cell); hover names the region (sheet, slot centre or edge, spacer wall, recess index and wall, pillar) with the position
in nm or µm and the coverage at every snapshot; tiles rise into place with the scroll progress. The WebGPU solver gained
a `gaaStack` builder (the same arithmetic as `gaa_stack`, periodic along x) and a "GAA, 3 sheets" preset whose live number
is the bottom slot centre (wp10's `gap_bottom_mid` region, read through a per-preset probe instead of the floor) against
the paper's 0.558 s; it runs with plain faces (area correction off) as the paper's run did, and at the reference budget the
Metal GPU gives 0.558 s.
