# Layout-to-recipe tool: partial for the live page

A section for `docs/index.html` that turns a zoned die map into a TMA pulse recipe with the supply rule, sized on the most
demanding zone, with the lateral-feed margin measured in the two non-uniform chamber runs (`paper/assumptions.md`, Night 3,
part 1; `figures/pattern_layout.csv`).

Files in this work package:

| File | What it is |
|---|---|
| `docs/partials/layout-tool.html` | the `<section class="story alt" id="layout-tool">` with its own scoped `<style>` block |
| `docs/assets/layout-tool.js` | classic script, no dependencies; `window.LayoutTool = { compute, parseCsv, toCsv, presets, chamberDefault, init }` |
| `tests/test_layout_tool.py` | the arithmetic in Python against the logged runs, and the JavaScript against the Python through Playwright (skipped without it) |
| `docs/partials/qa_{dark,light}_{1280,375}.png` | section screenshots from the headless QA run below |

## Inserting it

1. Paste the whole of `layout-tool.html` into `docs/index.html` after the `#run` section ("09 · Reproduce it") and before
   `<footer>`. It is numbered "10 · Layout to recipe". If a section without the `alt` background is wanted at that spot,
   drop `alt` from the opening tag; nothing else depends on it. The `<style>` at the top of the partial is scoped to
   `#layout-tool` and may be moved into the inline `<style>` block in the head instead.
2. Add one script tag after `assets/site.js` in the head:
   `<script defer src="assets/layout-tool.js"></script>`
   The script initialises on `siteready`; since `site.js` is an earlier deferred script and dispatches `siteready`
   synchronously, `layout-tool.js` finds `window.SITE.ready` already true and initialises at once, with `DOMContentLoaded`
   as the fallback when it is loaded some other way. It re-renders on the page's `themechange` event and uses only the
   CSS variables from `site.css`.
3. Add the nav link `<a href="#layout-tool" class="hide-sm">Layouts tool</a>` to `nav .links` (after "Try it" keeps the page order).
4. Nothing in `data.js`, `site.js`, `site.css` or the other scripts changes. The partial reuses the `.calc`, `.calc label`,
   `.calc .out`, `.tablewrap`, `table`, `.note`, `.tag`, `tr.hl` and `.warn` styles that `index.html` and `site.css`
   already define, so it must sit inside `index.html` (or another page carrying the same head styles).

`CLAIMS.md` entries the integrator may want to add: the recipe pulse is the supply rule of `#try` (`aldconf rule`); the
margin band 0.82–0.87 and 0.98 are the ratios 11.85/14.39, 12.50/14.36 and 14.05/14.39 from the Night 3 table; the bare-field
0.40–0.50 s and the 39 % depth spread are from the same log entry.

## What the tool computes

Per zone (concentric rings `r_in`–`r_out`, hole diameter `d`, aspect ratio `AR`, hole density per m²):
`A_wall/A_wafer = density × π d × AR d`, local rule `t = 1.2 (1 + A_wall/A_wafer) N₀ / (U c₀)` with
`c₀ = x p / (k_B T)` and `U = Q_sccm × 1.6667e-8 × (101325/p) × (T/273.15) / (π 0.15²)`, the arithmetic of the `#try`
calculator unchanged. The recipe pulse is the largest local rule; the wafer-average rule uses the area-weighted wall ratio
(area outside any zone counts as bare); the under-dose factor is their ratio. The expected band multiplies the recipe by
0.82–0.87 (non-uniform layout) or 0.98 (uniform: every zone patterned, same wall ratio within 1 %, whole wafer covered),
and the UI states that this is a calibration on two runs. Bare zones show the measured flat saturation next to an array
(0.40–0.50 s) rather than the rule, which does not apply to a blanket surface. Warnings: AR outside 20–100, TMA outside
0.5–10 %, overlapping zones, zones past 150 mm, uncovered wafer area.

Presets. "Dense centre" is the Night 3 layout (a): AR 100 and 1.78e13 /m² for r < 75 mm, the 75–100 mm linear ramp to zero
as one zone at its area-weighted mean density 8.48e12 /m² (0.476 ×, not 0.5 ×, because the outer, emptier radii carry more
area; this is what reproduces the logged wafer-average ratio 19.2 and 5.1 s), bare beyond. "Uniform AR-100 wafer" is the
reference run. "AR 100 → 60 ramp" is layout (b) as five 30 mm zones at AR 100, 90, 80, 70, 60.

## QA

Harness: a throwaway page in the session scratch area loading the real `site.css`, the head `<style>` of `index.html`,
the real `site.js`, the partial and `layout-tool.js`; headless Chromium through Playwright for Python
(`.venv/bin/python`), 1280 × 900 and 375 × 812 (device scale 2), dark and light, with a `layout-shift`
PerformanceObserver and a `scrollWidth` check. Run on 2026-11-10.

| Check | Result |
|---|---|
| Console errors / page errors, 4 contexts plus preset, CSV, add/remove, theme-toggle and bad-input interactions | 0 |
| Cumulative layout shift, 4 contexts | 0 |
| Horizontal overflow at 375 px | none (scrollWidth 375) |
| Dense centre: zone 1 wall ratio / local rule | 55.9 / 14.38 s (log: 14.39 s with U rounded to 0.0501 m/s) |
| Dense centre: wafer-average wall ratio / rule | 19.2 / 5.09 s (log: 19.2 / 5.10 s) |
| Dense centre: expected band | 11.8–12.5 s (measured 11.85 s) |
| Dense centre: under-dose factor | 2.82× |
| Uniform AR 100: rule / expected | 14.38 s / 14.09 s (measured 14.05 s) |
| AR ramp: rule / expected / wafer-average | 14.38 s / 11.8–12.5 s (measured 12.50 s) / 10.42 s, wall 40.3 (log: 10.62 s, 41.0 for the continuous ramp; the five-zone staircase is 1.7 % lower) |
| AR ramp: AR 60 edge zone local rule | 8.73 s (log: 8.76 s at AR 60.2); depth spread 39 % (log: 39 %) |
| CSV loader | the pasted dense-centre map gives 14.38 / 5.09 s; Table → CSV of the ramp preset loads back to the same five rows |
| `tests/test_layout_tool.py` | 6 passed (3 Python-only, 3 JavaScript against Python to 1e-9 relative) |

Screenshots: `qa_dark_1280.png`, `qa_light_1280.png`, `qa_dark_375.png`, `qa_light_375.png` (section only, dense-centre
preset). Not yet run: the full-page `docs/tools/page_qa.py` and Lighthouse on the integrated `index.html`; rerun both
after insertion (`python docs/build_assets.py --qa --lighthouse`).
