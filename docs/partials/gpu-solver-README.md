# 3D voxel Monte Carlo in the browser: partial for the live page

A section for `docs/index.html` that runs the feature-scale solver of `lib/aldmc` (3D voxel test-particle Monte
Carlo, molecular flow) on the visitor's own graphics card through WebGPU, with a plain-JavaScript fallback, and draws
the feature as its exposed voxel faces coloured by coverage. The six presets are the paper's voxel runs
(`paper/wp/wp9_voxel_mc.md`): straight 100 nm holes of AR 20 / 50 / 100, the tapered hole with a 60 % bottom, the
AR-50 slit and the word-line slit w100_Lr365; the reference t90 of each is printed beside the live number.

Files in this work package (nothing else in the repository changed):

| File | What it is |
|---|---|
| `docs/partials/gpu-solver.html` | the `<section class="story alt" id="gpu-solver">` ("12 · Run the 3D solver in your browser"): scoped `<style>`, markup, and a classic inline `<script>` that wires the controls and imports the modules on Run |
| `docs/assets/aldmc-gpu.js` | ES module: geometry builders (`hole`, `taperedHole`, `slit`, `recessedSlit`, `fromArray`), the shared grid preprocessing (`buildGrid`), the WebGPU backend, `createSolver({geometry, params, backend})`, `requestGpu`, `compileCheck`, `timeToLevel`, `numberDensity`, `meanSpeed` |
| `docs/assets/aldmc-gpu.wgsl` | the compute shader: `trace_main` (one invocation per particle) and `update_main` (one per face) |
| `docs/assets/aldmc-cpu.js` | the CPU fallback: the same tracer in plain JavaScript on typed arrays, plus the xoshiro128** / splitmix helpers |
| `docs/partials/qa_gpu_{dark,light}_{1280,375}.png` | section screenshots from the headless QA run below (Metal WebGPU, AR-20 hole paused mid-run) |

## Inserting it

1. Paste the whole of `gpu-solver.html` into `docs/index.html` after the `#layout-tool` section ("10 · Layout to
   recipe"; section 11 is `#night3`, which sits before it in the file) and before `<footer>`. It is numbered "12". The
   `alt` background matches the layout tool's; drop `alt` if the integrator wants it to alternate.
2. No script tag in the head. The partial's own button does `import('./assets/aldmc-gpu.js')`, `import('three')` and
   `import('three/addons/controls/OrbitControls.js')` when Run is clicked; the last two resolve through the import
   map already in the head (three.js 0.186.1), and the WGSL is fetched by the module relative to its own URL
   (`new URL('./aldmc-gpu.wgsl', import.meta.url)`). Before the click the section costs its markup and about 14 KB of
   inline script and style; no module, shader, buffer or WebGL context exists, and nothing runs on load (checked:
   no `three` or `aldmc` resource entry before the click in every QA context).
3. Nav link: `<a href="#gpu-solver" class="hide-sm">3D solver</a>` in `nav .links` after "Layouts tool".
4. The partial reuses `.calc`, `.calc label`, `.calc input`, `.calc .out` from the head `<style>` of `index.html` and
   `.story`, `.wrap`, `.eyebrow`, `.lead`, `.vizhead`, `.legend`, `.vizfoot`, `.note`, `.btn`, `.btn.ghost` and the
   CSS variables from `site.css`, so it must sit inside `index.html`. It listens to the page's `themechange` event to
   recolour the 3D view and reads `window.SITE.softGL` (if present) to refresh the view less often on a software
   WebGL renderer. It does not touch `data.js`, `site.js`, `body*.js`, `scene.js`, `viz.js` or `three-kit.js`.
5. GitHub Pages serves `.wgsl` as `application/octet-stream`; the module reads it with `response.text()`, which does
   not care. A local server that refuses unknown extensions would need `.wgsl` mapped to `text/plain`.
6. After insertion, rerun the page QA and Lighthouse (`python docs/build_assets.py --qa --lighthouse`). The partial
   adds nothing to the load path, so the Lighthouse scores should hold; the WebGPU work only starts on a click.

## Architecture

- **Geometry** (`aldmc-gpu.js`): the builders mirror `src/aldconf/voxel.py` voxel for voxel. `hole` gives every slice
  exactly round(π r²) gas voxels (the voxel centres nearest the axis, ties broken by numpy's ravel order), with the
  frustum radius falling linearly for the tapered hole; `slit` is periodic over 4 voxels along x; `recessedSlit`
  rounds g_r and p_r to whole voxels and records the actual values in `meta`. Occupancy is `Uint8Array` with the
  library's layout `occ[(k ny + j) nx + i]`.
- **Preprocessing** (`buildGrid`, both backends): Chebyshev distance field by a 26-neighbour BFS (0 solid, the mouth
  and far planes and any non-periodic boundary count as walls), exposed faces as (gas voxel, direction 0 +x … 5 −z)
  with a per-voxel block of six face ids, projected-area face weights w = |n · e_face| from the first moment of the
  occupancy in a cone kernel (1 − r/R, R = 5 voxels), stick and loss probabilities s w and r w per face, the mouth
  voxel list, mouth area and weighted site area. The weighted inventory comes out at 0.99995 (AR-20 cylinder), 1.0016
  (AR-100) and 1.0026 (frustum) of the smooth areas, the figures `lib/aldmc/README.md` and wp9 quote.
- **WebGPU backend**: six storage buffers and one uniform (the default limit is eight storage buffers per stage):
  `dist` (u8 packed four to a u32), `topo` (i32: face base per voxel, then face ids per block), `mouth` (u32),
  `theta` (f32 per face), `prob` (f32: stick then loss per face) and `acc` (atomic u32: eight counters, then hits
  and sticks per face). `trace_main` runs one invocation per particle, workgroup size 64, up to 65 535 × 64 =
  4.19 M particles per dispatch: it seeds xoshiro128** from (seed, step, invocation), picks a mouth voxel and a
  position in it, draws a cosine direction, and flies: free flights of D − 1 voxel lengths through the distance
  field, an exact DDA for the last voxels, a wall hit incrementing the face's atomic counter, then stick / loss /
  diffuse re-emission about the facet normal with the coverage frozen for the step. Positions are kept as (voxel
  index, fraction in [0, 1]) so f32 precision is uniform down a 1600-voxel hole. The fates and the flights (64-bit
  through a low word plus an overflow count) go to the counter block. `update_main`, one invocation per face,
  applies θ ← 1 − (1 − θ) exp(−s hits w_mol / (N₀ h²)) and clears the counters. One command buffer per step: zero
  the counters, trace, copy the counters to a staging buffer, update, copy θ to the staging buffer, submit, `mapAsync`;
  `step()` resolves to the counts with the coverage mirror refreshed. Integer atomics make the hit counts, and
  therefore the coverage, independent of scheduling for a given seed.
- **CPU backend** (`aldmc-cpu.js`): a line-for-line port of `aldmc_t::trace` and `aldmc_t::step` with doubles in
  voxel units, blocks of 64 particles with one xoshiro128** stream per (seed, step, block), and a cooperative yield
  to the event loop every 40 ms of work so the page stays responsive. It is chosen when `navigator.gpu` is missing,
  the adapter or device request fails, or the shader does not compile; `solver.fallbackReason` says why.
- **Flight cap**: `maxSegments` = 2 000 000 voxel segments per particle (the library's default is 1e8); a particle
  reaching it is counted "left". None did in any QA run (the AR-100 hole averages about 2 000 segments per particle
  when saturated; the word-line slit 2 400 at 3 s).
- **Random numbers**: xoshiro128** (32-bit words, which JavaScript and WGSL both handle at full speed) seeded through
  a splitmix-style finaliser from (seed, step, block or invocation), the shape of the library's `stream_seed`. The
  library uses xoshiro256**, so the streams are not bitwise the library's; the statistics are the same.
- **View** (inline in the partial): `THREE.InstancedMesh` of unit planes, one per exposed face, positioned at the
  face centre with the normal pointing into the gas and front-face culling, so the camera always looks at the inside
  of the cavity from any angle without a cut; depth is compressed when the feature is more than four times deeper
  than wide (stated in the caption); colours are the page's data palette, grey (θ 0) → cool (0.5) → green (1), with a
  per-direction shade baked in, re-read from the CSS variables on `themechange`; OrbitControls, render on demand, a
  ResizeObserver, `data-lenis-prevent` so the wheel zooms the view. A slice-mean profile is drawn as an inline SVG.
- **Driver** (inline in the partial): dt defaults to the preset's reference t90 / 200, scaled with N₀, pressure and
  temperature (not with s); n = p / (k_B T) and v̄ = √(8 k_B T / π m) with m = 72.09 g/mol; the loop steps until the
  floor reaches 99 % or four reference saturation times, and records t90 by linear interpolation
  (`time_to_level`'s rule). Changing s, N₀, temperature, voxels across or the preset rebuilds the solver on the next
  Run; pressure, M and dt are read every step. On the CPU fallback M is lowered to 5 000 once, with a message.

## QA

Harness: a page in the session scratch area (`page_test.html`, served over `http://localhost`, which counts as a
secure context) loading the real `site.css`, the import map and head `<style>` of `index.html`, the CDN GSAP / Lenis
builds, the real `site.js` and the partial. Playwright for Python 1.65 (`.venv/bin/python`), Chromium 158.0.8127.18
(`~/Library/Caches/ms-playwright/chromium-1254` and `chromium_headless_shell-1254`). Run on 2026-11-11 on an Apple M4
(10 cores).

WebGPU headlessly, what worked:

| Launch | `navigator.gpu` | Adapter |
|---|---|---|
| headless shell (Playwright default), any flags, `about:blank` | absent (not a secure context) | – |
| headless shell, `--enable-unsafe-webgpu`, page over `http://localhost` | present | `google / swiftshader` (GPU API emulated on the CPU) |
| full Chromium (`chrome-mac-arm64/Google Chrome for Testing.app`), `--headless=new --enable-unsafe-webgpu --ignore-gpu-blocklist`, page over `http://localhost` | present | `apple / metal-3` (the real GPU) |

The `--enable-features=Vulkan,UseSkiaRenderer --use-angle=swiftshader` variants add nothing on this Mac. Limits
reported by the Metal device: workgroups per dimension 65 535, invocations per workgroup 256, storage binding and
buffer size 4 GiB, storage buffers per stage 8. `compileCheck` (`createShaderModule` + `getCompilationInfo`) returns
no messages on either adapter.

Verification, AR-20 hole at the paper's conditions (s 5·10⁻³, N₀ 4·10¹⁸ m⁻², 573 K, 0.3 Pa TMA: n 3.792·10¹⁹ m⁻³,
v̄ 410.2 m/s; 16 voxels across, h 6.25 nm, grid 18 × 18 × 321, 20 681 faces, 201 mouth voxels; seed 1). Reference:
`paper/wp/wp9_voxel_mc.md` table (a), t90 1.148 s / t95 1.316 s at dt 5 ms, M 25 000 (1.153 / 1.320 s at dt 20 ms);
Machball converged 1.158 / 1.334 s.

| Backend, budget | floor t90 | t95 | vs 1.148 s | particles | throughput | wall |
|---|---|---|---|---|---|---|
| WebGPU, Apple Metal, dt 5 ms, M 25 000 (the reference budget), module harness | 1.144 s | 1.312 s | −0.3 % | 6.6 M (263 steps) | 461 k particles/s, 327 M flights/s | 14.3 s |
| WebGPU, Apple Metal, same, through the page's Run button (view updating) | 1.144 s | – | −0.3 % | 5.9 M | 500 k particles/s, 338 M flights/s | 12.2 s to t90 |
| WebGPU, SwiftShader (headless shell), dt 20 ms, M 5 000 | 1.165 s | 1.337 s | +1.5 % | 335 k | 29.8 k particles/s, 21 M flights/s | 11.2 s |
| CPU JavaScript, dt 20 ms, M 4 000, module harness | 1.128 s | 1.294 s | −1.8 % | 260 k | 29.3 k particles/s, 20.9 M flights/s | 8.9 s |
| CPU JavaScript, dt 5 ms, M 5 000 (the fallback's automatic M), through the page's Run button, software WebGL view updating | 1.142 s | – | −0.5 % | 1.4 M | 42 k particles/s, 31 M flights/s | 25 s to t90 |

Both backends are within the 10 % asked for; the GPU at the reference budget is within the seed-to-seed scatter of
the library (0.2–0.4 %). Flights per particle: 710 (GPU) / 713 (CPU) against the library's 720 for this hole. Every
run conserved the particle count (escaped + stuck + lost + left = launched) with zero "left".

Other presets on the Metal GPU (short runs, `run_presets_gpu.py`): slit AR 50 at dt 12.5 ms, M 25 000 reaches t90 at
1.0829 s against the paper's 1.083 s (grid 4 × 18 × 801, 6 464 faces, 672 k particles/s); the word-line slit at h 5 nm
builds the paper's 4 × 168 × 1001 grid with 124 880 faces (100 recesses per wall by the builder's arithmetic; wp9
quotes 98 with the same face count), 415 k particles/s at 2 355 flights per particle, slit floor 27 % at 3 s; the
tapered AR-100 hole 82 919 faces, 745 k particles/s; the straight AR-100 hole 102 601 faces, 1.14 M particles/s on the
empty hole (567 flights per particle at 1 s). Grid preprocessing on the main thread: 0.3 s (AR 20), 1.0–1.9 s (AR
100), 1.3–2.6 s (word line); the face-weight loop yields to the event loop and reports progress.

Page checks (Metal WebGPU unless stated; dark and light, 1280 × 900 and 375 × 812 at device scale 2; Run, 5 s, Pause,
orbit drag, two theme toggles, preset switch, Reset):

| Check | Result |
|---|---|
| Console errors / page errors, 4 contexts + the CPU context | 0 |
| Resources matching `three` or `aldmc` before the click | none in every context |
| Cumulative layout shift | 0 at 1280 (both themes); 0.0076 at 375 (both themes), from the readout cards' first numbers |
| Horizontal overflow at 375 px | none (scrollWidth 375) |
| 3D view drawn (canvas sized, placeholder hidden) after Run | yes in all contexts, including the CPU context under software WebGL |
| Backend line | "WebGPU · apple metal-3" with the device limits; CPU context: "CPU · JavaScript, one thread (WebGPU: navigator.gpu is not available …)" with M lowered to 5 000 and the message shown |
| Preset switch while paused | dt follows the preset (slit → 12.5 ms), title and reference line update, solver marked for rebuild |

Screenshots: `qa_gpu_dark_1280.png`, `qa_gpu_light_1280.png`, `qa_gpu_dark_375.png`, `qa_gpu_light_375.png` (section
only, AR-20 hole paused at 0.6–0.7 s of simulated time on the Metal GPU).

## Known limits

- WebGPU availability is the browser's: Chrome / Edge 113+ and Safari 26 on macOS and iOS have it; Firefox ships it on
  Windows and keeps it behind `dom.webgpu.enabled` elsewhere; it needs a secure context (https or localhost). Without
  it the CPU path runs at about 30 k particles/s on this machine against 0.45–1.1 M on the Metal GPU, so the AR-100
  cases take many minutes on the CPU.
- No film narrowing and no gas-phase collisions, by construction (as in `lib/aldmc`). Only the library's defaults
  are offered: area correction on, facet re-emission, impingement update, r = 0.
- Statistics: the floor t90 scatters by a few tenths of a per cent between seeds at the reference budgets; at the
  CPU fallback's M = 5 000 the scatter is larger (the 1.128 s above is one seed). The time step enters linearly
  through the coverage frozen within a step (0.6 % at the reference dt; the auto dt scales with N₀, p and T but not
  with s). Per-face coverages are noisy at a few hits per face per step.
- GPU: one dispatch per step, so M is capped at 4.19 M particles; the flight cap is 2 M segments per particle. On
  SwiftShader (or any CPU-emulated adapter) the "GPU" path is no faster than the JavaScript one. A device loss
  surfaces as an error message in the status line; Reset and Run start over.
- The 3D view is an unlit coverage map with front-face culling, which hides the near wall so the cavity is visible;
  looking from directly behind a slit shows nothing but the far wall's edge. Depth is compressed for deep features.
- Not yet run: the full-page `docs/tools/page_qa.py` and Lighthouse on the integrated `index.html`.
