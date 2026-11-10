/* Night-three scenes for the ALD page, drawn by the shared Viz renderer (assets/viz.js) from window.ALD.n3 (assets/data.js):
   sticking   saturating pulse against the TMA sticking coefficient: blanket, flat surface, hole bottoms, supply limits, literature band;
   stepcov    step coverage and thickness after 50 cycles against pulse time for AR 20, 50, 100 with the rule marked;
   cycle      the Al2O3 cycle as stacked segments for three recipes; wph: wafers per hour against wall-area ratio from the rule;
   exposure   saturating TMA exposure against temperature: this model and published values;
   gakis      the Fiji F200 reproduction: final coverage across the 200 mm wafer against the measured growth per cycle;
   wafer3d    the 3D chamber's hole-bottom saturation time over the wafer as a height map (deviation from the wedge) with the pump port;
   azimuth    saturation time against azimuth at three radii;
   taper      four cut-open holes with narrowing bottoms coloured by coverage at 10 s;
   recess     a word-line stack: slit and lateral recesses coloured by coverage at 9 s;
   chemistry  the rule's pulse per chemistry with the TDMAH chamber run. */
import { Viz } from './viz.js';
import { THREE } from './three-kit.js';
const D = window.ALD, N = D.n3; const lg = Math.log10; const f = (v, d = 0) => (+v).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d });
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const sphere = (v, A, x, y, z, col, r = .06, data = null) => { const m = new THREE.Mesh(new THREE.SphereGeometry(r, 14, 12), Viz.mat(col, { emissive: col, emissiveIntensity: .55 })); m.position.copy(A.map(x, y, z)); m.userData = data || {}; v.rig.add(m); v.pickables.push(m); return m; };
const fadeIn = (objs, p, from = .6) => objs.forEach(o => { o.material.transparent = true; o.material.opacity = clamp((p - from) / (1 - from), 0, 1); });
const viridis = t => { const s = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]]; t = clamp(t, 0, 1) * 4; const i = Math.floor(t), fr = t - i, a = s[i], b = s[Math.min(i + 1, 4)]; return [(a[0] + (b[0] - a[0]) * fr) / 255, (a[1] + (b[1] - a[1]) * fr) / 255, (a[2] + (b[2] - a[2]) * fr) / 255]; };
if (!N) console.warn('night-three data missing');

/* ---------- sticking: the rule against s ---------- */
N && Viz.view('[data-viz="sticking"]', v => {
  const S = N.sticking, R = N.sticking_refs; const bl = S.filter(q => q.kind === 'blanket').sort((a, b) => a.s - b.s), pa = S.filter(q => q.kind === 'patterned').sort((a, b) => a.s - b.s);
  const A = Viz.axes(v, { size: [4.4, 2.4, 2.4], x: { lo: lg(4e-4), hi: lg(1.1e-2), ticks: [lg(5e-4), lg(1e-3), lg(2e-3), lg(5e-3), lg(8.5e-3)], fmt: t => Math.pow(10, t).toExponential(0).replace('e-', 'e-'), label: 'TMA sticking coefficient s (log)' }, y: { lo: lg(.15), hi: lg(25), ticks: [lg(.2), lg(1), lg(10)], fmt: t => `${+Math.pow(10, t).toFixed(1)} s`, label: 'time to 90 % coverage (log)' }, z: { lo: -.5, hi: 2.5, ticks: [] } });
  const series = [['blanket wafer (worst face)', 'grey', bl.map(q => lg(q.s)), bl.map(q => lg(q.t90_worst)), 0, bl], ['patterned wafer: flat surface between holes', 'warm', pa.map(q => lg(q.s)), pa.map(q => lg(q.t90_flat)), 1, pa], ['patterned wafer: hole bottoms (worst face)', 'green', pa.map(q => lg(q.s)), pa.map(q => lg(q.t90_worst)), 2, pa]];
  const rbs = series.map(([name, col, X, Y, z, rows]) => { const rb = Viz.ribbon(v, A, X, Y, { color: col, width: .22, z, name }); rb.mesh.userData = { ...rb.mesh.userData, name, rows }; const sph = rows.map((q, i) => sphere(v, A, X[i], Y[i], z, col, .055, { q, name })); const lab = Viz.label(v, name, [A.map(lg(4e-4), 0, z).x - .2, A.map(0, Y[0]).y + .25, A.map(0, 0, z).z], { px: 18, scale: .08, color: col, weight: 700 }); return { rb, sph, lab }; });
  // supply-limit times as thin ribbons
  [[pa[0].t_supply, 2, 'patterned supply limit 11.99 s'], [bl[0].t_supply, 0, 'blanket supply limit 0.21 s']].forEach(([ts, z, name]) => { const rb = Viz.ribbon(v, A, [lg(4e-4), lg(1.1e-2)], [lg(ts), lg(ts)], { color: 'cool', width: .07, z, skirt: false }); rb.draw(1); rb.mesh.userData = { ...rb.mesh.userData, name, ts }; Viz.label(v, name, [A.map(lg(1.1e-2), 0, z).x + .9, A.map(0, lg(ts)).y, A.map(0, 0, z).z], { px: 16, scale: .07, color: 'cool' }); });
  // literature band and reference values
  const bw = A.map(lg(R.arts_band[1]), 0, 0).x - A.map(lg(R.arts_band[0]), 0, 0).x; const band = new THREE.Mesh(new THREE.BoxGeometry(bw, 2.4, 2.4), Viz.mat('green', { opacity: .1 })); band.position.set((A.map(lg(R.arts_band[1]), 0, 0).x + A.map(lg(R.arts_band[0]), 0, 0).x) / 2, 0, 0); v.rig.add(band);
  [[R.ylilammi, 'Ylilammi 2018: 5.7e-3'], [R.vandalon, 'Vandalon 2016: 4e-3'], [R.pillarhall, 'PillarHall fit: 8.5e-3'], [R.runs, 'chamber runs: 5e-3']].forEach(([s, name], i) => { Viz.line(v, [[A.map(lg(s), 0, -.5).x, A.map(0, lg(.15)).y, A.map(0, 0, -.5).z], [A.map(lg(s), 0, -.5).x, A.map(0, lg(25)).y, A.map(0, 0, -.5).z]], 'muted'); Viz.label(v, name, [A.map(lg(s), 0, -.5).x, A.map(0, lg(25)).y + .14 + (i % 2) * .16, A.map(0, 0, -.5).z], { px: 14, scale: .065, bg: false }); });
  Viz.label(v, 'Arts et al. 2019: 0.5 to 2e-3', [A.map(lg(1e-3), 0, -.5).x, A.map(0, lg(.15)).y - .3, A.map(0, 0, -.5).z + .3], { px: 15, scale: .07, color: 'green' });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; fadeIn(r.sph, p, .5); }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.q) return `<b>${u.name}, s = ${u.q.s}</b>t90 centre ${f(u.q.t90_centre, 2)} s · worst ${f(u.q.t90_worst, 2)} s${u.q.t90_flat ? ` · flat ${f(u.q.t90_flat, 2)} s` : ''} · ${u.q.kind === 'patterned' ? `${f(u.q.t90_worst / u.q.t_supply, 2)}× supply limit` : `${f(u.q.t90_worst / u.q.t_supply, 2)}× supply limit`} · utilization ${f(u.q.util * 100, 0)} %`; if (u.ts) return `<b>${u.name}</b>`; return null; };
  v.rot = -.35; v.tilt = .12;
}, { camera: [0, 2.2, 7.8], target: [0, -.1, 0], fov: 32 });

/* ---------- stepcov: step coverage and thickness after 50 cycles ---------- */
N && Viz.view('[data-viz="stepcov"]', v => {
  const C = N.stepcov.curves; const cols = { 20: 'green', 50: 'warm', 100: 'rose' }; const rules = {}; D.map.rows.filter(r => r.kind === 'patterned' && Math.abs(r.x - .01) < 1e-9).forEach(r => rules[r.AR] = 1.2 * r.t_supply);
  const A = Viz.axes(v, { size: [4.4, 2.2, 2.6], x: { lo: lg(.08), hi: lg(25), ticks: [lg(.1), lg(1), lg(10)], fmt: t => `${+Math.pow(10, t).toFixed(1)} s`, label: 'TMA pulse (log)' }, y: { lo: 0, hi: 1.05, ticks: [0, .5, .9, 1], fmt: t => `${Math.round(t * 100)} %`, label: 'step coverage · thickness / 5.15 nm' }, z: { lo: -.5, hi: 2.5, ticks: [] } });
  const rbs = C.map((c, i) => { const X = c.t.map(lg); const rb = Viz.ribbon(v, A, X, c.sc, { color: cols[c.AR], width: .22, z: i, name: `AR ${c.AR}` }); rb.mesh.userData = { ...rb.mesh.userData, c, kind: 'sc' }; const th = Viz.ribbon(v, A, X, c.thick.map(t => t / 5.15), { color: cols[c.AR], width: .08, z: i, skirt: false }); th.mesh.userData = { ...th.mesh.userData, c, kind: 'th' }; th.mesh.material.opacity = .7; th.mesh.material.transparent = true;
    const ru = rules[c.AR]; if (ru) { Viz.line(v, [[A.map(lg(ru), 0, i).x, A.map(0, 0).y, A.map(0, 0, i).z], [A.map(lg(ru), 1.05, i).x, A.map(0, 1.05).y, A.map(0, 0, i).z]], 'muted'); }
    const lab = Viz.label(v, `AR ${c.AR} · rule ${f(ru, 1)} s`, [A.map(lg(25), 0, i).x + .6, A.map(0, 1).y + .1, A.map(0, 0, i).z], { px: 18, scale: .08, color: cols[c.AR], weight: 700 }); return { rb, th, lab }; });
  Viz.line(v, [[A.map(lg(.08), .9, -.5).x, A.map(0, .9).y, A.map(0, 0, -.5).z], [A.map(lg(25), .9, 2.5).x, A.map(0, .9).y, A.map(0, 0, 2.5).z]], 'muted');
  Viz.label(v, 'thin ribbons: Al₂O₃ at the hole bottom after 50 cycles, / 5.15 nm', [A.map(lg(.3), 0, -.5).x, A.map(0, 1.05).y + .35, A.map(0, 0, -.5).z], { px: 15, scale: .07, color: 'muted' });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); r.rb.draw(p); r.th.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.c) return null; const i = Viz.ribbonIndex(hit, A); const c = u.c; return `<b>AR ${c.AR}, pulse ${f(c.t[i], 2)} s</b>step coverage ${f(c.sc[i] * 100, 1)} % · bottom ${f(c.thick[i], 2)} nm, top ${f(c.top[i], 2)} nm after 50 cycles`; };
  v.rot = -.35; v.tilt = .12;
}, { camera: [0, 2.1, 7.6], target: [0, -.1, 0], fov: 32 });

/* ---------- cycle: stacked segments ---------- */
N && Viz.view('[data-viz="cycle"]', v => {
  const rows = N.cycle.filter(q => q.nm === 5 && /kdes 0\.1/.test(q.basis)); const names = ['blanket\nbaseline', 'AR 100\nbaseline', 'AR 100\nsame-supply optimum'];
  const A = Viz.axes(v, { size: [4.2, 2.4, 1.4], x: { lo: -.6, hi: rows.length - .4, ticks: [] }, y: { lo: 0, hi: 42, ticks: [0, 10, 20, 30, 40], fmt: t => `${t} s`, label: 'time per ALD cycle' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const segCol = ['blue', 'grey', 'warm', 'steel'], segName = ['TMA pulse', 'purge (chamber to 1 %)', 'H₂O pulse', 'purge (k_des 0.1 /s)']; const boxes = [];
  rows.forEach((q, i) => { let y0 = 0; q.seg.forEach((h, k) => { const g = new THREE.BoxGeometry(.6, 1, .6); g.translate(0, .5, 0); const m = new THREE.Mesh(g, Viz.mat(segCol[k], { emissive: segCol[k], emissiveIntensity: .15 })); const b = A.map(i, y0, 0), t = A.map(i, y0 + h, 0); m.position.copy(b); m.userData = { q, k, h, y0, hgt: t.y - b.y }; m.scale.y = 1e-4; v.rig.add(m); v.pickables.push(m); boxes.push(m); y0 += h; });
    Viz.label(v, `${f(q.cycle, 1)} s per cycle\n${f(q.wph60, 1)} WPH at 5 nm`, [A.map(i, q.cycle, 0).x, A.map(0, q.cycle).y + .45, 0], { px: 18, scale: .08, color: 'ink', weight: 700 }); Viz.label(v, names[i], [A.map(i, 0, 0).x, A.map(0, 0).y - .42, .4], { px: 16, scale: .072, bg: false }); });
  segName.forEach((n, k) => Viz.label(v, n, [A.map(-.6, 0, 0).x - .3, A.map(0, 40 - k * 4).y, .6], { px: 15, scale: .068, color: segCol[k] }));
  v.update = vv => { boxes.forEach((m, i) => { const p = clamp(vv.p * 1.4 - (i % 4) * .12, 0, 1); m.scale.y = Math.max(1e-4, m.userData.hgt * p); }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.q) return null; return `<b>${u.q.recipe}: ${segName[u.k]}</b>${f(u.h, 2)} s of a ${f(u.q.cycle, 1)} s cycle · ${u.q.cycles} cycles for ${u.q.nm} nm · ${f(u.q.wph, 2)} WPH (${f(u.q.wph60, 2)} with 60 s handling)`; };
  v.rot = -.25; v.tilt = .12;
}, { camera: [0, 2.0, 7.6], target: [0, -.1, 0], fov: 32 });

/* ---------- wph: the rule as throughput ---------- */
N && Viz.view('[data-viz="wph"]', v => {
  const W = N.wph.rows; const xs = [...new Set(W.map(q => q.x))].sort((a, b) => a - b);
  const A = Viz.axes(v, { size: [4.4, 2.2, 1.8], x: { lo: 0, hi: 2.25, ticks: [0, 1, 2], fmt: t => `${Math.round(Math.pow(10, t))}×`, label: 'wall area / wafer area (log)' }, y: { lo: lg(.7), hi: lg(9), ticks: [0, lg(2), lg(5)], fmt: t => `${+Math.pow(10, t).toFixed(0)}`, label: 'wafers per hour, 5 nm Al₂O₃ (log)' }, z: { lo: -.5, hi: 1.5, ticks: [] } });
  const cols = { 0.01: 'grey', 0.03: 'warm' }; const rbs = xs.map((x, i) => { const rows = W.filter(q => q.x === x).sort((a, b) => a.wall - b.wall); const rb = Viz.ribbon(v, A, rows.map(q => lg(q.wall)), rows.map(q => lg(q.wph60)), { color: cols[x] || 'cyan', width: .22, z: i, name: `${x * 100} % TMA, 60 s handling` }); rb.mesh.userData = { ...rb.mesh.userData, rows, name: `${x * 100} % TMA and H₂O, 60 s handling per wafer` }; const r0 = Viz.ribbon(v, A, rows.map(q => lg(q.wall)), rows.map(q => lg(q.wph0)), { color: cols[x] || 'cyan', width: .07, z: i, skirt: false }); r0.mesh.userData = { ...r0.mesh.userData, rows, name: `${x * 100} %, no handling overhead`, no: true }; const lab = Viz.label(v, `${x * 100} % TMA`, [A.map(2.3, 0, i).x + .4, A.map(0, lg(rows[rows.length - 1].wph60)).y, A.map(0, 0, i).z], { px: 18, scale: .08, color: cols[x] || 'cyan', weight: 700 }); return { rb, r0, lab }; });
  const run = N.cycle.find(q => /AR-100 baseline|AR 100 baseline|patterned.*baseline/i.test(q.recipe) && q.nm === 5 && /kdes 0\.1/.test(q.basis)); if (run) { const m = sphere(v, A, lg(55.9), lg(run.wph60), 0, 'green', .1, { html: `<b>the chamber run: AR-100 wafer, 1 % TMA</b>cycle ${f(run.cycle, 1)} s (14.05 + 13.4 s pulses), ${f(run.wph60, 2)} wafers per hour with 60 s handling` }); Viz.label(v, 'chamber run', [m.position.x + .5, m.position.y + .3, .3], { px: 16, scale: .07, color: 'green', weight: 700 }); }
  [[11.2, 'AR 20'], [28, 'AR 50'], [55.9, 'AR 100']].forEach(([w, t]) => { Viz.line(v, [[A.map(lg(w), 0, -.5).x, A.map(0, lg(.7)).y, A.map(0, 0, -.5).z], [A.map(lg(w), 0, -.5).x, A.map(0, lg(9)).y, A.map(0, 0, -.5).z]], 'muted'); Viz.label(v, t, [A.map(lg(w), 0, -.5).x, A.map(0, lg(9)).y + .15, A.map(0, 0, -.5).z], { px: 15, scale: .068, bg: false }); });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); r.rb.draw(p); r.r0.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.html) return u.html; if (!u.rows) return null; const i = Viz.ribbonIndex(hit, A); const q = u.rows[i]; return `<b>${u.name}</b>wall ratio ${f(q.wall, 1)}×: cycle ${f(q.cycle, 1)} s, ${f(u.no ? q.wph0 : q.wph60, 2)} wafers per hour`; };
  v.rot = -.3; v.tilt = .1;
}, { camera: [0, 1.9, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- exposure: this model among published saturation exposures ---------- */
N && Viz.view('[data-viz="exposure"]', v => {
  const E = N.exposure; const A = Viz.axes(v, { size: [4.6, 2.4, 1.6], x: { lo: 0, hi: 360, ticks: [0, 100, 200, 300], fmt: t => `${t} °C`, label: 'temperature' }, y: { lo: -1.5, hi: 3, ticks: [-1, 0, 1, 2, 3], fmt: t => `${Math.pow(10, t) >= 1 ? Math.round(Math.pow(10, t)) : Math.pow(10, t).toFixed(1)} Pa·s`, label: 'saturating TMA exposure (log)' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const kindOf = q => /this work/.test(q.source) ? 'this' : /applied|feed|reactor p/.test(q.kind) ? 'applied' : /transport/.test(q.kind) ? 'excluded' : /CFD|impingement/.test(q.kind) ? 'cfd' : 'kinetic';
  const colOf = { this: 'green', applied: 'cool', kinetic: 'warm', cfd: 'cyan', excluded: 'grey' }; const pts = [];
  E.forEach(q => { if (!(q.Pa_s > 0)) return; const k = kindOf(q); pts.push(sphere(v, A, q.T, lg(q.Pa_s), 0, colOf[k], k === 'this' ? .12 : .075, { html: `<b>${q.source}</b>${q.T} °C · ${q.kind} · ${f(q.Pa_s, 3)} Pa·s (${f(q.L, 0)} L)${q.s ? ` · s ${q.s}` : ''}<br><span style="color:var(--muted)">${q.note}</span>` })); });
  // Langmuir 90 % exposure lines for s = 1e-3, 5e-3, 1e-2: E90 = ln(10) N0 sqrt(2 pi m kB T) / s
  const m = 72.09 * 1.66054e-27, kB = 1.380649e-23, N0 = 4e18;
  [[1e-3, 's = 1e-3'], [5e-3, 's = 5e-3'], [1e-2, 's = 1e-2']].forEach(([s, name], i) => { const T = [], Y = []; for (let tc = 0; tc <= 360; tc += 20) { T.push(tc); Y.push(lg(Math.log(10) * N0 * Math.sqrt(2 * Math.PI * m * kB * (tc + 273.15)) / s)); } const rb = Viz.ribbon(v, A, T, Y, { color: 'muted', width: .05, z: 0, skirt: false }); rb.draw(1); rb.mesh.userData = { ...rb.mesh.userData, name: `Langmuir 90 % exposure at N₀ = 4e18 /m², ${name}` }; Viz.label(v, `Langmuir E90, ${name}`, [A.map(370, 0, 0).x + .3, A.map(0, Y[Y.length - 1]).y, 0], { px: 14, scale: .065, bg: false, color: 'muted' }); });
  const tags = [Viz.label(v, 'applied doses (upper bounds)', [A.map(150, 0, 0).x, A.map(0, 2.3).y, .3], { px: 16, scale: .07, color: 'cool' }), Viz.label(v, 'kinetic values implied by a reported s', [A.map(120, 0, 0).x, A.map(0, -.9).y, .3], { px: 16, scale: .07, color: 'warm' }), Viz.label(v, 'this model: feed 1.03 Pa·s · surface 0.146 Pa·s', [A.map(300, 0, 0).x - 1.2, A.map(0, .4).y, .5], { px: 16, scale: .07, color: 'green', weight: 700 })];
  v.update = vv => { pts.forEach((p, i) => p.scale.setScalar(.001 + clamp(vv.p * 1.4 - i / pts.length * .4, 0, 1))); tags.forEach(t => t.material.opacity = clamp(vv.p * 2 - 1, 0, 1)); return false; };
  v.hover = hit => hit.object.userData.html || (hit.object.userData.name ? `<b>${hit.object.userData.name}</b>` : null); v.rot = -.15; v.tilt = .08;
}, { camera: [0, 1.7, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- gakis: coverage across the wafer against the measured profile ---------- */
N && Viz.view('[data-viz="gakis"]', v => {
  const G = N.gakis; const A = Viz.axes(v, { size: [4.4, 2.2, 2.2], x: { lo: 0, hi: 100, ticks: [0, 25, 50, 75, 100], fmt: t => `${t} mm`, label: 'wafer radius' }, y: { lo: .965, hi: 1.004, ticks: [.97, .98, .99, 1.0], fmt: t => t.toFixed(2), label: 'coverage after the cycle · GPC / GPC(60 ms)' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const cols = ['rose', 'warm', 'green', 'cyan']; const rbs = G.cases.map((c, i) => { const rb = Viz.ribbon(v, A, c.r_mm, c.theta_end.map(x => Math.max(.965, x)), { color: cols[i], width: .2, z: i, name: c.name, skirt: false }); rb.mesh.userData = { ...rb.mesh.userData, c }; const lab = Viz.label(v, `${c.name}: ${f(c.theta_end[0], 4)}`, [A.map(12, 0, i).x, A.map(0, 1.004).y + .18 + i * .17, A.map(0, 0, i).z], { px: 16, scale: .072, color: cols[i], weight: 700 }); return { rb, lab }; });
  const meas = G.meas_pos_mm.map((p, i) => sphere(v, A, Math.abs(p), G.meas_gpc[i], p < 0 ? -.25 : .25, p < 0 ? 'ink' : 'steel', .07, { html: `<b>measured, 25 ms, ${p < 0 ? 'side away from the loading door' : 'loading-door side'}</b>r = ${f(Math.abs(p), 1)} mm: GPC / GPC(60 ms) = ${f(G.meas_gpc[i], 4)} ± ${G.meas_err}` }));
  G.meas_pos_mm.forEach((p, i) => { const lo = A.map(Math.abs(p), G.meas_gpc[i] - G.meas_err, p < 0 ? -.25 : .25), hi = A.map(Math.abs(p), G.meas_gpc[i] + G.meas_err, p < 0 ? -.25 : .25); Viz.line(v, [[lo.x, lo.y, lo.z], [hi.x, hi.y, hi.z]], 'muted'); });
  const tag = Viz.label(v, 'measured 25 ms (Gakis 2019, Fig. 10): 2.4 % non-uniform\nmodel: 0.0007 %; their 3D CFD exposure 0.55–0.60 Pa·s, this 0.67', [A.map(50, 0, -.5).x, A.map(0, .972).y, A.map(0, 0, -.5).z + .3], { px: 15, scale: .07, color: 'ink' });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .08, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); fadeIn(meas, vv.p, .3); tag.material.opacity = clamp(vv.p * 2 - 1, 0, 1); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.html) return u.html; if (!u.c) return null; const i = Viz.ribbonIndex(hit, A); return `<b>${u.c.name}</b>r = ${f(u.c.r_mm[i], 0)} mm: coverage ${f(u.c.theta_end[i], 5)} · peak surface TMA ${f(Math.max(...u.c.p_centre_Pa), 1)} Pa at the centre`; };
  v.rot = -.3; v.tilt = .12;
}, { camera: [0, 2.0, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- wafer3d: the 3D chamber's saturation-time map ---------- */
N && Viz.view('[data-viz="wafer3d"]', v => {
  const C = N.chamber3d, n = C.n, SC = 1 / 60, HS = .9; const pos = new Float32Array(n * 3), col = new Float32Array(n * 3);
  let lo = Infinity, hi = -Infinity; for (let i = 0; i < n; i++) { lo = Math.min(lo, C.t90[i]); hi = Math.max(hi, C.t90[i]); }
  const dlo = Math.min(...C.dev), dhi = Math.max(...C.dev);
  for (let i = 0; i < n; i++) { pos[3 * i] = C.y[i] * SC; pos[3 * i + 2] = -C.z[i] * SC; pos[3 * i + 1] = (C.dev[i] - dlo) / ((dhi - dlo) || 1) * HS; const c = viridis((C.dev[i] - dlo) / ((dhi - dlo) || 1)); col[3 * i] = c[0]; col[3 * i + 1] = c[1]; col[3 * i + 2] = c[2]; }
  const geo = new THREE.BufferGeometry(); geo.setAttribute('position', new THREE.BufferAttribute(pos, 3)); geo.setAttribute('color', new THREE.BufferAttribute(col, 3)); geo.setIndex(C.tri); geo.computeVertexNormals();
  const mesh = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: .5, metalness: .1 })); mesh.userData.kind = 'map'; v.rig.add(mesh); v.pickables.push(mesh);
  const base = new THREE.Mesh(new THREE.CylinderGeometry(2.5, 2.5, .06, 96), Viz.mat('steel', { metalness: .5, roughness: .4 })); base.position.y = -.05; v.rig.add(base);
  const ring = new THREE.Mesh(new THREE.TorusGeometry(200 * SC, .025, 8, 96), Viz.mat('grey')); ring.rotation.x = Math.PI / 2; ring.position.y = .02; v.rig.add(ring);
  const port = new THREE.Mesh(new THREE.TorusGeometry(200 * SC, .08, 8, 32, Math.PI / 3), Viz.mat('rose', { emissive: 'rose', emissiveIntensity: .6 })); port.rotation.x = Math.PI / 2; port.rotation.z = -Math.PI / 6; port.position.y = .04; v.rig.add(port);
  Viz.label(v, 'pump port, 60°', [3.6, .3, 0], { px: 16, scale: .075, color: 'rose', weight: 700 });
  Viz.label(v, `hole-bottom t90 over the wafer: ${f(lo, 2)} → ${f(hi, 2)} s · colour and height = deviation from the wedge, ${f(dlo, 2)} → ${f(dhi, 2)} %`, [0, HS + 1.1, -2.6], { px: 16, scale: .075, color: 'ink' });
  Viz.label(v, `mean ${f(C.mean, 2)} s · worst ${f(C.worst, 2)} s on the far side · wedge worst ${C.wedge_worst} s · ${f(C.cells, 0)} cells`, [0, -.5, 2.9], { px: 15, scale: .07, color: 'muted' });
  const base0 = pos.slice();
  v.update = vv => { const e = vv.p; for (let i = 0; i < n; i++) pos[3 * i + 1] = base0[3 * i + 1] * e; geo.attributes.position.needsUpdate = true; return false; };
  v.hover = hit => { if (hit.object.userData.kind !== 'map' || !hit.face) return null; const i = hit.face.a; return `<b>r = ${f(Math.hypot(C.y[i], C.z[i]), 0)} mm, azimuth ${f(Math.atan2(C.z[i], C.y[i]) * 180 / Math.PI, 0)}° from the port</b>hole-bottom t90 ${f(C.t90[i], 3)} s · ${f(C.dev[i], 2)} % from the wedge`; };
  v.spin = .12; v.tilt = .35;
}, { camera: [0, 4.2, 7.4], target: [0, .4, 0], fov: 34 });

/* ---------- azimuth: t90 against azimuth at three radii ---------- */
N && Viz.view('[data-viz="azimuth"]', v => {
  const Z = N.chamber3d.azimuthal; const keys = Object.keys(Z); const cols = ['cyan', 'green', 'warm'];
  const A = Viz.axes(v, { size: [4.4, 2.0, 2.0], x: { lo: -180, hi: 180, ticks: [-180, -90, 0, 90, 180], fmt: t => `${t}°`, label: 'azimuth from the pump port' }, y: { lo: 13.95, hi: 14.2, ticks: [14, 14.05, 14.1, 14.15, 14.2], fmt: t => `${t.toFixed(2)} s`, label: 'hole-bottom t90' }, z: { lo: -.5, hi: 2.5, ticks: [] } });
  const rbs = keys.map((k, i) => { const q = Z[k]; const rb = Viz.ribbon(v, A, q.az, q.t90, { color: cols[i], width: .2, z: i, name: `r = ${k} mm` }); rb.mesh.userData = { ...rb.mesh.userData, q, k }; const mn = Math.min(...q.t90), mx = Math.max(...q.t90); const lab = Viz.label(v, `r = ${k} mm · spread ${f((mx - mn) / mn * 100, 2)} %`, [A.map(190, 0, i).x + .6, A.map(0, q.t90[q.t90.length - 1]).y + .1, A.map(0, 0, i).z], { px: 17, scale: .075, color: cols[i], weight: 700 }); return { rb, lab }; });
  Viz.line(v, [[A.map(-180, 14.05, -.5).x, A.map(0, 14.05).y, A.map(0, 0, -.5).z], [A.map(180, 14.05, 2.5).x, A.map(0, 14.05).y, A.map(0, 0, 2.5).z]], 'muted'); Viz.label(v, 'wedge, centre: 14.05 s', [A.map(-150, 14.05, -.5).x, A.map(0, 14.05).y + .12, A.map(0, 0, -.5).z], { px: 14, scale: .065, bg: false });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.q) return null; const i = Viz.ribbonIndex(hit, A); return `<b>r = ${u.k} mm, azimuth ${f(u.q.az[i], 0)}°</b>hole-bottom t90 ${f(u.q.t90[i], 3)} s`; };
  v.rot = -.3; v.tilt = .1;
}, { camera: [0, 1.9, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- taper: four cut-open holes ---------- */
N && Viz.view('[data-viz="taper"]', v => {
  const T = N.taper, HH = 2.6, NA = 24, ND = T.depths.length; const cBare = new THREE.Color(), cCov = new THREE.Color(), c = new THREE.Color(); const meshes = [];
  T.holes.forEach((h, k) => { const x0 = -2.4 + k * 1.6, R0 = .42; const geo = new THREE.BufferGeometry(); const pos = new Float32Array(ND * (NA + 1) * 3), col = new Float32Array(ND * (NA + 1) * 3), idx = [];
    for (let i = 0; i < ND; i++) { const fr = i / (ND - 1), R = R0 * (1 - (1 - h.ratio) * fr); for (let j = 0; j <= NA; j++) { const q = i * (NA + 1) + j, th = Math.PI * .5 + j / NA * Math.PI; pos[3 * q] = x0 + Math.cos(th) * R; pos[3 * q + 1] = HH / 2 - fr * HH; pos[3 * q + 2] = Math.sin(th) * R; if (i < ND - 1 && j < NA) { const a = q, b = q + 1, cc = q + NA + 1, d = q + NA + 2; idx.push(a, cc, b, b, cc, d); } } }
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3)); geo.setAttribute('color', new THREE.BufferAttribute(col, 3)); geo.setIndex(idx); geo.computeVertexNormals();
    const m = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: .5, metalness: .2 })); m.userData = { h, col, geo, ND, NA }; v.rig.add(m); v.pickables.push(m); meshes.push(m);
    const shell = new THREE.Mesh(new THREE.CylinderGeometry(R0 + .15, R0 * h.ratio + .15, HH, 24, 1, false, Math.PI * .5, Math.PI), Viz.mat('grey', { side: THREE.DoubleSide, opacity: .45 })); shell.position.x = x0; v.rig.add(shell);
    Viz.label(v, `bottom ${Math.round(h.ratio * 100)} %\nt90 ${f(h.t90_bottom, 2)} s`, [x0, HH / 2 + .45, 0], { px: 18, scale: .08, color: 'ink', weight: 700 }); });
  const base = new THREE.Mesh(new THREE.BoxGeometry(7.2, .12, 1.6), Viz.mat('grey', { roughness: .8 })); base.position.y = -HH / 2 - .06; v.rig.add(base);
  const tag = Viz.label(v, 'coverage at 10 s of the 1 % pulse, chamber-driven · straight hole 14.05 s', [0, -HH / 2 - .5, 1.0], { px: 15, scale: .07, color: 'muted' });
  const paint = (m, p) => { const { h, col, geo, ND, NA } = m.userData; cBare.set(v.dark ? 0x5b6676 : 0x8a95a5); cCov.set(v.pal.green); const snap = h.snap['10']; for (let i = 0; i < ND; i++) { c.copy(cBare).lerp(cCov, snap[i] * p); for (let j = 0; j <= NA; j++) { const q = i * (NA + 1) + j; col[3 * q] = c.r; col[3 * q + 1] = c.g; col[3 * q + 2] = c.b; } } geo.attributes.color.needsUpdate = true; };
  v.update = vv => { meshes.forEach(m => paint(m, vv.p)); return false; }; v.onTheme = vv => meshes.forEach(m => paint(m, vv.p));
  v.hover = hit => { const u = hit.object.userData; if (!u.h) return null; const i = clamp(Math.round((HH / 2 - hit.point.y) / HH * (u.ND - 1)), 0, u.ND - 1); return `<b>bottom at ${Math.round(u.h.ratio * 100)} % of the mouth, ${Math.round(T.depths[i] * 100)} % depth</b>coverage ${f(u.h.snap['5'][i] * 100, 0)} % at 5 s · ${f(u.h.snap['10'][i] * 100, 0)} % at 10 s · ${f(u.h.snap['13'][i] * 100, 0)} % at 13 s · bottom t90 ${f(u.h.t90_bottom, 2)} s`; };
  v.rot = .2; v.tilt = .1;
}, { camera: [0, .8, 8.0], target: [0, 0, 0], fov: 34 });

/* ---------- recess: the word-line stack ---------- */
N && Viz.view('[data-viz="recess"]', v => {
  const R = N.recess, snapKeys = Object.keys(R.snaps), key = '9.0'; const S = R.snaps[key] || R.snaps[snapKeys[snapKeys.length - 1]]; const nx = R.x_um.length, ny = R.y_nm.length;
  const HH = 3.0, GAP = .5, RL = 1.6, cBare = new THREE.Color(), cCov = new THREE.Color(), c = new THREE.Color();
  // the two slit walls as vertical strips coloured by theta_slit; the recesses as instanced slabs on both sides, coloured by theta_recess
  const walls = []; [-1, 1].forEach(side => { const g = new THREE.PlaneGeometry(.001, HH, 1, nx - 1); const pos = g.attributes.position; for (let i = 0; i < pos.count; i++) pos.setX(i, side * GAP / 2); const col = new Float32Array(pos.count * 3); g.setAttribute('color', new THREE.BufferAttribute(col, 3)); const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: .6 })); m.userData = { side, col, g }; v.rig.add(m); walls.push(m); });
  const cnt = nx * ny * 2; const slab = new THREE.BoxGeometry(RL / ny, HH / nx * .55, .9); const im = new THREE.InstancedMesh(slab, new THREE.MeshStandardMaterial({ roughness: .55 }), cnt); v.rig.add(im); v.pickables.push(im);
  const M = new THREE.Matrix4(), cc = new THREE.Color(); let k = 0; const meta = [];
  for (let i = 0; i < nx; i++) for (let side of [-1, 1]) for (let j = 0; j < ny; j++) { M.makeTranslation(side * (GAP / 2 + (j + .5) * RL / ny), HH / 2 - (i + .5) / nx * HH, 0); im.setMatrixAt(k, M); meta.push({ i, j, side }); k++; }
  im.instanceMatrix.needsUpdate = true;
  const stackL = new THREE.Mesh(new THREE.BoxGeometry(RL + .4, HH, 1.0), Viz.mat('steel', { opacity: .18, side: THREE.DoubleSide })); stackL.position.set(-(GAP / 2 + RL / 2 + .2), 0, 0); v.rig.add(stackL); const stackR = stackL.clone(); stackR.position.x = GAP / 2 + RL / 2 + .2; v.rig.add(stackR);
  Viz.label(v, `slit ${R.geom.w_nm} nm wide, ${R.geom.L_um} µm deep · recesses ${R.geom.g_nm} nm × ${R.geom.L_r_nm} nm at ${R.geom.pitch_nm} nm pitch · coverage at ${key} s`, [0, HH / 2 + .5, 0], { px: 16, scale: .075, color: 'ink' });
  Viz.label(v, 'slit', [0, -HH / 2 - .35, .6], { px: 15, scale: .07, bg: false }); Viz.label(v, 'word-line recesses →', [GAP / 2 + RL / 2, -HH / 2 - .35, .6], { px: 15, scale: .07, bg: false }); Viz.label(v, '← word-line recesses', [-(GAP / 2 + RL / 2), -HH / 2 - .35, .6], { px: 15, scale: .07, bg: false });
  const paint = p => { cBare.set(v.dark ? 0x5b6676 : 0x8a95a5); cCov.set(v.pal.green); walls.forEach(w => { const { col, g } = w.userData; const pos = g.attributes.position; for (let q = 0; q < pos.count; q++) { const row = Math.floor(q / 2); c.copy(cBare).lerp(cCov, S.slit[Math.min(nx - 1, row)] * p); col[3 * q] = c.r; col[3 * q + 1] = c.g; col[3 * q + 2] = c.b; } g.attributes.color.needsUpdate = true; }); meta.forEach((m, q) => { cc.copy(cBare).lerp(cCov, S.recess[m.i][m.j] * p); im.setColorAt(q, cc); }); im.instanceColor.needsUpdate = true; };
  paint(0);
  v.update = vv => { paint(vv.p); return false; }; v.onTheme = vv => paint(vv.p);
  v.hover = hit => { if (hit.object !== im) return null; const m = meta[hit.instanceId]; if (!m) return null; const parts = snapKeys.map(kk => `${f(R.snaps[kk].recess[m.i][m.j] * 100, 0)} % at ${kk} s`).join(' · '); return `<b>recess at ${f(R.x_um[m.i], 2)} µm down the slit, ${f(R.y_nm[m.j], 0)} nm into the recess</b>${parts}`; };
  v.rot = -.35; v.tilt = .1;
}, { camera: [0, .8, 8.6], target: [0, 0, 0], fov: 34 });

/* ---------- chemistry: the rule per chemistry ---------- */
N && Viz.view('[data-viz="chemistry"]', v => {
  const rows = N.chemistry.rows; const td = N.chemistry.tdmah_run || {};
  const A = Viz.axes(v, { size: [4.4, 2.2, 1.2], x: { lo: -.5, hi: rows.length - .5, ticks: [] }, y: { lo: 0, hi: 62, ticks: [0, 15, 30, 45, 60], fmt: t => `${t} s`, label: 'rule: saturating pulse, AR-100 wafer at 1 %' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const colOf = { Al2O3: 'green', HfO2: 'cyan', TiN: 'warm', W: 'rose', ZrO2: 'violet' };
  const bars = rows.map((q, i) => ({ x: i, z: 0, y: q.rule_AR100, color: colOf[q.film] || 'grey', q }));
  const B = Viz.bars(v, A, bars, { width: .5, depth: .5 });
  bars.forEach((b, i) => { const p = A.map(i, b.y, 0); b.lab = Viz.label(v, `${f(b.y, 1)} s`, [p.x, p.y + .18, p.z], { px: 18, scale: .08 }); Viz.label(v, `${b.q.film}\n${b.q.precursor}`, [p.x, A.map(0, 0).y - .36, .45], { px: 16, scale: .072, bg: false }); });
  const hf = rows.findIndex(q => q.film === 'HfO2'); const runT = td.TDMAH && (td.TDMAH.t90_bottom_centre || td.TDMAH.bottom_t90_centre || td.TDMAH.t90_bottom); let runS = null;
  if (hf >= 0 && runT) { runS = sphere(v, A, hf, runT, 0, 'ink', .1, { html: `<b>HfO₂ from TDMAH: the chamber run</b>hole bottoms at 90 % in ${f(runT, 2)} s against the rule's ${f(rows[hf].rule_AR100, 2)} s (${f(runT / rows[hf].rule_AR100, 2)}×); the excess is mostly the 2.3× lower gas diffusivity` }); Viz.label(v, `chamber run ${f(runT, 1)} s`, [runS.position.x + .5, runS.position.y + .25, .3], { px: 16, scale: .072, color: 'ink', weight: 700 }); }
  v.update = vv => { B.draw(vv.p); bars.forEach(b => b.lab.material.opacity = clamp(vv.p * 2 - 1, 0, 1)); return false; };
  v.hover = hit => { if (hit.object.userData.html) return hit.object.userData.html; const b = bars[hit.instanceId]; return b ? `<b>${b.q.film} from ${b.q.precursor}</b>s ${b.q.s} · N₀ ${Number(b.q.N0).toExponential(2)} /m² · GPC ${b.q.gpc} nm · v̄ ${f(b.q.vbar, 0)} m/s · rule ${f(b.q.rule_AR100, 2)} s (blanket ${f(b.q.rule_blanket, 3)} s)<br><span style="color:var(--muted)">${b.q.note}</span>` : null; };
  v.rot = -.2; v.tilt = .1;
}, { camera: [0, 1.8, 7.6], target: [0, -.2, 0], fov: 32 });
