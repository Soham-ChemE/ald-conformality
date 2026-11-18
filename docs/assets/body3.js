/* Re-analysis scenes (section 03b), drawn by the shared Viz renderer from window.ALD.lhar (assets/data-lhar.js, built by
   docs/build_lhar_assets.py from the committed result files):
   lharslope  the slope-relation constant A = s / slope^2 against front position: ideal-slit ballistic curves per sticking
              coefficient, the continuum constants and the corrected continuum model;
   lhars      the TMA sticking coefficient each data set gives under each transport description, with the optical band;
   lharzen    Zenodo 3974438: measured profile and the three fitted single-exposure models;
   lharyli    Ylilammi 2018 Fig. 6: the same. */
import { Viz } from './viz.js';
import { THREE } from './three-kit.js';
const L = window.ALD && window.ALD.lhar; const lg = Math.log10; const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const f = (v, d = 0) => (+v).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d });
const sci = v => { const e = Math.floor(lg(v)); return `${(v / Math.pow(10, e)).toFixed(1)}e${e}`; };
const fadeIn = (objs, p, from = .6) => objs.forEach(o => { o.material.transparent = true; o.material.opacity = clamp((p - from) / (1 - from), 0, 1); });
if (!L) console.warn('re-analysis data missing');

/* ---------- lharslope: A against front position ---------- */
L && Viz.view('[data-viz="lharslope"]', v => {
  const S = L.slope, C = S.curves, n = C.length;
  const A = Viz.axes(v, { size: [4.4, 2.4, 2.6], x: { lo: lg(60), hi: lg(1200), ticks: [lg(100), lg(200), lg(500), lg(1000)], fmt: t => `${Math.round(Math.pow(10, t))}`, label: 'front position x₅₀ / H (log)' }, y: { lo: 0, hi: 85, ticks: [0, 20, 40, 60, 80], fmt: t => `${t}`, label: 'A = s / slope²' }, z: { lo: -.6, hi: n - .4, ticks: [] } });
  const cols = ['violet', 'violet', 'blue', 'blue', 'cyan', 'cyan', 'green', 'warm', 'hot'];
  const rbs = C.map((c, i) => { const rb = Viz.ribbon(v, A, c.z50.map(lg), c.A, { color: cols[i % cols.length], width: .2, z: i, name: `ballistic, s = ${c.s}`, skirtOpacity: .06 }); rb.mesh.userData = { ...rb.mesh.userData, c, kind: 'ball' }; const last = c.z50.length - 1; const lab = Viz.label(v, `s = ${c.s}`, [A.map(lg(c.z50[last]), 0, i).x + .34, A.map(0, c.A[last]).y, A.map(0, 0, i).z], { px: 15, scale: .065, color: cols[i % cols.length], bg: false }); return { rb, lab }; });
  // continuum constants as flat sheets across every s, and the corrected continuum model as a ribbon on the front layer
  const sheet = (y, key, op) => { const m = new THREE.Mesh(new THREE.PlaneGeometry(4.4, 2.6), Viz.mat(key, { opacity: op, side: THREE.DoubleSide, depthWrite: false })); m.rotation.x = -Math.PI / 2; m.position.set(0, A.map(0, y).y, 0); m.userData = { kind: 'cont', y }; v.rig.add(m); v.pickables.push(m); return m; };
  const s1 = sheet(S.arts, 'rose', .13), s2 = sheet(S.gonsalves, 'grey', .1);
  const zz = []; for (let q = lg(60); q <= lg(1200) + 1e-9; q += .05) zz.push(q);
  const corr = Viz.ribbon(v, A, zz, zz.map(q => S.arts * .375 * Math.log(Math.pow(10, q))), { color: 'rose', width: .1, z: -.3, name: 'corrected continuum', skirt: false }); corr.mesh.userData = { ...corr.mesh.userData, kind: 'corr' };
  const hyd = S.hydraulic.z50.map((z, i) => { const m = new THREE.Mesh(new THREE.SphereGeometry(.04, 12, 10), Viz.mat('ink', { emissive: 'ink', emissiveIntensity: .4 })); m.position.copy(A.map(lg(z), S.hydraulic.A[i], -.3)); m.userData = { kind: 'hyd', z, A: S.hydraulic.A[i] }; v.rig.add(m); v.pickables.push(m); return m; });
  const labs = [Viz.label(v, `continuum, hydraulic diffusivity: ${S.arts} (Arts 2019) and ${S.gonsalves} (Gonsalves 2024)`, [A.map(lg(260), 0, 0).x, A.map(0, S.arts).y - .2, A.map(0, 0, n - .4).z + .3], { px: 17, scale: .075, color: 'rose', weight: 700 }),
    Viz.label(v, `ideal-slit ballistic (Machball): ${f(S.A_range[0], 0)} to ${f(S.A_range[1], 0)}`, [A.map(lg(180), 0, 0).x, A.map(0, 82).y, A.map(0, 0, n / 2).z], { px: 18, scale: .08, color: 'cyan', weight: 700 }),
    Viz.label(v, '13.9 × (3/8) ln(x₅₀/H)', [A.map(lg(700), 0, -.3).x, A.map(0, 30).y - .22, A.map(0, 0, -.3).z], { px: 15, scale: .068, color: 'rose' })];
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.5 - i * .06, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .97 ? 1 : 0; }); corr.draw(clamp(vv.p * 1.4, 0, 1)); fadeIn(hyd, vv.p, .2); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.kind === 'ball') { const i = Viz.ribbonIndex(hit, A); return `<b>ideal-slit ballistic, s = ${u.c.s}</b>front at ${f(u.c.z50[i], 0)} channel heights: A = ${f(u.c.A[i], 1)}, ${f(u.c.A[i] / S.arts, 1)}× the continuum constant`; } if (u.kind === 'cont') return `<b>continuum model, hydraulic diffusivity</b>A = ${u.y}, independent of front position and of s`; if (u.kind === 'corr') return `<b>corrected continuum, cut-off at the front</b>captures the trend, 1.35 to 1.6× below the ballistic curve at s = 8e-3`; if (u.kind === 'hyd') return `<b>this repository's 1D slit model</b>front at ${f(u.z, 0)}: A = ${f(u.A, 2)}`; return null; };
  v.rot = -.3; v.tilt = .05;
}, { camera: [0, 1.5, 8.4], target: [0, -.05, 0], fov: 32 });

/* ---------- lhars: sticking coefficient by data set and transport description ---------- */
L && Viz.view('[data-viz="lhars"]', v => {
  const M = L.methods, R = M.rows, nm = M.names.length, nr = R.length; const lo = lg(5e-4), hi = lg(4e-2);
  const A = Viz.axes(v, { size: [4.4, 2.5, 3.0], x: { lo: -.6, hi: nm - .4, ticks: [] }, y: { lo, hi, ticks: [lg(1e-3), lg(3e-3), lg(1e-2), lg(3e-2)], fmt: t => sci(Math.pow(10, t)), label: 'TMA sticking coefficient (log)' }, z: { lo: -.6, hi: nr - .4, ticks: [] } });
  const cols = { yim: 'blue', zen: 'rose', yli: 'violet', zenm: 'hot', ylim: 'warm', arts: 'green' };
  const items = []; R.forEach((row, j) => row.s.forEach((s, i) => { if (s != null) items.push({ x: i, y: lg(s), z: j, color: cols[row.key], s, row, method: M.names[i] }); }));
  const bars = Viz.bars(v, A, items, { width: .34, depth: .26 });
  // the optical value as a slab through the whole field, and the reactor fit as a line on the back wall
  const y0 = A.map(0, lg(M.sfg[0])).y, y1 = A.map(0, lg(M.sfg[2])).y; const slab = new THREE.Mesh(new THREE.BoxGeometry(4.4, y1 - y0, 3.0), Viz.mat('cyan', { opacity: .16, depthWrite: false })); slab.position.set(0, (y0 + y1) / 2, 0); slab.userData = { kind: 'sfg' }; v.rig.add(slab); v.pickables.push(slab);
  Viz.line(v, [[-2.2, A.map(0, lg(M.gakis)).y, -1.5], [2.2, A.map(0, lg(M.gakis)).y, -1.5]], 'muted');
  const mlab = M.names.map((t, i) => Viz.label(v, t.replace(', cut-off', ',\ncut-off').replace(' (Arts', '\n(Arts'), [A.map(i, 0, 0).x, A.map(0, lo).y - .32, A.map(0, 0, nr - .4).z + .35], { px: 16, scale: .076, color: 'ink', bg: false }));
  const rlab = R.map((row, j) => Viz.label(v, row.name, [A.map(nm - .4, 0, 0).x + 1.25, A.map(0, lo).y + .06, A.map(0, 0, j).z], { px: 16, scale: .074, color: cols[row.key], bg: false }));
  const slab1 = Viz.label(v, 'optical (SFG): 4 ± 1 × 10⁻³', [A.map(-.6, 0, 0).x - .1, (y0 + y1) / 2 + .2, A.map(0, 0, nr - .4).z], { px: 16, scale: .072, color: 'cyan', weight: 700 });
  v.update = vv => { bars.draw(vv.p); [...mlab, ...rlab, slab1].forEach(l => l.material.opacity = vv.p > .85 ? 1 : 0); return false; };
  v.hover = hit => { if (hit.object === bars.mesh && hit.instanceId != null) { const it = items[hit.instanceId]; const extra = it.row.samples ? ` · geometric mean of ${it.row.samples.length} samples` : ''; return `<b>${it.row.name}</b>${it.method}: s = ${sci(it.s)}${extra} · ${f(it.s / M.sfg[1], 1)}× the optical value`; } if (hit.object.userData.kind === 'sfg') return `<b>optical (sum-frequency generation), Vandalon and Kessels 2016</b>4 ± 1 × 10⁻³ on a planar surface at 300 °C; reactor-scale fit (Gakis 2019): ${sci(M.gakis)}`; return null; };
  v.rot = -.5; v.tilt = .2;
}, { camera: [0, 2.6, 8.8], target: [.35, -.15, 0], fov: 32 });

/* ---------- lharzen / lharyli: a measured profile and the three single-exposure fits ---------- */
const profile = (key, sel) => L && Viz.view(sel, v => {
  const P = L.fits[key], MD = P.models;
  const A = Viz.axes(v, { size: [4.4, 2.2, 2.4], x: { lo: 0, hi: 400, ticks: [0, 100, 200, 300, 400], label: 'depth / channel height' }, y: { lo: -.02, hi: 1.12, ticks: [0, .5, 1], fmt: t => `${t}`, label: 'thickness / plateau' }, z: { lo: -.5, hi: 4.5, ticks: [] } });
  const PM = L.pillar_mc && L.pillar_mc.cases.walls90_pillars49 && L.pillar_mc.cases.walls90_pillars49.fits[key];
  if (PM) MD.real = { y: PM.curve.y, s: PM.s, exposure: PM.exposure_Pa_s, rms: PM.rms, rms_front: PM.rms_front, rms_tail: PM.rms_tail };
  const defs = [['meas', 'measured', 'grey', P.measured.AR, P.measured.y, 0, .3], ['1d_hydraulic', 'continuum, hydraulic diffusivity', 'blue', P.AR, MD['1d_hydraulic'].y, 1, .2], ['1d_corr_1mm', 'continuum, corrected diffusivity', 'green', P.AR, MD['1d_corr_1mm'].y, 2, .2], ['ballistic', 'ideal-slit ballistic', 'rose', P.AR, MD.ballistic.y, 3, .2]].concat(PM ? [['real', 'real cavity (pillars at 49 µm, 90 µm side walls), Monte Carlo', 'warm', PM.curve.AR, PM.curve.y, 4, .2]] : []);
  const rbs = defs.map(([k, name, col, X, Y, z, w]) => { const rb = Viz.ribbon(v, A, X, Y, { color: col, width: w, z, name }); rb.mesh.userData = { ...rb.mesh.userData, k, name, X, Y }; return rb; });
  // the measured points repeated as small spheres on the ballistic layer, where the tails part company
  const dots = []; for (let i = 0; i < P.measured.AR.length; i += 2) { const m = new THREE.Mesh(new THREE.SphereGeometry(.028, 10, 8), Viz.mat('grey', { emissive: 'grey', emissiveIntensity: .5 })); m.position.copy(A.map(P.measured.AR[i], P.measured.y[i], PM ? 4 : 3)); v.rig.add(m); dots.push(m); }
  if (P.fit_to) { const x0 = A.map(P.fit_to, 0, 0).x, x1 = A.map(400, 0, 0).x; const box = new THREE.Mesh(new THREE.BoxGeometry(x1 - x0, 2.2, 2.4), Viz.mat('grey', { opacity: .1, depthWrite: false })); box.position.set((x0 + x1) / 2, 0, 0); v.rig.add(box); Viz.label(v, 'not fitted: the scans rise again here', [(x0 + x1) / 2, A.map(0, .62).y, A.map(0, 0, 0).z], { px: 15, scale: .066, color: 'muted', bg: false }); }
  const txt = k => `s = ${sci(MD[k].s)} · ${f(MD[k].exposure, 1)} Pa s · rms ${f(MD[k].rms, 3)}`;
  const labs = [Viz.label(v, 'measured', [A.map(60, 1, 0).x, A.map(0, 1.06).y + .22, A.map(0, 0, 0).z], { px: 18, scale: .078, color: 'grey', weight: 700 }), Viz.label(v, `hydraulic: ${txt('1d_hydraulic')}`, [A.map(330, .5, 1).x, A.map(0, .62).y, A.map(0, 0, 1).z], { px: 17, scale: .074, color: 'blue', weight: 700 }), Viz.label(v, `corrected: ${txt('1d_corr_1mm')}`, [A.map(330, .5, 2).x, A.map(0, .46).y, A.map(0, 0, 2).z], { px: 17, scale: .074, color: 'green', weight: 700 }), Viz.label(v, `ideal-slit ballistic: ${txt('ballistic')}`, [A.map(320, .5, 3).x, A.map(0, .3).y, A.map(0, 0, 3).z], { px: 17, scale: .074, color: 'rose', weight: 700 })].concat(PM ? [Viz.label(v, `real cavity, Monte Carlo: ${txt('real')}\nthe tail stays; the data have none`, [A.map(310, .5, 4).x, A.map(0, .14).y, A.map(0, 0, 4).z], { px: 17, scale: .074, color: 'warm', weight: 700 })] : []);
  v.update = vv => { rbs.forEach((r, i) => r.draw(clamp(vv.p * 1.3 - i * .08, 0, 1))); fadeIn(dots, vv.p, .5); labs.forEach(l => l.material.opacity = vv.p > .95 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.k) return null; const i = Viz.ribbonIndex(hit, A); if (u.k === 'meas') return `<b>${P.title}</b>depth ${f(u.X[i], 0)} heights: thickness / plateau ${f(u.Y[i], 3)}`; const m = MD[u.k]; return `<b>${u.name}</b>depth ${f(u.X[i], 0)}: ${f(u.Y[i], 3)} · fitted s = ${sci(m.s)}, residual ${f(m.rms, 3)} (front ${f(m.rms_front, 3)}, tail ${f(m.rms_tail, 3)})`; };
  v.rot = -.35; v.tilt = .12;
}, { camera: [0, 2.0, 7.6], target: [0, -.1, 0], fov: 32 });
profile('zenodo', '[data-viz="lharzen"]'); profile('ylilammi_al2o3', '[data-viz="lharyli"]');
