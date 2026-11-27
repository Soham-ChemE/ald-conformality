/* Scenes of section 03d (tested blind, and against an independent measurement), drawn by the shared Viz renderer from
   window.ALD_MECH.scenes (assets/data-mech.js, built by docs/build_mech_assets.py):
   edsfront   EDS line scans (left channel) of the 250, 500 and 1000-cycle samples against the model with and without channel narrowing
              (docs/blindtest/v2/posthoc/eds_cycle_series_curves.json);
   slitbench  diffusion along a slit: kinetic theory (Siewert and Valougeorgis 2004, Ne-Ar) against the Monte Carlo with
              hard-sphere and with relaxation collisions and against the hydraulic formula (figures/benchmark_slit_diffusion.csv). */
import { Viz } from './viz.js';
import { THREE } from './three-kit.js';
const M = window.ALD_MECH, SC = M && M.scenes; const lg = Math.log10; const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const f = (v, d = 0) => (+v).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d });
const fade = (objs, p, from = .5) => objs.forEach(o => { o.material.transparent = true; o.material.opacity = (o.userData.op ?? 1) * clamp((p - from) / (1 - from), 0, 1); });
const ball = (v, pos, key, r = .05, ud = {}, op = 1) => { const m = new THREE.Mesh(new THREE.SphereGeometry(r, 12, 8), Viz.mat(key, { emissive: key, emissiveIntensity: .35 })); m.position.copy(pos); m.userData = { ...ud, op }; v.rig.add(m); v.pickables.push(m); return m; };

/* ---------- edsfront: measured EDS fronts against the model with and without narrowing ---------- */
SC && SC.eds && Viz.view('[data-viz="edsfront"]', v => {
  const D = SC.eds, cols = ['blue', 'cyan', 'green'];
  const A = Viz.axes(v, { size: [4.6, 2.2, 3.0], x: { lo: 0, hi: 1.7, ticks: [0, .5, 1, 1.5], fmt: t => `${t}`, label: 'depth / half-thickness depth' }, y: { lo: 0, hi: 1.15, ticks: [0, .5, 1], fmt: t => `${t}`, label: 'Al signal / plateau' }, z: { lo: -.6, hi: 2.6, ticks: [] } });
  const dots = [], rbs = [], labs = [];
  D.forEach((d, i) => { const X = d.grid.map(z => z / d.z50);
    const a = Viz.ribbon(v, A, X, d.with_narrowing, { color: 'rose', width: .2, z: i + .16, name: 'with narrowing', skirtOpacity: .06 }); a.mesh.userData = { ...a.mesh.userData, d, kind: 'with channel narrowing', y: d.with_narrowing, X };
    const b = Viz.ribbon(v, A, X, d.without, { color: 'grey', width: .2, z: i - .16, name: 'without', skirtOpacity: .06 }); b.mesh.userData = { ...b.mesh.userData, d, kind: 'without narrowing', y: d.without, X }; rbs.push(a, b);
    for (let k = 0; k < d.z.length; k += 3) { const x = d.z[k] / d.z50; if (x > 1.7) break; dots.push(ball(v, A.map(x, clamp(d.y[k], 0, 1.15), i), cols[i], .032, { d, k })); }
    labs.push(Viz.label(v, `${d.cycles} cycles, ${Math.round(d.fill * 100)} % filled · rms ${f(d.rms_with, 3)} with narrowing, ${f(d.rms_without, 3)} without`, [A.map(1.7, 0, i).x + .2, A.map(0, .1).y, A.map(0, 0, i).z], { px: 15, scale: .064, color: cols[i], bg: false })); });
  labs.push(Viz.label(v, 'points: electron-microscope line scans · rose: model with narrowing · grey: without', [A.map(.75, 0, 0).x, A.map(0, 0).y - .32, A.map(0, 0, 2.6).z + .35], { px: 15, scale: .066, color: 'ink' }));
  v.update = vv => { rbs.forEach((r, i) => r.draw(clamp(vv.p * 1.4 - Math.floor(i / 2) * .1, 0, 1))); fade(dots, vv.p, .2); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.d) return null;
    if (u.k != null) return `<b>${u.d.cycles} cycles, measured (EDS)</b>depth ${f(u.d.z[u.k] / 2, 0)} µm (${f(u.d.z[u.k] / u.d.z50, 2)} half depths): signal ${f(u.d.y[u.k], 2)} of the plateau`;
    const i = Viz.ribbonIndex(hit, A); return `<b>${u.d.cycles} cycles, model ${u.kind}</b>at ${f(u.X[i], 2)} half depths: ${f(u.y[i], 2)} · rms against the scan ${f(u.kind.startsWith('with ') ? u.d.rms_with : u.d.rms_without, 3)} · hydraulic kernel, s = 0.04, 1 µm spot`; };
  v.rot = -.4; v.tilt = .12;
}, { camera: [0, 2.0, 8.0], target: [0, -.1, 0], fov: 32 });

/* ---------- slitbench: diffusion along a slit against kinetic theory ---------- */
SC && SC.bench && Viz.view('[data-viz="slitbench"]', v => {
  const B = SC.bench, X = B.width.map(lg), xl = -2.15, xh = 2.15, yl = lg(.008), yh = lg(4);
  const A = Viz.axes(v, { size: [4.4, 2.6, 2.8], x: { lo: xl, hi: xh, ticks: [-2, -1, 0, 1, 2], fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'slit width / mean free path (log)' }, y: { lo: yl, hi: yh, ticks: [.01, .1, 1].map(lg), fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'diffusivity along the slit, reduced (log)' }, z: { lo: -.6, hi: 2.6, ticks: [] } });
  const lanes = [['Monte Carlo, hard-sphere collisions', 'rose', B.mc], ['Monte Carlo, relaxation collisions', 'green', B.relax], ['hydraulic diameter + gas diffusion', 'grey', B.hydraulic]];
  const rbs = [], dots = [], labs = [];
  lanes.forEach(([name, key, Y], i) => { const r = Viz.ribbon(v, A, X, B.theory.map(lg), { color: 'cyan', width: .5, z: i, name: 'kinetic theory', skirtOpacity: .06 }); r.mesh.userData = { ...r.mesh.userData, th: 1 }; rbs.push(r);
    const dev = []; Y.forEach((y, k) => { if (y == null) return; dots.push(ball(v, A.map(X[k], lg(y), i), key, .07, { name, k, y, th: B.theory[k] })); dev.push(y / B.theory[k]); });
    labs.push(Viz.label(v, `${name} · ${f(Math.min(...dev), 2)} to ${f(Math.max(...dev), 2)} × theory`, [A.map(xh, 0, i).x + .2, A.map(0, yl).y + .12, A.map(0, 0, i).z], { px: 15, scale: .064, color: key, bg: false })); });
  labs.push(Viz.label(v, 'ribbon: kinetic theory for Ne–Ar between plates (2004)', [A.map(-.4, 0, 0).x, A.map(0, yh).y + .25, A.map(0, 0, 1).z], { px: 16, scale: .07, color: 'cyan', weight: 700 }), Viz.label(v, '← low pressure, narrow gaps', [A.map(-1.5, 0, 0).x, A.map(0, yl).y - .3, A.map(0, 0, 2.6).z + .35], { px: 15, scale: .066, color: 'muted' }));
  v.update = vv => { rbs.forEach((r, i) => r.draw(clamp(vv.p * 1.4 - i * .1, 0, 1))); fade(dots, vv.p, .4); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.name) return `<b>${u.name}</b>slit ${B.width[u.k]} mean free paths wide: ${f(u.y, 3)} against ${f(u.th, 3)} from kinetic theory (${f(u.y / u.th, 2)} ×)`; if (u.th) return '<b>kinetic theory</b>Siewert and Valougeorgis 2004: diffusion of a Ne–Ar mixture between diffuse plates, heavy species, in units of H v / 2'; return null; };
  v.rot = -.45; v.tilt = .1;
}, { camera: [0, 1.9, 8.2], target: [0, -.05, 0], fov: 32 });
