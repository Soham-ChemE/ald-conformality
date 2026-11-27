/* Scenes of section 03c (what sets the measured front), drawn by the shared Viz renderer from window.ALD_MECH.scenes
   (assets/data-mech.js, built by docs/build_mech_assets.py from figures/mechanism_*.csv/json and the Monte Carlo frames):
   mechheld    predicted against measured front slope for the held-out samples, one lane per model;
   mechcycles  the thickness profile after 1 to 1000 cycles from the narrowing law (Monte Carlo with N2 at 300 Pa);
   mechlimit   measured slope x depth against the filled fraction, under the narrowing limit;
   mechsens    front slope against sticking coefficient for five samples in the complete gas model, with the measured slopes;
   mechtail    coverage profiles for the wall-scattering laws and carrier-gas pressures at one front depth. */
import { Viz } from './viz.js';
import { THREE } from './three-kit.js';
const M = window.ALD_MECH, SC = M && M.scenes; const lg = Math.log10; const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const f = (v, d = 0) => (+v).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d });
const fade = (objs, p, from = .5) => objs.forEach(o => { o.material.transparent = true; o.material.opacity = (o.userData.op ?? 1) * clamp((p - from) / (1 - from), 0, 1); });
const ball = (v, pos, key, r = .05, ud = {}, op = 1) => { const m = new THREE.Mesh(new THREE.SphereGeometry(r, 14, 10), Viz.mat(key, { emissive: key, emissiveIntensity: .35 })); m.position.copy(pos); m.userData = { ...ud, op }; v.rig.add(m); v.pickables.push(m); return m; };
const nm = s => s.replace('Yim', 'Yim 2020 sample ').replace('HaimiA', 'Haimi 2021 specimen A (EDS)');
if (!SC) console.warn('mechanism scene data missing');

/* ---------- mechheld: predicted against measured slope, six lanes ---------- */
SC && Viz.view('[data-viz="mechheld"]', v => {
  const lanes = [['hydraulic', 'hydraulic', 'grey'], ['ballistic', 'ballistic', 'blue'], ['gas', 'N₂ 300 Pa', 'rose'], ['hydraulic+narrowing+spot', 'hydraulic', 'grey'], ['ballistic+narrowing+spot', 'ballistic', 'blue'], ['gas+narrowing+spot', 'N₂ 300 Pa', 'rose']];
  const lo = lg(.004), hi = lg(.065), tk = [.005, .01, .02, .05];
  const A = Viz.axes(v, { size: [3.2, 3.2, 4.2], x: { lo, hi, ticks: tk.map(lg), fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'measured front slope (log)' }, y: { lo, hi, ticks: tk.map(lg), fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'predicted (log)' }, z: { lo: -.6, hi: 5.6, ticks: [] } });
  const dots = [], stems = [], labs = [];
  lanes.forEach(([m, name, key], i) => { const z = i < 3 ? i : i + .0; Viz.line(v, [A.map(lo, lo, z).toArray(), A.map(hi, hi, z).toArray()], 'ink', { opacity: .55 });
    SC.heldout[m].forEach(d => { const a = A.map(lg(d.measured), lg(clamp(d.predicted, .004, .065)), z), b = A.map(lg(d.measured), lg(d.measured), z); stems.push(Viz.line(v, [a.toArray(), b.toArray()], key, { opacity: .5 })); dots.push(ball(v, a, key, .055, { d, m }, i < 3 ? .55 : 1)); });
    labs.push(Viz.label(v, `${name}${i < 3 ? '' : ' + narrowing + spot'} · rms ${SC.heldout_rms[m]}`, [A.map(hi, 0, z).x + .15, A.map(0, lo).y + .12, A.map(0, 0, z).z], { px: 15, scale: .066, color: key, bg: false })); });
  labs.push(Viz.label(v, 'transport model alone: nearly one slope for every sample', [A.map(lg(.012), 0, 0).x, A.map(0, hi).y + .28, A.map(0, 0, 1).z], { px: 17, scale: .075, color: 'muted', weight: 700 }), Viz.label(v, 'with channel narrowing and spot size: on the diagonal', [A.map(lg(.009), 0, 0).x, A.map(0, lo).y - .3, A.map(0, 0, 5.6).z + .35], { px: 17, scale: .075, color: 'green', weight: 700 }));
  v.update = vv => { fade(dots, vv.p, .15); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.d) return null; return `<b>${nm(u.d.sample)}</b>gap ${u.d.H_nm} nm, ${u.d.cycles} cycles · measured ${u.d.measured}, predicted ${f(u.d.predicted, 4)} (${u.m.replace(/\+/g, ' + ')}, s = ${SC.heldout_s[u.m]})`; };
  v.rot = -.62; v.tilt = .1;
}, { camera: [0, 1.9, 9.2], target: [0, 0, 0], fov: 32 });

/* ---------- mechcycles: the narrowing law applied to one simulated exposure ---------- */
SC && Viz.view('[data-viz="mechcycles"]', v => {
  const K = M.kernels.gas300, si = 3, Z = M.z, NZ = Z.length, E0 = 3.82, Es = M.E; let j = 1; while (j < Es.length - 1 && Es[j] < E0) j++; const w = (E0 - Es[j - 1]) / (Es[j] - Es[j - 1]);
  const th = Z.map((_, i) => K.theta[si][j - 1][i] + (K.theta[si][j][i] - K.theta[si][j - 1][i]) * w); const samp = z => { const u = z / 4; if (u >= NZ - 1) return 0; const i = Math.floor(u); return th[i] + (th[i + 1] - th[i]) * (u - i); };
  const cum = fill => { const o = Z.map(z => { if (fill <= 0) return samp(z); let a = 0; for (let k = 0; k <= 10; k++) a += samp(z / Math.pow(Math.max(.02, 1 - fill + fill * k / 10), K.gap_exponent)); return a / 11; }); const n0 = o[1] || 1; return o.map(x => Math.min(1, x / n0)); };
  const cases = [[1, 0, 'grey'], [250, .1, 'blue'], [500, .2, 'cyan'], [1000, .4, 'green'], [500, .96, 'rose']];
  const A = Viz.axes(v, { size: [4.6, 2.2, 2.8], x: { lo: 0, hi: 320, ticks: [0, 100, 200, 300], fmt: t => `${t}`, label: 'distance into the channel, gap heights' }, y: { lo: 0, hi: 1.05, ticks: [0, .5, 1], fmt: t => `${t}`, label: 'thickness / plateau' }, z: { lo: -.5, hi: 4.5, ticks: [] } });
  const X = Z.filter(z => z <= 320); const cross = (y, l) => { for (let i = 1; i < y.length; i++) if (y[i] < l) return Z[i - 1] + (y[i - 1] - l) / ((y[i - 1] - y[i]) || 1) * 4; return NaN; };
  const rbs = cases.map(([cyc, fill, key], i) => { const y = cum(fill), sl = .2 / (cross(y, .4) - cross(y, .6)), z50 = cross(y, .5); const rb = Viz.ribbon(v, A, X, y.slice(0, X.length), { color: key, width: .32, z: i, name: `${cyc} cycles`, skirtOpacity: .08 }); rb.mesh.userData = { ...rb.mesh.userData, cyc, fill, sl, z50, y };
    const lab = Viz.label(v, i === 4 ? `100 nm gap, 500 cycles: closed · slope ${f(sl, 4)}` : `${cyc} cycle${cyc > 1 ? 's' : ''}, ${Math.round(fill * 100)} % filled · slope ${f(sl, 4)}`, [A.map(322, 0, i).x + .2, A.map(0, .12).y, A.map(0, 0, i).z], { px: 15, scale: .066, color: key, bg: false }); return { rb, lab }; });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.4 - i * .08, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .97 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.cyc == null) return null; const i = Viz.ribbonIndex(hit, A); return `<b>${u.fill > .9 ? '100 nm gap, 500 cycles (closed)' : u.cyc + ' cycles, ' + Math.round(u.fill * 100) + ' % of the gap filled'}</b>depth ${f(X[i], 0)} gaps: thickness ${f(u.y[i], 3)} · half-thickness at ${f(u.z50, 0)}, slope ${f(u.sl, 4)} · same exposure (3.82 Pa s), s = 0.0158, N₂ at 300 Pa`; };
  v.rot = -.38; v.tilt = .12;
}, { camera: [0, 2.0, 7.8], target: [0, -.1, 0], fov: 32 });

/* ---------- mechlimit: measured slope x depth under the narrowing limit ---------- */
SC && Viz.view('[data-viz="mechlimit"]', v => {
  const kinds = ['reflectometry, 5-6 um spot', '25 um spot', '50 um spot', 'EDS (Haimi 2021)', 'Arts 2019'], cols = ['blue', 'warm', 'green', 'rose', 'violet'];
  const xl = lg(.035), xh = 0, yl = lg(.35), yh = lg(30);
  const A = Viz.axes(v, { size: [4.2, 2.6, 2.6], x: { lo: xl, hi: xh, ticks: [.05, .1, .2, .5, 1].map(lg), fmt: t => `${Math.round(Math.pow(10, t) * 100)} %`, label: 'fraction of the gap filled by the film (log)' }, y: { lo: yl, hi: yh, ticks: [.5, 1, 2, 5, 10, 20].map(lg), fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'slope × depth (log)' }, z: { lo: -.5, hi: 4.5, ticks: [] } });
  const ff = []; for (let q = xl; q <= 1e-9; q += .03) ff.push(q); const lim = q => { const fr = Math.pow(10, q); return lg((1 - fr / 2) / fr); };
  const sheet = Viz.ribbon(v, A, ff, ff.map(lim), { color: 'ink', width: 2.6, z: 2, name: 'narrowing limit', skirt: false }); sheet.mesh.material.transparent = true; sheet.mesh.material.opacity = .22; sheet.mesh.material.depthWrite = false; sheet.mesh.userData = { ...sheet.mesh.userData, lim: 1 };
  const dots = SC.limit.map(d => { const k = kinds.indexOf(d.kind); return ball(v, A.map(lg(d.fill), lg(d.sd), k), cols[k], .065, { d }); });
  const labs = kinds.map((k, i) => Viz.label(v, k.replace('um', 'µm'), [A.map(xh, 0, i).x + .2, A.map(0, yl).y + .1, A.map(0, 0, i).z], { px: 15, scale: .066, color: cols[i], bg: false }));
  labs.push(Viz.label(v, 'steepest model front: slope × depth = (1 − f/2) / f', [A.map(lg(.09), 0, 0).x, A.map(0, lg(16)).y, A.map(0, 0, 2).z], { px: 17, scale: .075, color: 'ink', weight: 700 }));
  v.update = vv => { sheet.draw(clamp(vv.p * 1.3, 0, 1)); fade(dots, vv.p, .4); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.d) { const L = (1 - u.d.fill / 2) / u.d.fill; return `<b>${u.d.sample === 'Arts' ? 'Arts 2019, TMA-limited front' : nm(u.d.sample)}</b>gap ${u.d.H_nm} nm, ${u.d.cycles} cycles, ${Math.round(u.d.fill * 100)} % filled · slope ${u.d.slope} × depth ${f(u.d.z50, 0)} = ${f(u.d.sd, 2)}, ${Math.round(100 * u.d.sd / L)} % of the limit (${f(L, 2)})`; } if (u.lim) return '<b>narrowing value</b>the steepest front the model produces when later cycles stop shorter in the narrowed gap (no carrier gas); with gas it is 1/n higher, about 10 % at 300 Pa'; return null; };
  v.rot = -.42; v.tilt = .1;
}, { camera: [0, 1.9, 8.0], target: [0, -.05, 0], fov: 32 });

/* ---------- mechsens: slope against s, complete gas model ---------- */
SC && Viz.view('[data-viz="mechsens"]', v => {
  const D = SC.sens, cols = ['violet', 'blue', 'warm', 'green', 'rose'], xl = lg(.002), xh = lg(.065), yl = lg(.005), yh = lg(.075);
  const A = Viz.axes(v, { size: [4.2, 2.6, 2.8], x: { lo: xl, hi: xh, ticks: [.003, .01, .03, .06].map(lg), fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'sticking coefficient s (log)' }, y: { lo: yl, hi: yh, ticks: [.005, .01, .02, .05].map(lg), fmt: t => `${+Math.pow(10, t).toPrecision(1)}`, label: 'front slope (log)' }, z: { lo: -1.5, hi: 4.5, ticks: [] } });
  const ss = []; for (let q = lg(.003); q <= xh; q += .04) ss.push(q); const ref = Viz.ribbon(v, A, ss, ss.map(q => lg(Math.sqrt(Math.pow(10, q) / 13.9))), { color: 'grey', width: .14, z: -1, name: 'slope method', skirt: false }); ref.mesh.userData = { ...ref.mesh.userData, ref: 1 };
  const rbs = D.map((d, i) => { const rb = Viz.ribbon(v, A, d.s.map(lg), d.slope.map(lg), { color: cols[i], width: .34, z: i, name: d.label, skirtOpacity: .07 }); rb.mesh.userData = { ...rb.mesh.userData, d };
    const y = A.map(0, lg(d.measured)).y, z = A.map(0, 0, i).z; const bar = new THREE.Mesh(new THREE.BoxGeometry(4.2, .012, .34), Viz.mat('ink', { opacity: .3, depthWrite: false })); bar.position.set(0, y, z); bar.userData = { meas: d, op: .3 }; v.rig.add(bar); v.pickables.push(bar);
    const lab = Viz.label(v, `${d.label} (${Math.round(d.fill * 100)} % filled)`, [A.map(xh, 0, i).x + .2, A.map(0, lg(d.slope[d.slope.length - 1])).y, z], { px: 15, scale: .066, color: cols[i], bg: false }); return { rb, bar, lab }; });
  const labs = [Viz.label(v, 'slope method: slope = √(s / 13.9)', [A.map(lg(.012), 0, 0).x, A.map(0, lg(.06)).y + .2, A.map(0, 0, -1).z], { px: 16, scale: .07, color: 'muted', weight: 700 }), Viz.label(v, 'flat bars: the measured slope of each sample', [A.map(lg(.006), 0, 0).x, A.map(0, yl).y - .32, A.map(0, 0, 4.5).z + .3], { px: 15, scale: .066, color: 'ink' })];
  v.update = vv => { ref.draw(clamp(vv.p * 1.3, 0, 1)); rbs.forEach((r, i) => { const p = clamp(vv.p * 1.4 - i * .08, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .97 ? 1 : 0; }); fade(rbs.map(r => r.bar), vv.p, .6); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.d) { const i = Viz.ribbonIndex(hit, A); return `<b>${u.d.label}</b>at s = ${u.d.s[i]} the model slope is ${f(u.d.slope[i], 4)} · measured ${u.d.measured}, ${Math.round(u.d.over_ceiling * 100)} % of the slope at the largest simulated s`; } if (u.meas) return `<b>measured: ${u.meas.label}</b>slope ${u.meas.measured}${u.meas.s_match ? ' · the model matches it at s = ' + u.meas.s_match : ' · above the model ceiling or not invertible'}`; if (u.ref) return '<b>slope method</b>slope = √(s / 13.9): what the slope would say if narrowing and the spot did not exist'; return null; };
  v.rot = -.4; v.tilt = .1;
}, { camera: [0, 1.9, 8.2], target: [0, -.05, 0], fov: 32 });

/* ---------- mechtail: coverage profiles for wall laws and carrier gas ---------- */
SC && Viz.view('[data-viz="mechtail"]', v => {
  const T = SC.transport, V = T.variants, n = V.length, col = { wall: 'blue', spec: 'warm', retro: 'green', lobe: 'violet', gas: 'rose', hyd: 'grey' };
  const A = Viz.axes(v, { size: [4.2, 2.2, 4.4], x: { lo: 120, hi: 360, ticks: [150, 215, 275, 350], fmt: t => `${t}`, label: 'distance, gap heights (front at 215)' }, y: { lo: 0, hi: 1.02, ticks: [0, .5, 1], fmt: t => `${t}`, label: 'coverage' }, z: { lo: -.6, hi: n - .4, ticks: [] } });
  const rbs = V.map((d, i) => { const rb = Viz.ribbon(v, A, T.z, d.theta, { color: col[d.family], width: .24, z: i, name: d.label, skirtOpacity: .07 }); rb.mesh.userData = { ...rb.mesh.userData, d };
    const lab = Viz.label(v, `${d.label} · ${f(d.tail, 3)}`, [A.map(360, 0, i).x + .2, A.map(0, .06).y, A.map(0, 0, i).z], { px: 15, scale: .064, color: col[d.family], bg: false }); return { rb, lab }; });
  const mark = new THREE.Mesh(new THREE.PlaneGeometry(.02, 2.2), Viz.mat('ink', { opacity: .0 })); v.rig.add(mark);
  Viz.line(v, [A.map(275, 0, -.6).toArray(), A.map(275, 0, n - .4).toArray()], 'ink', { opacity: .6 });
  const labs = [Viz.label(v, 'coverage 60 gaps ahead of the front →', [A.map(275, 0, 0).x, A.map(0, 0).y - .3, A.map(0, 0, n - .4).z + .35], { px: 16, scale: .07, color: 'ink', weight: 700 })];
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.5 - i * .05, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .97 ? 1 : 0; }); labs.forEach(l => l.material.opacity = vv.p > .9 ? 1 : 0); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.d) return null; const i = Viz.ribbonIndex(hit, A); return `<b>${u.d.label}</b>distance ${f(T.z[i], 0)}: coverage ${f(u.d.theta[i], 3)} · 60 gaps ahead of the front ${f(u.d.tail, 3)} · exposure for this depth ${f(u.d.exposure, 2)} Pa s (s = ${T.s})`; };
  v.rot = -.55; v.tilt = .12;
}, { camera: [0, 2.3, 9.4], target: [0, -.1, 0], fov: 32 });
