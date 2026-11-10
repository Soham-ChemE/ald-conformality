/* Body scenes for the ALD page, drawn by the shared Viz renderer (assets/viz.js) from window.ALD (assets/data.js):
   stack    the chamber as an exploded cross-section with a magnified patch of the hole array;
   map      the pulse-time map: saturating pulse against TMA fraction for blanket and AR 20/50/100 wafers, with supply limits;
   wiwnu    hole-bottom coverage across the wafer for a pulse stopped at 8, 10, 12 and 14 s;
   rule     precursor utilization against total reactive area: chamber runs and two measured powder reactors;
   pillar   the PillarHall saturation profile and the fitted hole model;
   hole     one hole cut open, coloured by coverage against depth at the time of the sticky step in view;
   doe      the JMP response surface over TMA fraction and carrier flow with the 25 runs;
   coo      precursor use against pulse time (cost of ownership) with the surrogate's Pareto front;
   water    water coverage and outlet signal for four desorption constants;
   peald    the recombination ceiling: largest reachable aspect ratio against wall loss probability;
   spatial  temporal against spatial cycle time from the rule;
   robust   the saturating pulse under every removed idealization as a bar field. */
import { Viz } from './viz.js';
import { THREE, redrawLabel } from './three-kit.js';
const D = window.ALD, H = D.hero; const lg = Math.log10; const f = (v, d = 0) => (+v).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d });
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const sphere = (v, A, x, y, z, col, r = .06, data = null) => { const m = new THREE.Mesh(new THREE.SphereGeometry(r, 14, 12), Viz.mat(col, { emissive: col, emissiveIntensity: .55 })); m.position.copy(A.map(x, y, z)); m.userData = data || {}; v.rig.add(m); v.pickables.push(m); return m; };
const fadeIn = (objs, p, from = .6) => objs.forEach(o => { o.material.transparent = true; o.material.opacity = clamp((p - from) / (1 - from), 0, 1); });

/* ---------- stack: chamber cross-section and the hole pattern ---------- */
Viz.view('[data-viz="stack"]', v => {
  const G = D.geometry, SC = 1 / 60, WR = G.wafer_mm / 2 * SC, CR = G.chamber_mm / 2 * SC, gap = G.gap_mm * SC * 3;
  const wall = new THREE.Mesh(new THREE.CylinderGeometry(CR, CR, gap + 1.2, 64, 1, true, Math.PI * .15, Math.PI * 1.7), Viz.mat('steel', { metalness: .6, roughness: .4, side: THREE.DoubleSide, opacity: .4 })); wall.position.y = .25;
  const wafer = new THREE.Mesh(new THREE.CylinderGeometry(WR, WR, .06, 96), Viz.mat('violet', { metalness: .5, roughness: .35 }));
  const stage = new THREE.Mesh(new THREE.CylinderGeometry(WR + .25, WR + .25, .22, 96), Viz.mat('grey', { roughness: .6 })); stage.position.y = -.14;
  const floor = new THREE.Mesh(new THREE.RingGeometry(WR + .02, CR, 96), Viz.mat('grey', { roughness: .9, side: THREE.DoubleSide, opacity: .8 })); floor.rotation.x = -Math.PI / 2; floor.position.y = -.01;
  const head = new THREE.Mesh(new THREE.CylinderGeometry(WR, WR, .18, 64), Viz.mat('blue', { metalness: .4, roughness: .4, emissive: 'blue', emissiveIntensity: .2 })); head.position.y = gap + .09;
  const plate = new THREE.Mesh(new THREE.CylinderGeometry(CR, CR, .1, 64), Viz.mat('steel', { metalness: .6, roughness: .4, opacity: .7 })); plate.position.y = gap + .25;
  const outlet = new THREE.Mesh(new THREE.TorusGeometry(CR, .07, 12, 64), Viz.mat('hot', { emissive: 'hot', emissiveIntensity: .35 })); outlet.rotation.x = Math.PI / 2; outlet.position.y = gap / 2;
  // magnified patch of the pattern: a slab with a hexagonal array of holes, drawn as dark cylinders standing in a green slab
  const patch = new THREE.Group(); const slab = new THREE.Mesh(new THREE.BoxGeometry(1.6, .5, 1.4), Viz.mat('green', { roughness: .5, opacity: .9 })); slab.position.y = .25; patch.add(slab);
  const hg = new THREE.CylinderGeometry(.065, .065, .52, 12), hm = Viz.mat('ink', { roughness: .9, metalness: 0 }); for (let i = -3; i <= 3; i++) for (let j = -2; j <= 2; j++) { const m = new THREE.Mesh(hg, hm); m.position.set(i * .22 + (j % 2 ? .11 : 0), .26, j * .25); patch.add(m); }
  patch.position.set(WR * .35, gap + 1.5, WR * .55); patch.scale.setScalar(.001);
  const ex = Viz.exploded(v, [{ mesh: stage, from: [0, -1.6, 0], delay: 0, span: .22 }, { mesh: wafer, from: [0, -1.2, 0], delay: .1, span: .22 }, { mesh: floor, from: [0, -1, 0], delay: .18, span: .22 }, { mesh: outlet, from: [2, 0, 0], delay: .26, span: .22 }, { mesh: wall, from: [0, 1.2, 0], delay: .34, span: .25 }, { mesh: plate, from: [0, 2.2, 0], delay: .44, span: .25 }, { mesh: head, from: [0, 1.6, 0], delay: .5, span: .25 }, { mesh: patch, from: [0, 1.4, 0], delay: .66, span: .3, fade: false }]);
  const labels = [[`300 mm wafer on its stage · 150 faces, each with its own hole model`, [0, -.6, WR + .6], .2], [`${G.gap_mm} mm gap (drawn 3×) · N₂ 300 Pa, 573 K`, [CR + 1.0, gap / 2, .4], .4], ['outlet, r = 200 mm', [-CR - 1.2, gap / 2, .8], .35], ['showerhead: total-flux inlet, 300 sccm', [0, gap + .75, .9], .6], ['the pattern: 100 nm holes, AR 100,\n180 nm hexagonal pitch, 14 % open area', [WR * .35, gap + 2.45, WR * .55], .85]];
  const labs = labels.map(([t, p, d]) => { const s = Viz.label(v, t, p, { px: 20, scale: .09 }); s.userData.d = d; return s; });
  v.update = vv => { ex.draw(vv.p); const ps = clamp((vv.p - .66) / .3, 0, 1); patch.scale.setScalar(.001 + ps * (1 - .001)); patch.rotation.y = vv.p * .6; labs.forEach(l => l.material.opacity = clamp((vv.p - l.userData.d - .12) * 4, 0, 1)); return false; };
  v.hover = () => null; v.spin = .1; v.tilt = .22;
}, { camera: [0, 3.8, 11], target: [0, .9, 0], fov: 32 });

/* ---------- map: pulse time against TMA fraction, with supply limits ---------- */
Viz.view('[data-viz="map"]', v => {
  const R = D.map.rows;
  const A = Viz.axes(v, { size: [4.2, 2.4, 2.6], x: { lo: lg(.004), hi: lg(.12), ticks: [lg(.005), lg(.01), lg(.03), lg(.1)], fmt: t => `${+(Math.pow(10, t) * 100).toFixed(1)} %`, label: 'TMA fraction (log)' }, y: { lo: lg(.05), hi: lg(60), ticks: [lg(.1), lg(1), lg(10)], fmt: t => `${+Math.pow(10, t).toFixed(1)} s`, label: 'saturating pulse (log)' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const series = [['blanket', 0, 'grey', 'blanket wafer'], ['patterned', 20, 'green', 'AR 20'], ['patterned', 50, 'warm', 'AR 50'], ['patterned', 100, 'rose', 'AR 100']];
  const rbs = series.map(([kind, ar, col, name], i) => { const pts = R.filter(r => r.kind === kind && (kind === 'blanket' || r.AR === ar)).sort((a, b) => a.x - b.x); const rb = Viz.ribbon(v, A, pts.map(r => lg(r.x)), pts.map(r => lg(r.t90_worst)), { color: col, width: .22, z: i, name }); rb.mesh.userData.pts = pts; rb.mesh.userData.name = name;
    const sup = Viz.ribbon(v, A, pts.map(r => lg(r.x)), pts.map(r => lg(r.t_supply)), { color: 'cool', width: .08, z: i, skirt: false }); sup.mesh.userData.pts = pts; sup.mesh.userData.name = name + ' supply limit'; sup.mesh.material.opacity = .85; sup.mesh.material.transparent = true;
    const sph = pts.map(r => sphere(v, A, lg(r.x), lg(r.t90_worst), i, col, .055, { r, name }));
    const lab = Viz.label(v, `${name} · ${kind === 'blanket' ? `${D.map.blanket_ratio[0]}–${D.map.blanket_ratio[1]}× supply` : `wall ${f(pts[0].wall, 0)}×`}`, [A.map(lg(.004), 0, i).x - .3, A.map(0, lg(pts[0].t90_worst)).y + .22, A.map(0, 0, i).z], { px: 19, scale: .08, color: col, weight: 700 }); return { rb, sup, sph, lab }; });
  v.update = vv => { rbs.forEach((s, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); s.rb.draw(p); s.sup.draw(p); s.lab.material.opacity = p > .95 ? 1 : 0; fadeIn(s.sph, p, .5); }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.r) return `<b>${u.name}, ${f(u.r.x * 100, 1)} % TMA</b>saturating pulse ${f(u.r.t90_worst, 2)} s · supply limit ${f(u.r.t_supply, 2)} s · ratio ${f(u.r.ratio, 2)}`; if (u.pts) { const i = Viz.ribbonIndex(hit, A); const r = u.pts[i]; return `<b>${u.name}</b>${f(r.x * 100, 1)} % TMA: ${u.name.includes('supply') ? f(r.t_supply, 2) : f(r.t90_worst, 2)} s`; } return null; };
  v.rot = -.4; v.tilt = .12;
}, { camera: [0, 2.2, 7.6], target: [0, -.1, 0], fov: 32 });

/* ---------- wiwnu: hole-bottom coverage across the wafer at four stop times ---------- */
Viz.view('[data-viz="wiwnu"]', v => {
  const cuts = [[8, 'rose'], [10, 'warm'], [12, 'cyan'], [14, 'green']]; const rr = H.radii_mm;
  const A = Viz.axes(v, { size: [4.6, 1.9, 2.2], x: { lo: 0, hi: 150, ticks: [0, 50, 100, 150], fmt: t => `${t} mm`, label: 'wafer radius' }, y: { lo: 0, hi: 1.02, ticks: [0, .5, .9, 1], fmt: t => `${Math.round(t * 100)} %`, label: 'hole-bottom coverage' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const idx = t => { let b = 0, bd = 1e9; H.times.forEach((tt, i) => { const d = Math.abs(tt - t); if (d < bd) { bd = d; b = i; } }); return b; };
  const rbs = cuts.map(([t, col], i) => { const row = H.cov.thetaBottom[idx(t)]; const rb = Viz.ribbon(v, A, rr, row, { color: col, width: .24, z: i, name: `pulse stopped at ${t} s` }); rb.mesh.userData.t = t; rb.mesh.userData.row = row; const lab = Viz.label(v, `${t} s · centre ${Math.round(row[0] * 100)} %, edge ${Math.round(row[rr.length - 1] * 100)} %`, [A.map(158, 0, i).x + .5, A.map(0, row[0]).y + .18, A.map(0, 0, i).z], { px: 19, scale: .08, color: col, weight: 700 }); return { rb, lab }; });
  Viz.line(v, [[A.map(0, .9, -.5).x, A.map(0, .9).y, A.map(0, 0, -.5).z], [A.map(150, .9, 3.5).x, A.map(0, .9).y, A.map(0, 0, 3.5).z]], 'muted');
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.row) return null; const i = Viz.ribbonIndex(hit, A); return `<b>pulse stopped at ${u.t} s</b>r = ${f(rr[i], 0)} mm: hole-bottom coverage ${f(u.row[i] * 100, 1)} %`; };
  v.rot = -.3; v.tilt = .1;
}, { camera: [0, 1.9, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- rule: utilization against total reactive area ---------- */
Viz.view('[data-viz="rule"]', v => {
  const V = D.validation; const A = Viz.axes(v, { size: [4.4, 2.3, 1.6], x: { lo: lg(.05), hi: lg(3000), ticks: [lg(.1), lg(1), lg(10), lg(100), lg(1000)], fmt: t => `${+Math.pow(10, t).toFixed(1)} m²`, label: 'total reactive surface area (log)' }, y: { lo: 0, hi: 1.05, ticks: [0, .25, .5, .75, 1], fmt: t => `${Math.round(t * 100)} %`, label: 'precursor utilization' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const pts = [];
  V.points.forEach(p => pts.push(sphere(v, A, lg(p.S_m2), p.util, 0, p.kind === 'blanket' ? 'grey' : 'green', .07, { html: `<b>${p.kind === 'blanket' ? 'blanket wafer' : `patterned wafer, AR ${p.AR}`}, ${f(p.x * 100, 1)} % TMA (chamber model)</b>${f(p.S_m2, 2)} m² · utilization ${f(p.util * 100, 0)} % · saturation time ${f(p.tau, 2)}× supply limit` })));
  const fb = V.fluidized; pts.push(sphere(v, A, lg(fb.S_m2), fb.util, 0, 'rose', .1, { html: `<b>measured: fluidized bed, 940 m² silica gel</b>${fb.source}<br>utilization ${f(fb.util * 100, 1)} % (+${f(fb.util_hi * 100, 0)} / −${f(fb.util_lo * 100, 0)}) · saturation time ${f(fb.tau, 2)}× supply limit` }));
  Viz.line(v, [[A.map(lg(fb.S_m2), fb.util - fb.util_lo, 0).x, A.map(0, fb.util - fb.util_lo).y, 0], [A.map(lg(fb.S_m2), fb.util + fb.util_hi, 0).x, A.map(0, Math.min(1, fb.util + fb.util_hi)).y, 0]], 'rose');
  V.drum.S_m2.forEach((s, i) => pts.push(sphere(v, A, lg(s), V.drum.util[i], 0, 'violet', .08, { html: `<b>measured: rotary drum, ${V.drum.mass_g[i]} g silica gel (${f(s, 0)} m²)</b>${V.drum.source}<br>utilization ${f(V.drum.util[i] * 100, 1)} %` })));
  const tags = [Viz.label(v, 'blanket 300 mm wafers:\nkinetics-limited', [A.map(lg(.08), .35, 0).x, A.map(0, .35).y, .4], { px: 18, scale: .08, color: 'grey' }), Viz.label(v, 'patterned wafers:\nsupply-limited', [A.map(lg(4), .98, 0).x, A.map(0, .98).y + .3, .3], { px: 18, scale: .08, color: 'green' }), Viz.label(v, 'fluidized bed', [A.map(lg(940), .92, 0).x, A.map(0, .92).y + .3, .3], { px: 18, scale: .08, color: 'rose' }), Viz.label(v, 'rotary drum', [A.map(lg(300), .5, 0).x, A.map(0, .5).y - .3, .3], { px: 18, scale: .08, color: 'violet' })];
  v.update = vv => { pts.forEach((m, i) => { const p = clamp(vv.p * 1.4 - i / pts.length * .4, 0, 1); m.scale.setScalar(.001 + p); }); tags.forEach(t => t.material.opacity = clamp(vv.p * 2 - 1, 0, 1)); return false; };
  v.hover = hit => hit.object.userData.html || null; v.rot = -.15; v.tilt = .08;
}, { camera: [0, 1.6, 7.2], target: [0, -.1, 0], fov: 32 });

/* ---------- pillar: the PillarHall profile and the fitted hole model ---------- */
Viz.view('[data-viz="pillar"]', v => {
  const P = D.pillarhall; const ok = P.measured.map(y => y != null); const xs = P.AR.filter((_, i) => ok[i]), ys = P.measured.filter((_, i) => ok[i]), sd = P.std.filter((_, i) => ok[i]);
  const A = Viz.axes(v, { size: [4.4, 2.2, 2.2], x: { lo: 0, hi: 360, ticks: [0, 100, 200, 300], label: 'depth / gap (aspect ratio)' }, y: { lo: -.02, hi: 1.12, ticks: [0, .5, 1], fmt: t => `${t}`, label: 'thickness / plateau' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const meas = Viz.ribbon(v, A, xs, ys, { color: 'grey', width: .3, z: 0, name: 'measured' }); meas.mesh.userData.kind = 'meas';
  const band = []; for (let i = 0; i < xs.length; i += 3) { const m = new THREE.Mesh(new THREE.BoxGeometry(.03, Math.max(.01, (A.map(0, sd[i]).y - A.map(0, 0).y) * 2), .03), Viz.mat('grey', { opacity: .5 })); m.position.copy(A.map(xs[i], ys[i], 0)); v.rig.add(m); band.push(m); }
  const fit = Viz.ribbon(v, A, P.AR, P.fit, { color: 'rose', width: .24, z: 1, name: 'fit' }); fit.mesh.userData.kind = 'fit';
  const a1 = Viz.ribbon(v, A, P.AR, P.alt['1e-3'], { color: 'cool', width: .14, z: 2, name: 's = 1e-3', skirt: false }); a1.mesh.userData.kind = 'a1';
  const a2 = Viz.ribbon(v, A, P.AR, P.alt['1e-2'], { color: 'cool', width: .14, z: 3, name: 's = 1e-2', skirt: false }); a2.mesh.userData.kind = 'a2';
  const labs = [Viz.label(v, `measured, ${P.n_sheets} profiles · plateau ${P.plateau_nm} nm`, [A.map(40, 1, 0).x, A.map(0, 1.06).y + .25, A.map(0, 0, 0).z], { px: 19, scale: .08, color: 'grey', weight: 700 }), Viz.label(v, `fit: s = ${P.s_fit}, ${P.p_fit_Pa} Pa TMA · rms ${P.rms}`, [A.map(250, .5, 1).x, A.map(0, .55).y, A.map(0, 0, 1).z], { px: 19, scale: .08, color: 'rose', weight: 700 }), Viz.label(v, 's = 1e-3 (too shallow)', [A.map(120, .5, 2).x, A.map(0, .6).y, A.map(0, 0, 2).z], { px: 17, scale: .075, color: 'cool' }), Viz.label(v, 's = 1e-2 (too sharp)', [A.map(300, .4, 3).x, A.map(0, .45).y, A.map(0, 0, 3).z], { px: 17, scale: .075, color: 'cool' })];
  const series = { meas: [xs, ys], fit: [P.AR, P.fit], a1: [P.AR, P.alt['1e-3']], a2: [P.AR, P.alt['1e-2']] };
  v.update = vv => { [meas, fit, a1, a2].forEach((r, i) => r.draw(clamp(vv.p * 1.3 - i * .08, 0, 1))); fadeIn(band, vv.p, .3); labs.forEach(l => l.material.opacity = vv.p > .95 ? 1 : 0); return false; };
  v.hover = hit => { const k = hit.object.userData.kind; if (!k) return null; const i = Viz.ribbonIndex(hit, A); const [X, Y] = series[k]; return `<b>${{ meas: 'Zenodo 3974438, mean of 8 profiles', fit: 'hole model, fitted', a1: 'hole model, s = 1e-3', a2: 'hole model, s = 1e-2' }[k]}</b>AR ${f(X[i], 0)}: thickness / plateau ${f(Y[i], 3)}`; };
  v.rot = -.35; v.tilt = .12;
}, { camera: [0, 2.0, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- hole: cut open, coloured by coverage at the step time ---------- */
const holeView = Viz.view('[data-viz="hole"]', v => {
  const HH = 2.8, HR = .5, ND = H.hole.depths.length, NA = 28; const geo = new THREE.BufferGeometry(); const pos = new Float32Array(ND * (NA + 1) * 3), col = new Float32Array(ND * (NA + 1) * 3), idx = [];
  for (let i = 0; i < ND; i++) for (let j = 0; j <= NA; j++) { const k = i * (NA + 1) + j, th = Math.PI * .5 + j / NA * Math.PI, y = HH / 2 - i / (ND - 1) * HH; pos[3 * k] = -1.2 + Math.cos(th) * HR; pos[3 * k + 1] = y; pos[3 * k + 2] = Math.sin(th) * HR; if (i < ND - 1 && j < NA) { const a = k, b = k + 1, c = k + NA + 1, d = k + NA + 2; idx.push(a, c, b, b, c, d); } }
  geo.setAttribute('position', new THREE.BufferAttribute(pos, 3)); geo.setAttribute('color', new THREE.BufferAttribute(col, 3)); geo.setIndex(idx); geo.computeVertexNormals();
  const hole = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: .5, metalness: .2 })); v.rig.add(hole);
  const shell = new THREE.Mesh(new THREE.CylinderGeometry(HR + .2, HR + .2, HH, 32, 1, false, Math.PI * .5, Math.PI), Viz.mat('grey', { side: THREE.DoubleSide, opacity: .5 })); shell.position.x = -1.2; v.rig.add(shell);
  const base = new THREE.Mesh(new THREE.BoxGeometry(2.2, .12, 1.6), Viz.mat('grey', { roughness: .8 })); base.position.set(-1.2, -HH / 2 - .06, 0); v.rig.add(base);
  const front = new THREE.Mesh(new THREE.TorusGeometry(HR + .22, .025, 8, 48, Math.PI), Viz.mat('rose', { emissive: 'rose', emissiveIntensity: .8 })); front.rotation.x = Math.PI / 2; front.rotation.z = Math.PI / 2; front.position.set(-1.2, HH / 2, 0); v.rig.add(front);
  // ribbon of coverage against depth beside it: x = coverage, y = depth (mouth at the top)
  const A = Viz.axes(v, { size: [1.6, HH, .8], x: { lo: 0, hi: 1, ticks: [0, .5, 1], fmt: t => `${Math.round(t * 100)} %`, label: 'coverage' }, y: { lo: 1, hi: 0, ticks: [0, 1], fmt: t => t === 0 ? 'mouth' : 'bottom', label: '' }, z: { lo: -.5, hi: .5, ticks: [] } }); A.group.position.x = 1.3;
  const rbGeo = new THREE.BufferGeometry(); const rpos = new Float32Array(ND * 2 * 3), ridx = []; for (let i = 1; i < ND; i++) ridx.push(i * 2 - 2, i * 2 - 1, i * 2, i * 2 - 1, i * 2 + 1, i * 2); rbGeo.setAttribute('position', new THREE.BufferAttribute(rpos, 3)); rbGeo.setIndex(ridx);
  const rb = new THREE.Mesh(rbGeo, Viz.mat('green', { side: THREE.DoubleSide, emissive: 'green', emissiveIntensity: .35 })); rb.position.x = 1.3; v.rig.add(rb);
  const tLab = Viz.label(v, 't = 0.5 s', [1.3, HH / 2 + .35, 0], { px: 24, scale: .1, color: 'ink', weight: 700 });
  const cBare = new THREE.Color(), cCov = new THREE.Color(), c = new THREE.Color();
  const thetaAt = (t, i) => { const T = H.hole.times; if (t <= T[0]) return H.hole.theta[0][i] * t / T[0]; let k = 1; while (k < T.length - 1 && T[k] < t) k++; const fr = clamp((t - T[k - 1]) / (T[k] - T[k - 1]), 0, 1); return H.hole.theta[k - 1][i] + (H.hole.theta[k][i] - H.hole.theta[k - 1][i]) * fr; };
  v.time = .5; v.shown = -1;
  const paint = t => { cBare.set(v.dark ? 0x5b6676 : 0x8a95a5); cCov.set(v.pal.green); let fr = 0, bot = 0;
    for (let i = 0; i < ND; i++) { const th = thetaAt(t, i); if (th >= .5) fr = i / (ND - 1); if (i === ND - 1) bot = th; c.copy(cBare).lerp(cCov, th); for (let j = 0; j <= NA; j++) { const k = i * (NA + 1) + j; col[3 * k] = c.r; col[3 * k + 1] = c.g; col[3 * k + 2] = c.b; } const q = A.map(th, H.hole.depths[i], 0); rpos[6 * i] = q.x; rpos[6 * i + 1] = q.y; rpos[6 * i + 2] = q.z - .12; rpos[6 * i + 3] = q.x; rpos[6 * i + 4] = q.y; rpos[6 * i + 5] = q.z + .12; }
    geo.attributes.color.needsUpdate = true; rbGeo.attributes.position.needsUpdate = true; rbGeo.computeVertexNormals(); front.position.y = HH / 2 - fr * HH;
    tLab.userData.label.text = `t = ${t} s`; redrawLabel(tLab, v.dark);
    const b = document.getElementById('sfBot'), fEl = document.getElementById('sfFront'); if (b) b.textContent = Math.round(bot * 100); if (fEl) fEl.textContent = Math.round(fr * 100); v.shown = t; v.dirty = true; };
  v.update = vv => { if (vv.shown !== vv.time) { paint(vv.time); return true; } return false; };
  v.onTheme = () => paint(v.time); v.hover = () => null; v.rot = .25; v.tilt = .08; paint(.5);
}, { camera: [0, .6, 6.6], target: [0, 0, 0], fov: 34 });
document.addEventListener('ald-step', e => { if (holeView) { holeView.time = e.detail; holeView.dirty = true; } });

/* ---------- doe: response surface over TMA fraction and carrier flow ---------- */
Viz.view('[data-viz="doe"]', v => {
  const Q = D.doe; const C = Q.rsm.log10_t90_worst.coef; const pred = (cx, cQ, cp = 0, cg = 0) => C.Intercept + C.c_x * cx + C.c_Q * cQ + C.c_p * cp + C.c_g * cg;
  const xs = [], zs = []; for (let i = 0; i <= 12; i++) { xs.push(-1 + i / 6); zs.push(-1 + i / 6); } const Z = zs.map(cQ => xs.map(cx => pred(cx, cQ)));
  const A = Viz.axes(v, { size: [3.6, 2.2, 3.0], x: { lo: -1, hi: 1, ticks: [-1, 0, 1], fmt: t => `${[0.5, 1.22, 3][t + 1]} %`, label: 'TMA fraction' }, y: { lo: lg(1.5), hi: lg(80), ticks: [lg(2), lg(5), lg(10), lg(20), lg(50)], fmt: t => `${+Math.pow(10, t).toFixed(0)} s`, label: 'saturating pulse (log)' }, z: { lo: -1, hi: 1, ticks: [-1, 0, 1], fmt: t => `${[150, 300, 600][t + 1]} sccm`, label: 'carrier flow' } });
  const surf = Viz.surface(v, A, xs, zs, Z, { stops: ['green', 'cyan', 'warm', 'rose'] });
  const pcol = p => p < 150 ? 'cyan' : p < 500 ? 'green' : 'warm';
  const sph = Q.runs.map(r => sphere(v, A, r.cx, lg(r.t90), r.cQ, pcol(r.p), .07, { r }));
  const tag = Viz.label(v, `JMP fit, R² ${Q.rsm.log10_t90_worst.r2.toFixed(3)}:\nlog₁₀ t₉₀ = ${C.Intercept.toFixed(2)} ${C.c_x.toFixed(2)} x̂ ${C.c_Q.toFixed(2)} Q̂ ${C.c_p.toFixed(2)} p̂ + ${C.c_g.toFixed(2)} ĝ`, [A.map(-1, lg(80), -1).x + .4, A.map(0, lg(80)).y + .5, A.map(0, 0, -1).z], { px: 18, scale: .08, color: 'ink', weight: 600 });
  v.update = vv => { surf.draw(vv.p); fadeIn(sph, vv.p, .55); tag.material.opacity = clamp(vv.p * 2 - 1, 0, 1); return false; };
  v.hover = hit => { const r = hit.object.userData.r; if (r) return `<b>run ${r.run}: ${f(r.x * 100, 2)} % TMA, ${f(r.Q, 0)} sccm, ${f(r.p, 0)} Pa, ${f(r.g, 1)} mm</b>saturating pulse ${f(r.t90, 2)} s · utilization ${f(r.util * 100, 0)} % · edge lag ${f(r.edge_lag, 1)} %`; if (hit.object.userData.Z) return `<b>fitted surface</b>at the centre pressure (316 Pa) and gap (17.5 mm)`; return null; };
  v.rot = -.6; v.tilt = .2;
}, { camera: [0, 2.6, 7.6], target: [0, -.1, 0], fov: 32 });

/* ---------- coo: precursor use against pulse time ---------- */
Viz.view('[data-viz="coo"]', v => {
  const Q = D.doe; const A = Viz.axes(v, { size: [4.4, 2.2, 1.4], x: { lo: lg(1.8), hi: lg(70), ticks: [lg(2), lg(5), lg(10), lg(20), lg(50)], fmt: t => `${+Math.pow(10, t).toFixed(0)} s`, label: 'saturating pulse (log)' }, y: { lo: 1, hi: 2.1, ticks: [1, 1.25, 1.5, 1.75, 2], fmt: t => `${t.toFixed(2)}×`, label: 'precursor delivered / stoichiometric minimum' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const sph = Q.runs.map(r => sphere(v, A, lg(r.t90), 1 / r.util, 0, 'grey', .07, { r }));
  let frontRb = null, flab = null;
  if (D.pareto && D.pareto.length) { const F = [...D.pareto].sort((a, b) => a.t90_s - b.t90_s); frontRb = Viz.ribbon(v, A, F.map(q => lg(q.t90_s)), F.map(q => q.precursor_over_minimum), { color: 'green', width: .16, z: 0, name: 'front' }); frontRb.mesh.userData.F = F; flab = Viz.label(v, 'surrogate front: least precursor for a given pulse, edge lag ≤ 5 %', [A.map(lg(10), 1.05, 0).x, A.map(0, 1.0).y - .35, .4], { px: 18, scale: .08, color: 'green', weight: 700 }); }
  const centre = Q.runs.find(r => r.cx === 0 && r.cQ === 0 && r.cp === 0 && r.cg === 0); if (centre) sphere(v, A, lg(centre.t90), 1 / centre.util, 0, 'rose', .1, { r: centre });
  const clab = centre ? Viz.label(v, `centre: ${f(centre.t90, 1)} s, ${f(centre.util * 100, 0)} % used`, [A.map(lg(centre.t90), 1 / centre.util, 0).x + .9, A.map(0, 1 / centre.util).y + .3, .3], { px: 18, scale: .08, color: 'rose', weight: 700 }) : null;
  const hi = Viz.label(v, 'high carrier flow, low pressure:\nfast gas, half the TMA wasted', [A.map(lg(5), 1.95, 0).x, A.map(0, 1.95).y + .25, .3], { px: 17, scale: .075, color: 'muted' });
  v.update = vv => { fadeIn(sph, vv.p, .3); if (frontRb) frontRb.draw(clamp(vv.p * 1.3 - .2, 0, 1)); [flab, clab, hi].forEach(l => l && (l.material.opacity = clamp(vv.p * 2 - 1, 0, 1))); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.r) return `<b>run ${u.r.run}: ${f(u.r.x * 100, 2)} % TMA, ${f(u.r.Q, 0)} sccm, ${f(u.r.p, 0)} Pa, ${f(u.r.g, 1)} mm</b>pulse ${f(u.r.t90, 2)} s · precursor ${f(1 / u.r.util, 2)}× the minimum (utilization ${f(u.r.util * 100, 0)} %)`; if (u.F) { const i = Viz.ribbonIndex(hit, A); const q = u.F[i]; return `<b>front: ${f(q.t90_s, 2)} s</b>${f(q.precursor_over_minimum, 2)}× the minimum · ${f(q.x * 100, 2)} % TMA, ${f(q.Q_sccm, 0)} sccm, ${f(q.p_Pa, 0)} Pa, ${f(q.g_mm, 1)} mm`; } return null; };
  v.rot = -.2; v.tilt = .1;
}, { camera: [0, 1.7, 7.2], target: [0, -.1, 0], fov: 32 });

/* ---------- water: coverage and outlet for four desorption constants ---------- */
Viz.view('[data-viz="water"]', v => {
  const W = D.water.kdes; const cols = ['cyan', 'green', 'warm', 'rose'];
  const A = Viz.axes(v, { size: [3.8, 2.0, 2.6], x: { lo: 0, hi: 8, ticks: [0, 2, 4, 6, 8], fmt: t => `${t} s`, label: 'time (5 s pulse, then purge)' }, y: { lo: 0, hi: 1.05, ticks: [0, .5, .9, 1], fmt: t => `${Math.round(t * 100)} %`, label: 'coverage · outlet / feed' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const rbs = W.map((w, i) => { const rb = Viz.ribbon(v, A, w.times, w.centre, { color: cols[i], width: .22, z: i, name: `kdes ${w.kdes}` }); rb.mesh.userData.w = w; rb.mesh.userData.kind = 'cov'; const ob = Viz.ribbon(v, A, w.times, w.outlet, { color: cols[i], width: .08, z: i, skirt: false }); ob.mesh.userData.w = w; ob.mesh.userData.kind = 'out'; ob.mesh.material.opacity = .7; ob.mesh.material.transparent = true;
    const lab = Viz.label(v, `k_des = ${w.kdes} /s · ${w.t90 != null ? `90 % at ${f(w.t90, 2)} s` : 'never reaches 90 %'}`, [A.map(8.3, 0, i).x + .5, A.map(0, w.centre[w.centre.length - 1]).y + .15, A.map(0, 0, i).z], { px: 18, scale: .078, color: cols[i], weight: 700 }); return { rb, ob, lab }; });
  Viz.line(v, [[A.map(5, 0, -.5).x, A.map(0, 0).y, A.map(0, 0, -.5).z], [A.map(5, 0, 3.5).x, A.map(0, 0).y, A.map(0, 0, 3.5).z], [A.map(5, 1.05, 3.5).x, A.map(0, 1.05).y, A.map(0, 0, 3.5).z]], 'muted');
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .08, 0, 1); r.rb.draw(p); r.ob.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.w) return null; const i = Viz.ribbonIndex(hit, A); return `<b>k_des = ${u.w.kdes} /s, ${u.kind === 'cov' ? 'coverage at the wafer centre' : 'outlet H₂O / feed'}</b>t = ${f(u.w.times[i], 2)} s: ${f((u.kind === 'cov' ? u.w.centre : u.w.outlet)[i] * 100, 1)} %`; };
  v.rot = -.45; v.tilt = .15;
}, { camera: [0, 2.2, 7.4], target: [0, -.1, 0], fov: 32 });

/* ---------- peald: the recombination ceiling ---------- */
Viz.view('[data-viz="peald"]', v => {
  const P = D.peald; const A = Viz.axes(v, { size: [4.0, 2.2, 1.4], x: { lo: -5.2, hi: -1.3, ticks: [-5, -4, -3, -2], fmt: t => `1e${t}`, label: 'wall recombination probability (log)' }, y: { lo: lg(10), hi: lg(120), ticks: [lg(10), lg(20), lg(50), lg(100)], fmt: t => `${Math.round(Math.pow(10, t))}`, label: 'largest AR reached in 5 s (log)' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const xs = P.ceiling_r.map(lg), ys = P.ceiling_AR.map(lg); const rb = Viz.ribbon(v, A, xs, ys, { color: 'rose', width: .24, z: 0, name: 'ceiling' }); rb.mesh.userData.kind = 'c';
  const sph = P.ceiling_r.map((r, i) => sphere(v, A, lg(r), lg(P.ceiling_AR[i]), 0, 'rose', .06, { html: `<b>r = ${r.toExponential(1)}</b>largest aspect ratio whose bottom reaches 90 % in a 5 s exposure: ${f(P.ceiling_AR[i], 0)}` }));
  Viz.line(v, [[A.map(-5.2, lg(100), 0).x, A.map(0, lg(100)).y, -.5], [A.map(-1.3, lg(100), 0).x, A.map(0, lg(100)).y, -.5]], 'muted');
  const bw = A.map(lg(8e-5), 0, 0).x - A.map(lg(5e-5), 0, 0).x; const band = new THREE.Mesh(new THREE.BoxGeometry(bw, 2.2, .9), Viz.mat('cyan', { opacity: .15 })); band.position.set((A.map(lg(8e-5), 0, 0).x + A.map(lg(5e-5), 0, 0).x) / 2, 0, 0); v.rig.add(band);
  const labs = [Viz.label(v, 'AR 100: never reached in 5 s\nat this radical density', [A.map(-2.2, lg(100), 0).x, A.map(0, lg(100)).y + .3, .3], { px: 17, scale: .075, color: 'muted' }), Viz.label(v, 'O on SiO₂ / TiO₂\n(Arts 2019)', [A.map(lg(6e-5), lg(12), 0).x, A.map(0, lg(12)).y - .1, .5], { px: 17, scale: .075, color: 'cyan' })];
  v.update = vv => { rb.draw(vv.p); fadeIn(sph, vv.p, .5); labs.forEach(l => l.material.opacity = clamp(vv.p * 2 - 1, 0, 1)); return false; };
  v.hover = hit => { const u = hit.object.userData; if (u.html) return u.html; if (u.kind === 'c') { const i = Viz.ribbonIndex(hit, A); return `<b>r = ${P.ceiling_r[i].toExponential(1)}</b>largest AR in 5 s: ${f(P.ceiling_AR[i], 0)}`; } return null; };
  v.rot = -.3; v.tilt = .1;
}, { camera: [0, 1.8, 7.0], target: [0, -.1, 0], fov: 32 });

/* ---------- spatial: temporal against spatial cycle time ---------- */
Viz.view('[data-viz="spatial"]', v => {
  const Sp = D.spatial; const J = Sp.U * Sp.c0; const wall = []; for (let i = 0; i <= 40; i++) wall.push(Math.pow(10, i / 40 * 2.2));
  const phi = w => 1.2 * (1 + w) * Sp.N0; const tempo = (w, xm = 1) => 2 * phi(w) / (xm * J) + 2 * Sp.t_purge; const spat = (w, ratio) => 2 * Math.ceil(phi(w) / (ratio * J * Sp.f_zone * Sp.tau_rot)) * Sp.tau_rot;
  const A = Viz.axes(v, { size: [4.0, 2.2, 2.6], x: { lo: 0, hi: 2.2, ticks: [0, 1, 2], fmt: t => `${Math.round(Math.pow(10, t))}×`, label: 'wall area / wafer area (log)' }, y: { lo: 0, hi: lg(400), ticks: [0, 1, 2], fmt: t => `${Math.round(Math.pow(10, t))} s`, label: 'cycle time, TMA + H₂O (log)' }, z: { lo: -.5, hi: 4.5, ticks: [] } });
  const series = [['temporal, 1 % precursor', 'grey', w => tempo(w, 1)], ['temporal, 3 %', 'grey', w => tempo(w, 3)], ['spatial, zone 10× showerhead', 'cyan', w => spat(w, 10)], ['spatial, 30×', 'green', w => spat(w, 30)], ['spatial, 100×', 'warm', w => spat(w, 100)]];
  const rbs = series.map(([name, col, fn], i) => { const ys = wall.map(w => lg(fn(w))); const rb = Viz.ribbon(v, A, wall.map(lg), ys, { color: col, width: .2, z: i, name, skirt: i < 2 }); rb.mesh.userData.fn = fn; rb.mesh.userData.name = name; const lab = Viz.label(v, name, [A.map(2.3, 0, i).x + .6, A.map(0, ys[ys.length - 1]).y + .1, A.map(0, 0, i).z], { px: 17, scale: .075, color: col, weight: 700 }); return { rb, lab }; });
  [[11, 'AR 20'], [56, 'AR 100']].forEach(([w, t]) => { Viz.line(v, [[A.map(lg(w), 0, -.5).x, A.map(0, 0).y, A.map(0, 0, -.5).z], [A.map(lg(w), lg(400), -.5).x, A.map(0, lg(400)).y, A.map(0, 0, -.5).z]], 'muted'); Viz.label(v, t, [A.map(lg(w), 0, -.5).x, A.map(0, lg(400)).y + .15, A.map(0, 0, -.5).z], { px: 16, scale: .07, bg: false }); });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .06, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.fn) return null; const i = Viz.ribbonIndex(hit, A); return `<b>${u.name}</b>wall area ${f(wall[i], 0)}× the wafer: cycle ${f(u.fn(wall[i]), 1)} s`; };
  v.rot = -.45; v.tilt = .14;
}, { camera: [0, 2.3, 7.6], target: [0, -.1, 0], fov: 32 });

/* ---------- robust: every removed idealization as a bar ---------- */
Viz.view('[data-viz="robust"]', v => {
  const R = D.robust, W = D.water; const items = [
    { name: 'bare wafer', y: R.edge_ring.patterned_bare_bottom.centre, w: R.edge_ring.patterned_bare_bottom.worst, col: 'green', src: 'reactor/showerhead_pulse_patterned_AR100_long' }, { name: 'edge ring added', y: R.edge_ring.patterned_ring_bottom.centre, w: R.edge_ring.patterned_ring_bottom.worst, col: 'green', src: 'reactor/edgering_pulse_patterned_AR100' },
    ...R.temperature.map(t => ({ name: `${t.T_C} °C`, y: t.bottom_centre, w: t.bottom_worst, col: 'green', src: 'reactor/temperature (fixed sccm and fraction)' })),
    { name: '50 nm holes', y: R.hole_size[0].centre, w: R.hole_size[0].worst, col: 'green', src: 'reactor/holesize_d50nm_AR100 (same wall area)' }, { name: '100 nm holes', y: R.edge_ring.patterned_bare_bottom.centre, w: R.edge_ring.patterned_bare_bottom.worst, col: 'green', src: 'reference case' }, { name: '200 nm holes', y: R.hole_size[1].centre, w: R.hole_size[1].worst, col: 'green', src: 'reactor/holesize_d200nm_AR100 (same wall area)' },
    { name: 'H₂O half-cycle', y: W.patterned_h2o_t90_bottom, w: W.patterned_h2o_t90_bottom, col: 'cyan', src: 'reactor/showerhead_h2o_patterned_AR100, s = 1e-4' }];
  const A = Viz.axes(v, { size: [4.8, 2.0, 1.2], x: { lo: -.5, hi: items.length - .5, ticks: [] }, y: { lo: 0, hi: 16, ticks: [0, 5, 10, 14.05], fmt: t => `${t} s`, label: 'saturating pulse, hole bottoms' }, z: { lo: -.5, hi: .5, ticks: [] } });
  const bars = items.map((it, i) => ({ x: i, z: 0, y: it.y, color: it.col, it }));
  const B = Viz.bars(v, A, bars, { width: .5, depth: .5 });
  Viz.line(v, [[A.map(-.5, 14.05, 0).x, A.map(0, 14.05).y, .35], [A.map(items.length - .5, 14.05, 0).x, A.map(0, 14.05).y, .35]], 'muted');
  bars.forEach((b, i) => { const p = A.map(i, b.y, 0); b.lab = Viz.label(v, `${f(b.y, 2)}`, [p.x, p.y + .18, p.z], { px: 18, scale: .08 }); b.nm = Viz.label(v, b.it.name, [p.x, A.map(0, 0).y - .3, .45], { px: 16, scale: .072, bg: false }); });
  v.update = vv => { B.draw(vv.p); bars.forEach(b => b.lab.material.opacity = clamp(vv.p * 2 - 1, 0, 1)); return false; };
  v.hover = hit => { const b = bars[hit.instanceId]; return b ? `<b>${b.it.name}</b>centre ${f(b.y, 2)} s · worst on wafer ${f(b.it.w, 2)} s<br><span style="color:var(--muted)">${b.it.src}</span>` : null; };
  v.rot = -.2; v.tilt = .1;
}, { camera: [0, 1.8, 7.8], target: [0, -.2, 0], fov: 32 });

/* ---------- layout: non-uniform patterns across the wafer ---------- */
Viz.view('[data-viz="layout"]', v => {
  const L = D.layout; if (!L) return; const rr = L.radius_mm;
  // uniform reference: per-radius hole-bottom t90 from the hero case
  const uni = H.radii_mm.map((_, j) => { const col = H.cov.thetaBottom.map(row => row[j]); const i = col.findIndex(x => x >= .9); return i < 0 ? null : H.times[i]; });
  const A = Viz.axes(v, { size: [4.6, 2.2, 2.4], x: { lo: 0, hi: 150, ticks: [0, 50, 100, 150], fmt: t => `${t} mm`, label: 'wafer radius' }, y: { lo: 0, hi: 16, ticks: [0, 5, 10, 14.05], fmt: t => `${t} s`, label: 'saturating pulse, hole bottoms' }, z: { lo: -.5, hi: 2.5, ticks: [] } });
  const clean = (xs, ys) => { const X = [], Y = []; xs.forEach((x, i) => { if (ys[i] != null) { X.push(x); Y.push(ys[i]); } }); return [X, Y]; };
  const sets = [['uniform AR 100 (reference)', 'grey', ...clean(H.radii_mm, uni), 0], ['dense centre, bare edge', 'rose', ...clean(rr, L.dense.t90_bottom), 1], ['AR 100 → 60 ramp', 'warm', ...clean(rr, L.ramp.t90_bottom), 2]];
  const rbs = sets.map(([name, col, X, Y, z]) => { const rb = Viz.ribbon(v, A, X, Y, { color: col, width: .24, z, name }); rb.mesh.userData = { ...rb.mesh.userData, X, Y, name }; const lab = Viz.label(v, `${name}: ${Y[0].toFixed(2)} s at the centre`, [A.map(38, 0, z).x, A.map(0, Y[0]).y + .3, A.map(0, 0, z).z + .25], { px: 18, scale: .08, color: col, weight: 700 }); return { rb, lab }; });
  // the patterns themselves as thin ribbons at the floor: hole density (dense centre) and AR (ramp), both scaled to the axis height
  const dens = Viz.ribbon(v, A, rr, L.dense.density.map(d => d / 1.78e13 * 3), { color: 'rose', width: .08, z: 1, skirt: false }); dens.mesh.userData = { ...dens.mesh.userData, pattern: 'hole density, dense centre', X: rr, Y: L.dense.density.map(d => d / 1.78e13) }; dens.mesh.material.opacity = .6; dens.mesh.material.transparent = true;
  const ars = Viz.ribbon(v, A, rr, L.ramp.AR.map(a => a / 100 * 3), { color: 'warm', width: .08, z: 2, skirt: false }); ars.mesh.userData = { ...ars.mesh.userData, pattern: 'aspect ratio, ramp', X: rr, Y: L.ramp.AR }; ars.mesh.material.opacity = .6; ars.mesh.material.transparent = true;
  const note = Viz.label(v, 'thin ribbons: the pattern itself (hole density, aspect ratio), scaled to 3 s', [A.map(75, 0, 2.5).x, A.map(0, 3.6).y, A.map(0, 0, 2.5).z + .3], { px: 16, scale: .07, color: 'muted' });
  v.update = vv => { rbs.forEach((r, i) => { const p = clamp(vv.p * 1.3 - i * .1, 0, 1); r.rb.draw(p); r.lab.material.opacity = p > .95 ? 1 : 0; }); dens.draw(vv.p); ars.draw(vv.p); note.material.opacity = clamp(vv.p * 2 - 1, 0, 1); return false; };
  v.hover = hit => { const u = hit.object.userData; if (!u.X) return null; const i = Viz.ribbonIndex(hit, A); if (u.pattern) return `<b>${u.pattern}</b>r = ${f(u.X[i], 0)} mm: ${u.pattern.startsWith('hole') ? `${f(u.Y[i] * 100, 0)} % of 1.78e13 /m²` : `AR ${f(u.Y[i], 1)}`}`; return `<b>${u.name}</b>r = ${f(u.X[i], 0)} mm: hole bottoms reach 90 % at ${f(u.Y[i], 2)} s`; };
  v.rot = -.35; v.tilt = .12;
}, { camera: [0, 2.1, 7.6], target: [0, -.1, 0], fov: 32 });
