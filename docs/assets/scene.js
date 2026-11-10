/* Hero scene for the ALD page. Data: window.ALD.hero from assets/data.js: the TMA field in the gap (40 radial x 10 vertical bins at
   29 times, from the OpenFOAM case showerhead_pulse_patterned_AR100_long), the per-face coverage of the flat surface and of the
   hole bottoms against time, and the coverage-against-depth history of one hole at the wafer centre (aldconf, driven by the
   chamber's face concentration). Scroll scrubs time from 0 to 16 s: gas particles brighten with the local TMA density, the wafer
   top colours with flat coverage, and a cut-open hole beside the chamber fills from the mouth down. */
import { THREE, createStage, cameraPath, makeLabel, redrawLabel } from './three-kit.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
const D = window.ALD, H = D.hero, S = window.SITE, reduced = S.reduced, phone = S.phone;
const T_END = 16, SC = 20, GAPX = 3;                       // 1 m = 20 units; the 15 mm gap drawn 3x
const WR = .15 * SC, CR = .2 * SC, GAP = .015 * SC * GAPX;
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const lerp = (a, b, t) => a + (b - a) * t;
const mulberry = a => () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
// time interpolation helpers over the exported series
function at(times, t) { let i = 1; while (i < times.length - 1 && times[i] < t) i++; const a = times[i - 1], b = times[i]; return [i - 1, i, clamp((t - a) / ((b - a) || 1), 0, 1)]; }
const fieldAt = (() => { const [NR, NZ] = H.field_shape; return (t, ir, iz) => { const [i, j, f] = at(H.field_times, t); const k = ir * NZ + iz; return lerp(H.field[i][k], H.field[j][k], f); }; })();
const covAt = (key, t, jr) => { const [i, j, f] = at(H.times, t); return lerp(H.cov[key][i][jr], H.cov[key][j][jr], f); };
const holeAt = (t, jd) => { if (t <= H.hole.times[0]) return H.hole.theta[0][jd] * clamp(t / H.hole.times[0], 0, 1); const [i, j, f] = at(H.hole.times, t); return lerp(H.hole.theta[i][jd], H.hole.theta[j][jd], f); };

(function () {
  const canvas = document.getElementById('heroCanvas'); if (!canvas) return;
  const st = createStage(canvas, { fov: 38, bloom: { strength: .32, radius: .5, threshold: .85 }, dof: { aperture: .0004, maxblur: .003 }, exposure: 1.0, envIntensity: .5 });
  if (st.failed) return; const { scene, camera } = st;
  const ctr = new OrbitControls(camera, canvas); ctr.enableDamping = true; ctr.enableZoom = false; ctr.enablePan = false; ctr.autoRotate = !reduced; ctr.autoRotateSpeed = .4; ctr.minPolarAngle = .8; ctr.maxPolarAngle = 1.45;
  const hemi = new THREE.HemisphereLight(0xffffff, 0x223355, .9); scene.add(hemi); const dl = new THREE.DirectionalLight(0xffffff, 1.3); dl.position.set(4, 8, 5); scene.add(dl);
  const g = new THREE.Group(); scene.add(g);
  const P = () => st.pal;
  // ---- chamber parts ----
  const wall = new THREE.Mesh(new THREE.CylinderGeometry(CR, CR, GAP + .9, 64, 1, true), new THREE.MeshStandardMaterial({ color: 0x8b98a9, transparent: true, opacity: .07, side: THREE.DoubleSide, depthWrite: false })); wall.position.y = (GAP + .9) / 2 - .35; g.add(wall);
  const stageM = new THREE.Mesh(new THREE.CylinderGeometry(WR + .25, WR + .25, .22, 96), new THREE.MeshStandardMaterial({ color: 0x2b3444, roughness: .8 })); stageM.position.y = -.17; g.add(stageM);
  const floor = new THREE.Mesh(new THREE.RingGeometry(WR + .02, CR, 96), new THREE.MeshStandardMaterial({ color: 0x3a4557, roughness: .9, side: THREE.DoubleSide })); floor.rotation.x = -Math.PI / 2; floor.position.y = -.01; g.add(floor);
  const outlet = new THREE.Mesh(new THREE.TorusGeometry(CR, .045, 8, 96), new THREE.MeshStandardMaterial({ color: 0xf97316, emissive: 0xf97316, emissiveIntensity: .45 })); outlet.rotation.x = Math.PI / 2; outlet.position.y = GAP / 2; g.add(outlet);
  const head = new THREE.Mesh(new THREE.CylinderGeometry(CR, CR, .14, 64), new THREE.MeshStandardMaterial({ color: 0x3a4557, roughness: .7, transparent: true, opacity: .35, depthWrite: false })); head.position.y = GAP + .07; g.add(head);
  const headRing = new THREE.Mesh(new THREE.TorusGeometry(CR, .03, 8, 96), new THREE.MeshStandardMaterial({ color: 0x8b98a9, metalness: .6 })); headRing.rotation.x = Math.PI / 2; headRing.position.y = GAP + .07; g.add(headRing);
  // showerhead holes: a disc of small emissive dots over the wafer radius (the inlet patch spans r < 150 mm)
  const NH = phone ? 160 : 331; const dotG = new THREE.CylinderGeometry(.035, .035, .03, 8), dotM = new THREE.MeshStandardMaterial({ color: 0x60a5fa, emissive: 0x60a5fa, emissiveIntensity: .9 }); const dots = new THREE.InstancedMesh(dotG, dotM, NH); g.add(dots);
  { const m4 = new THREE.Matrix4(); let k = 0; const rings = phone ? 7 : 10; for (let i = 0; i <= rings && k < NH; i++) { const rr = i / rings * (WR - .15); const n = i === 0 ? 1 : Math.round(6 * i); for (let j = 0; j < n && k < NH; j++) { const th = j / n * Math.PI * 2; m4.makeTranslation(Math.cos(th) * rr, GAP - .005, Math.sin(th) * rr); dots.setMatrixAt(k++, m4); } } dots.count = k; dots.instanceMatrix.needsUpdate = true; }
  // ---- wafer top: polar grid coloured by flat coverage ----
  const NR = H.radii_mm.length, NT = phone ? 64 : 96; const wg = new THREE.BufferGeometry(); const wpos = new Float32Array((NR + 1) * (NT + 1) * 3), wcol = new Float32Array((NR + 1) * (NT + 1) * 3), widx = [];
  for (let i = 0; i <= NR; i++) for (let j = 0; j <= NT; j++) { const k = i * (NT + 1) + j, r = i / NR * WR, th = j / NT * Math.PI * 2; wpos[3 * k] = Math.cos(th) * r; wpos[3 * k + 1] = 0; wpos[3 * k + 2] = Math.sin(th) * r; if (i < NR && j < NT) { const a = k, b = k + 1, c = k + NT + 1, d = k + NT + 2; widx.push(a, c, b, b, c, d); } }
  wg.setAttribute('position', new THREE.BufferAttribute(wpos, 3)); wg.setAttribute('color', new THREE.BufferAttribute(wcol, 3)); wg.setIndex(widx); wg.computeVertexNormals();
  const wafer = new THREE.Mesh(wg, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: .35, metalness: .45 })); g.add(wafer);
  const waferEdge = new THREE.Mesh(new THREE.CylinderGeometry(WR, WR, .06, 96, 1, true), new THREE.MeshStandardMaterial({ color: 0x5b6676, metalness: .6, roughness: .3 })); waferEdge.position.y = -.03; g.add(waferEdge);
  // ---- gas: instanced particles in the gap, coloured and sized by the local TMA density ----
  const NP = phone ? 900 : 2600; const pg = new THREE.SphereGeometry(.024, 6, 5), pm = new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: .85, depthWrite: false }); const parts = new THREE.InstancedMesh(pg, pm, NP); g.add(parts);
  const rnd = mulberry(7), PR = [], PZ = [], PTH = []; for (let n = 0; n < NP; n++) { PR.push(Math.sqrt(rnd()) * CR * .995); PZ.push(rnd()); PTH.push(rnd() * Math.PI * 2); }
  const m4 = new THREE.Matrix4(), col = new THREE.Color(), cLow = new THREE.Color(), cHigh = new THREE.Color();
  function paintGas(t) {
    cLow.set(st.dark ? 0x1e2a44 : 0xb9c4d3); cHigh.set(P().warm);
    for (let n = 0; n < NP; n++) { const r = PR[n], z = PZ[n]; const ir = clamp(Math.floor(r / CR * H.field_shape[0]), 0, H.field_shape[0] - 1), iz = clamp(Math.floor(z * H.field_shape[1]), 0, H.field_shape[1] - 1); const c = t <= 0 ? 0 : fieldAt(t, ir, iz); const u = clamp(Math.pow(c, .6) * 1.15, 0, 1); const sc = .35 + 1.3 * u; m4.makeScale(sc, sc, sc); m4.setPosition(Math.cos(PTH[n]) * r, .03 + z * (GAP - .06), Math.sin(PTH[n]) * r); parts.setMatrixAt(n, m4); col.copy(cLow).lerp(cHigh, u); parts.setColorAt(n, col); }
    parts.instanceMatrix.needsUpdate = true; parts.instanceColor.needsUpdate = true;
  }
  const cBare = new THREE.Color(), cCov = new THREE.Color();
  function paintWafer(t) {
    cBare.set(st.dark ? 0x5b6676 : 0x8a95a5); cCov.set(P().green);
    for (let i = 0; i <= NR; i++) { const th = t <= 0 ? 0 : covAt('thetaFlat', t, Math.min(NR - 1, i)); col.copy(cBare).lerp(cCov, th); for (let j = 0; j <= NT; j++) { const k = i * (NT + 1) + j; wcol[3 * k] = col.r; wcol[3 * k + 1] = col.g; wcol[3 * k + 2] = col.b; } }
    wg.attributes.color.needsUpdate = true;
  }
  // ---- the cut-open hole beside the chamber: a half cylinder whose wall is coloured by coverage against depth ----
  const HX = CR + 1.2, HR = .42, HH = 2.4, ND = H.hole.depths.length, NA = 24; const hg = new THREE.BufferGeometry(); const hpos = new Float32Array(ND * (NA + 1) * 3), hcol = new Float32Array(ND * (NA + 1) * 3), hidx = [];
  for (let i = 0; i < ND; i++) for (let j = 0; j <= NA; j++) { const k = i * (NA + 1) + j, th = Math.PI * .5 + j / NA * Math.PI, y = HH - i / (ND - 1) * HH; hpos[3 * k] = HX + Math.cos(th) * HR; hpos[3 * k + 1] = y; hpos[3 * k + 2] = Math.sin(th) * HR; if (i < ND - 1 && j < NA) { const a = k, b = k + 1, c = k + NA + 1, d = k + NA + 2; hidx.push(a, c, b, b, c, d); } }
  hg.setAttribute('position', new THREE.BufferAttribute(hpos, 3)); hg.setAttribute('color', new THREE.BufferAttribute(hcol, 3)); hg.setIndex(hidx); hg.computeVertexNormals();
  const hole = new THREE.Mesh(hg, new THREE.MeshStandardMaterial({ vertexColors: true, side: THREE.DoubleSide, roughness: .5, metalness: .2 })); g.add(hole);
  const holeShell = new THREE.Mesh(new THREE.CylinderGeometry(HR + .18, HR + .18, HH, 32, 1, false, Math.PI * .5, Math.PI), new THREE.MeshStandardMaterial({ color: 0x3a4557, roughness: .8, side: THREE.DoubleSide, transparent: true, opacity: .55 })); holeShell.position.set(HX, HH / 2, 0); g.add(holeShell);
  const holeBase = new THREE.Mesh(new THREE.CylinderGeometry(HR + .5, HR + .5, .12, 48), new THREE.MeshStandardMaterial({ color: 0x2b3444, roughness: .8 })); holeBase.position.set(HX, -.06, 0); g.add(holeBase);
  const front = new THREE.Mesh(new THREE.TorusGeometry(HR + .2, .025, 8, 48, Math.PI), new THREE.MeshStandardMaterial({ color: 0xfb7185, emissive: 0xfb7185, emissiveIntensity: .8 })); front.rotation.x = Math.PI / 2; front.rotation.z = Math.PI / 2; front.position.set(HX, HH, 0); g.add(front);
  function paintHole(t) {
    cBare.set(st.dark ? 0x5b6676 : 0x8a95a5); cCov.set(P().green); let fr = 0;
    for (let i = 0; i < ND; i++) { const th = t <= 0 ? 0 : holeAt(t, i); if (th >= .5) fr = i / (ND - 1); col.copy(cBare).lerp(cCov, th); for (let j = 0; j <= NA; j++) { const k = i * (NA + 1) + j; hcol[3 * k] = col.r; hcol[3 * k + 1] = col.g; hcol[3 * k + 2] = col.b; } }
    hg.attributes.color.needsUpdate = true; front.position.y = HH - fr * HH; return fr;
  }
  // ---- labels ----
  const labs = []; const lab = (t, x, y, z, o = {}) => { const s = makeLabel(t, { px: 24, scale: .085, color: 'muted', ...o }); s.position.set(x, y, z); g.add(s); labs.push(s); return s; };
  lab('showerhead, 300 mm · 1 % TMA in N₂, 300 sccm', 0, GAP + .55, 0); lab('wafer: 100 nm holes, AR 100, 56× wall area', 0, -.55, WR + .9); lab('outlet', -CR - .6, GAP / 2, 0); lab('one hole, cut open · 100 nm × 10 µm, drawn 1:200,000', HX, HH + .45, 0, { color: 'ink', weight: 700 }); lab('coverage front', HX + HR + .9, HH, 0, { color: 'rose' }); const frontLab = labs[labs.length - 1];
  st.onTheme = dark => { hemi.intensity = .9 * st.lights; dl.intensity = 1.3 * st.lights; labs.forEach(l => redrawLabel(l, dark)); paintGas(cur); paintWafer(cur); paintHole(cur); st.dirty = true; };
  // ---- readout ----
  const ro = { t: document.getElementById('roT'), bar: document.getElementById('roBar'), flat: document.getElementById('roFlat'), bot: document.getElementById('roBot'), ml: document.getElementById('roML'), hint: document.getElementById('roHint') };
  let cur = -1;
  function setTime(t) {
    paintGas(t); paintWafer(t); const fr = paintHole(t); frontLab.position.y = front.position.y;
    if (ro.t) { ro.t.textContent = t.toFixed(2); ro.bar.style.width = (t / T_END * 100) + '%'; ro.flat.textContent = Math.round(100 * (t <= 0 ? 0 : covAt('thetaFlat', t, 0))) + ' %'; ro.bot.textContent = Math.round(100 * (t <= 0 ? 0 : covAt('thetaBottom', t, 0))) + ' %'; ro.ml.textContent = (t / H.t_supply).toFixed(2);
      ro.hint.textContent = t < .02 ? 'scroll to run the pulse' : t < H.blanket.t90_centre ? 'a blanket wafer would be saturated at 0.35 s…' : t < 2 ? 'the flat surface between the holes saturates (1.95 s); the holes are draining the gas' : t < 13.9 ? `front at ${Math.round(fr * 100)} % of the depth · ${(t / H.t_supply).toFixed(2)} monolayers delivered` : t < 14.3 ? 'hole bottoms reach 90 % at 14.05 s: the rule said 14.4 s' : 'saturated · 1.2 monolayers of supply, as the rule says'; }
    cur = t; st.dirty = true;
  }
  // ---- camera path and layout ----
  const wide = () => st.wide !== false;
  function layout() { g.position.set(wide() ? .6 : 0, wide() ? -.5 : .6, 0); g.scale.setScalar(wide() ? 1 : .68); }
  const pathWide = cameraPath([[8.5, 3.4, 9.5], [9.6, 1.8, 5.0], [10.2, 2.6, 3.8], [9.4, 3.2, 5.2]], [[1.0, .4, 0], [1.8, .3, 0], [3.6, .9, 0], [3.4, 1.1, 0]]);
  const pathNarrow = cameraPath([[0, 4.0, 11.5], [5.0, 2.6, 7.5], [6.5, 3.0, 4.5], [5.5, 3.6, 4.0]], [[0, .8, 0], [.5, .8, 0], [2.5, 1.4, 0], [3.5, 1.6, 0]]);
  const pos3 = new THREE.Vector3(), tgt3 = new THREE.Vector3(); let dragging = false, lastT = -1;
  canvas.addEventListener('pointerdown', () => dragging = true); addEventListener('pointerup', () => setTimeout(() => dragging = false, 2500)); ctr.addEventListener('change', () => st.dirty = true);
  layout(); setTime(reduced ? T_END : 0); (wide() ? pathWide : pathNarrow).at(reduced ? 1 : 0, pos3, tgt3); camera.position.copy(pos3); ctr.target.copy(tgt3); ctr.update();
  window.__hero = { setTime, T_END };
  st.loop((t, dtMs) => {
    const T = S.hero.t; let changed = false; layout();
    const time = reduced ? T_END : clamp((T - .04) / .9, 0, 1) * T_END;        // a short hold at t = 0, then 16 s over the track
    if (Math.abs(time - cur) > 1e-3) { setTime(time); changed = true; }
    (wide() ? pathWide : pathNarrow).at(reduced ? 1 : T, pos3, tgt3);
    if (!dragging && camera.position.distanceToSquared(pos3) > 1e-6) { camera.position.lerp(pos3, reduced ? 1 : .08); changed = true; }
    ctr.target.lerp(tgt3, .08); ctr.autoRotate = !reduced && T < .04 && !dragging; if (ctr.autoRotate || dragging) changed = true; ctr.update();
    st.setBloom(.26 + .18 * S.hero.glow); if (st.bokeh) st.bokeh.uniforms.focus.value = camera.position.distanceTo(ctr.target);
    if (Math.abs(T - lastT) > 1e-4) { lastT = T; changed = true; }
    return changed;
  });
})();
