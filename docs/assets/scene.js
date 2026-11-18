/* Hero scene for the ALD page. Data: window.ALD.hero from assets/data.js: the TMA field in the gap (40 radial x 10 vertical bins at
   29 times, from the OpenFOAM case showerhead_pulse_patterned_AR100_long), the per-face coverage of the flat surface and of the
   hole bottoms against time, and the coverage-against-depth history of one hole at the wafer centre (aldconf, driven by the
   chamber's face concentration). The first half of the scroll scrubs time from 0 to 16 s: gas particles brighten with the local
   TMA density, the wafer top colours with flat coverage, and a cut-open hole beside the chamber fills from the mouth down.
   The second half moves to a PillarHall lateral channel (window.ALD_MECH from assets/data-mech.js, the pillar-cavity Monte Carlo):
   one TMA exposure with the coverage front advancing, purge and water, then hundreds of cycles with the film narrowing the gap.
   In the channel the wall coverage, the film profile, the front position and every readout are data; the molecule paths are an
   illustration (gap drawn 28x taller, flights stretched with it, sticking probability exaggerated). */
import { THREE, createStage, cameraPath, makeLabel, redrawLabel } from './three-kit.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
const D = window.ALD, H = D.hero, S = window.SITE, reduced = S.reduced, phone = S.phone, M = window.ALD_MECH || null;
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

  // ---- the PillarHall channel (steps 3 and 4) ----
  const SX0 = HX + 7.5, LX = 9, G0 = .8, DZ = 2.4, ZMAX = 400, SI = 3;                 // 9 units = 400 gap heights; s = M.s[SI] = 0.0158
  const slit = { on: !!M, mode: 'gas300', E: 0, cycles: 1, fill: 0, z50: 0, tail: 0, phase: 'dose', userCyc: 0, userMode: 0 };
  const sg = new THREE.Group(); sg.position.set(SX0, 0, 0); g.add(sg);
  let paintSlit = () => {}, stepSlit = () => false, slitLabs = [];
  if (M) {
    const NZ = M.z.length, th1 = new Float32Array(NZ), thk = new Float32Array(NZ);
    const single = (key, e) => { const K = M.kernels[key].theta[SI], Es = M.E; if (e <= Es[0]) { const f = clamp(e / Es[0], 0, 1); for (let i = 0; i < NZ; i++) th1[i] = K[0][i] * f; return; } let j = 1; while (j < Es.length - 1 && Es[j] < e) j++; const f = clamp((e - Es[j - 1]) / (Es[j] - Es[j - 1]), 0, 1); for (let i = 0; i < NZ; i++) th1[i] = lerp(K[j - 1][i], K[j][i], f); };
    const samp = (arr, z) => { const u = z / 4; if (u >= NZ - 1) return 0; const i = Math.floor(u); return lerp(arr[i], arr[i + 1], u - i); };
    // channel narrowing: the N-cycle thickness is the single-exposure coverage averaged over the gap stretch f^n, f uniform on [1 - fill, 1]
    const cumulative = (fill, n) => { if (fill <= 1e-4) { thk.set(th1); return; } for (let i = 0; i < NZ; i++) { let a = 0; for (let k = 0; k <= 10; k++) { const f = Math.pow(Math.max(.02, 1 - fill + fill * k / 10), n); a += samp(th1, M.z[i] / f); } thk[i] = a / 11; } const n0 = thk[1] || 1; for (let i = 0; i < NZ; i++) thk[i] = Math.min(1, thk[i] / n0); };
    const cross = (arr, l) => { for (let i = 1; i < NZ; i++) if (arr[i] < l) return Math.max(0, M.z[i - 1] + (arr[i - 1] - l) / ((arr[i - 1] - arr[i]) || 1) * 4); return ZMAX; };
    const plateM = new THREE.MeshStandardMaterial({ color: 0x2b3444, roughness: .8 });
    const base = new THREE.Mesh(new THREE.BoxGeometry(LX + .4, .14, DZ + .3), plateM); base.position.set(LX / 2, -.07, 0); sg.add(base);
    const lid = new THREE.Mesh(new THREE.BoxGeometry(LX + .4, .1, DZ + .3), new THREE.MeshStandardMaterial({ color: 0x8b98a9, transparent: true, opacity: .16, roughness: .4, depthWrite: false })); lid.position.set(LX / 2, G0 + .05, 0); sg.add(lid);
    const endw = new THREE.Mesh(new THREE.BoxGeometry(.1, G0, DZ + .3), plateM); endw.position.set(LX + .15, G0 / 2, 0); sg.add(endw);
    const mk = (opacity) => { const geo = new THREE.BufferGeometry(), pos = new Float32Array(NZ * 2 * 3), colr = new Float32Array(NZ * 2 * 3), idx = []; for (let i = 0; i < NZ; i++) for (let j = 0; j < 2; j++) { const k = i * 2 + j; pos[3 * k] = M.z[i] / ZMAX * LX; pos[3 * k + 2] = (j - .5) * DZ; if (i < NZ - 1 && j === 0) idx.push(k, k + 1, k + 2, k + 1, k + 3, k + 2); }
      geo.setAttribute('position', new THREE.BufferAttribute(pos, 3)); geo.setAttribute('color', new THREE.BufferAttribute(colr, 3)); geo.setIndex(idx); const m = new THREE.Mesh(geo, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: .45, metalness: .25, side: THREE.DoubleSide, transparent: opacity < 1, opacity, depthWrite: opacity >= 1 })); sg.add(m); return { geo, pos, colr }; };
    const fb = mk(1), ft = mk(.5);
    const pilM = new THREE.MeshStandardMaterial({ color: 0x8b98a9, roughness: .5, metalness: .3, transparent: true, opacity: .7 });
    [[1.2, .62], [1.2, -.62], [3.32, 0], [5.44, .62], [5.44, -.62], [7.56, 0]].forEach(([x, z]) => { const c = new THREE.Mesh(new THREE.CylinderGeometry(.1, .1, G0, 16), pilM); c.position.set(x, G0 / 2, z); sg.add(c); });
    const frontS = new THREE.Mesh(new THREE.BoxGeometry(.035, G0 * 1.15, DZ + .1), new THREE.MeshStandardMaterial({ color: 0xfb7185, emissive: 0xfb7185, emissiveIntensity: .9, transparent: true, opacity: .75 })); frontS.position.set(0, G0 / 2, 0); sg.add(frontS);
    const film = x => { const z = clamp(x / LX, 0, 1) * ZMAX; return .004 + slit.fill / 2 * G0 * samp(thk, z); };
    paintSlit = () => {
      const K = M.kernels[slit.mode]; single(slit.mode, slit.E); cumulative(slit.fill, K.gap_exponent);
      cBare.set(st.dark ? 0x5b6676 : 0x8a95a5); cCov.set(P().green);
      for (let i = 0; i < NZ; i++) { const t = thk[i], h = .004 + slit.fill / 2 * G0 * t; col.copy(cBare).lerp(cCov, t); for (let j = 0; j < 2; j++) { const k = i * 2 + j; fb.pos[3 * k + 1] = h; ft.pos[3 * k + 1] = G0 - h; fb.colr[3 * k] = ft.colr[3 * k] = col.r; fb.colr[3 * k + 1] = ft.colr[3 * k + 1] = col.g; fb.colr[3 * k + 2] = ft.colr[3 * k + 2] = col.b; } }
      [fb, ft].forEach(f => { f.geo.attributes.position.needsUpdate = true; f.geo.attributes.color.needsUpdate = true; f.geo.computeVertexNormals(); });
      slit.z50 = cross(thk, .5); slit.tail = samp(thk, slit.z50 + 60); frontS.position.x = slit.z50 / ZMAX * LX; frontS.visible = slit.E > .02;
    };
    // molecules: an illustration, not the Monte Carlo. They enter at the mouth, re-emit with the cosine law and stick where the wall is still bare
    const NT = phone ? 130 : 300, NB = phone ? 90 : 220, V = 4.2, LAM = 1.0;
    const tg = new THREE.SphereGeometry(.035, 8, 6), tmM = new THREE.MeshBasicMaterial({ color: 0xffffff }); const tma = new THREE.InstancedMesh(tg, tmM, NT); sg.add(tma);
    const n2 = new THREE.InstancedMesh(new THREE.SphereGeometry(.02, 6, 5), new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: .7 }), NB); sg.add(n2);
    const r2 = mulberry(11), X = new Float32Array(NT), Y = new Float32Array(NT), Zp = new Float32Array(NT), UX = new Float32Array(NT), UY = new Float32Array(NT), UZ = new Float32Array(NT), ON = new Uint8Array(NT);
    const BX = new Float32Array(NB), BY = new Float32Array(NB), BZ = new Float32Array(NB), BV = new Float32Array(NB * 3);
    for (let i = 0; i < NB; i++) { BX[i] = r2() * LX; BY[i] = .05 + r2() * (G0 - .1); BZ[i] = (r2() - .5) * DZ; for (let k = 0; k < 3; k++) BV[3 * i + k] = (r2() - .5) * 1.6; }
    const cosine = (i, ny) => { const r = r2(), ph = 2 * Math.PI * r2(), ct = Math.sqrt(r), sn = Math.sqrt(1 - r); UY[i] = ny * ct; UX[i] = sn * Math.cos(ph); UZ[i] = sn * Math.sin(ph); };
    const spawn = i => { const xf = slit.z50 / ZMAX * LX; X[i] = r2() < .3 ? 0 : r2() * Math.max(0, xf * .9); Y[i] = film(X[i]) + .03 + r2() * Math.max(.05, G0 - 2 * film(X[i]) - .06); Zp[i] = (r2() - .5) * DZ; const r = r2(), ph = 2 * Math.PI * r2(), ct = Math.sqrt(r), sn = Math.sqrt(1 - r); UX[i] = ct; UY[i] = sn * Math.cos(ph); UZ[i] = sn * Math.sin(ph); ON[i] = 1; };
    const cT = new THREE.Color(), cW = new THREE.Color(), cN = new THREE.Color(), zero = new THREE.Matrix4().makeScale(0, 0, 0);
    for (let i = 0; i < NT; i++) { tma.setMatrixAt(i, zero); tma.setColorAt(i, col.set(0xffffff)); } for (let i = 0; i < NB; i++) { n2.setMatrixAt(i, zero); n2.setColorAt(i, col.set(0xffffff)); }
    stepSlit = (dt) => {
      const flying = slit.phase === 'dose' || slit.phase === 'water' || slit.phase === 'cycles', gas = slit.mode !== 'free', water = slit.phase === 'water';
      cT.set(P().warm); cW.set(P().cyan); cN.set(P().blue);
      for (let i = 0; i < NT; i++) {
        if (!flying) { if (ON[i]) { ON[i] = 0; tma.setMatrixAt(i, zero); } continue; }
        if (!ON[i]) { if (r2() < dt * 2.5) spawn(i); else continue; }
        let left = V * dt, guard = 0;
        while (left > 1e-5 && guard++ < 8) {
          const lo = film(X[i]), hi = G0 - lo; let t = left;
          if (UY[i] > 1e-6) t = Math.min(t, (hi - Y[i]) / UY[i]); else if (UY[i] < -1e-6) t = Math.min(t, (lo - Y[i]) / UY[i]);
          t = Math.max(t, 0); X[i] += UX[i] * t; Y[i] += UY[i] * t; Zp[i] += UZ[i] * t; left -= t;
          if (Zp[i] > DZ / 2) Zp[i] -= DZ; else if (Zp[i] < -DZ / 2) Zp[i] += DZ;
          if (X[i] <= 0) { spawn(i); break; } if (X[i] >= LX) { X[i] = LX - 1e-3; UX[i] = -Math.abs(UX[i]); }
          if (left > 1e-5) {                                              // reached a wall
            const cov = samp(thk, X[i] / LX * ZMAX), bare = water ? cov * .25 : (1 - (slit.fill > 0 ? samp(th1, X[i] / LX * ZMAX / Math.pow(Math.max(.02, 1 - slit.fill), M.kernels[slit.mode].gap_exponent)) : cov));
            if (r2() < .55 * bare) { spawn(i); break; }
            const top = Y[i] > G0 / 2; Y[i] = top ? hi - 1e-3 : lo + 1e-3; cosine(i, top ? -1 : 1);
          }
        }
        if (gas && r2() < V * dt / LAM) { const c = 2 * r2() - 1, sn = Math.sqrt(1 - c * c), ph = 2 * Math.PI * r2(); UX[i] = sn * Math.cos(ph); UY[i] = c; UZ[i] = sn * Math.sin(ph); }
        m4.makeTranslation(X[i], Y[i], Zp[i]); tma.setMatrixAt(i, m4); tma.setColorAt(i, water ? cW : cT);
      }
      for (let i = 0; i < NB; i++) {
        if (!gas) { n2.setMatrixAt(i, zero); continue; }
        BX[i] += BV[3 * i] * dt; BY[i] += BV[3 * i + 1] * dt * .4; BZ[i] += BV[3 * i + 2] * dt;
        if (BX[i] < 0 || BX[i] > LX) BV[3 * i] *= -1; if (BY[i] < .05 || BY[i] > G0 - .05) BV[3 * i + 1] *= -1; if (BZ[i] < -DZ / 2 || BZ[i] > DZ / 2) BV[3 * i + 2] *= -1;
        BX[i] = clamp(BX[i], 0, LX); BY[i] = clamp(BY[i], .05, G0 - .05); BZ[i] = clamp(BZ[i], -DZ / 2, DZ / 2);
        m4.makeTranslation(BX[i], BY[i], BZ[i]); n2.setMatrixAt(i, m4); n2.setColorAt(i, cN);
      }
      tma.instanceMatrix.needsUpdate = true; n2.instanceMatrix.needsUpdate = true; if (tma.instanceColor) tma.instanceColor.needsUpdate = true; if (n2.instanceColor) n2.instanceColor.needsUpdate = true;
      return true;
    };
    const sl = (t, x, y, z, o = {}) => { const sp = makeLabel(t, { px: 24, scale: .085, color: 'muted', ...o }); sp.position.set(x, y, z); sg.add(sp); slitLabs.push(sp); return sp; };
    sl('PillarHall lateral channel · gap 500 nm, drawn 28× taller', LX * .3, G0 + .6, -DZ / 2, { color: 'ink', weight: 700 }); sl('mouth', -.55, G0 / 2, 0); sl('support pillars, 4 µm at 49 µm', 3.32, -.42, DZ / 2 + .5); sl('closed end', LX + .8, G0 / 2, 0);
    paintSlit();
  }
  // ---- labels ----
  const labs = []; const lab = (t, x, y, z, o = {}) => { const s = makeLabel(t, { px: 24, scale: .085, color: 'muted', ...o }); s.position.set(x, y, z); g.add(s); labs.push(s); return s; };
  lab('showerhead, 300 mm · 1 % TMA in N₂, 300 sccm', 0, GAP + .55, 0); lab('wafer: 100 nm holes, AR 100, 56× wall area', 0, -.55, WR + .9); lab('outlet', -CR - .6, GAP / 2, 0); lab('one hole, cut open · 100 nm × 10 µm, drawn 1:200,000', HX, HH + .45, 0, { color: 'ink', weight: 700 }); lab('coverage front', HX + HR + .9, HH, 0, { color: 'rose' }); const frontLab = labs[labs.length - 1];
  st.onTheme = dark => { hemi.intensity = .9 * st.lights; dl.intensity = 1.3 * st.lights; labs.concat(slitLabs).forEach(l => redrawLabel(l, dark)); paintGas(cur); paintWafer(cur); paintHole(cur); paintSlit(); st.dirty = true; };
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
  // ---- the story: chamber and hole for the first half of the track, the channel for the second ----
  const SPLIT = slit.on ? .46 : 1, E_SHOW = 4.5, CYC_MAX = 1000;
  const el = id => document.getElementById(id), roC = el('roChamber'), roS = el('roSlit'), stepEls = [...document.querySelectorAll('#heroSteps li')], capEl = el('heroCap');
  const rs = { E: el('rsE'), z: el('rsZ'), tail: el('rsTail'), cyc: el('rsCyc'), fill: el('rsFill'), hint: el('rsHint'), range: el('rsRange'), gasOn: el('rsGasOn'), gasOff: el('rsGasOff') };
  const setMode = m => { slit.mode = m; if (rs.gasOn) { rs.gasOn.classList.toggle('on', m !== 'free'); rs.gasOff.classList.toggle('on', m === 'free'); rs.gasOn.setAttribute('aria-pressed', m !== 'free'); rs.gasOff.setAttribute('aria-pressed', m === 'free'); } slitKey = ''; };
  if (rs.gasOn) { rs.gasOn.addEventListener('click', () => setMode('gas300')); rs.gasOff.addEventListener('click', () => setMode('free')); rs.range.addEventListener('input', () => { slit.userCyc = performance.now(); slit.cycles = +rs.range.value; slitKey = ''; }); }
  let slitKey = '', stepNow = -1;
  function slitState(v) {                                              // v: 0 to 1 over the channel half of the track
    const e = clamp((v - .12) / .40, 0, 1) * E_SHOW; let phase = v < .54 ? 'dose' : v < .61 ? 'purge' : v < .70 ? 'water' : 'cycles';
    const user = performance.now() - slit.userCyc < 6000 && slit.userCyc > 0;
    const cyc = user ? slit.cycles : phase === 'cycles' ? Math.max(1, Math.round(clamp((v - .70) / .27, 0, 1) * CYC_MAX)) : 1; if (user && cyc > 1) phase = 'cycles';
    slit.E = reduced ? E_SHOW : e; slit.cycles = reduced ? CYC_MAX : cyc; slit.phase = reduced ? 'idle' : phase; slit.fill = slit.cycles <= 1 ? 0 : Math.min(.95, 2 * slit.cycles * M.gpc_nm / M.gap_nm);
    const key = slit.mode + '|' + slit.E.toFixed(3) + '|' + slit.cycles + '|' + st.dark;
    if (key !== slitKey) { slitKey = key; paintSlit();
      if (rs.E) { rs.E.textContent = slit.E.toFixed(2); rs.z.textContent = slit.E > .02 ? Math.round(slit.z50) + ' gaps · ' + Math.round(slit.z50 * M.gap_nm / 1000) + ' µm' : '0'; rs.tail.textContent = slit.E > .02 ? slit.tail.toFixed(3) : '0'; rs.cyc.textContent = slit.cycles; rs.fill.textContent = Math.round(slit.fill * 100) + ' %'; if (!user) rs.range.value = slit.cycles;
        rs.hint.textContent = phase === 'dose' ? (slit.mode === 'free' ? 'no carrier gas: long flights carry TMA far ahead of the front' : 'N₂ at 300 Pa: flights end after about 39 gaps, the tail ahead of the front is four times smaller') : phase === 'purge' ? 'purge: the channel empties' : phase === 'water' ? 'water pulse: the surface is ready for the next TMA dose' : 'each cycle leaves film on both walls; the narrower gap lets later cycles penetrate less, so the measured front smears'; }
      return true; }
    return false;
  }
  function setStep(k) { if (k === stepNow) return; stepNow = k; stepEls.forEach((li, i) => li.classList.toggle('on', i === k)); if (roC && roS) { roC.hidden = k >= 2; roS.hidden = k < 2; }
    if (capEl) capEl.textContent = k < 2 ? 'Gas field, wafer and hole coverage: OpenFOAM chamber run with the feature model on every wafer face (assets/data.js).' : 'Wall coverage, film profile, front and readouts: pillar-cavity Monte Carlo at s = 0.0158 (assets/data-mech.js). Molecule paths are an illustration: gap drawn 28× taller, sticking exaggerated.'; }
  const slitC = new THREE.Vector3(), off = new THREE.Vector3();
  function slitCam(v, outPos, outTgt) { const w = wide(), a = clamp(v / .7, 0, 1); slitC.set(SX0 + LX * lerp(.36, .5, a), G0 / 2, 0); g.localToWorld(slitC); outTgt.copy(slitC);
    if (w) off.set(lerp(-2.2, .4, a), lerp(2.6, 1.5, a), lerp(8.6, 6.4, a)); else off.set(lerp(-1, 0, a), lerp(5.5, 4.4, a), lerp(23, 20.5, a)); outPos.copy(slitC).addScaledVector(off, g.scale.x); }
  layout(); g.updateMatrixWorld(true); setTime(reduced && !slit.on ? T_END : 0); (wide() ? pathWide : pathNarrow).at(0, pos3, tgt3);
  if (reduced && slit.on) { slitState(1); setStep(3); slitCam(1, pos3, tgt3); setTime(T_END); } else setStep(0);
  camera.position.copy(pos3); ctr.target.copy(tgt3); ctr.update();
  window.__hero = { setTime, T_END, slit, SPLIT };
  st.loop((t, dtMs) => {
    const T = S.hero.t; let changed = false; layout(); g.updateMatrixWorld(true);
    const u = clamp(T / SPLIT, 0, 1), v = slit.on ? clamp((T - SPLIT) / (1 - SPLIT), 0, 1) : 0, inSlit = slit.on && (T > SPLIT || reduced);
    const time = reduced ? T_END : clamp((u - .05) / .9, 0, 1) * T_END;      // a short hold at t = 0, then 16 s over the chamber half
    if (Math.abs(time - cur) > 1e-3) { setTime(time); changed = true; }
    if (inSlit) { if (slitState(reduced ? 1 : v)) changed = true; if (!reduced && stepSlit(Math.min(.05, dtMs / 1000))) changed = true; }
    setStep(!inSlit ? (time < 2 ? 0 : 1) : slit.phase === 'cycles' || reduced ? 3 : 2);
    if (inSlit) slitCam(reduced ? 1 : v, pos3, tgt3); else (wide() ? pathWide : pathNarrow).at(u, pos3, tgt3);
    const kf = reduced ? 1 : 1 - Math.exp(-Math.min(200, dtMs) / 170);
    if (!dragging && camera.position.distanceToSquared(pos3) > 1e-6) { camera.position.lerp(pos3, kf); changed = true; }
    ctr.target.lerp(tgt3, kf); ctr.autoRotate = !reduced && T < .02 && !dragging; if (ctr.autoRotate || dragging) changed = true; ctr.update();
    st.setBloom(.26 + .18 * S.hero.glow); if (st.bokeh) st.bokeh.uniforms.focus.value = camera.position.distanceTo(ctr.target);
    if (Math.abs(T - lastT) > 1e-4) { lastT = T; changed = true; }
    return changed;
  });
})();
