/* Body scenes: every graph, schematic and process diagram on the page is a WebGL scene drawn by ONE renderer into
   a fixed, transparent, pointer-transparent canvas, scissor-clipped to the rectangle of its placeholder element.
   This is the three.js "webgl_multiple_elements" pattern (one context for any number of views, no per-view canvas,
   no browser context cap): https://github.com/mrdoob/three.js/blob/master/examples/webgl_multiple_elements.html
   Each view assembles or draws progressively as it scrolls into the viewport (progress 0 -> 1), renders only while on
   screen and only when something changed, supports hover readouts (raycast), drag-to-orbit, and both themes.
   Reduced motion or no WebGL: the views are not built and each placeholder shows its static PNG (docs/build_assets.py).
   Helpers: axes, ribbon, bars, surface, tiles, exploded, label, mat. ES module; three.js via the page's import map. */
import { THREE, palette, THEME, roomEnvironment, makeLabel, redrawLabel } from './three-kit.js';

const S = () => window.SITE || { reduced: matchMedia('(prefers-reduced-motion:reduce)').matches, phone: matchMedia('(max-width:640px)').matches, isDark: () => matchMedia('(prefers-color-scheme:dark)').matches, onFrame: cb => requestAnimationFrame(function l(t) { cb(t / 1000, 16); requestAnimationFrame(l) }) };
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const ease = t => t < .5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
const NAV_H = 52;

export const Viz = {
  views: [], renderer: null, canvas: null, ready: false, static: false, dark: true, pal: palette(true),

  /* Build the shared renderer; called once by the first view(). Returns false in static mode. */
  init() {
    if (this.renderer || this.static) return !this.static;
    const s = S(); this.S = s; this.dark = s.isDark(); this.pal = palette(this.dark);
    const toStatic = () => { this.static = true; document.documentElement.classList.add('viz-static'); this.applyStatic(); return false; };
    if (s.reduced || (window.VIZ_FORCE_STATIC)) return toStatic();
    try { const p = document.createElement('canvas'); if (!(p.getContext('webgl2') || p.getContext('webgl'))) return toStatic(); } catch (e) { return toStatic(); }
    const canvas = document.createElement('canvas'); canvas.className = 'viz-canvas'; canvas.setAttribute('aria-hidden', 'true'); document.body.appendChild(canvas); this.canvas = canvas;
    let r; try { r = new THREE.WebGLRenderer({ canvas, alpha: true, antialias: true, powerPreference: 'high-performance' }); } catch (e) { canvas.remove(); return toStatic(); }
    this.renderer = r; r.setPixelRatio(Math.min(devicePixelRatio || 1, s.phone ? 1.5 : 2)); r.setClearColor(0x000000, 0); r.toneMapping = THREE.ACESFilmicToneMapping; r.autoClear = false;
    this.env = roomEnvironment(r);
    const tip = document.createElement('div'); tip.className = 'viz-tip'; tip.hidden = true; document.body.appendChild(tip); this.tip = tip;
    const fit = () => { const w = innerWidth, h = innerHeight; r.setSize(w, h, false); this.w = w; this.h = h; this.views.forEach(v => v.dirty = true); };
    fit(); addEventListener('resize', fit);
    this.applyTheme(); document.addEventListener('themechange', () => this.applyTheme());
    let lastY = -1, lastW = 0;
    s.onFrame((t, dt) => {
      if (this.disposed) return; const scrolled = scrollY !== lastY || innerWidth !== lastW; lastY = scrollY; lastW = innerWidth;
      let any = false, anyVisible = false; const vh = innerHeight, vw = innerWidth;
      for (const v of this.views) {
        const rc = v.el.getBoundingClientRect(); v.rect = rc;
        v.visible = rc.bottom > NAV_H && rc.top < vh && rc.right > 0 && rc.left < vw && rc.width > 2 && rc.height > 2;
        if (!v.visible) continue; anyVisible = true;
        // scroll progress: starts when the top enters at 92 % of the viewport, completes when the centre reaches 58 %
        const start = vh * .92, end = Math.min(vh * .58 - rc.height / 2, start - vh * .25);
        const target = v.pinned != null ? v.pinned : clamp((start - rc.top) / Math.max(1, start - end), 0, 1);
        if (Math.abs(target - v.p) > 1e-3) { v.p += (target - v.p) * (v.pinned != null ? 1 : .14); if (Math.abs(target - v.p) < 1e-3) v.p = target; v.dirty = true; }
        const changed = v.update ? v.update(v, t, dt) : false; if (changed || v.dirty || scrolled) any = true;
      }
      if (any && anyVisible) this.renderAll();
    });
    this.ready = true; return true;
  },

  applyTheme() {
    const dark = this.S.isDark(), T = THEME[dark ? 'dark' : 'light']; this.dark = dark; this.pal = palette(dark);
    if (this.renderer) this.renderer.toneMappingExposure = T.exposure;
    for (const v of this.views) {
      v.dark = dark; v.pal = this.pal; if (v.scene.environment) v.scene.environmentIntensity = v.envIntensity * T.env;
      v.lights.forEach(l => l.intensity = l.userData.base * T.lights);
      v.scene.traverse(o => { if (o.userData.label) redrawLabel(o, dark); if (o.material && o.material.userData && o.material.userData.pal) { const m = o.material, k = m.userData.pal; m.color.set(this.pal[k] || k); if (m.emissive && m.userData.emissive) m.emissive.set(this.pal[m.userData.emissive] || m.userData.emissive); } });
      v.onTheme && v.onTheme(v, dark); v.dirty = true;
    }
    if (this.static) this.applyStatic();
  },
  applyStatic() { const t = this.S ? (this.S.isDark() ? 'dark' : 'light') : 'dark'; document.querySelectorAll('.viz img.fallback[data-viz-fig]').forEach(i => { i.src = (i.dataset.dir || 'assets/') + i.dataset.vizFig + '_' + t + '.png'; }); },

  renderAll() {
    const r = this.renderer, h = this.h; r.setScissorTest(false); r.clear();
    r.setScissorTest(true);
    for (const v of this.views) {
      if (!v.visible) continue; const rc = v.rect;
      const top = Math.max(rc.top, NAV_H), bottom = Math.min(rc.bottom, h), left = Math.max(rc.left, 0), right = Math.min(rc.right, this.w);
      if (bottom - top < 2 || right - left < 2) continue;
      // the camera keeps the whole element's aspect; the scissor clips to the visible part
      const asp = rc.width / rc.height; if (Math.abs(v.camera.aspect - asp) > 1e-3) { v.camera.aspect = asp; v.camera.updateProjectionMatrix(); }
      r.setViewport(rc.left, h - rc.bottom, rc.width, rc.height); r.setScissor(left, h - bottom, right - left, bottom - top);
      r.render(v.scene, v.camera); v.dirty = false;
    }
    r.setScissorTest(false);
  },

  /* Register a view on element el. build(v) adds objects to v.scene and may set v.update(v,t,dt) -> changed, v.hover(hit) -> html, v.onTheme(v,dark). */
  view(el, build, o = {}) {
    if (typeof el === 'string') el = document.querySelector(el); if (!el) return null;
    if (!this.init()) return null;
    const scene = new THREE.Scene(); scene.environment = this.env; const envIntensity = o.envIntensity ?? .5; scene.environmentIntensity = envIntensity * THEME[this.dark ? 'dark' : 'light'].env;
    const camera = new THREE.PerspectiveCamera(o.fov || 34, 1.6, .05, 200); camera.position.set(...(o.camera || [4, 3, 6])); const target = new THREE.Vector3(...(o.target || [0, 0, 0])); camera.lookAt(target);
    const rig = new THREE.Group(); scene.add(rig);
    const v = { el, scene, camera, rig, target, p: 0, pinned: null, dirty: true, visible: false, dark: this.dark, pal: this.pal, lights: [], pickables: [], envIntensity, name: el.dataset.viz || '', spin: o.spin || 0, orbit: o.orbit !== false, rot: 0, tilt: 0 };
    const hemi = new THREE.HemisphereLight(0xffffff, 0x334466, o.hemi ?? .9), dl = new THREE.DirectionalLight(0xffffff, o.sun ?? 1.3); dl.position.set(4, 8, 6); hemi.userData.base = hemi.intensity; dl.userData.base = dl.intensity; scene.add(hemi, dl); v.lights.push(hemi, dl);
    v.lights.forEach(l => l.intensity = l.userData.base * THEME[this.dark ? 'dark' : 'light'].lights);
    el.classList.add('viz-live');
    // pointer: drag to orbit the rig, hover readouts through a raycaster on v.pickables
    const ray = new THREE.Raycaster(), ndc = new THREE.Vector2(); let drag = null;
    const ndcOf = e => { const rc = v.rect || el.getBoundingClientRect(); ndc.set(((e.clientX - rc.left) / rc.width) * 2 - 1, -((e.clientY - rc.top) / rc.height) * 2 + 1); return ndc; };
    el.addEventListener('pointerdown', e => { if (!v.orbit) return; drag = { x: e.clientX, y: e.clientY, rot: v.rot, tilt: v.tilt, moved: false }; el.setPointerCapture && el.setPointerCapture(e.pointerId); });
    addEventListener('pointerup', () => drag = null);
    el.addEventListener('pointermove', e => {
      if (drag) { const dx = e.clientX - drag.x, dy = e.clientY - drag.y; if (Math.abs(dx) + Math.abs(dy) > 3) drag.moved = true; v.rot = drag.rot + dx * .008; v.tilt = clamp(drag.tilt + dy * .005, -.6, .6); v.dirty = true; this.tip.hidden = true; return; }
      if (!v.hover || !v.pickables.length || this.S.phone) return;
      ray.setFromCamera(ndcOf(e), camera); const hits = ray.intersectObjects(v.pickables, true); const html = hits.length ? v.hover(hits[0], v) : null;
      if (html) { this.tip.innerHTML = html; this.tip.hidden = false; const tw = this.tip.offsetWidth; this.tip.style.left = Math.min(innerWidth - tw - 8, e.clientX + 14) + 'px'; this.tip.style.top = (e.clientY + 14) + 'px'; } else this.tip.hidden = true;
    });
    el.addEventListener('pointerleave', () => { this.tip.hidden = true; });
    this.views.push(v);
    build(v, this);
    this.halos(v);
    const userUpdate = v.update;
    v.update = (vv, t, dt) => { let ch = false; if (vv.halo) this.syncHalos(vv); if (vv.fx) for (const fx of vv.fx) ch = fx.step(Math.min(.05, dt / 1000)) || ch; if (vv.spin && !drag) { vv.rot += vv.spin * (dt / 1000); ch = true; } rig.rotation.y = vv.rot; rig.rotation.x = vv.tilt; if (userUpdate) ch = userUpdate(vv, t, dt) || ch; return ch; };
    return v;
  },

  /* Force every view to its final state and render (used by build_assets.py's snapshot pass and by QA). */
  force(p = 1) { for (const v of this.views) { v.pinned = p; v.p = p; v.dirty = true; } },

  /* ---------- materials and labels ---------- */
  mat(key, o = {}) { const pal = this.pal; const m = new THREE.MeshStandardMaterial({ color: pal[key] || key, roughness: o.roughness ?? .45, metalness: o.metalness ?? .15, transparent: !!o.opacity && o.opacity < 1, opacity: o.opacity ?? 1, side: o.side || THREE.FrontSide, emissive: o.emissive ? (pal[o.emissive] || o.emissive) : 0x000000, emissiveIntensity: o.emissiveIntensity ?? .25, flatShading: !!o.flat, depthWrite: o.depthWrite ?? true }); m.userData.pal = key; if (o.emissive) m.userData.emissive = o.emissive; return m; },
  label(v, text, pos, o = {}) { const sp = makeLabel(text, { ...o, dark: v.dark }); sp.position.set(...pos); (o.parent || v.rig).add(sp); return sp; },
  line(v, pts, key = 'line', o = {}) { const g = new THREE.BufferGeometry().setFromPoints(pts.map(p => new THREE.Vector3(...p))); const m = new THREE.LineBasicMaterial({ color: this.pal[key] || key, transparent: o.opacity != null, opacity: o.opacity ?? 1 }); m.userData.pal = key; const l = new THREE.Line(g, m); (o.parent || v.rig).add(l); return l; },

  /* ---------- glow and sparks (the body scenes have no bloom pass, so glow is an additive sprite with a soft radial falloff) ---------- */
  glowTex() {
    if (this._glow) return this._glow; const c = document.createElement('canvas'); c.width = c.height = 64; const x = c.getContext('2d'), g = x.createRadialGradient(32, 32, 0, 32, 32, 32);
    g.addColorStop(0, 'rgba(255,255,255,1)'); g.addColorStop(.25, 'rgba(255,255,255,.55)'); g.addColorStop(.6, 'rgba(255,255,255,.12)'); g.addColorStop(1, 'rgba(255,255,255,0)'); x.fillStyle = g; x.fillRect(0, 0, 64, 64);
    this._glow = new THREE.CanvasTexture(c); return this._glow;
  },
  /* A soft halo behind every emissive data marker (sphere) of a view, dark theme only. */
  halos(v) {
    const list = []; v.rig.traverse(o => { if (o.isMesh && o.geometry && o.geometry.type === 'SphereGeometry' && o.material && o.material.userData && o.material.userData.emissive) list.push(o); });
    if (!list.length || list.length > 400) return;
    v.halo = list.map(m => { const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: this.glowTex(), transparent: true, depthWrite: false, blending: THREE.AdditiveBlending, opacity: 0 })); const r = m.geometry.parameters.radius || .06; sp.scale.setScalar(r * 6.5); sp.raycast = () => {}; sp.renderOrder = 3; m.add(sp); return { m, sp }; });
    this.syncHalos(v);
  },
  syncHalos(v) { for (const { m, sp } of v.halo) { sp.visible = v.dark; if (!v.dark) continue; sp.material.color.copy(m.material.emissive); sp.material.opacity = .5 * (m.material.transparent ? m.material.opacity : 1); } },
  /* Streaks that travel along a path and spark where they land. o.sample(R) returns null or { path: [[x,y,z], ...] (the last point is where
     it lands; one point = no streak), p: probability of a spark, key: palette key } in the coordinates of o.parent (default v.rig).
     o.active() gates the effect; o.overlay draws it over the scene (for gas inside a closed part). A spark is a glow that swells and fades plus two embers. Returns null under reduced motion. */
  fx(v, o = {}) {
    const S = this.S; if (!S || S.reduced || this.static) return null;
    const n = Math.max(4, Math.round((o.n ?? 26) * (S.phone ? .5 : 1))), NS = Math.max(8, Math.round((o.sparks ?? 36) * (S.phone ? .6 : 1))), size = o.size ?? 1, speed = o.speed ?? 2.2, gap = o.gap ?? .9, active = o.active || (() => v.p > .95);
    const G = new THREE.Group(); (o.parent || v.rig).add(G); const zero = new THREE.Matrix4().makeScale(0, 0, 0), q = new THREE.Quaternion(), up = new THREE.Vector3(0, 1, 0), a = new THREE.Vector3(), b = new THREE.Vector3(), d = new THREE.Vector3(), sc = new THREE.Vector3(), M = new THREE.Matrix4(), col = new THREE.Color();
    let seed = 97 + this.views.length * 31; const R = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
    const dt0 = !o.overlay, sm = new THREE.MeshBasicMaterial({ transparent: true, depthWrite: false, depthTest: dt0 }), streak = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 5, 1), sm, n); streak.frustumCulled = false; streak.renderOrder = 5; for (let i = 0; i < n; i++) { streak.setMatrixAt(i, zero); streak.setColorAt(i, col.set(0xffffff)); } G.add(streak);
    const mk = () => { const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: this.glowTex(), transparent: true, depthWrite: false, depthTest: dt0, opacity: 0 })); sp.visible = false; sp.raycast = () => {}; sp.renderOrder = 6; G.add(sp); return sp; };
    const SK = Array.from({ length: NS }, () => ({ sp: mk(), life: 0 })), EM = Array.from({ length: NS * 2 }, () => ({ sp: mk(), life: 0, v: new THREE.Vector3() })), ST = Array.from({ length: n }, () => ({ live: false, wait: R() * gap * 1.5, path: null, seg: 0, u: 0, p: 0, key: 'warm' }));
    let sk = 0, em = 0, dark = null, shown = false;
    const paint = () => { dark = v.dark; const bl = dark ? THREE.AdditiveBlending : THREE.NormalBlending; sm.blending = bl; sm.opacity = dark ? .95 : .8; sm.needsUpdate = true; SK.concat(EM).forEach(e => { e.sp.material.blending = bl; e.sp.material.needsUpdate = true; }); };
    const emit = (pt, key) => { const s = SK[sk++ % NS]; s.life = 1; s.sp.position.set(pt[0], pt[1], pt[2]); s.sp.material.color.set(v.pal[key] || key); s.sp.visible = true;
      for (let j = 0; j < 2; j++) { const e = EM[em++ % (NS * 2)]; e.life = .8 + .2 * R(); e.sp.position.copy(s.sp.position); e.v.set(R() - .5, R() - .2, R() - .5).normalize().multiplyScalar(size * (.35 + .4 * R())); e.sp.material.color.copy(s.sp.material.color); e.sp.visible = true; } };
    const hide = () => { for (let i = 0; i < n; i++) { streak.setMatrixAt(i, zero); ST[i].live = false; } streak.instanceMatrix.needsUpdate = true; SK.concat(EM).forEach(e => { e.life = 0; e.sp.visible = false; }); shown = false; };
    const step = dt => {
      if (!active()) { if (shown) { hide(); return true; } return false; }
      if (dark !== v.dark) paint(); shown = true;
      for (let i = 0; i < n; i++) { const s = ST[i];
        if (!s.live) { streak.setMatrixAt(i, zero); s.wait -= dt; if (s.wait > 0) continue; const r = o.sample(R); if (!r || !r.path || !r.path.length) { s.wait = .3 + R() * gap; continue; }
          if (r.path.length < 2) { if (R() < r.p) emit(r.path[0], r.key || 'warm'); s.wait = .15 + R() * gap; continue; }
          s.live = true; s.path = r.path; s.seg = 0; s.u = 0; s.p = r.p; s.key = r.key || 'warm'; streak.setColorAt(i, col.set(v.pal[s.key] || s.key)); streak.instanceColor.needsUpdate = true; }
        a.fromArray(s.path[s.seg]); b.fromArray(s.path[s.seg + 1]); d.subVectors(b, a); const len = d.length() || 1e-6; s.u += speed * size * dt / len;
        if (s.u >= 1) { s.seg++; s.u = 0; if (s.seg >= s.path.length - 1) { if (R() < s.p) emit(s.path[s.path.length - 1], s.key); s.live = false; s.wait = R() * gap; streak.setMatrixAt(i, zero); } continue; }
        const L = Math.min(.16 * size, len * s.u + .02); d.multiplyScalar(1 / len); q.setFromUnitVectors(up, d); a.addScaledVector(d, len * s.u - L / 2); sc.set(.011 * size, L, .011 * size); M.compose(a, q, sc); streak.setMatrixAt(i, M); }
      streak.instanceMatrix.needsUpdate = true;
      for (const s of SK) { if (s.life <= 0) continue; s.life -= dt * 2.3; if (s.life <= 0) { s.sp.visible = false; continue; } const k = s.life, w = size * (dark ? .24 + .6 * (1 - k) : .1 + .1 * k); s.sp.scale.set(w, w, 1); s.sp.material.opacity = dark ? Math.pow(k, 1.3) : .85 * k; }
      for (const e of EM) { if (e.life <= 0) continue; e.life -= dt * 1.8; if (e.life <= 0) { e.sp.visible = false; continue; } e.sp.position.addScaledVector(e.v, dt); e.v.multiplyScalar(Math.max(0, 1 - 2.4 * dt)); const w = size * (dark ? .1 : .05) * e.life; e.sp.scale.set(w, w, 1); e.sp.material.opacity = (dark ? .9 : .8) * e.life; }
      return true;
    };
    const fx = { step, emit, group: G }; (v.fx = v.fx || []).push(fx); return fx;
  },

  /* ---------- 3D axes box: maps data ranges to a world box of size [sx,sy,sz] centred on the origin ---------- */
  axes(v, spec) {
    const [sx, sy, sz] = spec.size || [4, 2.2, 2.4]; const ax = { x: spec.x, y: spec.y, z: spec.z || { lo: 0, hi: 1 } };
    const map = (x, y, z = 0) => new THREE.Vector3((x - ax.x.lo) / (ax.x.hi - ax.x.lo) * sx - sx / 2, (y - ax.y.lo) / (ax.y.hi - ax.y.lo) * sy - sy / 2, (z - ax.z.lo) / (ax.z.hi - ax.z.lo) * sz - sz / 2);
    const g = new THREE.Group(); v.rig.add(g); const pts = [];
    const x0 = -sx / 2, y0 = -sy / 2, z0 = -sz / 2, z1 = sz / 2, x1 = sx / 2, y1 = sy / 2;
    // floor grid and back wall grid, in the theme line colour
    const nx = (ax.x.ticks || []).length, ny = (ax.y.ticks || []).length;
    for (const tx of ax.x.ticks || []) { const px = map(tx, ax.y.lo, ax.z.lo).x; pts.push([px, y0, z0], [px, y0, z1], [px, y0, z0], [px, y1, z0]); }
    for (const ty of ax.y.ticks || []) { const py = map(ax.x.lo, ty, ax.z.lo).y; pts.push([x0, py, z0], [x1, py, z0], [x0, py, z0], [x0, py, z1]); }
    for (const tz of ax.z.ticks || []) { const pz = map(ax.x.lo, ax.y.lo, tz).z; pts.push([x0, y0, pz], [x1, y0, pz], [x0, y0, pz], [x0, y1, pz]); }
    pts.push([x0, y0, z0], [x1, y0, z0], [x0, y0, z0], [x0, y1, z0], [x0, y0, z0], [x0, y0, z1], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1], [x1, y0, z1], [x0, y1, z0], [x1, y1, z0], [x0, y1, z0], [x0, y1, z1]);
    const geo = new THREE.BufferGeometry().setFromPoints(pts.map(p => new THREE.Vector3(...p))); const lm = new THREE.LineBasicMaterial({ color: this.pal.line, transparent: true, opacity: .9 }); lm.userData.pal = 'line'; g.add(new THREE.LineSegments(geo, lm));
    const px = spec.px || 22, lab = (t, p, o) => { const s = makeLabel(t, { px, bg: false, color: 'muted', dark: v.dark, scale: spec.labelScale || .09, ...o }); s.position.copy(p); g.add(s); return s; };
    const fx = ax.x.fmt || (t => t), fy = ax.y.fmt || (t => t), fz = ax.z.fmt || (t => t);
    for (const tx of ax.x.ticks || []) lab(fx(tx), map(tx, ax.y.lo, ax.z.lo).add(new THREE.Vector3(0, -.14, .18 + sz * .02)));
    for (const ty of ax.y.ticks || []) lab(fy(ty), map(ax.x.lo, ty, ax.z.lo).add(new THREE.Vector3(-.22 - sx * .02, 0, 0)));
    for (const tz of ax.z.ticks || []) lab(fz(tz), map(ax.x.lo, ax.y.lo, tz).add(new THREE.Vector3(.2 + sx * .03, -.1, 0)), { color: 'ink' });
    if (ax.x.label) lab(ax.x.label, new THREE.Vector3(0, y0 - .42, z1 + .12), { color: 'ink', px: px + 2, weight: 700 });
    if (ax.y.label) lab(ax.y.label, new THREE.Vector3(x0 - .55 - sx * .03, 0, z0), { color: 'ink', px: px + 2, weight: 700 });
    if (ax.z.label) lab(ax.z.label, new THREE.Vector3(x1 + .45 + sx * .04, y0 - .1, 0), { color: 'ink', px: px + 2, weight: 700 });
    return { map, group: g, size: [sx, sy, sz], ax };
  },

  /* ---------- a data series as a 3D ribbon: a band of width w (in z) following (x, y), drawn progressively ---------- */
  ribbon(v, A, xs, ys, o = {}) {
    const z = o.z ?? 0, w = o.width ?? .12, n = xs.length, key = o.color || 'cyan';
    const pos = new Float32Array(n * 2 * 3), idx = []; const p = new THREE.Vector3();
    for (let i = 0; i < n; i++) { A.map(xs[i], ys[i], z).toArray(pos, i * 6); pos[i * 6 + 2] -= w / 2; pos[i * 6 + 3] = pos[i * 6]; pos[i * 6 + 4] = pos[i * 6 + 1]; pos[i * 6 + 5] = pos[i * 6 + 2] + w; if (i) idx.push(i * 2 - 2, i * 2 - 1, i * 2, i * 2 - 1, i * 2 + 1, i * 2); }
    const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.BufferAttribute(pos, 3)); g.setIndex(idx); g.computeVertexNormals();
    const m = this.mat(key, { side: THREE.DoubleSide, roughness: .35, metalness: .25, emissive: key, emissiveIntensity: .35 }); const mesh = new THREE.Mesh(g, m); mesh.userData = { xs, ys, series: o.name || '', unit: o.unit || '' }; v.rig.add(mesh); v.pickables.push(mesh);
    let skirt = null;
    if (o.skirt !== false) { const sp = new Float32Array(n * 2 * 3), si = []; for (let i = 0; i < n; i++) { const q = A.map(xs[i], ys[i], z); q.toArray(sp, i * 6); const b = A.map(xs[i], A.ax.y.lo, z); b.toArray(sp, i * 6 + 3); if (i) si.push(i * 2 - 2, i * 2 - 1, i * 2, i * 2 - 1, i * 2 + 1, i * 2); } const sg = new THREE.BufferGeometry(); sg.setAttribute('position', new THREE.BufferAttribute(sp, 3)); sg.setIndex(si); const sm = this.mat(key, { side: THREE.DoubleSide, opacity: o.skirtOpacity ?? .14, depthWrite: false }); skirt = new THREE.Mesh(sg, sm); v.rig.add(skirt); }
    const draw = p => { const k = Math.floor(clamp(p, 0, 1) * (n - 1)); g.setDrawRange(0, k * 6); if (skirt) skirt.geometry.setDrawRange(0, k * 6); }; draw(0);
    return { mesh, skirt, draw, n };
  },
  /* nearest data index for a hover hit on a ribbon */
  ribbonIndex(hit, A) { const xs = hit.object.userData.xs; const lx = (hit.point.x + A.size[0] / 2) / A.size[0] * (A.ax.x.hi - A.ax.x.lo) + A.ax.x.lo; let best = 0, bd = 1e30; for (let i = 0; i < xs.length; i++) { const d = Math.abs(xs[i] - lx); if (d < bd) { bd = d; best = i; } } return best; },

  /* ---------- bar field: instanced boxes at (x, z) data positions, heights y, growing with progress ---------- */
  bars(v, A, items, o = {}) {
    const n = items.length, w = o.width ?? .3, d = o.depth ?? w; const geo = new THREE.BoxGeometry(1, 1, 1); geo.translate(0, .5, 0);
    const m = new THREE.MeshStandardMaterial({ roughness: .4, metalness: .2, vertexColors: false }); const im = new THREE.InstancedMesh(geo, m, n); im.userData.items = items; im.userData.pal = true;
    const M = new THREE.Matrix4(), col = new THREE.Color(); v.rig.add(im); v.pickables.push(im);
    const paint = () => { items.forEach((it, i) => { col.set(this.pal[it.color] || it.color || this.pal.cyan); im.setColorAt(i, col); }); im.instanceColor.needsUpdate = true; };
    const draw = p => { items.forEach((it, i) => { const pi = ease(clamp(p * (1 + .5) - (i / n) * .5, 0, 1)); const b = A.map(it.x, A.ax.y.lo, it.z ?? 0), top = A.map(it.x, it.y, it.z ?? 0); M.makeScale(w, Math.max(1e-4, (top.y - b.y) * pi), d); M.setPosition(b); im.setMatrixAt(i, M); }); im.instanceMatrix.needsUpdate = true; };
    paint(); draw(0); im.userData.paint = paint; v.__barFields = (v.__barFields || []).concat([paint]);
    const oldTheme = v.onTheme; v.onTheme = (vv, dark) => { paint(); oldTheme && oldTheme(vv, dark); };
    return { mesh: im, draw };
  },

  /* ---------- height-field surface over a grid, colours by height, rising with progress ---------- */
  surface(v, A, xs, zs, Z, o = {}) {
    const nx = xs.length, nz = zs.length; const g = new THREE.PlaneGeometry(1, 1, nx - 1, nz - 1); g.rotateX(-Math.PI / 2);
    const pos = g.attributes.position; const base = new Float32Array(pos.count), colors = new Float32Array(pos.count * 3); const stops = (o.stops || ['cool', 'cyan', 'warm', 'hot']).map(k => new THREE.Color(this.pal[k] || k));
    let lo = Infinity, hi = -Infinity; for (const row of Z) for (const z of row) { lo = Math.min(lo, z); hi = Math.max(hi, z); } const vlo = o.lo ?? lo, vhi = o.hi ?? hi;
    const tmp = new THREE.Color();
    const paint = () => { const st = (o.stops || ['cool', 'cyan', 'warm', 'hot']).map(k => new THREE.Color(this.pal[k] || k)); for (let j = 0; j < nz; j++) for (let i = 0; i < nx; i++) { const k = j * nx + i, u = clamp((Z[j][i] - vlo) / ((vhi - vlo) || 1), 0, 1) * (st.length - 1), a = Math.min(st.length - 2, Math.floor(u)); tmp.copy(st[a]).lerp(st[a + 1], u - a); colors[k * 3] = tmp.r; colors[k * 3 + 1] = tmp.g; colors[k * 3 + 2] = tmp.b; } g.attributes.color.needsUpdate = true; };
    for (let j = 0; j < nz; j++) for (let i = 0; i < nx; i++) { const k = j * nx + i, q = A.map(xs[i], Z[j][i], zs[j]); pos.setXYZ(k, q.x, q.y, q.z); base[k] = q.y; }
    g.setAttribute('color', new THREE.BufferAttribute(colors, 3)); paint();
    const m = new THREE.MeshStandardMaterial({ vertexColors: true, roughness: .5, metalness: .1, side: THREE.DoubleSide, flatShading: !!o.flat }); const mesh = new THREE.Mesh(g, m); mesh.userData = { xs, zs, Z }; v.rig.add(mesh); v.pickables.push(mesh);
    const wire = o.wire === false ? null : new THREE.LineSegments(new THREE.WireframeGeometry(g), new THREE.LineBasicMaterial({ color: this.pal.ink, transparent: true, opacity: .12 })); if (wire) { wire.material.userData.pal = 'ink'; v.rig.add(wire); }
    const y0 = A.map(0, A.ax.y.lo, 0).y;
    const draw = p => { const e = ease(clamp(p, 0, 1)); for (let k = 0; k < pos.count; k++) pos.setY(k, y0 + (base[k] - y0) * e); pos.needsUpdate = true; g.computeVertexNormals(); if (wire) { wire.geometry.dispose(); wire.geometry = new THREE.WireframeGeometry(g); } };
    draw(0); const oldTheme = v.onTheme; v.onTheme = (vv, dark) => { paint(); oldTheme && oldTheme(vv, dark); };
    return { mesh, draw, wire };
  },

  /* ---------- tile field: a categorical matrix as instanced slabs that drop into place ---------- */
  tiles(v, cells, o = {}) {
    const n = cells.length, w = o.width ?? .9, d = o.depth ?? .9, h = o.height ?? .16; const geo = new THREE.BoxGeometry(w, h, d); geo.translate(0, h / 2, 0);
    const im = new THREE.InstancedMesh(geo, new THREE.MeshStandardMaterial({ roughness: .45, metalness: .15 }), n); im.userData.cells = cells; v.rig.add(im); v.pickables.push(im);
    const M = new THREE.Matrix4(), col = new THREE.Color();
    const paint = () => { cells.forEach((c, i) => { col.set(this.pal[c.color] || c.color); im.setColorAt(i, col); }); im.instanceColor.needsUpdate = true; };
    const draw = p => { cells.forEach((c, i) => { const pi = ease(clamp(p * 1.6 - (c.delay ?? i / n) * .6, 0, 1)); M.makeScale(1, Math.max(1e-3, pi * (c.scale ?? 1)), 1); M.setPosition(c.x, (1 - pi) * 1.2, c.z); im.setMatrixAt(i, M); }); im.instanceMatrix.needsUpdate = true; };
    paint(); draw(0); const oldTheme = v.onTheme; v.onTheme = (vv, dark) => { paint(); oldTheme && oldTheme(vv, dark); };
    return { mesh: im, draw };
  },

  /* ---------- exploded assembly: parts fly in from an offset as progress rises, staggered ---------- */
  exploded(v, parts) {
    parts.forEach(pt => { pt.rest = pt.mesh.position.clone(); pt.from = new THREE.Vector3(...(pt.from || [0, 2, 0])); pt.delay = pt.delay ?? 0; pt.span = pt.span ?? .6; if (pt.mesh.parent == null) v.rig.add(pt.mesh); });
    const draw = p => { parts.forEach(pt => { const e = ease(clamp((p - pt.delay) / pt.span, 0, 1)); pt.mesh.position.copy(pt.rest).addScaledVector(pt.from, 1 - e); if (pt.mesh.material && pt.fade !== false) { const ms = [].concat(pt.mesh.material); ms.forEach(m => { if (m.userData.baseOpacity == null) m.userData.baseOpacity = m.opacity; m.opacity = m.userData.baseOpacity * (.15 + .85 * e); m.transparent = true; }); } if (pt.onDraw) pt.onDraw(e); }); };
    draw(0); return { draw, parts };
  },

  /* ---------- count-up number bound to a view's progress ---------- */
  counter(el, target, p, fmt) { if (!el) return; el.textContent = fmt ? fmt(target * ease(p)) : Math.round(target * ease(p)).toLocaleString('en-US'); },

  dispose() { this.disposed = true; for (const v of this.views) v.scene.traverse(o => { o.geometry && o.geometry.dispose(); if (o.material) [].concat(o.material).forEach(m => { for (const k in m) if (m[k] && m[k].isTexture) m[k].dispose(); m.dispose(); }); }); if (this.renderer) { if (this.renderer.__env) this.renderer.__env.dispose(); this.renderer.dispose(); this.renderer.forceContextLoss && this.renderer.forceContextLoss(); } }
};
addEventListener('pagehide', () => Viz.dispose(), { once: true });
window.__viz = Viz;
export default Viz;
