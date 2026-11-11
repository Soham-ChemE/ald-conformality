/* aldmc-cpu.js: the CPU fallback of the browser voxel Monte Carlo solver (see aldmc-gpu.js for the entry point, the
   geometry builders and the WebGPU backend). Plain JavaScript on typed arrays, a line-for-line port of the tracer in
   lib/aldmc/src/aldmc.cpp (aldmc_t::trace and aldmc_t::step): molecular flow, cosine entry over the open mouth voxels,
   straight flights through the Chebyshev distance field (a voxel at distance D allows D - 1 voxel lengths without a
   check) with an exact DDA for the last voxels, sticking s w (1 - theta_face), loss r w, diffuse re-emission about the
   facet normal, escape through the mouth plane z = 0, "left" through any other boundary or at the flight cap. The M
   particles of a step are processed in blocks of 64; block b of step k has its own xoshiro128** stream seeded from
   (seed, k, b), so a run is reproducible for a given seed. Positions are doubles in voxel units as in the C++ code.
   Differences from the library: xoshiro128** (32-bit words, which JavaScript handles at full speed) instead of
   xoshiro256**, so the random streams are not bitwise those of libaldmc; the statistics are the same. Only the
   defaults of the library are implemented: projected-area face weights, facet re-emission, the impingement
   (exponential) coverage update. */

const FATE_ESCAPED = 0, FATE_STUCK = 1, FATE_LOST = 2, FATE_LEFT = 3;

/* splitmix32-style finaliser (the 32-bit murmur3 fmix). The WGSL shader uses the same function. */
export function fmix32(x) {
  let z = x >>> 0;
  z = Math.imul(z ^ (z >>> 16), 0x85EBCA6B) >>> 0;
  z = Math.imul(z ^ (z >>> 13), 0xC2B2AE35) >>> 0;
  return (z ^ (z >>> 16)) >>> 0;
}

/* Stream seed from (seed, step, block): the shape of stream_seed() in aldmc.cpp with 32-bit constants. */
export function streamSeed(seed, step, block) {
  let x = (seed ^ 0xD1B54A32) >>> 0;
  const a = fmix32((x + 0x9E3779B9) >>> 0);
  x = (a ^ Math.imul(step >>> 0, 0x9E3779B9)) >>> 0;
  const b = fmix32((x + 0x9E3779B9) >>> 0);
  x = (b ^ Math.imul(block >>> 0, 0x85EBCA6B)) >>> 0;
  return fmix32((x + 0x9E3779B9) >>> 0);
}

/* xoshiro128** with the 4-word state in a Uint32Array; seeded by a splitmix sequence from the stream seed. */
export class Xoshiro128 {
  constructor() { this.s = new Uint32Array(4); }
  seed(x) {
    let base = x >>> 0;
    for (let i = 0; i < 4; i++) { base = (base + 0x9E3779B9) >>> 0; this.s[i] = fmix32(base); }
    if (!(this.s[0] | this.s[1] | this.s[2] | this.s[3])) this.s[0] = 1;
  }
  next() {
    const s = this.s;
    const r5 = Math.imul(s[1], 5) >>> 0;
    const result = Math.imul(((r5 << 7) | (r5 >>> 25)) >>> 0, 9) >>> 0;
    const t = (s[1] << 9) >>> 0;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = ((s[3] << 11) | (s[3] >>> 21)) >>> 0;
    return result;
  }
  uniform() { return this.next() * 2.3283064365386963e-10; }   /* [0, 1) */
}

/**
 * CPU backend on a preprocessed grid (aldmc-gpu.js buildGrid). Same interface as the WebGPU backend:
 * step(nMouth, dt, M, opts) -> Promise<counts>, coverage(), floorCoverage(), sliceMean(), reset(), time, stepIndex.
 */
export function createCpuBackend(grid, params) {
  const { nx, ny, nz, dist, faceBase, faceOfDir, stickProb, lossProb, nFaces, mouth, periodicX, periodicY, h } = grid;
  const nMouthVox = mouth.length;
  const theta = new Float64Array(nFaces);
  const hits = new Uint32Array(nFaces), sticks = new Uint32Array(nFaces);
  const maxSeg = params.maxSegments;
  const rng = new Xoshiro128();
  const cum = { escaped: 0, stuck: 0, lost: 0, left: 0, flights: 0 };
  let last = { escaped: 0, stuck: 0, lost: 0, left: 0, flights: 0 };
  let time = 0, stepIndex = 0, busy = false;

  /* direction and its inverse live in module-level scalars (no allocation per particle) */
  let dx = 0, dy = 0, dz = 0, ix = 0, iy = 0, iz = 0;
  function setInv() { ix = dx !== 0 ? 1 / dx : Infinity; iy = dy !== 0 ? 1 / dy : Infinity; iz = dz !== 0 ? 1 / dz : Infinity; }
  /* Lambertian direction about `axis` (0 x, 1 y, 2 z) with sign: the tangent projection of a cosine-distributed unit
     vector is uniform over the unit disk, sampled by rejection, as cosine_direction() in aldmc.cpp. */
  function cosineDir(axis, sign) {
    let x, y, r2;
    do { x = 2 * rng.uniform() - 1; y = 2 * rng.uniform() - 1; r2 = x * x + y * y; } while (r2 >= 1);
    const nrm = Math.sqrt(Math.max(0, 1 - r2));
    const a = sign * (nrm > 0 ? nrm : 1e-9);
    if (axis === 0) { dx = a; dy = x; dz = y; } else if (axis === 1) { dy = a; dz = x; dx = y; } else { dz = a; dx = x; dy = y; }
  }

  let segs = 0;   /* flights (voxel segments) of the current step */

  function trace(i, j, k, px, py, pz) {
    setInv();
    let v = (k * ny + j) * nx + i;
    let n = 0;
    for (;;) {
      if (++n > maxSeg) { segs += n; return FATE_LEFT; }
      const D = dist[v];
      if (D >= 2) {
        const L = D - 1;
        px += dx * L; py += dy * L; pz += dz * L;
        i = Math.floor(px); j = Math.floor(py); k = Math.floor(pz);
        if (periodicX) { while (i < 0) { i += nx; px += nx; } while (i >= nx) { i -= nx; px -= nx; } }
        if (periodicY) { while (j < 0) { j += ny; py += ny; } while (j >= ny) { j -= ny; py -= ny; } }
        if (i < 0) i = 0; else if (i >= nx) i = nx - 1;
        if (j < 0) j = 0; else if (j >= ny) j = ny - 1;
        if (k < 0) k = 0; else if (k >= nz) k = nz - 1;
        v = (k * ny + j) * nx + i;
        continue;
      }
      const tx = dx > 0 ? (i + 1 - px) * ix : (dx < 0 ? (i - px) * ix : Infinity);
      const ty = dy > 0 ? (j + 1 - py) * iy : (dy < 0 ? (j - py) * iy : Infinity);
      const tz = dz > 0 ? (k + 1 - pz) * iz : (dz < 0 ? (k - pz) * iz : Infinity);
      let axis, t;
      if (tx <= ty && tx <= tz) { axis = 0; t = tx; } else if (ty <= tz) { axis = 1; t = ty; } else { axis = 2; t = tz; }
      if (!(t >= 0)) t = 0;
      px += dx * t; py += dy * t; pz += dz * t;
      const pi = i, pj = j, pk = k, pv = v;
      let dirIndex;
      if (axis === 0) {
        if (dx > 0) { ++i; px = i; dirIndex = 0; } else { px = i; --i; dirIndex = 1; }
        if (i < 0) { if (!periodicX) { segs += n; return FATE_LEFT; } i += nx; px += nx; }
        else if (i >= nx) { if (!periodicX) { segs += n; return FATE_LEFT; } i -= nx; px -= nx; }
      } else if (axis === 1) {
        if (dy > 0) { ++j; py = j; dirIndex = 2; } else { py = j; --j; dirIndex = 3; }
        if (j < 0) { if (!periodicY) { segs += n; return FATE_LEFT; } j += ny; py += ny; }
        else if (j >= ny) { if (!periodicY) { segs += n; return FATE_LEFT; } j -= ny; py -= ny; }
      } else {
        if (dz > 0) { ++k; pz = k; dirIndex = 4; } else { pz = k; --k; dirIndex = 5; }
        if (k < 0) { segs += n; return FATE_ESCAPED; }
        if (k >= nz) { segs += n; return FATE_LEFT; }
      }
      v = (k * ny + j) * nx + i;
      if (dist[v] !== 0) continue;
      /* wall collision on face (pv, dirIndex) */
      const base = faceBase[pv];
      const f = base >= 0 ? faceOfDir[base * 6 + dirIndex] : -1;
      if (f < 0) { segs += n; return FATE_LEFT; }
      hits[f]++;
      const u = rng.uniform();
      const ps = stickProb[f] * (1 - theta[f]);
      if (u < ps) { sticks[f]++; segs += n; return FATE_STUCK; }
      if (u < ps + lossProb[f]) { segs += n; return FATE_LOST; }
      i = pi; j = pj; k = pk; v = pv;
      const inward = (axis === 0 ? dx : axis === 1 ? dy : dz) > 0 ? -1 : 1;
      cosineDir(axis, inward);
      setInv();
    }
  }

  const yieldNow = () => new Promise(r => setTimeout(r, 0));

  /**
   * Launch M particles for mouth density nMouth (1/m3) over dt (s), then apply the coverage update.
   * opts.yieldMs: give the event loop a turn after this much work (default 40 ms; 0 = never).
   * Resolves to the step's counts {escaped, stuck, lost, left, flights, particles, wallMs}.
   */
  async function step(nMouth, dt, M, opts = {}) {
    if (busy) throw new Error('step already running');
    if (!(dt >= 0) || !(nMouth >= 0) || !(M >= 0)) throw new Error('dt, nMouth must be >= 0 and M >= 0');
    busy = true;
    const t0 = performance.now();
    const yieldMs = opts.yieldMs === undefined ? 40 : opts.yieldMs;
    hits.fill(0); sticks.fill(0); segs = 0;
    let escaped = 0, stuck = 0, lost = 0, left = 0;
    M = Math.floor(M);
    if (M > 0 && nMouth > 0 && dt > 0 && nMouthVox > 0) {
      const wmol = nMouth * params.vbar / 4 * grid.mouthArea * dt / M;
      const nBlocks = Math.ceil(M / 64);
      let lastYield = performance.now();
      for (let b = 0; b < nBlocks; b++) {
        rng.seed(streamSeed(params.seed, stepIndex, b));
        const p1 = Math.min(M, (b + 1) * 64);
        for (let q = b * 64; q < p1; q++) {
          const m = Math.min(nMouthVox - 1, Math.floor(rng.uniform() * nMouthVox));
          const v = mouth[m];
          const i = v % nx, j = ((v - i) / nx) % ny;
          const px = i + rng.uniform(), py = j + rng.uniform();
          cosineDir(2, 1);
          const fate = trace(i, j, 0, px, py, 0);
          if (fate === FATE_ESCAPED) escaped++; else if (fate === FATE_STUCK) stuck++; else if (fate === FATE_LOST) lost++; else left++;
        }
        if (yieldMs > 0 && (b & 7) === 7) { const now = performance.now(); if (now - lastYield > yieldMs) { await yieldNow(); lastYield = performance.now(); } }
      }
      /* impingement update: theta <- 1 - (1 - theta) exp(-s hits w / (N0 h^2)) */
      const perHit = params.s * wmol / (params.N0 * h * h);
      for (let f = 0; f < nFaces; f++) if (hits[f]) theta[f] = 1 - (1 - theta[f]) * Math.exp(-perHit * hits[f]);
    }
    time += dt; stepIndex++;
    last = { escaped, stuck, lost, left, flights: segs, particles: escaped + stuck + lost + left, wallMs: performance.now() - t0 };
    cum.escaped += escaped; cum.stuck += stuck; cum.lost += lost; cum.left += left; cum.flights += segs;
    busy = false;
    return last;
  }

  function reset() { theta.fill(0); hits.fill(0); sticks.fill(0); time = 0; stepIndex = 0; cum.escaped = cum.stuck = cum.lost = cum.left = cum.flights = 0; }

  return {
    backend: 'cpu', limits: null, info: { description: 'JavaScript on the main thread' },
    step, reset, coverage: () => theta, lastHits: () => hits, lastSticks: () => sticks,
    counts: () => ({ ...cum }), lastCounts: () => ({ ...last }),
    get time() { return time; }, get stepIndex() { return stepIndex; },
    destroy() {}
  };
}
