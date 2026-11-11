/* aldmc-gpu.js: the 3D voxel test-particle Monte Carlo feature solver of lib/aldmc, running in the browser.
   ES module. Nothing here touches the DOM; docs/partials/gpu-solver.html imports it when the visitor clicks Run.

     import { createSolver, hole, taperedHole, slit, recessedSlit, numberDensity, meanSpeed } from './aldmc-gpu.js';
     const g = hole(100e-9, 2e-6, 100e-9 / 16);                      // AR-20 hole, 16 voxels across
     const mc = await createSolver({ geometry: g, params: { s: 5e-3, N0: 4e18, vbar: 410, seed: 1 }, backend: 'webgpu' });
     const counts = await mc.step(3.79e19, 0.005, 25000);           // n_mouth (1/m3), dt (s), M particles
     mc.floorCoverage(); mc.sliceMean(); mc.coverage(); mc.backend; mc.limits;

   Physics and discretisation are those of lib/aldmc (lib/aldmc/README.md, include/aldmc.h): molecular flow; particles
   enter the open mouth plane z = 0 with a cosine angular distribution and a uniform position over the open mouth
   voxels; straight flights through a Chebyshev distance field with an exact DDA for the last voxels; on a solid voxel
   face they stick with probability s w (1 - theta_face), are lost with probability r w, and otherwise re-emit
   diffusely about the facet normal; particles crossing z = 0 outward escape. A face carries N0 h^2 w sites with the
   projected-area weight w = |n . e_face| of the smooth normal estimated from the first moment of the occupancy
   (cone kernel 1 - r/R, R = 5 voxels), the staircase correction. A step of M particles represents n vbar/4 A_mouth dt
   molecules; the coverage is frozen during the step and then relaxed per face,
       theta <- 1 - (1 - theta) exp(-s hits w_mol / (N0 h^2)),   w_mol = n vbar/4 A_mouth dt / M.
   Random numbers: xoshiro128** streams, one per block of 64 particles on the CPU (seeded from seed, step, block) and
   one per invocation on the GPU (seed, step, invocation); the library uses xoshiro256** per block, so the streams are
   not bitwise the library's; the statistics are. Only the library's defaults are implemented (area correction on,
   facet re-emission, impingement update). Film growth does not narrow the feature; there are no gas-phase collisions.

   Backends. "webgpu": occupancy/distance grid as a storage buffer (u8 packed into u32), per-face hit and stick
   counters as atomic<u32>, one compute invocation per particle (aldmc-gpu.wgsl, trace_main), then a second pass over
   the faces for the coverage update (update_main); the per-step fates and the coverage are read back through one
   staging buffer. "cpu": the same algorithm in plain JavaScript on typed arrays (aldmc-cpu.js), used when
   navigator.gpu is absent or the adapter or device request fails. Both resolve step() to the same count object.
   Flight cap: a particle that has traversed maxSegments voxel segments (default 2,000,000; the library's default is
   1e8) is counted "left"; at 16 voxels across an AR-100 hole averages about 2,000 segments per particle. */

import { createCpuBackend } from './aldmc-cpu.js';

export const VERSION = '1.0.0';
export const K_B = 1.380649e-23;
export const TMA_MOLAR_MASS = 72.09;          /* g/mol, Al(CH3)3 */
export const DEFAULT_MAX_SEGMENTS = 2000000;

/** Number density (1/m3) at partial pressure p (Pa) and temperature T (K). */
export const numberDensity = (p, T) => p / (K_B * T);
/** Mean thermal speed (m/s) of a molecule of molar mass M (g/mol) at T (K): sqrt(8 kT / (pi m)). */
export const N_A = 6.02214076e23;
export const meanSpeed = (T, M = TMA_MOLAR_MASS) => Math.sqrt(8 * K_B * T / (Math.PI * (M * 1e-3 / N_A)));

export function defaultParams(p = {}) {
  return {
    s: 5e-3, r: 0, N0: 4e18, vbar: 410, seed: 1, areaCorrection: true, normalRadius: 5,
    maxSegments: DEFAULT_MAX_SEGMENTS, ...p
  };
}

/* ---------------------------------------------------------------------------------------------------------------
   Geometry builders (src/aldconf/voxel.py). Occupancy is a Uint8Array over (nx, ny, nz) stored as
   occ[(k ny + j) nx + i] (the layout aldmc_create takes), 1 solid, 0 gas; z increases into the feature, the plane
   z = 0 is the open mouth; lengths in metres. */
const roundVoxels = (length, voxel, minimum = 1) => Math.max(minimum, Math.round(length / voxel));

export function fromArray(occupancy, nx, ny, nz, voxel, periodic = [false, false], meta = {}) {
  if (occupancy.length !== nx * ny * nz) throw new Error('occupancy length must be nx*ny*nz');
  return { occupancy, nx, ny, nz, voxel, periodicX: !!periodic[0], periodicY: !!periodic[1], meta };
}

/** Cylindrical hole of diameter d and depth L (frustum when bottomRatio < 1), closed by a floor layer, `pad` solid
    voxels around. Every slice has exactly round(pi r^2) gas voxels: the nearest voxel centres to the axis, ties broken
    by the ravel order ix * n + iy as numpy's stable argsort does, so the digital cross-section matches voxel.hole. */
export function hole(d, L, voxel, pad = 1, bottomRatio = 1.0) {
  const nd = roundVoxels(d, voxel, 2), nL = roundVoxels(L, voxel, 1);
  const n = nd + 2 * pad, c = n / 2;
  const r2 = new Float64Array(n * n), order = new Array(n * n);
  for (let ix = 0; ix < n; ix++) for (let iy = 0; iy < n; iy++) {
    const x = ix + 0.5 - c, y = iy + 0.5 - c; r2[ix * n + iy] = x * x + y * y; order[ix * n + iy] = ix * n + iy;
  }
  order.sort((a, b) => r2[a] - r2[b] || a - b);
  const nz = nL + 1;
  const occ = new Uint8Array(n * n * nz).fill(1);
  for (let k = 0; k < nL; k++) {
    const zc = (k + 0.5) / nL;
    const radius = 0.5 * nd * (1 - (1 - bottomRatio) * zc);
    const count = Math.round(Math.PI * radius * radius);
    for (let q = 0; q < count; q++) {
      const ix = Math.floor(order[q] / n), iy = order[q] % n;          /* numpy ravel: index = ix * n + iy -> occ[ix, iy] */
      occ[(k * n + iy) * n + ix] = 0;
    }
  }
  const meta = { kind: bottomRatio === 1 ? 'hole' : 'tapered_hole', d, L, bottomRatio, voxelsAcross: nd, voxelsDeep: nL, dActual: nd * voxel, LActual: nL * voxel };
  return fromArray(occ, n, n, nz, voxel, [false, false], meta);
}

export function taperedHole(d, L, bottomRatio, voxel, pad = 1) {
  if (!(bottomRatio > 0 && bottomRatio <= 1)) throw new Error('bottomRatio must be in (0, 1]');
  return hole(d, L, voxel, pad, bottomRatio);
}

/** Slit of gap w (along y) and depth L, infinite along x (periodic over nx voxels), closed by a floor. */
export function slit(w, L, voxel, nx = 4, pad = 1) {
  const nw = roundVoxels(w, voxel, 1), nL = roundVoxels(L, voxel, 1);
  const ny = nw + 2 * pad, nz = nL + 1;
  const occ = new Uint8Array(nx * ny * nz).fill(1);
  for (let k = 0; k < nL; k++) for (let j = pad; j < pad + nw; j++) for (let i = 0; i < nx; i++) occ[(k * ny + j) * nx + i] = 0;
  const meta = { kind: 'slit', w, L, voxelsAcross: nw, voxelsDeep: nL, wActual: nw * voxel, LActual: nL * voxel };
  return fromArray(occ, nx, ny, nz, voxel, [true, false], meta);
}

/** The word-line structure (aldconf.recess.RecessedSlit): a slit of gap w and depth L whose two walls carry lateral
    recesses of height gR (along z), lateral length LR (along y) and vertical pitch pR, periodic along the slit (x).
    The stack starts with a solid (oxide) layer of thickness pR - gR at the mouth when oxideFirst. gR and pR are
    rounded to whole voxels; meta records the actual values. */
export function recessedSlit(w, L, gR, LR, pR, voxel, nx = 4, pad = 1, oxideFirst = true) {
  const nw = roundVoxels(w, voxel, 1), nL = roundVoxels(L, voxel, 1), ng = roundVoxels(gR, voxel, 1);
  const npitch = roundVoxels(pR, voxel, ng + 1), nLr = roundVoxels(LR, voxel, 1);
  const ny = nw + 2 * nLr + 2 * pad, nz = nL + 1;
  const occ = new Uint8Array(nx * ny * nz).fill(1);
  const y0 = pad + nLr;
  const clear = (j0, j1, k0, k1) => { for (let k = k0; k < k1; k++) for (let j = j0; j < j1; j++) for (let i = 0; i < nx; i++) occ[(k * ny + j) * nx + i] = 0; };
  clear(y0, y0 + nw, 0, nL);
  const ox = npitch - ng, recesses = [];
  for (let k = 0; k < nL; k += npitch) {
    const start = k + (oxideFirst ? ox : 0), end = Math.min(start + ng, nL);
    if (start >= nL) break;
    if (end - start === ng) { clear(pad, y0, start, end); clear(y0 + nw, y0 + nw + nLr, start, end); recesses.push([start, end]); }
  }
  const meta = { kind: 'recessed_slit', w, L, gR, LR, pR, voxelsAcross: nw, voxelsDeep: nL, wActual: nw * voxel, LActual: nL * voxel,
    gRActual: ng * voxel, pRActual: npitch * voxel, LRActual: nLr * voxel, nRecesses: recesses.length, recessK: recesses, pad, ySlit: [y0, y0 + nw] };
  return fromArray(occ, nx, ny, nz, voxel, [true, false], meta);
}

/* ---------------------------------------------------------------------------------------------------------------
   Grid preprocessing shared by both backends (aldmc_t::build): distance field, exposed faces, face weights, mouth. */
const yieldNow = () => new Promise(r => setTimeout(r, 0));

export async function buildGrid(geom, params, onProgress) {
  const { nx, ny, nz, occupancy: occ, periodicX, periodicY, voxel: h } = geom;
  const n = nx * ny * nz, nxy = nx * ny;
  const idx = (i, j, k) => (k * ny + j) * nx + i;
  /* Chebyshev distance to the nearest solid voxel or non-periodic boundary (mouth and far plane included), by a
     multi-source BFS over the 26-neighbourhood; 0 solid, capped at 255. */
  const dist = new Uint8Array(n);
  for (let v = 0; v < n; v++) dist[v] = occ[v] ? 0 : 255;
  let frontier = [];
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v = idx(i, j, k);
    if (dist[v] === 0) continue;
    let boundary = (k === 0 || k === nz - 1) || (!periodicX && (i === 0 || i === nx - 1)) || (!periodicY && (j === 0 || j === ny - 1));
    if (!boundary) {
      outer: for (let dk = -1; dk <= 1; dk++) for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) {
        if (!di && !dj && !dk) continue;
        let ii = i + di, jj = j + dj, kk = k + dk;
        if (periodicX) ii = (ii + nx) % nx;
        if (periodicY) jj = (jj + ny) % ny;
        if (dist[idx(ii, jj, kk)] === 0) { boundary = true; break outer; }
      }
    }
    if (boundary) { dist[v] = 1; frontier.push(v); }
  }
  let level = 1;
  while (frontier.length && level < 255) {
    const next = [];
    for (const v of frontier) {
      const i = v % nx, j = ((v - i) / nx) % ny, k = Math.floor(v / nxy);
      for (let dk = -1; dk <= 1; dk++) {
        const kk = k + dk; if (kk < 0 || kk >= nz) continue;
        for (let dj = -1; dj <= 1; dj++) {
          let jj = j + dj; if (periodicY) jj = (jj + ny) % ny; else if (jj < 0 || jj >= ny) continue;
          for (let di = -1; di <= 1; di++) {
            let ii = i + di; if (periodicX) ii = (ii + nx) % nx; else if (ii < 0 || ii >= nx) continue;
            const u = idx(ii, jj, kk);
            if (dist[u] === 255) { dist[u] = level + 1; next.push(u); }
          }
        }
      }
    }
    frontier = next; level++;
  }
  if (onProgress) { onProgress('faces'); await yieldNow(); }
  /* exposed faces: (gas voxel, direction 0 +x, 1 -x, 2 +y, 3 -y, 4 +z, 5 -z toward the solid neighbour) */
  const faceBase = new Int32Array(n).fill(-1);
  const faceOfDir = [], faceVoxel = [], faceDir = [], faceK = [];
  const neighbour = (i, j, k, dir) => {
    let ii = i, jj = j, kk = k;
    switch (dir) { case 0: ii++; break; case 1: ii--; break; case 2: jj++; break; case 3: jj--; break; case 4: kk++; break; default: kk--; }
    if (kk < 0 || kk >= nz) return -1;
    if (ii < 0 || ii >= nx) { if (!periodicX) return -1; ii = (ii + nx) % nx; }
    if (jj < 0 || jj >= ny) { if (!periodicY) return -1; jj = (jj + ny) % ny; }
    return idx(ii, jj, kk);
  };
  let nBlocks = 0;
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const v = idx(i, j, k);
    if (dist[v] === 0) continue;
    let any = false; const ids = [-1, -1, -1, -1, -1, -1];
    for (let d = 0; d < 6; d++) {
      const u = neighbour(i, j, k, d);
      if (u >= 0 && dist[u] === 0) { ids[d] = faceVoxel.length; faceVoxel.push(v); faceDir.push(d); faceK.push(k); any = true; }
    }
    if (any) { faceBase[v] = nBlocks++; for (let d = 0; d < 6; d++) faceOfDir.push(ids[d]); }
  }
  const nFaces = faceVoxel.length;
  /* projected-area weights from the smooth normal (first moment of the occupancy, cone kernel 1 - r/R) */
  const faceW = new Float64Array(nFaces).fill(1);
  if (params.areaCorrection) {
    const R = params.normalRadius > 0 ? params.normalRadius : 5, R2 = R * R + 1e-9;
    const occClamped = (i, j, k) => {
      if (k < 0) k = 0; else if (k >= nz) k = nz - 1;
      if (periodicX) i = ((i % nx) + nx) % nx; else if (i < 0) i = 0; else if (i >= nx) i = nx - 1;
      if (periodicY) j = ((j % ny) + ny) % ny; else if (j < 0) j = 0; else if (j >= ny) j = ny - 1;
      return dist[idx(i, j, k)] === 0;
    };
    /* the kernel offsets inside the ball, once */
    const offs = [];
    for (let dk = -R - 1; dk <= R + 1; dk++) for (let dj = -R - 1; dj <= R + 1; dj++) for (let di = -R - 1; di <= R + 1; di++) offs.push([di, dj, dk]);
    let lastYield = performance.now();
    for (let f = 0; f < nFaces; f++) {
      const v = faceVoxel[f], i = v % nx, j = ((v - i) / nx) % ny, k = Math.floor(v / nxy), d = faceDir[f], axis = d >> 1;
      const c = [i + 0.5, j + 0.5, k + 0.5]; c[axis] += (d % 2 === 0) ? 0.5 : -0.5;
      let mx = 0, my = 0, mz = 0;
      for (let q = 0; q < offs.length; q++) {
        const o = offs[q], ii = i + o[0], jj = j + o[1], kk = k + o[2];
        const rx = ii + 0.5 - c[0], ry = jj + 0.5 - c[1], rz = kk + 0.5 - c[2], r2 = rx * rx + ry * ry + rz * rz;
        if (r2 > R2) continue;
        if (!occClamped(ii, jj, kk)) continue;
        const wk = 1 - Math.sqrt(r2) / R;
        mx += wk * rx; my += wk * ry; mz += wk * rz;
      }
      const norm = Math.sqrt(mx * mx + my * my + mz * mz);
      let wf = 1;
      if (norm > 0) wf = Math.abs(axis === 0 ? mx : axis === 1 ? my : mz) / norm;
      faceW[f] = Math.min(1, Math.max(1e-3, wf));
      if (onProgress && (f & 1023) === 1023) { const now = performance.now(); if (now - lastYield > 60) { onProgress('weights', f / nFaces); await yieldNow(); lastYield = performance.now(); } }
    }
  }
  const stickProb = new Float64Array(nFaces), lossProb = new Float64Array(nFaces);
  for (let f = 0; f < nFaces; f++) { stickProb[f] = params.s * faceW[f]; lossProb[f] = params.r * faceW[f]; }
  const mouthList = [];
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if (dist[idx(i, j, 0)] !== 0) mouthList.push(idx(i, j, 0));
  const mouth = Uint32Array.from(mouthList);
  let siteArea = 0; for (let f = 0; f < nFaces; f++) siteArea += faceW[f] * h * h;
  return {
    nx, ny, nz, h, periodicX, periodicY, dist, faceBase, faceOfDir: Int32Array.from(faceOfDir), nBlocks,
    faceVoxel: Int32Array.from(faceVoxel), faceDir: Uint8Array.from(faceDir), faceK: Int32Array.from(faceK),
    faceW, stickProb, lossProb, nFaces, mouth, mouthArea: mouth.length * h * h, siteArea
  };
}

/* ---------------------------------------------------------------------------------------------------------------
   WebGPU backend */
let wgslSource = null;
async function loadWgsl() {
  if (wgslSource) return wgslSource;
  const res = await fetch(new URL('./aldmc-gpu.wgsl', import.meta.url));
  if (!res.ok) throw new Error('could not load aldmc-gpu.wgsl: ' + res.status);
  wgslSource = await res.text();
  return wgslSource;
}

/** Request an adapter and device; returns {adapter, device, info} or throws. */
export async function requestGpu(opts = {}) {
  if (!navigator.gpu) throw new Error('navigator.gpu is not available (WebGPU off or an insecure context)');
  const adapter = await navigator.gpu.requestAdapter({ powerPreference: opts.powerPreference || 'high-performance' });
  if (!adapter) throw new Error('no WebGPU adapter');
  const info = adapter.info || (adapter.requestAdapterInfo ? await adapter.requestAdapterInfo() : {});
  const want = {};
  for (const k of ['maxStorageBufferBindingSize', 'maxBufferSize', 'maxComputeWorkgroupsPerDimension']) if (adapter.limits[k] !== undefined) want[k] = adapter.limits[k];
  const device = await adapter.requestDevice({ requiredLimits: want });
  return { adapter, device, info: { vendor: info.vendor || '', architecture: info.architecture || '', device: info.device || '', description: info.description || '' } };
}

/** Compile the shader and return the compilation messages (QA: an empty error list means the WGSL compiles). */
export async function compileCheck(device) {
  const module = device.createShaderModule({ code: await loadWgsl() });
  const info = await module.getCompilationInfo();
  return info.messages.map(m => ({ type: m.type, line: m.lineNum, col: m.linePos, message: m.message }));
}

async function createGpuBackend(grid, params, gpu) {
  const { device } = gpu;
  const { nx, ny, nz, dist, faceBase, faceOfDir, stickProb, lossProb, nFaces, mouth, periodicX, periodicY, h } = grid;
  const nVox = nx * ny * nz;
  if (nFaces === 0) throw new Error('geometry has no exposed faces');
  if (mouth.length === 0) throw new Error('geometry has no open mouth voxels');
  const maxWg = device.limits.maxComputeWorkgroupsPerDimension;
  const maxM = Math.min(maxWg * 64, 1 << 22);

  const code = await loadWgsl();
  const module = device.createShaderModule({ code, label: 'aldmc' });
  const ci = await module.getCompilationInfo();
  const errs = ci.messages.filter(m => m.type === 'error');
  if (errs.length) throw new Error('WGSL: ' + errs.map(m => `${m.lineNum}:${m.linePos} ${m.message}`).join('; '));

  const mk = (size, usage, label) => device.createBuffer({ size: Math.max(16, Math.ceil(size / 4) * 4), usage, label });
  const ST = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC;
  /* dist as u8 packed in u32 */
  const packed = new Uint32Array(Math.ceil(nVox / 4));
  for (let v = 0; v < nVox; v++) packed[v >> 2] |= dist[v] << ((v & 3) * 8);
  const bDist = mk(packed.byteLength, ST, 'dist'); device.queue.writeBuffer(bDist, 0, packed);
  const topo = new Int32Array(nVox + faceOfDir.length); topo.set(faceBase, 0); topo.set(faceOfDir, nVox);
  const bTopo = mk(topo.byteLength, ST, 'topo'); device.queue.writeBuffer(bTopo, 0, topo);
  const bMouth = mk(mouth.byteLength, ST, 'mouth'); device.queue.writeBuffer(bMouth, 0, mouth);
  const bTheta = mk(nFaces * 4, ST, 'theta');
  const prob = new Float32Array(2 * nFaces); prob.set(stickProb, 0); prob.set(lossProb, nFaces);
  const bProb = mk(prob.byteLength, ST, 'prob'); device.queue.writeBuffer(bProb, 0, prob);
  const accWords = 8 + 2 * nFaces;
  const bAcc = mk(accWords * 4, ST, 'acc');
  const bParams = device.createBuffer({ size: 64, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST, label: 'params' });
  const stagingBytes = 32 + Math.ceil(nFaces * 4 / 4) * 4;
  const bStaging = device.createBuffer({ size: stagingBytes, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST, label: 'staging' });

  const layout = device.createBindGroupLayout({ entries: [
    { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
    { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
    { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
    { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
    { binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
    { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
    { binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } }
  ] });
  const bind = device.createBindGroup({ layout, entries: [
    { binding: 0, resource: { buffer: bParams } }, { binding: 1, resource: { buffer: bDist } }, { binding: 2, resource: { buffer: bTopo } },
    { binding: 3, resource: { buffer: bMouth } }, { binding: 4, resource: { buffer: bTheta } }, { binding: 5, resource: { buffer: bProb } },
    { binding: 6, resource: { buffer: bAcc } }
  ] });
  const pl = device.createPipelineLayout({ bindGroupLayouts: [layout] });
  const [pTrace, pUpdate] = await Promise.all([
    device.createComputePipelineAsync({ layout: pl, compute: { module, entryPoint: 'trace_main' }, label: 'trace' }),
    device.createComputePipelineAsync({ layout: pl, compute: { module, entryPoint: 'update_main' }, label: 'update' })
  ]);

  const theta = new Float32Array(nFaces);
  const paramsU32 = new Uint32Array(16), paramsF32 = new Float32Array(paramsU32.buffer);
  const zeros8 = new Uint32Array(8);
  const cum = { escaped: 0, stuck: 0, lost: 0, left: 0, flights: 0 };
  let last = { escaped: 0, stuck: 0, lost: 0, left: 0, flights: 0 };
  let time = 0, stepIndex = 0, busy = false, lost = false;
  device.lost.then(() => { lost = true; });

  async function step(nMouth, dt, M) {
    if (busy) throw new Error('step already running');
    if (lost) throw new Error('GPU device lost');
    if (!(dt >= 0) || !(nMouth >= 0) || !(M >= 0)) throw new Error('dt, nMouth must be >= 0 and M >= 0');
    busy = true;
    const t0 = performance.now();
    M = Math.min(maxM, Math.floor(M));
    let escaped = 0, stuck = 0, lostN = 0, left = 0, flights = 0;
    try {
      if (M > 0 && nMouth > 0 && dt > 0) {
        const wmol = nMouth * params.vbar / 4 * grid.mouthArea * dt / M;
        const perHit = params.s * wmol / (params.N0 * h * h);
        paramsU32.set([nx, ny, nz, nVox, mouth.length, M, nFaces, stepIndex >>> 0, params.seed >>> 0, periodicX ? 1 : 0, periodicY ? 1 : 0, params.maxSegments >>> 0]);
        paramsF32[12] = perHit;
        device.queue.writeBuffer(bParams, 0, paramsU32);
        device.queue.writeBuffer(bAcc, 0, zeros8);
        const enc = device.createCommandEncoder();
        let pass = enc.beginComputePass({ label: 'trace' });
        pass.setPipeline(pTrace); pass.setBindGroup(0, bind); pass.dispatchWorkgroups(Math.ceil(M / 64)); pass.end();
        enc.copyBufferToBuffer(bAcc, 0, bStaging, 0, 32);
        pass = enc.beginComputePass({ label: 'update' });
        pass.setPipeline(pUpdate); pass.setBindGroup(0, bind); pass.dispatchWorkgroups(Math.ceil(nFaces / 64)); pass.end();
        enc.copyBufferToBuffer(bTheta, 0, bStaging, 32, nFaces * 4);
        device.queue.submit([enc.finish()]);
        await bStaging.mapAsync(GPUMapMode.READ);
        const range = bStaging.getMappedRange();
        const counts = new Uint32Array(range, 0, 8);
        escaped = counts[0]; stuck = counts[1]; lostN = counts[2]; left = counts[3]; flights = counts[4] + counts[5] * 4294967296;
        theta.set(new Float32Array(range, 32, nFaces));
        bStaging.unmap();
      }
      time += dt; stepIndex++;
      last = { escaped, stuck, lost: lostN, left, flights, particles: escaped + stuck + lostN + left, wallMs: performance.now() - t0 };
      cum.escaped += escaped; cum.stuck += stuck; cum.lost += lostN; cum.left += left; cum.flights += flights;
      return last;
    } finally { busy = false; }
  }

  function reset() {
    theta.fill(0); device.queue.writeBuffer(bTheta, 0, theta);
    device.queue.writeBuffer(bAcc, 0, new Uint32Array(accWords));
    time = 0; stepIndex = 0; cum.escaped = cum.stuck = cum.lost = cum.left = cum.flights = 0;
  }

  const limits = {};
  for (const k of ['maxComputeWorkgroupsPerDimension', 'maxComputeInvocationsPerWorkgroup', 'maxStorageBufferBindingSize', 'maxBufferSize', 'maxStorageBuffersPerShaderStage']) limits[k] = device.limits[k];

  return {
    backend: 'webgpu', limits, info: gpu.info, maxParticlesPerStep: maxM, device,
    step, reset, coverage: () => theta,
    counts: () => ({ ...cum }), lastCounts: () => ({ ...last }),
    get time() { return time; }, get stepIndex() { return stepIndex; },
    destroy() { for (const b of [bDist, bTopo, bMouth, bTheta, bProb, bAcc, bParams, bStaging]) b.destroy(); }
  };
}

/* ---------------------------------------------------------------------------------------------------------------
   createSolver: geometry + params -> solver (webgpu or cpu) with the common queries. */
/**
 * @param {object} o
 * @param {object} o.geometry  from hole/taperedHole/slit/recessedSlit/fromArray
 * @param {object} o.params    {s, r, N0, vbar, seed, areaCorrection, normalRadius, maxSegments}
 * @param {string} o.backend   "webgpu" (falls back to cpu when unavailable unless o.strict) | "cpu"
 * @param {function} o.onProgress optional (stage, fraction) during preprocessing
 */
export async function createSolver(o) {
  const params = defaultParams(o.params || {});
  if (!(params.s >= 0 && params.s <= 1)) throw new Error('s must be in [0, 1]');
  if (!(params.r >= 0 && params.r <= 1) || params.s + params.r > 1) throw new Error('r must be in [0, 1] and s + r <= 1');
  if (!(params.N0 > 0) || !(params.vbar > 0)) throw new Error('N0 and vbar must be > 0');
  const geom = o.geometry;
  const grid = await buildGrid(geom, params, o.onProgress);
  let backend = null, fallbackReason = null;
  const want = o.backend || 'webgpu';
  if (want === 'webgpu') {
    try {
      const gpu = o.gpu || await requestGpu(o);
      backend = await createGpuBackend(grid, params, gpu);
    } catch (e) {
      if (o.strict) throw e;
      fallbackReason = String(e && e.message || e);
    }
  }
  if (!backend) backend = createCpuBackend(grid, params);

  const { nFaces, faceDir, faceK, faceW, nz } = grid;
  let kFloor = -1;
  for (let f = 0; f < nFaces; f++) if (faceDir[f] === 4 && faceK[f] > kFloor) kFloor = faceK[f];
  const solver = {
    geometry: geom, params, grid, backend: backend.backend, limits: backend.limits, info: backend.info, fallbackReason,
    maxParticlesPerStep: backend.maxParticlesPerStep || Infinity,
    nFaces, mouthArea: grid.mouthArea, siteArea: grid.siteArea, nMouthVoxels: grid.mouth.length,
    step: (nMouth, dt, M, opts) => backend.step(nMouth, dt, M, opts),
    reset: () => backend.reset(),
    coverage: () => backend.coverage(),
    counts: () => backend.counts(), lastCounts: () => backend.lastCounts(),
    get time() { return backend.time; }, get stepIndex() { return backend.stepIndex; },
    /** weighted mean coverage of the floor faces (+z) in the deepest slice that has any (aldmc_bottom_coverage) */
    floorCoverage() {
      if (kFloor < 0) return 0;
      const th = backend.coverage(); let num = 0, den = 0;
      for (let f = 0; f < nFaces; f++) if (faceDir[f] === 4 && faceK[f] === kFloor) { num += faceW[f] * th[f]; den += faceW[f]; }
      return den > 0 ? num / den : 0;
    },
    /** weighted mean coverage per z slice (aldmc_slice_coverage): {z (m, centres), mean, area (m2)} */
    sliceMean() {
      const th = backend.coverage(), num = new Float64Array(nz), den = new Float64Array(nz);
      for (let f = 0; f < nFaces; f++) { num[faceK[f]] += faceW[f] * th[f]; den[faceK[f]] += faceW[f]; }
      const mean = new Float64Array(nz), area = new Float64Array(nz), z = new Float64Array(nz);
      for (let k = 0; k < nz; k++) { mean[k] = den[k] > 0 ? num[k] / den[k] : 0; area[k] = den[k] * grid.h * grid.h; z[k] = (k + 0.5) * grid.h; }
      return { z, mean, area };
    },
    /** per-face geometry for drawing: gas voxel index, direction, slice, weight */
    faces() { return { voxel: grid.faceVoxel, direction: faceDir, k: faceK, weight: faceW, nx: grid.nx, ny: grid.ny, nz, h: grid.h }; },
    destroy() { backend.destroy(); }
  };
  return solver;
}

/** First crossing of `level` by linear interpolation between (t[i-1], y[i-1]) and (t[i], y[i]); NaN if never. */
export function timeToLevel(t, y, level = 0.9) {
  for (let i = 0; i < y.length; i++) if (y[i] >= level) {
    if (i === 0) return t[0];
    return t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1]);
  }
  return NaN;
}
