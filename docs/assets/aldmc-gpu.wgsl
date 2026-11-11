// aldmc-gpu.wgsl: the voxel test-particle Monte Carlo of lib/aldmc as two compute passes.
//
// trace_main: one invocation traces one particle from the open mouth plane z = 0 to its fate (escaped through the
// mouth, stuck, lost, left through another boundary or at the flight cap). Molecular flow: straight flights through
// the Chebyshev distance field (a voxel at distance D allows D - 1 voxel lengths with no check), an exact DDA for the
// last voxels, sticking s w (1 - theta_face) on the face's own coverage frozen for the step, loss r w, diffuse
// (cosine) re-emission about the facet normal. Wall collisions and sticking events are counted per face with atomics,
// the particle fates and the flights (voxel segments) in a small counter block. Every invocation has its own
// xoshiro128** stream seeded from (seed, step, invocation). Positions are kept as (voxel index, fraction in [0, 1])
// so f32 precision is the same in every voxel of a 1600-voxel-deep hole.
//
// update_main: one invocation per face applies the impingement (Langmuir) update
//     theta <- 1 - (1 - theta) exp(-s hits w_mol / (N0 h^2))      (per_hit = s w_mol / (N0 h^2) is passed in)
// and clears the face counters for the next step.
//
// Buffers (6 storage + 1 uniform; the default limit is 8 storage buffers per stage):
//   dist   u32, four u8 per word: 0 solid, else Chebyshev distance of the gas voxel to the nearest solid or boundary
//   topo   i32: [face_base per voxel (-1 or a block index) | face_of_dir, 6 per block (-1 if none)]
//   mouth  u32: linear indices of the gas voxels in slice k = 0
//   theta  f32 per face
//   prob   f32: [stick probability s w per face | loss probability r w per face]
//   acc    atomic<u32>: [0 escaped, 1 stuck, 2 lost, 3 left, 4 flights low word, 5 flights high word, 6-7 unused |
//                        hits per face | sticks per face]

struct Params {
  nx: u32, ny: u32, nz: u32, nvox: u32,
  nmouth: u32, m: u32, nfaces: u32, step: u32,
  seed: u32, periodic_x: u32, periodic_y: u32, max_seg: u32,
  per_hit: f32, pad0: u32, pad1: u32, pad2: u32,
}

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> dist: array<u32>;
@group(0) @binding(2) var<storage, read> topo: array<i32>;
@group(0) @binding(3) var<storage, read> mouth: array<u32>;
@group(0) @binding(4) var<storage, read_write> theta: array<f32>;
@group(0) @binding(5) var<storage, read> prob: array<f32>;
@group(0) @binding(6) var<storage, read_write> acc: array<atomic<u32>>;

const C_ESCAPED: u32 = 0u;
const C_STUCK: u32 = 1u;
const C_LOST: u32 = 2u;
const C_LEFT: u32 = 3u;
const C_SEG_LO: u32 = 4u;
const C_SEG_HI: u32 = 5u;
const HITS0: u32 = 8u;

// ---- xoshiro128** (per-invocation state) ----------------------------------------------------------------------
var<private> s0: u32;
var<private> s1: u32;
var<private> s2: u32;
var<private> s3: u32;

fn rotl(x: u32, k: u32) -> u32 { return (x << k) | (x >> (32u - k)); }

fn fmix(x: u32) -> u32 {
  var z = x;
  z = (z ^ (z >> 16u)) * 0x85EBCA6Bu;
  z = (z ^ (z >> 13u)) * 0xC2B2AE35u;
  return z ^ (z >> 16u);
}

fn seed_rng(seed: u32, step: u32, block: u32) {
  var x = seed ^ 0xD1B54A32u;
  let a = fmix(x + 0x9E3779B9u);
  x = a ^ (step * 0x9E3779B9u);
  let b = fmix(x + 0x9E3779B9u);
  x = b ^ (block * 0x85EBCA6Bu);
  var base = fmix(x + 0x9E3779B9u);
  base += 0x9E3779B9u; s0 = fmix(base);
  base += 0x9E3779B9u; s1 = fmix(base);
  base += 0x9E3779B9u; s2 = fmix(base);
  base += 0x9E3779B9u; s3 = fmix(base);
  if ((s0 | s1 | s2 | s3) == 0u) { s0 = 1u; }
}

fn next_u32() -> u32 {
  let result = rotl(s1 * 5u, 7u) * 9u;
  let t = s1 << 9u;
  s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3;
  s2 ^= t;
  s3 = rotl(s3, 11u);
  return result;
}

// uniform in [0, 1) with 24 random bits (exact in f32)
fn rnd01() -> f32 { return f32(next_u32() >> 8u) * (1.0 / 16777216.0); }

// Lambertian direction about `axis` with sign: tangent projection uniform over the unit disk, by rejection.
fn cosine_dir(axis: u32, sign: f32) -> vec3f {
  var x = 0.0; var y = 0.0; var r2 = 2.0;
  for (var tries = 0u; tries < 64u; tries++) {
    x = 2.0 * rnd01() - 1.0; y = 2.0 * rnd01() - 1.0; r2 = x * x + y * y;
    if (r2 < 1.0) { break; }
  }
  if (r2 >= 1.0) { x = 0.0; y = 0.0; r2 = 0.0; }
  let nrm = sqrt(max(0.0, 1.0 - r2));
  var d = vec3f(0.0, 0.0, 0.0);
  d[axis] = sign * max(nrm, 1e-6);
  d[(axis + 1u) % 3u] = x;
  d[(axis + 2u) % 3u] = y;
  return d;
}

fn dist_at(v: u32) -> u32 { return (dist[v >> 2u] >> ((v & 3u) * 8u)) & 0xFFu; }

fn wrap(i: i32, n: i32) -> i32 { return ((i % n) + n) % n; }

// One particle from gas voxel (i0, j0, k0), fraction f0 in it, direction d0. Returns the fate; *segs the flights.
fn trace(i0: i32, j0: i32, k0: i32, f0: vec3f, d0: vec3f, segs: ptr<function, u32>) -> u32 {
  var i = i0; var j = j0; var k = k0;
  var f = f0; var d = d0;
  let nx = i32(P.nx); let ny = i32(P.ny); let nz = i32(P.nz);
  var n = 0u;
  var v = u32((k * ny + j) * nx + i);
  loop {
    n += 1u;
    if (n > P.max_seg) { *segs = n; return C_LEFT; }
    let D = dist_at(v);
    if (D >= 2u) {
      // free flight of D - 1 voxel lengths (minus a rounding margin) reaches neither a solid nor a boundary
      let L = f32(D - 1u) - 1e-4;
      f += d * L;
      let fl = floor(f);
      i += i32(fl.x); j += i32(fl.y); k += i32(fl.z);
      f -= fl;
      if (P.periodic_x == 1u) { i = wrap(i, nx); }
      if (P.periodic_y == 1u) { j = wrap(j, ny); }
      i = clamp(i, 0, nx - 1); j = clamp(j, 0, ny - 1); k = clamp(k, 0, nz - 1);
      v = u32((k * ny + j) * nx + i);
      continue;
    }
    // one DDA step to the nearest voxel boundary
    let tx = select(select(1e30, -f.x / d.x, d.x < 0.0), (1.0 - f.x) / d.x, d.x > 0.0);
    let ty = select(select(1e30, -f.y / d.y, d.y < 0.0), (1.0 - f.y) / d.y, d.y > 0.0);
    let tz = select(select(1e30, -f.z / d.z, d.z < 0.0), (1.0 - f.z) / d.z, d.z > 0.0);
    var axis = 0u; var t = tx;
    if (ty < t) { axis = 1u; t = ty; }
    if (tz < t) { axis = 2u; t = tz; }
    t = max(t, 0.0);
    f += d * t;
    let pi = i; let pj = j; let pk = k; let pv = v;
    var dir_index = 0u;
    if (axis == 0u) {
      if (d.x > 0.0) { i += 1; f.x = 0.0; dir_index = 0u; } else { i -= 1; f.x = 1.0; dir_index = 1u; }
      if (i < 0) { if (P.periodic_x == 0u) { *segs = n; return C_LEFT; } i += nx; }
      else if (i >= nx) { if (P.periodic_x == 0u) { *segs = n; return C_LEFT; } i -= nx; }
    } else if (axis == 1u) {
      if (d.y > 0.0) { j += 1; f.y = 0.0; dir_index = 2u; } else { j -= 1; f.y = 1.0; dir_index = 3u; }
      if (j < 0) { if (P.periodic_y == 0u) { *segs = n; return C_LEFT; } j += ny; }
      else if (j >= ny) { if (P.periodic_y == 0u) { *segs = n; return C_LEFT; } j -= ny; }
    } else {
      if (d.z > 0.0) { k += 1; f.z = 0.0; dir_index = 4u; } else { k -= 1; f.z = 1.0; dir_index = 5u; }
      if (k < 0) { *segs = n; return C_ESCAPED; }
      if (k >= nz) { *segs = n; return C_LEFT; }
    }
    v = u32((k * ny + j) * nx + i);
    if (dist_at(v) != 0u) { continue; }
    // wall collision on face (pv, dir_index)
    let base = topo[pv];
    if (base < 0) { *segs = n; return C_LEFT; }
    let fidx = topo[P.nvox + u32(base) * 6u + dir_index];
    if (fidx < 0) { *segs = n; return C_LEFT; }
    let fi = u32(fidx);
    atomicAdd(&acc[HITS0 + fi], 1u);
    let u = rnd01();
    let ps = prob[fi] * (1.0 - theta[fi]);
    if (u < ps) { atomicAdd(&acc[HITS0 + P.nfaces + fi], 1u); *segs = n; return C_STUCK; }
    if (u < ps + prob[P.nfaces + fi]) { *segs = n; return C_LOST; }
    // diffuse re-emission back into the gas voxel pv from the face it hit
    i = pi; j = pj; k = pk; v = pv;
    let forward = d[axis] > 0.0;
    f[axis] = select(0.0, 1.0, forward);
    d = cosine_dir(axis, select(1.0, -1.0, forward));
  }
}

@compute @workgroup_size(64)
fn trace_main(@builtin(global_invocation_id) gid: vec3u) {
  let q = gid.x;
  if (q >= P.m) { return; }
  seed_rng(P.seed, P.step, q);
  // uniform mouth voxel, uniform position in it, cosine direction about +z
  let mi = min(P.nmouth - 1u, u32(rnd01() * f32(P.nmouth)));
  let v = mouth[mi];
  let i = i32(v % P.nx);
  let j = i32((v / P.nx) % P.ny);
  let f = vec3f(rnd01(), rnd01(), 0.0);
  let d = cosine_dir(2u, 1.0);
  var segs = 0u;
  let fate = trace(i, j, 0, f, d, &segs);
  atomicAdd(&acc[fate], 1u);
  let old = atomicAdd(&acc[C_SEG_LO], segs);
  if (old > 0xFFFFFFFFu - segs) { atomicAdd(&acc[C_SEG_HI], 1u); }
}

@compute @workgroup_size(64)
fn update_main(@builtin(global_invocation_id) gid: vec3u) {
  let fi = gid.x;
  if (fi >= P.nfaces) { return; }
  let hcount = atomicExchange(&acc[HITS0 + fi], 0u);
  atomicStore(&acc[HITS0 + P.nfaces + fi], 0u);
  if (hcount > 0u) { theta[fi] = 1.0 - (1.0 - theta[fi]) * exp(-P.per_hit * f32(hcount)); }
}
