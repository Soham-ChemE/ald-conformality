// aldmc: 3D voxel test-particle Monte Carlo for ALD in the molecular-flow regime. See include/aldmc.h.
#include "aldmc.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <thread>
#include <vector>

namespace {

thread_local std::string g_create_error;

// ---------------------------------------------------------------------------------------------------------------
// Random numbers: xoshiro256** seeded through splitmix64. Own implementation so that the streams are identical on
// every platform and compiler (the standard library's distributions are implementation-defined).
struct Rng {
    uint64_t s[4];
    static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    static uint64_t splitmix(uint64_t& x) {
        x += 0x9E3779B97F4A7C15ULL;
        uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    void seed(uint64_t x) {
        for (int i = 0; i < 4; ++i) s[i] = splitmix(x);
        if (s[0] == 0 && s[1] == 0 && s[2] == 0 && s[3] == 0) s[0] = 1;
    }
    uint64_t next() {
        const uint64_t result = rotl(s[1] * 5, 7) * 9;
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t; s[3] = rotl(s[3], 45);
        return result;
    }
    double uniform() { return (double)(next() >> 11) * 0x1.0p-53; }   // [0, 1)
};

inline uint64_t stream_seed(uint64_t seed, uint64_t step, uint64_t block) {
    uint64_t x = seed ^ 0xD1B54A32D192ED03ULL;
    uint64_t a = Rng::splitmix(x);
    x = a ^ (step * 0x9E3779B97F4A7C15ULL);
    uint64_t b = Rng::splitmix(x);
    x = b ^ (block * 0xBF58476D1CE4E5B9ULL);
    return Rng::splitmix(x);
}

// Lambertian (cosine) direction about axis `axis` with sign `sign`: the projection of a cosine-distributed unit
// vector on the tangent plane is uniform over the unit disk, so sample the disk by rejection (no trigonometry).
inline void cosine_direction(Rng& rng, int axis, double sign, double d[3]) {
    double x, y, r2;
    do {
        x = 2.0 * rng.uniform() - 1.0;
        y = 2.0 * rng.uniform() - 1.0;
        r2 = x * x + y * y;
    } while (r2 >= 1.0);
    const double nrm = std::sqrt(std::max(0.0, 1.0 - r2));
    const int b = (axis + 1) % 3, c = (axis + 2) % 3;
    d[axis] = sign * (nrm > 0.0 ? nrm : 1e-9);
    d[b] = x;
    d[c] = y;
}

enum Fate { FATE_ESCAPED = 0, FATE_STUCK = 1, FATE_LOST = 2, FATE_LEFT = 3 };

struct Counts {
    uint64_t escaped = 0, stuck = 0, lost = 0, left = 0, segments = 0;
    void add(const Counts& o) { escaped += o.escaped; stuck += o.stuck; lost += o.lost; left += o.left; segments += o.segments; }
};

struct ThreadBuf {
    std::vector<uint32_t> hits, sticks;
    Counts counts;
};

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
struct aldmc_t {
    int nx = 0, ny = 0, nz = 0;
    double h = 0.0;
    aldmc_params p{};
    int threads = 1;
    long long maxSegments = 100000000LL;

    std::vector<uint8_t> dist;          // 0 solid; >= 1 gas, Chebyshev distance to the nearest solid (capped 255)
    std::vector<int32_t> faceBase;      // per voxel: block index into faceOfDir (6 entries) or -1
    std::vector<int32_t> faceOfDir;     // face id per (boundary gas voxel, direction), -1 if none

    int nFaces = 0;
    std::vector<int64_t> faceVoxel;
    std::vector<uint8_t> faceDir;
    std::vector<int32_t> faceK;
    std::vector<double> faceW, theta, stickProb, lossProb;
    std::vector<double> faceN;          // 3 per face: estimated smooth normal pointing into the gas (unit)
    bool smoothEmission = false;

    std::vector<int64_t> mouthVoxels;   // linear indices of the gas voxels with k = 0
    double mouthArea = 0.0, siteArea = 0.0, time = 0.0;
    long long stepIndex = 0;

    Counts cum, last;
    std::vector<uint64_t> lastHits, lastSticks;
    std::vector<ThreadBuf> bufs;

    std::string err;

    int64_t idx(int i, int j, int k) const { return ((int64_t)k * ny + j) * nx + i; }
    bool solid(int i, int j, int k) const { return dist[idx(i, j, k)] == 0; }

    // --- construction -------------------------------------------------------------------------------------------
    void build(const unsigned char* occ) {
        const int64_t n = (int64_t)nx * ny * nz;
        dist.assign(n, 0);
        for (int64_t v = 0; v < n; ++v) dist[v] = occ[v] ? 0 : 255;
        distance_transform();
        enumerate_faces();
        face_weights();
        mouthVoxels.clear();
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                if (!solid(i, j, 0)) mouthVoxels.push_back(idx(i, j, 0));
        mouthArea = (double)mouthVoxels.size() * h * h;
        siteArea = 0.0;
        for (int f = 0; f < nFaces; ++f) siteArea += faceW[f] * h * h;
        theta.assign(nFaces, 0.0);
        lastHits.assign(nFaces, 0);
        lastSticks.assign(nFaces, 0);
        bufs.resize(threads);
        for (auto& b : bufs) { b.hits.assign(nFaces, 0); b.sticks.assign(nFaces, 0); }
    }

    // Chebyshev (L-infinity) distance of every gas voxel to the nearest solid voxel or to the domain boundary in a
    // non-periodic direction (the mouth plane and the far plane count as boundaries), by a multi-source BFS over
    // the 26-neighbourhood. A particle in a voxel of distance D may fly D - 1 voxel lengths without any check.
    void distance_transform() {
        const int64_t n = (int64_t)nx * ny * nz;
        std::vector<int64_t> frontier, next;
        frontier.reserve(n / 8 + 1);
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int64_t v = idx(i, j, k);
                    if (dist[v] == 0) continue;
                    bool boundary = (k == 0 || k == nz - 1) || (!p.periodicX && (i == 0 || i == nx - 1)) ||
                                    (!p.periodicY && (j == 0 || j == ny - 1));
                    if (!boundary) {
                        for (int dk = -1; dk <= 1 && !boundary; ++dk)
                            for (int dj = -1; dj <= 1 && !boundary; ++dj)
                                for (int di = -1; di <= 1 && !boundary; ++di) {
                                    if (!di && !dj && !dk) continue;
                                    int ii = i + di, jj = j + dj, kk = k + dk;
                                    if (p.periodicX) ii = (ii + nx) % nx;
                                    if (p.periodicY) jj = (jj + ny) % ny;
                                    if (dist[idx(ii, jj, kk)] == 0) boundary = true;
                                }
                    }
                    if (boundary) { dist[v] = 1; frontier.push_back(v); }
                }
        uint8_t level = 1;
        while (!frontier.empty() && level < 255) {
            next.clear();
            for (int64_t v : frontier) {
                const int i = (int)(v % nx), j = (int)((v / nx) % ny), k = (int)(v / ((int64_t)nx * ny));
                for (int dk = -1; dk <= 1; ++dk)
                    for (int dj = -1; dj <= 1; ++dj)
                        for (int di = -1; di <= 1; ++di) {
                            int ii = i + di, jj = j + dj, kk = k + dk;
                            if (kk < 0 || kk >= nz) continue;
                            if (p.periodicX) ii = (ii + nx) % nx; else if (ii < 0 || ii >= nx) continue;
                            if (p.periodicY) jj = (jj + ny) % ny; else if (jj < 0 || jj >= ny) continue;
                            const int64_t u = idx(ii, jj, kk);
                            if (dist[u] == 255) { dist[u] = (uint8_t)(level + 1); next.push_back(u); }
                        }
            }
            frontier.swap(next);
            ++level;
        }
        // voxels still at 255 are farther than 254 from any solid: keep 255 (a safe lower bound is used anyway)
    }

    // neighbour with wrap; returns false when outside a non-periodic boundary
    bool neighbour(int i, int j, int k, int dir, int& ii, int& jj, int& kk) const {
        ii = i; jj = j; kk = k;
        switch (dir) {
            case 0: ++ii; break; case 1: --ii; break;
            case 2: ++jj; break; case 3: --jj; break;
            case 4: ++kk; break; default: --kk; break;
        }
        if (kk < 0 || kk >= nz) return false;
        if (ii < 0 || ii >= nx) { if (!p.periodicX) return false; ii = (ii + nx) % nx; }
        if (jj < 0 || jj >= ny) { if (!p.periodicY) return false; jj = (jj + ny) % ny; }
        return true;
    }

    void enumerate_faces() {
        const int64_t n = (int64_t)nx * ny * nz;
        faceBase.assign(n, -1);
        faceOfDir.clear();
        faceVoxel.clear(); faceDir.clear(); faceK.clear();
        int32_t nBlocks = 0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const int64_t v = idx(i, j, k);
                    if (dist[v] == 0) continue;
                    int32_t ids[6];
                    bool any = false;
                    for (int d = 0; d < 6; ++d) {
                        int ii, jj, kk;
                        ids[d] = -1;
                        if (neighbour(i, j, k, d, ii, jj, kk) && solid(ii, jj, kk)) {
                            ids[d] = (int32_t)faceVoxel.size();
                            faceVoxel.push_back(v); faceDir.push_back((uint8_t)d); faceK.push_back(k);
                            any = true;
                        }
                    }
                    if (any) {
                        faceBase[v] = nBlocks++;
                        for (int d = 0; d < 6; ++d) faceOfDir.push_back(ids[d]);
                    }
                }
        nFaces = (int)faceVoxel.size();
    }

    // occupancy with clamping outside non-periodic boundaries (and above the mouth / below the far plane)
    int occ_clamped(int i, int j, int k) const {
        if (k < 0) k = 0; else if (k >= nz) k = nz - 1;
        if (p.periodicX) i = ((i % nx) + nx) % nx; else if (i < 0) i = 0; else if (i >= nx) i = nx - 1;
        if (p.periodicY) j = ((j % ny) + ny) % ny; else if (j < 0) j = 0; else if (j >= ny) j = ny - 1;
        return dist[idx(i, j, k)] == 0 ? 1 : 0;
    }

    // Projected-area weight of every face from the smooth normal estimated as the (negative) first moment of the
    // occupancy around the face centre with the cone kernel (1 - r/R), R voxels (default 5). The cone kernel's
    // staircase bias is below 2 % of the inventory for cylinders of 12 to 32 voxels across; a flat-topped (ball)
    // kernel over-counts the stair corners by 3 to 7 % at the same radius.
    void face_weights() {
        faceW.assign(nFaces, 1.0);
        stickProb.assign(nFaces, p.s);
        lossProb.assign(nFaces, p.r);
        faceN.assign((size_t)nFaces * 3, 0.0);
        for (int f = 0; f < nFaces; ++f) faceN[(size_t)f * 3 + faceDir[f] / 2] = (faceDir[f] % 2 == 0) ? -1.0 : 1.0;
        smoothEmission = p.areaCorrection && p.emission == ALDMC_EMISSION_SMOOTH;
        if (!p.areaCorrection) return;
        const int R = p.normalRadius > 0 ? p.normalRadius : 5;
        const double R2 = (double)R * R + 1e-9;
        std::vector<double> w(nFaces, 1.0);
        const int T = threads;
        std::vector<std::thread> pool;
        auto work = [&](int tid) {
            for (int f = tid; f < nFaces; f += T) {
                const int64_t v = faceVoxel[f];
                const int i = (int)(v % nx), j = (int)((v / nx) % ny), k = (int)(v / ((int64_t)nx * ny));
                const int d = faceDir[f];
                // face centre in voxel units relative to the gas voxel's centre
                double c[3] = {i + 0.5, j + 0.5, k + 0.5};
                const int axis = d / 2;
                c[axis] += (d % 2 == 0) ? 0.5 : -0.5;
                double m[3] = {0.0, 0.0, 0.0};
                for (int dk = -R - 1; dk <= R + 1; ++dk)
                    for (int dj = -R - 1; dj <= R + 1; ++dj)
                        for (int di = -R - 1; di <= R + 1; ++di) {
                            const int ii = i + di, jj = j + dj, kk = k + dk;
                            const double rx = ii + 0.5 - c[0], ry = jj + 0.5 - c[1], rz = kk + 0.5 - c[2];
                            const double r2 = rx * rx + ry * ry + rz * rz;
                            if (r2 > R2) continue;
                            if (!occ_clamped(ii, jj, kk)) continue;
                            const double wk = 1.0 - std::sqrt(r2) / (double)R;
                            m[0] += wk * rx; m[1] += wk * ry; m[2] += wk * rz;
                        }
                const double norm = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
                double wf = 1.0;
                if (norm > 0.0) {
                    wf = std::fabs(m[axis]) / norm;
                    // the moment points into the solid; the normal into the gas is its negative, and it must point
                    // out of this facet (same side as the facet's inward normal)
                    double nn[3] = {-m[0] / norm, -m[1] / norm, -m[2] / norm};
                    const double inward = (d % 2 == 0) ? -1.0 : 1.0;
                    if (nn[axis] * inward <= 1e-3) { nn[0] = nn[1] = nn[2] = 0.0; nn[axis] = inward; }
                    for (int a = 0; a < 3; ++a) faceN[(size_t)f * 3 + a] = nn[a];
                }
                w[f] = std::min(1.0, std::max(1e-3, wf));
            }
        };
        for (int t = 1; t < T; ++t) pool.emplace_back(work, t);
        work(0);
        for (auto& th : pool) th.join();
        faceW.swap(w);
        for (int f = 0; f < nFaces; ++f) { stickProb[f] = p.s * faceW[f]; lossProb[f] = p.r * faceW[f]; }
    }

    // --- one particle -----------------------------------------------------------------------------------------
    // Returns the fate; accumulates face hits and sticks in buf.
    int trace(Rng& rng, ThreadBuf& buf, int i, int j, int k, double px, double py, double pz, const double dir0[3]) {
        double d[3] = {dir0[0], dir0[1], dir0[2]};
        double inv[3];
        auto set_inv = [&]() {
            for (int a = 0; a < 3; ++a) inv[a] = (d[a] != 0.0) ? 1.0 / d[a] : std::numeric_limits<double>::infinity();
        };
        set_inv();
        const double INF = std::numeric_limits<double>::infinity();
        int64_t v = idx(i, j, k);
        uint64_t segments = 0;
        while (true) {
            if (++segments > (uint64_t)maxSegments) { buf.counts.segments += segments; return FATE_LEFT; }
            const uint8_t D = dist[v];
            if (D >= 2) {
                // free flight of D - 1 voxel lengths cannot reach a solid or a boundary
                const double L = (double)(D - 1);
                px += d[0] * L; py += d[1] * L; pz += d[2] * L;
                i = (int)std::floor(px); j = (int)std::floor(py); k = (int)std::floor(pz);
                if (p.periodicX) { while (i < 0) { i += nx; px += nx; } while (i >= nx) { i -= nx; px -= nx; } }
                if (p.periodicY) { while (j < 0) { j += ny; py += ny; } while (j >= ny) { j -= ny; py -= ny; } }
                // (non-periodic indices stay inside the domain by construction of the distance field)
                i = std::min(std::max(i, 0), nx - 1); j = std::min(std::max(j, 0), ny - 1); k = std::min(std::max(k, 0), nz - 1);
                v = idx(i, j, k);
                continue;
            }
            // one DDA step to the nearest voxel boundary
            const double tx = d[0] > 0.0 ? (i + 1 - px) * inv[0] : (d[0] < 0.0 ? (i - px) * inv[0] : INF);
            const double ty = d[1] > 0.0 ? (j + 1 - py) * inv[1] : (d[1] < 0.0 ? (j - py) * inv[1] : INF);
            const double tz = d[2] > 0.0 ? (k + 1 - pz) * inv[2] : (d[2] < 0.0 ? (k - pz) * inv[2] : INF);
            int axis; double t;
            if (tx <= ty && tx <= tz) { axis = 0; t = tx; } else if (ty <= tz) { axis = 1; t = ty; } else { axis = 2; t = tz; }
            if (!(t >= 0.0)) t = 0.0;
            px += d[0] * t; py += d[1] * t; pz += d[2] * t;
            const int pi = i, pj = j, pk = k;
            const int64_t pv = v;
            int dirIndex;
            if (axis == 0) {
                if (d[0] > 0.0) { ++i; px = (double)i; dirIndex = 0; } else { px = (double)i; --i; dirIndex = 1; }
                if (i < 0) { if (!p.periodicX) { buf.counts.segments += segments; return FATE_LEFT; } i += nx; px += nx; }
                else if (i >= nx) { if (!p.periodicX) { buf.counts.segments += segments; return FATE_LEFT; } i -= nx; px -= nx; }
            } else if (axis == 1) {
                if (d[1] > 0.0) { ++j; py = (double)j; dirIndex = 2; } else { py = (double)j; --j; dirIndex = 3; }
                if (j < 0) { if (!p.periodicY) { buf.counts.segments += segments; return FATE_LEFT; } j += ny; py += ny; }
                else if (j >= ny) { if (!p.periodicY) { buf.counts.segments += segments; return FATE_LEFT; } j -= ny; py -= ny; }
            } else {
                if (d[2] > 0.0) { ++k; pz = (double)k; dirIndex = 4; } else { pz = (double)k; --k; dirIndex = 5; }
                if (k < 0) { buf.counts.segments += segments; return FATE_ESCAPED; }
                if (k >= nz) { buf.counts.segments += segments; return FATE_LEFT; }
            }
            v = idx(i, j, k);
            if (dist[v] != 0) continue;   // gas: keep flying
            // wall collision on face (pv, dirIndex)
            const int32_t base = faceBase[pv];
            const int32_t f = base >= 0 ? faceOfDir[(int64_t)base * 6 + dirIndex] : -1;
            if (f < 0) { buf.counts.segments += segments; return FATE_LEFT; }   // cannot happen for a consistent grid
            ++buf.hits[f];
            const double u = rng.uniform();
            const double ps = stickProb[f] * (1.0 - theta[f]);
            if (u < ps) { ++buf.sticks[f]; buf.counts.segments += segments; return FATE_STUCK; }
            if (u < ps + lossProb[f]) { buf.counts.segments += segments; return FATE_LOST; }
            // diffuse re-emission into the gas voxel pv: about the smooth normal (directions into the facet
            // resampled, facet normal after 16 tries) or about the facet normal (axis, sign opposite to the step)
            i = pi; j = pj; k = pk; v = pv;
            const double inward = d[axis] > 0.0 ? -1.0 : 1.0;
            if (smoothEmission) {
                const double* nn = &faceN[(size_t)f * 3];
                bool ok = false;
                for (int tries = 0; tries < 16 && !ok; ++tries) {
                    double loc[3];
                    cosine_direction(rng, 2, 1.0, loc);               // about +z in a local frame
                    // local frame: e3 = nn, e1 = normalised (facet axis unit vector minus its nn component) or any
                    double e1[3] = {0.0, 0.0, 0.0};
                    const int ta = (axis + 1) % 3;
                    e1[ta] = 1.0;
                    double dp = e1[0] * nn[0] + e1[1] * nn[1] + e1[2] * nn[2];
                    for (int a = 0; a < 3; ++a) e1[a] -= dp * nn[a];
                    double l1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
                    if (l1 < 1e-9) { e1[0] = e1[1] = e1[2] = 0.0; e1[(axis + 2) % 3] = 1.0; dp = e1[0] * nn[0] + e1[1] * nn[1] + e1[2] * nn[2]; for (int a = 0; a < 3; ++a) e1[a] -= dp * nn[a]; l1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]); }
                    for (int a = 0; a < 3; ++a) e1[a] /= l1;
                    const double e2[3] = {nn[1] * e1[2] - nn[2] * e1[1], nn[2] * e1[0] - nn[0] * e1[2], nn[0] * e1[1] - nn[1] * e1[0]};
                    double cand[3];
                    for (int a = 0; a < 3; ++a) cand[a] = loc[0] * e1[a] + loc[1] * e2[a] + loc[2] * nn[a];
                    if (cand[axis] * inward > 1e-9) { d[0] = cand[0]; d[1] = cand[1]; d[2] = cand[2]; ok = true; }
                }
                if (!ok) cosine_direction(rng, axis, inward, d);
            } else {
                cosine_direction(rng, axis, inward, d);
            }
            set_inv();
        }
    }

    // --- one step ---------------------------------------------------------------------------------------------
    int step(double n_mouth, double dt, long long M) {
        if (!(dt >= 0.0) || !(n_mouth >= 0.0) || M < 0) { err = "aldmc_step: dt, n_mouth must be >= 0 and M >= 0"; return ALDMC_ERR_ARG; }
        last = Counts();
        std::fill(lastHits.begin(), lastHits.end(), 0);
        std::fill(lastSticks.begin(), lastSticks.end(), 0);
        if (M == 0 || n_mouth == 0.0 || dt == 0.0 || mouthVoxels.empty()) {
            time += dt; ++stepIndex; return ALDMC_OK;
        }
        const double wmol = n_mouth * p.vbar / 4.0 * mouthArea * dt / (double)M;
        const long long blockSize = 64;
        const long long nBlocks = (M + blockSize - 1) / blockSize;
        std::atomic<long long> nextBlock{0};
        const int T = (int)std::min<long long>(threads, nBlocks);
        for (int t = 0; t < threads; ++t) {
            std::fill(bufs[t].hits.begin(), bufs[t].hits.end(), 0u);
            std::fill(bufs[t].sticks.begin(), bufs[t].sticks.end(), 0u);
            bufs[t].counts = Counts();
        }
        const uint64_t stepSeedIndex = (uint64_t)stepIndex;
        const size_t nMouth = mouthVoxels.size();
        auto work = [&](int tid) {
            ThreadBuf& buf = bufs[tid];
            Rng rng;
            long long b;
            while ((b = nextBlock.fetch_add(1)) < nBlocks) {
                rng.seed(stream_seed(p.seed, stepSeedIndex, (uint64_t)b));
                const long long p0 = b * blockSize, p1 = std::min(M, p0 + blockSize);
                for (long long q = p0; q < p1; ++q) {
                    // uniform mouth voxel, uniform position in it, cosine direction about +z
                    const size_t m = std::min(nMouth - 1, (size_t)(rng.uniform() * (double)nMouth));
                    const int64_t v = mouthVoxels[m];
                    const int i = (int)(v % nx), j = (int)((v / nx) % ny);
                    const double px = i + rng.uniform(), py = j + rng.uniform(), pz = 0.0;
                    double d[3];
                    cosine_direction(rng, 2, 1.0, d);
                    const int fate = trace(rng, buf, i, j, 0, px, py, pz, d);
                    switch (fate) {
                        case FATE_ESCAPED: ++buf.counts.escaped; break;
                        case FATE_STUCK: ++buf.counts.stuck; break;
                        case FATE_LOST: ++buf.counts.lost; break;
                        default: ++buf.counts.left; break;
                    }
                }
            }
        };
        std::vector<std::thread> pool;
        for (int t = 1; t < T; ++t) pool.emplace_back(work, t);
        work(0);
        for (auto& th : pool) th.join();
        // merge (integer sums: independent of the thread count and of scheduling)
        for (int t = 0; t < T; ++t) {
            const ThreadBuf& buf = bufs[t];
            for (int f = 0; f < nFaces; ++f) { lastHits[f] += buf.hits[f]; lastSticks[f] += buf.sticks[f]; }
            last.add(buf.counts);
        }
        cum.add(last);
        // coverage update
        const double perHit = p.s * wmol / (p.N0 * h * h);
        if (p.update == ALDMC_UPDATE_DEPOSIT) {
            for (int f = 0; f < nFaces; ++f)
                if (lastSticks[f]) theta[f] = std::min(1.0, theta[f] + (double)lastSticks[f] * wmol / (p.N0 * h * h * faceW[f]));
        } else {
            for (int f = 0; f < nFaces; ++f)
                if (lastHits[f]) theta[f] = 1.0 - (1.0 - theta[f]) * std::exp(-perHit * (double)lastHits[f]);
        }
        time += dt;
        ++stepIndex;
        return ALDMC_OK;
    }
};

// ---------------------------------------------------------------------------------------------------------------
namespace {
int check_params(const aldmc_params& p, std::string& msg) {
    if (!(p.s >= 0.0 && p.s <= 1.0)) { msg = "s must be in [0, 1]"; return ALDMC_ERR_PARAM; }
    if (!(p.r >= 0.0 && p.r <= 1.0)) { msg = "r must be in [0, 1]"; return ALDMC_ERR_PARAM; }
    if (p.s + p.r > 1.0) { msg = "s + r must not exceed 1"; return ALDMC_ERR_PARAM; }
    if (!(p.N0 > 0.0)) { msg = "N0 must be > 0"; return ALDMC_ERR_PARAM; }
    if (!(p.vbar > 0.0)) { msg = "vbar must be > 0"; return ALDMC_ERR_PARAM; }
    if (p.threads < 0) { msg = "threads must be >= 0"; return ALDMC_ERR_PARAM; }
    if (p.update != ALDMC_UPDATE_EXPONENTIAL && p.update != ALDMC_UPDATE_DEPOSIT) { msg = "update must be 0 or 1"; return ALDMC_ERR_PARAM; }
    if (p.normalRadius < 0 || p.normalRadius > 8) { msg = "normalRadius must be in [0, 8]"; return ALDMC_ERR_PARAM; }
    if (p.emission != ALDMC_EMISSION_SMOOTH && p.emission != ALDMC_EMISSION_FACET) { msg = "emission must be 0 or 1"; return ALDMC_ERR_PARAM; }
    if (p.maxSegments < 0) { msg = "maxSegments must be >= 0"; return ALDMC_ERR_PARAM; }
    return ALDMC_OK;
}
}  // namespace

extern "C" {

ALDMC_API const char* aldmc_version(void) { return "1.0.0"; }

ALDMC_API void aldmc_params_default(aldmc_params* p) {
    if (!p) return;
    p->s = 1.0; p->r = 0.0; p->N0 = 4e18; p->vbar = 410.0; p->seed = 1ULL; p->threads = 0;
    p->periodicX = 0; p->periodicY = 0; p->areaCorrection = 1; p->update = ALDMC_UPDATE_EXPONENTIAL;
    p->normalRadius = 0; p->emission = ALDMC_EMISSION_FACET; p->maxSegments = 0;
}

ALDMC_API aldmc_t* aldmc_create(int nx, int ny, int nz, double h, const unsigned char* occupancy, const aldmc_params* p) {
    try {
        if (!occupancy || !p) { g_create_error = "aldmc_create: null pointer"; return nullptr; }
        if (nx < 1 || ny < 1 || nz < 1) { g_create_error = "aldmc_create: grid dimensions must be >= 1"; return nullptr; }
        if ((long long)nx * ny * nz > 2000000000LL) { g_create_error = "aldmc_create: grid too large"; return nullptr; }
        if (!(h > 0.0)) { g_create_error = "aldmc_create: voxel size must be > 0"; return nullptr; }
        std::string msg;
        if (check_params(*p, msg) != ALDMC_OK) { g_create_error = "aldmc_create: " + msg; return nullptr; }
        aldmc_t* s = new aldmc_t();
        s->nx = nx; s->ny = ny; s->nz = nz; s->h = h; s->p = *p;
        s->threads = p->threads > 0 ? p->threads : (int)std::max(1u, std::thread::hardware_concurrency());
        s->maxSegments = p->maxSegments > 0 ? p->maxSegments : 100000000LL;
        s->build(occupancy);
        return s;
    } catch (const std::bad_alloc&) {
        g_create_error = "aldmc_create: allocation failed"; return nullptr;
    } catch (const std::exception& e) {
        g_create_error = std::string("aldmc_create: ") + e.what(); return nullptr;
    } catch (...) {
        g_create_error = "aldmc_create: unknown error"; return nullptr;
    }
}

ALDMC_API void aldmc_destroy(aldmc_t* h) { delete h; }

ALDMC_API const char* aldmc_last_error(const aldmc_t* h) { return h ? h->err.c_str() : g_create_error.c_str(); }

ALDMC_API int aldmc_get_params(const aldmc_t* h, aldmc_params* p) {
    if (!h || !p) return ALDMC_ERR_ARG;
    *p = h->p; p->threads = h->threads; return ALDMC_OK;
}

ALDMC_API int aldmc_step(aldmc_t* h, double n_mouth, double dt, long long M) {
    if (!h) return ALDMC_ERR_ARG;
    try { return h->step(n_mouth, dt, M); }
    catch (const std::bad_alloc&) { h->err = "aldmc_step: allocation failed"; return ALDMC_ERR_ALLOC; }
    catch (const std::exception& e) { h->err = std::string("aldmc_step: ") + e.what(); return ALDMC_ERR_INTERNAL; }
    catch (...) { h->err = "aldmc_step: unknown error"; return ALDMC_ERR_INTERNAL; }
}

ALDMC_API int aldmc_reset(aldmc_t* h) {
    if (!h) return ALDMC_ERR_ARG;
    std::fill(h->theta.begin(), h->theta.end(), 0.0);
    std::fill(h->lastHits.begin(), h->lastHits.end(), 0);
    std::fill(h->lastSticks.begin(), h->lastSticks.end(), 0);
    h->cum = Counts(); h->last = Counts(); h->time = 0.0; h->stepIndex = 0;
    return ALDMC_OK;
}

ALDMC_API int aldmc_n_faces(const aldmc_t* h) { return h ? h->nFaces : -1; }
ALDMC_API int aldmc_n_mouth_voxels(const aldmc_t* h) { return h ? (int)h->mouthVoxels.size() : -1; }
ALDMC_API int aldmc_mouth_area(const aldmc_t* h, double* out) { if (!h || !out) return ALDMC_ERR_ARG; *out = h->mouthArea; return ALDMC_OK; }
ALDMC_API int aldmc_site_area(const aldmc_t* h, double* out) { if (!h || !out) return ALDMC_ERR_ARG; *out = h->siteArea; return ALDMC_OK; }
ALDMC_API int aldmc_time(const aldmc_t* h, double* out) { if (!h || !out) return ALDMC_ERR_ARG; *out = h->time; return ALDMC_OK; }

ALDMC_API int aldmc_face_info(const aldmc_t* h, long long* voxel, int* direction, double* weight, double* cx, double* cy, double* cz) {
    if (!h) return ALDMC_ERR_ARG;
    for (int f = 0; f < h->nFaces; ++f) {
        const int64_t v = h->faceVoxel[f];
        const int i = (int)(v % h->nx), j = (int)((v / h->nx) % h->ny), k = (int)(v / ((int64_t)h->nx * h->ny));
        const int d = h->faceDir[f];
        if (voxel) voxel[f] = (long long)v;
        if (direction) direction[f] = d;
        if (weight) weight[f] = h->faceW[f];
        double c[3] = {i + 0.5, j + 0.5, k + 0.5};
        c[d / 2] += (d % 2 == 0) ? 0.5 : -0.5;
        if (cx) cx[f] = c[0] * h->h;
        if (cy) cy[f] = c[1] * h->h;
        if (cz) cz[f] = c[2] * h->h;
    }
    return ALDMC_OK;
}

ALDMC_API int aldmc_face_coverage(const aldmc_t* h, double* out) {
    if (!h || !out) return ALDMC_ERR_ARG;
    std::copy(h->theta.begin(), h->theta.end(), out); return ALDMC_OK;
}

ALDMC_API int aldmc_set_face_coverage(aldmc_t* h, const double* in) {
    if (!h || !in) return ALDMC_ERR_ARG;
    for (int f = 0; f < h->nFaces; ++f) {
        if (!(in[f] >= 0.0 && in[f] <= 1.0)) { h->err = "aldmc_set_face_coverage: values must be in [0, 1]"; return ALDMC_ERR_ARG; }
    }
    std::copy(in, in + h->nFaces, h->theta.begin()); return ALDMC_OK;
}

ALDMC_API int aldmc_face_hits(const aldmc_t* h, unsigned long long* hits, unsigned long long* sticks) {
    if (!h) return ALDMC_ERR_ARG;
    if (hits) std::copy(h->lastHits.begin(), h->lastHits.end(), hits);
    if (sticks) std::copy(h->lastSticks.begin(), h->lastSticks.end(), sticks);
    return ALDMC_OK;
}

ALDMC_API int aldmc_slice_coverage(const aldmc_t* h, double* out, double* weight_out) {
    if (!h || !out) return ALDMC_ERR_ARG;
    std::vector<double> num(h->nz, 0.0), den(h->nz, 0.0);
    for (int f = 0; f < h->nFaces; ++f) { num[h->faceK[f]] += h->faceW[f] * h->theta[f]; den[h->faceK[f]] += h->faceW[f]; }
    for (int k = 0; k < h->nz; ++k) {
        out[k] = den[k] > 0.0 ? num[k] / den[k] : 0.0;
        if (weight_out) weight_out[k] = den[k] * h->h * h->h;
    }
    return ALDMC_OK;
}

ALDMC_API int aldmc_bottom_coverage(const aldmc_t* h, double* out) {
    if (!h || !out) return ALDMC_ERR_ARG;
    int kmax = -1;
    for (int f = 0; f < h->nFaces; ++f) if (h->faceDir[f] == ALDMC_DIR_PZ && h->faceK[f] > kmax) kmax = h->faceK[f];
    if (kmax < 0) { *out = 0.0; return ALDMC_OK; }
    double num = 0.0, den = 0.0;
    for (int f = 0; f < h->nFaces; ++f)
        if (h->faceDir[f] == ALDMC_DIR_PZ && h->faceK[f] == kmax) { num += h->faceW[f] * h->theta[f]; den += h->faceW[f]; }
    *out = den > 0.0 ? num / den : 0.0;
    return ALDMC_OK;
}

ALDMC_API int aldmc_counts(const aldmc_t* h, unsigned long long* escaped, unsigned long long* stuck, unsigned long long* lost, unsigned long long* left) {
    if (!h) return ALDMC_ERR_ARG;
    if (escaped) *escaped = h->cum.escaped; if (stuck) *stuck = h->cum.stuck; if (lost) *lost = h->cum.lost; if (left) *left = h->cum.left;
    return ALDMC_OK;
}

ALDMC_API int aldmc_last_counts(const aldmc_t* h, unsigned long long* escaped, unsigned long long* stuck, unsigned long long* lost, unsigned long long* left) {
    if (!h) return ALDMC_ERR_ARG;
    if (escaped) *escaped = h->last.escaped; if (stuck) *stuck = h->last.stuck; if (lost) *lost = h->last.lost; if (left) *left = h->last.left;
    return ALDMC_OK;
}

ALDMC_API int aldmc_last_segments(const aldmc_t* h, unsigned long long* out) {
    if (!h || !out) return ALDMC_ERR_ARG; *out = h->last.segments; return ALDMC_OK;
}

}  // extern "C"
