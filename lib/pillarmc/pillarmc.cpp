// pillarmc: test-particle Monte Carlo of a single ALD exposure in a lateral high-aspect-ratio cavity with support pillars.
//
// Geometry (exact, no voxels): two parallel plates at z = 0 and z = H; the cavity opens at x = 0 and is closed at x = L;
// optional diffuse side walls at y = 0 and y = W (otherwise the cavity is unbounded in y, which with the hexagonal pillar
// lattice is the periodic cavity); optional cylindrical pillars of radius R spanning the gap on a hexagonal lattice of
// inter-distance a (rows a*sqrt(3)/2 apart in x, alternate rows shifted by a/2 in y).
// Physics, as in lib/aldmc: molecular flow; particles enter the mouth with a cosine distribution and uniform position; on a
// wall they stick with probability s (1 - theta) of that wall element, otherwise re-emit with the cosine law about the local
// normal; particles that return to the mouth are lost. Coverage is frozen within an exposure step and advanced by the
// Langmuir relaxation of the collisions each element received. Plates are binned in x (both plates together); every pillar
// row and each side-wall bin carries its own coverage.
// Exposure steps are spaced so that the front advances roughly uniformly: E_k = E_end (k/N)^2.
//
// usage: pillarmc s E_end_Pa_s N_steps M_per_step H_um L_um  a_um(0 = no pillars) R_um  W_um(0 = unbounded)  seed threads out.csv
// output: one row per step: exposure_Pa_s, then the plate coverage per x-bin (bin width = H).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <thread>
#include <algorithm>
struct Rng { uint64_t s[4];
  static uint64_t sm(uint64_t &x) { uint64_t z = (x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
  explicit Rng(uint64_t seed) { for (auto &v : s) v = sm(seed); }
  static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  inline uint64_t next() { uint64_t r = rotl(s[0] + s[3], 23) + s[0], t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
  inline double u() { return (next() >> 11) * (1.0 / 9007199254740992.0); } };
struct Geo { double H, L, a, R, W, rowdx, x1; int nrows; bool pillars, walls; };
static const double PI = 3.14159265358979323846, INF = 1e300;
// cosine-law direction about unit normal n with tangents t1, t2
static inline void cosdir(Rng &g, const double n[3], const double t1[3], const double t2[3], double u[3]) {
  double r = g.u(), ph = 2 * PI * g.u(), ct = std::sqrt(r), st = std::sqrt(1 - r), a = st * std::cos(ph), b = st * std::sin(ph);
  for (int i = 0; i < 3; i++) u[i] = ct * n[i] + a * t1[i] + b * t2[i]; }
// first pillar hit along the xy projection within 3D path length tmax; returns t (3D) and the pillar centre and row
static inline double pillar_hit(const Geo &G, double x, double y, const double u[3], double tmax, double &cx, double &cy, int &row) {
  if (!G.pillars) return INF;
  double best = INF; double xe = x + u[0] * tmax; double xlo = std::min(x, xe) - G.R, xhi = std::max(x, xe) + G.R;
  int k0 = (int)std::ceil((xlo - G.x1) / G.rowdx), k1 = (int)std::floor((xhi - G.x1) / G.rowdx); if (k0 < 0) k0 = 0; if (k1 >= G.nrows) k1 = G.nrows - 1;
  if (k0 > k1) return INF;
  double uxy2 = u[0] * u[0] + u[1] * u[1]; if (uxy2 < 1e-300) return INF;
  int step = u[0] >= 0 ? 1 : -1; int ka = step > 0 ? k0 : k1, kb = step > 0 ? k1 : k0;
  for (int k = ka; step > 0 ? k <= kb : k >= kb; k += step) {
    double xk = G.x1 + k * G.rowdx, off = (k & 1) ? 0.5 * G.a : 0.0;
    // y range of the path while |x - xk| <= R
    double ta, tb; if (std::fabs(u[0]) > 1e-14) { ta = (xk - G.R - x) / u[0]; tb = (xk + G.R - x) / u[0]; if (ta > tb) std::swap(ta, tb); } else { ta = 0; tb = tmax; }
    if (ta < 0) ta = 0; if (tb > tmax) tb = tmax; if (ta > tb) continue;
    double ya = y + u[1] * ta, yb = y + u[1] * tb; if (ya > yb) std::swap(ya, yb);
    long j0 = (long)std::floor((ya - G.R - off) / G.a), j1 = (long)std::ceil((yb + G.R - off) / G.a);
    if (G.walls) { j0 = std::max(j0, 0L); j1 = std::min(j1, (long)std::floor((G.W - off) / G.a)); }
    if (j1 - j0 > 4000) j1 = j0 + 4000;
    for (long j = j0; j <= j1; j++) {
      double pcx = xk, pcy = off + j * G.a; if (G.walls && (pcy < G.R * 1.5 || pcy > G.W - G.R * 1.5)) continue;
      double ox = x - pcx, oy = y - pcy, b = ox * u[0] + oy * u[1], c = ox * ox + oy * oy - G.R * G.R, disc = b * b - uxy2 * c; if (disc <= 0) continue;
      double sq = std::sqrt(disc), t = (-b - sq) / uxy2; if (t <= 1e-9) continue;   // outside-in hit only (flights start outside or on a pillar moving outward)
      if (t < best && t <= tmax) { best = t; cx = pcx; cy = pcy; row = k; } }
    if (best < INF) break; }
  return best; }
struct Tally { std::vector<double> plate, wall, pil; };
static void worker(const Geo &G, double s, long M, uint64_t seed, const std::vector<double> &thP, const std::vector<double> &thW, const std::vector<double> &thPil, double &thEnd_hits, Tally &T, double ymouth) {
  Rng g(seed); int nb = (int)thP.size(); double dxb = G.H; thEnd_hits = 0;
  static const double nxp[3] = {1, 0, 0}, ey[3] = {0, 1, 0}, ez[3] = {0, 0, 1}, nxm[3] = {-1, 0, 0}, nzp[3] = {0, 0, 1}, nzm[3] = {0, 0, -1}, ex[3] = {1, 0, 0}, nyp[3] = {0, 1, 0}, nym[3] = {0, -1, 0};
  for (long m = 0; m < M; m++) {
    double p[3] = {0, ymouth * g.u(), G.H * g.u()}, u[3]; cosdir(g, nxp, ey, ez, u);
    for (long bounce = 0; bounce < 200000000L; bounce++) {
      double tz = u[2] > 0 ? (G.H - p[2]) / u[2] : (u[2] < 0 ? -p[2] / u[2] : INF);
      double tx0 = u[0] < 0 ? -p[0] / u[0] : INF, txL = u[0] > 0 ? (G.L - p[0]) / u[0] : INF;
      double ty = INF; if (G.walls) ty = u[1] > 0 ? (G.W - p[1]) / u[1] : (u[1] < 0 ? -p[1] / u[1] : INF);
      double tmax = std::min(std::min(tz, tx0), std::min(txL, ty));
      double cx = 0, cy = 0; int row = -1; double tp = pillar_hit(G, p[0], p[1], u, tmax, cx, cy, row);
      int surf; double t;                                   // 0 plate, 1 mouth, 2 end wall, 3 side wall, 4 pillar
      if (tp < tmax) { surf = 4; t = tp; } else { t = tmax; surf = (tmax == tz) ? 0 : (tmax == tx0) ? 1 : (tmax == txL) ? 2 : 3; }
      for (int i = 0; i < 3; i++) p[i] += u[i] * t;
      if (surf == 1) break;                                  // escaped through the mouth
      double theta; int ib = std::min(nb - 1, std::max(0, (int)(p[0] / dxb)));
      if (surf == 0) { theta = thP[ib]; T.plate[ib] += 1; }
      else if (surf == 3) { theta = thW[ib]; T.wall[ib] += 1; }
      else if (surf == 4) { theta = thPil[row]; T.pil[row] += 1; }
      else { theta = 1.0; thEnd_hits += 1; }                 // end wall: inert
      if (g.u() < s * (1 - theta)) break;                    // stuck
      if (surf == 0) { if (p[2] > 0.5 * G.H) { p[2] = G.H; cosdir(g, nzm, ex, ey, u); } else { p[2] = 0; cosdir(g, nzp, ex, ey, u); } }
      else if (surf == 2) { p[0] = G.L; cosdir(g, nxm, ey, ez, u); }
      else if (surf == 3) { if (p[1] > 0.5 * G.W) { p[1] = G.W; cosdir(g, nym, ex, ez, u); } else { p[1] = 0; cosdir(g, nyp, ex, ez, u); } }
      else { double n[3] = {(p[0] - cx) / G.R, (p[1] - cy) / G.R, 0}, t1[3] = {-n[1], n[0], 0}; p[0] = cx + n[0] * G.R * (1 + 1e-12); p[1] = cy + n[1] * G.R * (1 + 1e-12); cosdir(g, n, t1, ez, u); }
    } } }
int main(int argc, char **argv) {
  if (argc < 13) { std::fprintf(stderr, "usage: pillarmc s E_end N M H_um L_um a_um R_um W_um seed threads out.csv\n"); return 2; }
  double s = atof(argv[1]), Eend = atof(argv[2]); int N = atoi(argv[3]); long M = atol(argv[4]);
  Geo G; G.H = atof(argv[5]) * 1e-6; G.L = atof(argv[6]) * 1e-6; G.a = atof(argv[7]) * 1e-6; G.R = atof(argv[8]) * 1e-6; G.W = atof(argv[9]) * 1e-6;
  uint64_t seed = strtoull(argv[10], 0, 10); int nt = atoi(argv[11]); const char *out = argv[12];
  G.pillars = G.a > 0; G.walls = G.W > 0; G.rowdx = G.a * std::sqrt(3.0) / 2; G.x1 = G.rowdx; G.nrows = G.pillars ? (int)std::floor((G.L - G.x1) / G.rowdx) + 1 : 0;
  const double N0 = 4e18, J1 = 1.0 / std::sqrt(2 * PI * 72.09 * 1.66053906660e-27 * 1.380649e-23 * 573.15);   // molecules / m2 / s / Pa
  int nb = (int)std::llround(G.L / G.H); std::vector<double> thP(nb, 0), thW(nb, 0), thPil(std::max(1, G.nrows), 0);
  double Wy = G.walls ? G.W : (G.pillars ? G.a : G.H);        // mouth width represented by the launched particles
  double Aplate = 2 * G.H * Wy, Awall = 2 * G.H * G.H;         // per x-bin: both plates; both side walls
  double npil_per_row = G.pillars ? (G.walls ? std::floor(G.W / G.a) : 1.0) : 0, Apil = npil_per_row * 2 * PI * G.R * G.H;
  FILE *f = std::fopen(out, "w"); if (!f) return 3;
  std::fprintf(f, "# s=%g H=%g L=%g a=%g R=%g W=%g M=%ld N=%d seed=%llu\n", s, G.H, G.L, G.a, G.R, G.W, M, N, (unsigned long long)seed);
  double Eprev = 0;
  for (int k = 1; k <= N; k++) {
    double E = Eend * (double(k) / N) * (double(k) / N), dE = E - Eprev; Eprev = E;
    double w = J1 * dE * G.H * Wy / double(M);                // real molecules per test particle
    std::vector<Tally> T(nt); std::vector<double> endh(nt, 0); std::vector<std::thread> th;
    for (int i = 0; i < nt; i++) { T[i].plate.assign(nb, 0); T[i].wall.assign(nb, 0); T[i].pil.assign(thPil.size(), 0);
      th.emplace_back(worker, std::cref(G), s, M / nt, seed * 1000003ULL + (uint64_t)k * 7919ULL + i, std::cref(thP), std::cref(thW), std::cref(thPil), std::ref(endh[i]), std::ref(T[i]), Wy); }
    for (auto &t : th) t.join();
    double scale = double(M) / double((M / nt) * nt);
    for (int b = 0; b < nb; b++) { double c = 0, cw = 0; for (int i = 0; i < nt; i++) { c += T[i].plate[b]; cw += T[i].wall[b]; }
      thP[b] = 1 - (1 - thP[b]) * std::exp(-s * w * scale * c / (N0 * Aplate)); if (G.walls) thW[b] = 1 - (1 - thW[b]) * std::exp(-s * w * scale * cw / (N0 * Awall)); }
    for (size_t r = 0; r < thPil.size() && G.pillars; r++) { double c = 0; for (int i = 0; i < nt; i++) c += T[i].pil[r]; thPil[r] = 1 - (1 - thPil[r]) * std::exp(-s * w * scale * c / (N0 * Apil)); }
    std::fprintf(f, "%.6g", E); for (int b = 0; b < nb; b++) std::fprintf(f, ",%.5f", thP[b]); std::fprintf(f, "\n");
    if (k % 20 == 0 || k == N) { int i50 = 0; while (i50 < nb && thP[i50] > 0.5) i50++; std::fprintf(stderr, "step %d/%d  E %.3f Pa s  z50 ~ %d H\n", k, N, E, i50); }
  }
  std::fclose(f); return 0; }
