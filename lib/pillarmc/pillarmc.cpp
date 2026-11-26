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
// Wall scattering (optional trailing arguments; defaults give the cosine law used everywhere else):
//   n_cos    re-emission lobe cos^n(theta) about the local normal (1 = cosine law). n > 1 is the phenomenological
//            'directional re-emission' of feature-scale simulators; it is not reciprocal and so not an equilibrium law.
//   f_spec   fraction of non-sticking collisions reflected specularly (Maxwell's model; reciprocal).
//   f_retro  fraction reflected straight back along the incoming direction (back-scattering by roughness; reciprocal).
// Gas-phase collisions (optional): lambda_um is the mean free path of the precursor in the carrier gas (0 = none, free
//   molecular). A flight ends after an exponentially distributed length with that mean, and the molecule leaves the collision
//   in a uniformly random direction (a test particle in an equilibrium background gas; velocity persistence is neglected, so
//   lambda is the momentum-transfer mean free path). The flux entering the mouth is unchanged (equilibrium reservoir).
//   p=<Pa> switches on the physical model instead: hard-sphere collisions of the precursor (mass m1, the test particle now
//   carries its speed) with a Maxwellian carrier gas (mass m2, total pressure p, temperature T, collision diameter d12).
//   The collision rate depends on the molecule's speed, the partner velocity is drawn from the collision-weighted
//   Maxwellian, and scattering is isotropic in the centre-of-mass frame, so the persistence of velocity of a heavy
//   molecule in a light gas is kept. Walls re-emit (and the mouth admits) with the flux-weighted Maxwell speed at T.
//   Defaults: T=573.15 K, m1=72.09 (TMA monomer), m2=28.0134 (N2), d12=4.825 A (TMA 5.91 A, N2 3.74 A; Ylilammi 2018).
// usage: pillarmc s E_end_Pa_s N_steps M_per_step H_um L_um  a_um(0 = no pillars) R_um  W_um(0 = unbounded)  seed threads out.csv [n_cos f_spec f_retro lambda_um] [key=value ...]
//   keys: ncos fspec fretro lambda(um) p(Pa) T(K) m1(amu) m2(amu) d12(Angstrom) tfile(path)
// Film already on the plates (channel narrowing): tfile=<path> gives one film thickness in um per x-bin (bin = one initial gap
//   height H); the lower plate surface is then z = t(x) and the upper z = H - t(x), linear between bin centres, with the cosine
//   law about the local (slightly tilted) normal. The mouth admits molecules over its open height H - 2 t(0). Pillars keep
//   their radius. Used by docs/paper/multicycle_mc.py to follow hundreds of cycles in groups.
// output: one row per step: exposure_Pa_s, then the plate coverage per x-bin (bin width = H).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <thread>
#include <algorithm>
#include <string>
struct Rng { uint64_t s[4];
  static uint64_t sm(uint64_t &x) { uint64_t z = (x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
  explicit Rng(uint64_t seed) { for (auto &v : s) v = sm(seed); }
  static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  inline uint64_t next() { uint64_t r = rotl(s[0] + s[3], 23) + s[0], t = s[1] << 17; s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return r; }
  inline double u() { return (next() >> 11) * (1.0 / 9007199254740992.0); } };
struct Geo { double H, L, a, R, W, rowdx, x1; int nrows; bool pillars, walls; };
static const double PI = 3.14159265358979323846, INF = 1e300;
static double NCOS = 1.0, FSPEC = 0.0, FRETRO = 0.0, LAMBDA = 0.0;
static double PGAS = 0.0, TEMP = 573.15, M1 = 72.09, M2 = 28.0134, D12 = 4.825;       // Pa, K, amu, amu, Angstrom
static double C1 = 0, C2 = 0, NSIGMA = 0;                                               // most probable speeds (m/s), n*sigma (1/m)
static const double KB = 1.380649e-23, AMU = 1.66053906660e-27;
static std::vector<double> TF;                                                        // film thickness per bin (m); empty = bare plates
// first hit of the ray p + u t on the film-covered plates within tlim; top = which plate, slope = d t_film / dx there
static inline double plate_hit(double H, const double p[3], const double u[3], double tlim, bool &top, double &slope) {
  const int n = (int)TF.size(); const double dx = H; double tau0 = 0; int j = (int)std::floor(p[0] / dx - 0.5);
  for (int it = 0; it < 6000; it++) {
    if (j < -1) j = -1; if (j > n - 1) j = n - 1;
    const bool flat = (j < 0 || j >= n - 1); const double xa = j < 0 ? -1e300 : (j + 0.5) * dx, xb = j >= n - 1 ? 1e300 : (j + 1.5) * dx;
    const double Tj = TF[std::max(0, std::min(n - 1, j))], m = flat ? 0.0 : (TF[j + 1] - TF[j]) / dx, xr = flat ? p[0] : (j + 0.5) * dx;
    double te = 1e300; if (u[0] > 1e-300) te = (xb - p[0]) / u[0]; else if (u[0] < -1e-300) te = (xa - p[0]) / u[0];
    double best = 1e300; bool bt = false; const double dl = u[2] - m * u[0], du = u[2] + m * u[0];
    if (dl < 0) { double t = (Tj + m * (p[0] - xr) - p[2]) / dl; if (t >= tau0 - 1e-14 && t <= te) { best = t; bt = false; } }
    if (du > 0) { double t = (H - Tj - m * (p[0] - xr) - p[2]) / du; if (t >= tau0 - 1e-14 && t <= te && t < best) { best = t; bt = true; } }
    if (best < 1e299) { top = bt; slope = m; return best > 0 ? best : 0.0; }
    if (te >= tlim || te > 1e299) return 1e300;
    tau0 = te; j += u[0] > 0 ? 1 : -1; }
  return 1e300; }
static inline double filmAt(double H, double x) { const int n = (int)TF.size(); double q = x / H - 0.5; if (q <= 0) return TF[0]; if (q >= n - 1) return TF[n - 1]; int j = (int)q; return TF[j] + (TF[j + 1] - TF[j]) * (q - j); }
static inline double gauss(Rng &g) { double a = std::sqrt(-2 * std::log(1.0 - g.u())); return a * std::cos(2 * PI * g.u()); }
// speed of a molecule leaving a wall or entering through the mouth: flux-weighted Maxwellian, pdf ~ v^3 exp(-v^2/c^2)
static inline double fluxspeed(Rng &g, double c) { return c * std::sqrt(-std::log((1.0 - g.u()) * (1.0 - g.u()))); }
// hard-sphere collision rate per unit path length for a molecule of speed v in the Maxwellian carrier gas
static inline double inv_path(double v) { double x = v / C2; return NSIGMA * (std::exp(-x * x) / std::sqrt(PI) + (x + 0.5 / x) * std::erf(x)) / x; }
// one collision: partner from the collision-weighted Maxwellian, isotropic scattering in the centre-of-mass frame
static inline void gas_collide(Rng &g, double u[3], double &v) {
  double w[3], g3[3], gr, gmax = v + 6 * C2;
  do { for (int i = 0; i < 3; i++) { w[i] = C2 * gauss(g) * 0.7071067811865476; g3[i] = v * u[i] - w[i]; } gr = std::sqrt(g3[0] * g3[0] + g3[1] * g3[1] + g3[2] * g3[2]); } while (g.u() * gmax > gr);
  double c = 2 * g.u() - 1, sn = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * PI * g.u(), e[3] = {sn * std::cos(ph), sn * std::sin(ph), c}, f1 = M1 / (M1 + M2), f2 = M2 / (M1 + M2), vn[3], v2 = 0;
  for (int i = 0; i < 3; i++) { vn[i] = f1 * v * u[i] + f2 * w[i] + f2 * gr * e[i]; v2 += vn[i] * vn[i]; }
  v = std::sqrt(v2); for (int i = 0; i < 3; i++) u[i] = vn[i] / v; }
// cosine-law direction about unit normal n with tangents t1, t2
static inline void cosdir(Rng &g, const double n[3], const double t1[3], const double t2[3], double u[3]) {
  double r = g.u(), ph = 2 * PI * g.u(), ct = std::sqrt(r), st = std::sqrt(1 - r), a = st * std::cos(ph), b = st * std::sin(ph);
  for (int i = 0; i < 3; i++) u[i] = ct * n[i] + a * t1[i] + b * t2[i]; }
// re-emission from a wall with unit normal n (pointing into the gas): specular, retro-reflected or a cos^n lobe
static inline void reemit(Rng &g, const double n[3], const double t1[3], const double t2[3], double u[3], double &v) {
  if (FSPEC > 0 || FRETRO > 0) { double q = g.u();
    if (q < FSPEC) { double d = u[0] * n[0] + u[1] * n[1] + u[2] * n[2]; for (int i = 0; i < 3; i++) u[i] -= 2 * d * n[i]; return; }
    if (q < FSPEC + FRETRO) { for (int i = 0; i < 3; i++) u[i] = -u[i]; return; } }
  if (PGAS > 0) v = fluxspeed(g, C1);                       // diffuse re-emission: fully accommodated speed
  if (NCOS == 1.0) { cosdir(g, n, t1, t2, u); return; }
  double ct = std::pow(g.u(), 1.0 / (NCOS + 1.0)), st = std::sqrt(std::max(0.0, 1 - ct * ct)), ph = 2 * PI * g.u(), a = st * std::cos(ph), b = st * std::sin(ph);
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
    const bool film = !TF.empty(); double p[3] = {0, ymouth * g.u(), film ? TF[0] + (G.H - 2 * TF[0]) * g.u() : G.H * g.u()}, u[3]; cosdir(g, nxp, ey, ez, u); double v = PGAS > 0 ? fluxspeed(g, C1) : 0;
    for (long bounce = 0; bounce < 200000000L; bounce++) {
      double tx0 = u[0] < 0 ? -p[0] / u[0] : INF, txL = u[0] > 0 ? (G.L - p[0]) / u[0] : INF;
      double ty = INF; if (G.walls) ty = u[1] > 0 ? (G.W - p[1]) / u[1] : (u[1] < 0 ? -p[1] / u[1] : INF);
      bool ptop = false; double pm = 0;
      double tz = film ? plate_hit(G.H, p, u, std::min(std::min(tx0, txL), ty), ptop, pm) : (u[2] > 0 ? (G.H - p[2]) / u[2] : (u[2] < 0 ? -p[2] / u[2] : INF));
      double tmax = std::min(std::min(tz, tx0), std::min(txL, ty));
      double cx = 0, cy = 0; int row = -1; double tp = pillar_hit(G, p[0], p[1], u, tmax, cx, cy, row);
      int surf; double t;                                   // 0 plate, 1 mouth, 2 end wall, 3 side wall, 4 pillar
      if (tp < tmax) { surf = 4; t = tp; } else { t = tmax; surf = (tmax == tz) ? 0 : (tmax == tx0) ? 1 : (tmax == txL) ? 2 : 3; }
      if (PGAS > 0) { double lf = -std::log(1.0 - g.u()) / inv_path(v);
        if (lf < t) { for (int i = 0; i < 3; i++) p[i] += u[i] * lf; gas_collide(g, u, v); continue; } }
      if (LAMBDA > 0) { double lf = -LAMBDA * std::log(1.0 - g.u());
        if (lf < t) { for (int i = 0; i < 3; i++) p[i] += u[i] * lf;
          double c = 2 * g.u() - 1, sn = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * PI * g.u(); u[0] = sn * std::cos(ph); u[1] = sn * std::sin(ph); u[2] = c; continue; } }
      for (int i = 0; i < 3; i++) p[i] += u[i] * t;
      if (surf == 1) break;                                  // escaped through the mouth
      double theta; int ib = std::min(nb - 1, std::max(0, (int)(p[0] / dxb)));
      if (surf == 0) { theta = thP[ib]; T.plate[ib] += 1; }
      else if (surf == 3) { theta = thW[ib]; T.wall[ib] += 1; }
      else if (surf == 4) { theta = thPil[row]; T.pil[row] += 1; }
      else { theta = 1.0; thEnd_hits += 1; }                 // end wall: inert
      if (g.u() < s * (1 - theta)) break;                    // stuck
      if (surf == 0 && film) { const double tf = filmAt(G.H, p[0]), nn = 1.0 / std::sqrt(1 + pm * pm);
        if (ptop) { p[2] = G.H - tf - 1e-10 * G.H; const double n[3] = {-pm * nn, 0, -nn}, t2[3] = {nn, 0, -pm * nn}; reemit(g, n, ey, t2, u, v); }
        else { p[2] = tf + 1e-10 * G.H; const double n[3] = {-pm * nn, 0, nn}, t2[3] = {nn, 0, pm * nn}; reemit(g, n, ey, t2, u, v); } }
      else if (surf == 0) { if (p[2] > 0.5 * G.H) { p[2] = G.H; reemit(g, nzm, ex, ey, u, v); } else { p[2] = 0; reemit(g, nzp, ex, ey, u, v); } }
      else if (surf == 2) { p[0] = G.L; reemit(g, nxm, ey, ez, u, v); }
      else if (surf == 3) { if (p[1] > 0.5 * G.W) { p[1] = G.W; reemit(g, nym, ex, ez, u, v); } else { p[1] = 0; reemit(g, nyp, ex, ez, u, v); } }
      else { double n[3] = {(p[0] - cx) / G.R, (p[1] - cy) / G.R, 0}, t1[3] = {-n[1], n[0], 0}; p[0] = cx + n[0] * G.R * (1 + 1e-12); p[1] = cy + n[1] * G.R * (1 + 1e-12); reemit(g, n, t1, ez, u, v); }
    } } }
// self-test of the gas model: tracer diffusion coefficient of the precursor in the unbounded carrier gas (m2/s), from the
// mean-square displacement of molecules that start with a Maxwellian velocity. usage: pillarmc diffusion p_Pa [n_mol] [T m1 m2 d12]
static int diffusion_test(int argc, char **argv) {
  PGAS = atof(argv[2]); long n = argc > 3 ? atol(argv[3]) : 20000; if (argc > 4) TEMP = atof(argv[4]); if (argc > 5) M1 = atof(argv[5]); if (argc > 6) M2 = atof(argv[6]); if (argc > 7) D12 = atof(argv[7]);
  C1 = std::sqrt(2 * KB * TEMP / (M1 * AMU)); C2 = std::sqrt(2 * KB * TEMP / (M2 * AMU)); NSIGMA = PGAS / (KB * TEMP) * PI * D12 * D12 * 1e-20;
  Rng g(12345); double vbar = 2 * C1 / std::sqrt(PI), tend = 400.0 / (inv_path(vbar) * vbar), r2 = 0, v2sum = 0;
  for (long m = 0; m < n; m++) { double w[3], v = 0; for (int i = 0; i < 3; i++) { w[i] = C1 * gauss(g) * 0.7071067811865476; v += w[i] * w[i]; } v = std::sqrt(v);
    double u[3] = {w[0] / v, w[1] / v, w[2] / v}, p[3] = {0, 0, 0}, t = 0;
    for (;;) { double lf = -std::log(1.0 - g.u()) / inv_path(v), dt = lf / v; if (t + dt >= tend) { for (int i = 0; i < 3; i++) p[i] += u[i] * v * (tend - t); break; }
      for (int i = 0; i < 3; i++) p[i] += u[i] * lf; t += dt; gas_collide(g, u, v); }
    r2 += p[0] * p[0] + p[1] * p[1] + p[2] * p[2]; v2sum += v * v; }
  double mu = M1 * M2 / (M1 + M2) * AMU, nd = PGAS / (KB * TEMP), ce = 3.0 / (8 * nd * D12 * D12 * 1e-20) * std::sqrt(KB * TEMP / (2 * PI * mu));
  std::printf("D_mc %.6g  D_chapman_enskog_1 %.6g  ratio %.4f  mean_v2_end/(1.5 c1^2) %.4f  mean_free_path_at_vbar_um %.4g\n", r2 / n / (6 * tend), ce, r2 / n / (6 * tend) / ce, v2sum / n / (1.5 * C1 * C1), 1e6 / inv_path(vbar));
  return 0; }
int main(int argc, char **argv) {
  if (argc > 2 && std::string(argv[1]) == "diffusion") return diffusion_test(argc, argv);
  if (argc < 13) { std::fprintf(stderr, "usage: pillarmc s E_end N M H_um L_um a_um R_um W_um seed threads out.csv\n"); return 2; }
  double s = atof(argv[1]), Eend = atof(argv[2]); int N = atoi(argv[3]); long M = atol(argv[4]);
  Geo G; G.H = atof(argv[5]) * 1e-6; G.L = atof(argv[6]) * 1e-6; G.a = atof(argv[7]) * 1e-6; G.R = atof(argv[8]) * 1e-6; G.W = atof(argv[9]) * 1e-6;
  uint64_t seed = strtoull(argv[10], 0, 10); int nt = atoi(argv[11]); const char *out = argv[12];
  { double *pos[4] = {&NCOS, &FSPEC, &FRETRO, &LAMBDA}; const char *keys[9] = {"ncos", "fspec", "fretro", "lambda", "p", "T", "m1", "m2", "d12"}; double *val[9] = {&NCOS, &FSPEC, &FRETRO, &LAMBDA, &PGAS, &TEMP, &M1, &M2, &D12}; int np = 0;
    for (int i = 13; i < argc; i++) { std::string a(argv[i]); size_t e = a.find('=');
      if (e == std::string::npos) { if (np < 4) *pos[np++] = atof(argv[i]); continue; }
      if (a.substr(0, e) == "tfile") { FILE *tf = std::fopen(a.c_str() + e + 1, "r"); if (!tf) { std::fprintf(stderr, "cannot read %s\n", a.c_str() + e + 1); return 2; } double x; while (std::fscanf(tf, "%lf", &x) == 1) TF.push_back(x * 1e-6); std::fclose(tf); continue; }
      bool ok = false; for (int k = 0; k < 9; k++) if (a.substr(0, e) == keys[k]) { *val[k] = atof(a.c_str() + e + 1); ok = true; }
      if (!ok) { std::fprintf(stderr, "unknown option %s\n", argv[i]); return 2; } }
    LAMBDA *= 1e-6; C1 = std::sqrt(2 * KB * TEMP / (M1 * AMU)); C2 = std::sqrt(2 * KB * TEMP / (M2 * AMU)); NSIGMA = PGAS / (KB * TEMP) * PI * D12 * D12 * 1e-20; }
  G.pillars = G.a > 0; G.walls = G.W > 0; G.rowdx = G.a * std::sqrt(3.0) / 2; G.x1 = G.rowdx; G.nrows = G.pillars ? (int)std::floor((G.L - G.x1) / G.rowdx) + 1 : 0;
  const double N0 = 4e18, J1 = 1.0 / std::sqrt(2 * PI * M1 * AMU * KB * TEMP);   // molecules / m2 / s / Pa
  int nb = (int)std::llround(G.L / G.H); std::vector<double> thP(nb, 0), thW(nb, 0), thPil(std::max(1, G.nrows), 0);
  if (!TF.empty() && (int)TF.size() != nb) { std::fprintf(stderr, "tfile has %d values, the cavity has %d bins\n", (int)TF.size(), nb); return 2; }
  double Wy = G.walls ? G.W : (G.pillars ? G.a : G.H);        // mouth width represented by the launched particles
  double Aplate = 2 * G.H * Wy, Awall = 2 * G.H * G.H;         // per x-bin: both plates; both side walls
  double npil_per_row = G.pillars ? (G.walls ? std::floor(G.W / G.a) : 1.0) : 0, Apil = npil_per_row * 2 * PI * G.R * G.H;
  FILE *f = std::fopen(out, "w"); if (!f) return 3;
  std::fprintf(f, "# s=%g H=%g L=%g a=%g R=%g W=%g M=%ld N=%d seed=%llu n_cos=%g f_spec=%g f_retro=%g lambda=%g p=%g T=%g m1=%g m2=%g d12=%g\n", s, G.H, G.L, G.a, G.R, G.W, M, N, (unsigned long long)seed, NCOS, FSPEC, FRETRO, LAMBDA, PGAS, TEMP, M1, M2, D12);
  double Eprev = 0;
  for (int k = 1; k <= N; k++) {
    double E = Eend * (double(k) / N) * (double(k) / N), dE = E - Eprev; Eprev = E;
    double w = J1 * dE * (TF.empty() ? G.H : G.H - 2 * TF[0]) * Wy / double(M);                // real molecules per test particle
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
