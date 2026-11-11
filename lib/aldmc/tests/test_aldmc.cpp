// ctest suite of libaldmc: one test per command-line name, no framework. Every test is deterministic (fixed seed).
#include "aldmc.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static int g_failures = 0;
#define CHECK(cond) do { if (!(cond)) { std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_failures; } } while (0)
#define CHECK_NEAR(a, b, tol) do { double _a = (a), _b = (b); if (!(std::fabs(_a - _b) <= (tol))) { \
    std::printf("FAIL %s:%d: %s = %.9g vs %s = %.9g (tol %.3g)\n", __FILE__, __LINE__, #a, _a, #b, _b, (double)(tol)); ++g_failures; } } while (0)

typedef unsigned long long u64;

struct Grid {
    int nx, ny, nz;
    std::vector<unsigned char> occ;
    Grid(int x, int y, int z, unsigned char fill = 0) : nx(x), ny(y), nz(z), occ((size_t)x * y * z, fill) {}
    unsigned char& at(int i, int j, int k) { return occ[((size_t)k * ny + j) * nx + i]; }
};

static aldmc_params defaults(double s, double r = 0.0) {
    aldmc_params p; aldmc_params_default(&p);
    p.s = s; p.r = r; p.N0 = 4e18; p.vbar = 400.0; p.seed = 7; p.threads = 1;
    return p;
}

// floor at k = nz - 1, periodic in x and y: every particle hits a floor face once and leaves through the mouth
static Grid flat_floor(int n, int nz) {
    Grid g(n, n, nz);
    for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) g.at(i, j, nz - 1) = 1;
    return g;
}

// slit of gap `gap` voxels in y, periodic in x, depth nz - 1 with a floor
static Grid slit(int nx, int gap, int depth) {
    Grid g(nx, gap + 2, depth + 1, 1);
    for (int k = 0; k < depth; ++k) for (int j = 1; j <= gap; ++j) for (int i = 0; i < nx; ++i) g.at(i, j, k) = 0;
    return g;
}

// cylinder of diameter d voxels, depth L voxels, one voxel of solid padding around and a floor
static Grid hole(int d, int L, int pad = 1) {
    const int n = d + 2 * pad;
    Grid g(n, n, L + 1, 1);
    const double c = n / 2.0, r2 = (d / 2.0) * (d / 2.0);
    for (int k = 0; k < L; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        const double x = i + 0.5 - c, y = j + 0.5 - c;
        if (x * x + y * y < r2) g.at(i, j, k) = 0;
    }
    return g;
}

static int run(const std::string& name) {
    const double h = 5e-9;
    if (name == "version_and_defaults") {
        CHECK(std::strcmp(aldmc_version(), "1.0.0") == 0);
        aldmc_params p; aldmc_params_default(&p);
        CHECK(p.s == 1.0 && p.r == 0.0 && p.N0 == 4e18 && p.vbar == 410.0 && p.areaCorrection == 1);
        CHECK(p.update == ALDMC_UPDATE_EXPONENTIAL && p.periodicX == 0 && p.periodicY == 0);
    } else if (name == "open_box_loses_all") {
        Grid g(8, 8, 8);
        aldmc_params p = defaults(0.5); p.periodicX = p.periodicY = 1;
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        CHECK(aldmc_n_faces(s) == 0);
        CHECK(aldmc_step(s, 1e20, 1e-3, 5000) == ALDMC_OK);
        u64 e, st, lo, le;
        CHECK(aldmc_counts(s, &e, &st, &lo, &le) == ALDMC_OK);
        CHECK(st == 0 && lo == 0 && e + le == 5000 && le > 0);
        double b; CHECK(aldmc_bottom_coverage(s, &b) == ALDMC_OK); CHECK(b == 0.0);
        aldmc_destroy(s);
    } else if (name == "flat_langmuir_analytic" || name == "flat_langmuir_deposit_scheme") {
        // theta(t) = 1 - exp(-s n vbar/4 t / N0) on a flat surface: rate 1/s here, 10 steps of 0.1 s
        const bool deposit = name == "flat_langmuir_deposit_scheme";
        Grid g = flat_floor(8, 3);
        aldmc_params p = defaults(0.5); p.periodicX = p.periodicY = 1; p.threads = 2;
        if (deposit) p.update = ALDMC_UPDATE_DEPOSIT;
        const double rate = 1.0;                                   // s n vbar/4 / N0
        const double n = rate * p.N0 * 4.0 / (p.s * p.vbar);
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        CHECK(aldmc_n_faces(s) == 64);
        double area; aldmc_mouth_area(s, &area); CHECK_NEAR(area, 64 * h * h, 1e-30);
        aldmc_site_area(s, &area); CHECK_NEAR(area, 64 * h * h, 1e-24);   // flat: every weight is 1
        for (int it = 0; it < 10; ++it) CHECK(aldmc_step(s, n, 0.1, 64000) == ALDMC_OK);
        std::vector<double> th(64); aldmc_face_coverage(s, th.data());
        double mean = 0.0; for (double v : th) mean += v / 64.0;
        if (deposit) CHECK_NEAR(mean, 1.0 - std::pow(0.9, 10), 0.01);   // forward Euler of the same rate
        else         CHECK_NEAR(mean, 1.0 - std::exp(-1.0), 0.005);
        u64 e, st, lo, le; aldmc_counts(s, &e, &st, &lo, &le);
        CHECK(le == 0 && lo == 0 && e + st == 640000);
        std::vector<double> sl(3), w(3); aldmc_slice_coverage(s, sl.data(), w.data());
        CHECK_NEAR(sl[1], mean, 1e-12); CHECK(w[0] == 0.0 && w[2] == 0.0);
        double b; aldmc_bottom_coverage(s, &b); CHECK_NEAR(b, mean, 1e-12);
        double t; aldmc_time(s, &t); CHECK_NEAR(t, 1.0, 1e-12);
        aldmc_destroy(s);
    } else if (name == "periodic_slit_conserves") {
        Grid g = slit(4, 10, 40);
        aldmc_params p = defaults(5e-3); p.periodicX = 1; p.threads = 3;
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        for (int it = 0; it < 5; ++it) CHECK(aldmc_step(s, 3e19, 0.05, 20000) == ALDMC_OK);
        u64 e, st, lo, le; aldmc_counts(s, &e, &st, &lo, &le);
        CHECK(e + st + lo + le == 100000);
        CHECK(le == 0 && lo == 0 && st > 0 && e > 0);
        aldmc_destroy(s);
    } else if (name == "recombination_loses_all") {
        Grid g = slit(4, 10, 40);
        aldmc_params p = defaults(0.0, 1.0); p.periodicX = 1; p.areaCorrection = 0;
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        CHECK(aldmc_step(s, 3e19, 0.05, 10000) == ALDMC_OK);
        u64 e, st, lo, le; aldmc_counts(s, &e, &st, &lo, &le);
        CHECK(lo == 10000 && e == 0 && st == 0 && le == 0);
        aldmc_destroy(s);
    } else if (name == "determinism_and_thread_invariance") {
        Grid g = hole(10, 40);
        std::vector<std::vector<double>> th;
        std::vector<u64> stuck;
        const int threadCounts[3] = {1, 1, 4};
        for (int c = 0; c < 3; ++c) {
            aldmc_params p = defaults(5e-3); p.threads = threadCounts[c]; p.seed = 12345;
            aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
            CHECK(s != nullptr);
            for (int it = 0; it < 6; ++it) CHECK(aldmc_step(s, 3.8e19, 0.05, 30000) == ALDMC_OK);
            std::vector<double> t(aldmc_n_faces(s)); aldmc_face_coverage(s, t.data());
            u64 e, st, lo, le; aldmc_counts(s, &e, &st, &lo, &le);
            th.push_back(t); stuck.push_back(st);
            aldmc_destroy(s);
        }
        CHECK(th[0] == th[1] && stuck[0] == stuck[1]);     // same seed, single thread: bitwise
        CHECK(th[0] == th[2] && stuck[0] == stuck[2]);     // 1 vs 4 threads: bitwise (integer merge, per-block streams)
        // a different seed gives a different sample
        aldmc_params p = defaults(5e-3); p.seed = 999;
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        for (int it = 0; it < 6; ++it) aldmc_step(s, 3.8e19, 0.05, 30000);
        std::vector<double> t(aldmc_n_faces(s)); aldmc_face_coverage(s, t.data());
        CHECK(t != th[0]);
        aldmc_destroy(s);
    } else if (name == "area_weights_cylinder") {
        const int d = 24, L = 48;
        Grid g = hole(d, L, 2);
        const double smooth = M_PI * d * L * h * h + M_PI * d * d / 4.0 * h * h;   // wall + floor
        aldmc_params p = defaults(5e-3);
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        double corrected; aldmc_site_area(s, &corrected);
        aldmc_destroy(s);
        p.areaCorrection = 0;
        s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        double plain; aldmc_site_area(s, &plain);
        int nf = aldmc_n_faces(s);
        std::vector<int> dir(nf); std::vector<double> w(nf);
        aldmc_face_info(s, nullptr, dir.data(), w.data(), nullptr, nullptr, nullptr);
        for (int f = 0; f < nf; ++f) CHECK(w[f] == 1.0);
        aldmc_destroy(s);
        std::printf("smooth %.4g m2, corrected %.4g (%.3f), plain voxel faces %.4g (%.3f; staircase 4/pi = 1.273 on the wall)\n",
                    smooth, corrected, corrected / smooth, plain, plain / smooth);
        CHECK(std::fabs(corrected / smooth - 1.0) < 0.03);
        CHECK(plain / smooth > 1.15);
    } else if (name == "hole_saturates") {
        Grid g = hole(12, 60);
        aldmc_params p = defaults(5e-3); p.threads = 4;
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        double prev = 0.0, b = 0.0, mid_profile_ok = 1.0;
        std::vector<double> sl(g.nz);
        for (int it = 0; it < 60; ++it) {
            CHECK(aldmc_step(s, 3.8e19, 0.05, 40000) == ALDMC_OK);
            aldmc_bottom_coverage(s, &b);
            CHECK(b >= prev - 1e-12 && b <= 1.0);
            prev = b;
            if (it == 10) {
                aldmc_slice_coverage(s, sl.data(), nullptr);
                if (!(sl[0] > sl[30] && sl[30] > sl[59])) mid_profile_ok = 0.0;   // coverage falls with depth
            }
        }
        CHECK(mid_profile_ok == 1.0);
        CHECK(b > 0.9);
        u64 e, st, lo, le; aldmc_counts(s, &e, &st, &lo, &le);
        CHECK(le == 0 && e + st == 60ULL * 40000ULL);
        aldmc_destroy(s);
    } else if (name == "error_codes") {
        Grid g = hole(8, 8);
        aldmc_params p = defaults(5e-3);
        CHECK(aldmc_create(0, 1, 1, h, g.occ.data(), &p) == nullptr);
        CHECK(std::strlen(aldmc_last_error(nullptr)) > 0);
        CHECK(aldmc_create(g.nx, g.ny, g.nz, h, nullptr, &p) == nullptr);
        CHECK(aldmc_create(g.nx, g.ny, g.nz, 0.0, g.occ.data(), &p) == nullptr);
        p.s = 2.0; CHECK(aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p) == nullptr); p.s = 5e-3;
        p.r = 0.999; CHECK(aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p) == nullptr); p.r = 0.0;
        p.update = 7; CHECK(aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p) == nullptr); p.update = 0;
        aldmc_t* s = aldmc_create(g.nx, g.ny, g.nz, h, g.occ.data(), &p);
        CHECK(s != nullptr);
        CHECK(aldmc_step(nullptr, 1.0, 1.0, 1) == ALDMC_ERR_ARG);
        CHECK(aldmc_step(s, -1.0, 1.0, 1) == ALDMC_ERR_ARG);
        CHECK(std::strlen(aldmc_last_error(s)) > 0);
        CHECK(aldmc_step(s, 1e19, 0.0, 100) == ALDMC_OK);        // zero dt: nothing launched, time unchanged
        CHECK(aldmc_step(s, 0.0, 0.1, 100) == ALDMC_OK);         // zero density: time advances
        double t; aldmc_time(s, &t); CHECK_NEAR(t, 0.1, 1e-15);
        std::vector<double> th(aldmc_n_faces(s), 2.0);
        CHECK(aldmc_set_face_coverage(s, th.data()) == ALDMC_ERR_ARG);
        std::fill(th.begin(), th.end(), 0.25);
        CHECK(aldmc_set_face_coverage(s, th.data()) == ALDMC_OK);
        double b; aldmc_bottom_coverage(s, &b); CHECK_NEAR(b, 0.25, 1e-15);
        CHECK(aldmc_reset(s) == ALDMC_OK);
        aldmc_bottom_coverage(s, &b); CHECK(b == 0.0);
        aldmc_time(s, &t); CHECK(t == 0.0);
        CHECK(aldmc_bottom_coverage(s, nullptr) == ALDMC_ERR_ARG);
        aldmc_destroy(s);
        aldmc_destroy(nullptr);
    } else {
        std::printf("unknown test %s\n", name.c_str());
        return 2;
    }
    if (g_failures) { std::printf("%s: %d failure(s)\n", name.c_str(), g_failures); return 1; }
    std::printf("%s: ok\n", name.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: test_aldmc <name>\n"); return 2; }
    return run(argv[1]);
}
