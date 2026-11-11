// ctest suite for aldkin. Each test is selected by name on the command line (see CMakeLists.txt); no test framework.
#include "aldkin.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

int failures = 0;

#define CHECK(cond) do { if (!(cond)) { std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++failures; } } while (0)
#define CHECK_CLOSE(a, b, tol) do { const double _a = (a), _b = (b); if (!(std::fabs(_a - _b) <= (tol))) { \
    std::printf("  FAIL %s:%d: %s = %.15g vs %s = %.15g (|diff| %.3e > %.1e)\n", __FILE__, __LINE__, #a, _a, #b, _b, std::fabs(_a - _b), (double)(tol)); ++failures; } } while (0)
#define CHECK_OK(call) do { const int _rc = (call); if (_rc != ALDKIN_OK) { std::printf("  FAIL %s:%d: %s -> %d\n", __FILE__, __LINE__, #call, _rc); ++failures; } } while (0)

// the chamber runs' TMA parameters (reactor/showerhead_pulse_patterned_AR100_long)
aldkin_params tma_params() {
    aldkin_params p;
    aldkin_params_default(&p);
    p.s = 5e-3; p.r = 0.0; p.kdes = 0.0; p.vbar = 410.0; p.N0 = 4e18; p.c0 = 3.79e20;
    p.d = 1e-7; p.AR = 100.0; p.holeDensity = 1.78e13; p.nCells = 200;
    p.flatOnly = 0; p.scheme = ALDKIN_SCHEME_EULER;
    return p;
}

void test_version_and_defaults() {
    CHECK(std::strcmp(aldkin_version(), "1.1.0") == 0);
    aldkin_params p; aldkin_params_default(&p);
    CHECK(p.s == 1.0 && p.r == 0.0 && p.kdes == 0.0 && p.vbar == 410.0 && p.N0 == 4e18 && p.c0 == 1.0);
    CHECK(p.d == 1e-7 && p.AR == 50.0 && p.holeDensity == 0.0 && p.nCells == 100 && p.flatOnly == 0 && p.scheme == ALDKIN_SCHEME_EULER);
    CHECK(p.nu == 0.0 && p.vbarB == 0.0 && p.sB == 0.0 && p.kB == 0.0 && p.kEtch == 0.0 && p.cB0 == 0.0);
    aldkin_t* h = aldkin_create(&p);
    CHECK(h != nullptr);
    aldkin_params q; CHECK_OK(aldkin_get_params(h, &q));
    CHECK(q.nCells == 100 && q.vbar == 410.0);
    CHECK(aldkin_n_faces(h) == 0 && aldkin_n_cells(h) == 100);
    aldkin_destroy(h);
}

// A flat Langmuir face at constant density: theta(t) = 1 - exp(-s (vbar/4) n t / N0). The exponential scheme is this
// solution step by step, so it must agree to round-off at any dt; checked to 1e-9 over 10 s at 1 ms steps.
void test_flat_langmuir_analytic() {
    aldkin_params p = tma_params();
    p.flatOnly = 1; p.scheme = ALDKIN_SCHEME_EXPONENTIAL; p.c0 = 1.0;
    aldkin_t* h = aldkin_create(&p);
    CHECK(h != nullptr);
    const double n = 2.0e18;                       // 1/m3 (c0 = 1: absolute density)
    CHECK_OK(aldkin_set_faces(h, 3, nullptr));
    const double T[3] = {n, n, n};
    const double dt = 1e-3;
    const double a = p.s * (p.vbar / 4.0) * n / p.N0;   // 0.256 1/s
    double theta[3], K[3], flux[3];
    double maxerr = 0.0;
    for (int k = 1; k <= 10000; ++k) {
        CHECK_OK(aldkin_step(h, T, dt, flux, K));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, theta));
        const double exact = 1.0 - std::exp(-a * k * dt);
        for (int f = 0; f < 3; ++f) maxerr = std::fmax(maxerr, std::fabs(theta[f] - exact));
    }
    std::printf("  exponential scheme vs analytic: max |diff| = %.3e over 10 s (a = %.4g 1/s)\n", maxerr, a);
    CHECK(maxerr <= 1e-9);
    // K_eff = (vbar/4) s (1 - theta) and flux = K n at the end
    CHECK_CLOSE(K[0], (p.vbar / 4.0) * p.s * (1.0 - theta[0]), 1e-15);
    CHECK_CLOSE(flux[0], K[0] * n, 1e-3);
    // MID/BOTTOM/HOLE read as 0 on a flat-only handle
    double z[3] = {1, 1, 1};
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, z));
    CHECK(z[0] == 0.0 && z[2] == 0.0);
    aldkin_destroy(h);
}

// The Euler scheme (the boundary conditions') is first order: its error against the analytic solution halves with dt.
void test_flat_langmuir_euler_convergence() {
    aldkin_params p = tma_params();
    p.flatOnly = 1; p.scheme = ALDKIN_SCHEME_EULER; p.c0 = 1.0;
    const double n = 2.0e18, tend = 4.0;
    const double a = p.s * (p.vbar / 4.0) * n / p.N0;
    double err[3];
    const double dts[3] = {4e-3, 2e-3, 1e-3};
    for (int j = 0; j < 3; ++j) {
        aldkin_t* h = aldkin_create(&p);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const int steps = static_cast<int>(std::lround(tend / dts[j]));
        double th;
        for (int k = 0; k < steps; ++k) CHECK_OK(aldkin_step(h, &n, dts[j], nullptr, nullptr));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &th));
        err[j] = std::fabs(th - (1.0 - std::exp(-a * tend)));
        aldkin_destroy(h);
    }
    std::printf("  Euler error at dt 4e-3, 2e-3, 1e-3 s: %.3e %.3e %.3e (ratios %.3f %.3f)\n",
                err[0], err[1], err[2], err[0] / err[1], err[1] / err[2]);
    CHECK(err[0] / err[1] > 1.9 && err[0] / err[1] < 2.1);
    CHECK(err[1] / err[2] > 1.9 && err[1] / err[2] < 2.1);
    CHECK(err[2] < 1e-4);
}

// Flat Langmuir with desorption (the langmuirSticking kdes term): the Euler step is the literal boundary-condition
// expression, the steady state is theta = a/(a + kdes), refValue = J_des/K in T units, the uptake flux is K c - J_des.
void test_flat_langmuir_desorption() {
    aldkin_params p = tma_params();
    p.flatOnly = 1; p.kdes = 0.5; p.c0 = 3.79e20; p.s = 1e-4;        // water-like
    for (int scheme = 0; scheme < 2; ++scheme) {
        p.scheme = scheme;
        aldkin_t* h = aldkin_create(&p);
        CHECK(h != nullptr);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const double Tf = 0.3, dt = 1e-3;
        const double a = (p.vbar / 4.0) * p.s * p.c0 * Tf / p.N0;          // 0.2914 1/s
        // one literal Euler step from theta = 0.2
        double th0 = 0.2;
        CHECK_OK(aldkin_set_theta(h, ALDKIN_THETA_FLAT, &th0));
        double K, flux, th;
        CHECK_OK(aldkin_step(h, &Tf, dt, &flux, &K));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &th));
        if (scheme == ALDKIN_SCHEME_EULER) {
            const double lit = th0 + dt * ((p.vbar / 4.0) * p.s * (1.0 - th0) * p.c0 * Tf / p.N0 - p.kdes * th0);
            CHECK(th == lit);                                                 // bit-identical to the BC line
        }
        CHECK_CLOSE(K, (p.vbar / 4.0) * (p.s * (1.0 - th) + p.r), 1e-16);
        CHECK_CLOSE(flux, K * p.c0 * Tf - p.kdes * th * p.N0, 1e3);          // 1e3 of ~1e19
        double refV, refG, vf, Dd = 0.0084 / 2.5e-4, Jd;
        CHECK_OK(aldkin_mixed_bc(h, &Dd, &refV, &refG, &vf));
        CHECK_OK(aldkin_desorption_flux(h, &Jd));
        CHECK_CLOSE(refV, (p.kdes * th * p.N0 / p.c0) / K, 1e-12);
        CHECK(refG == 0.0);
        CHECK_CLOSE(vf, K / (K + Dd), 1e-16);
        CHECK_CLOSE(Jd, p.kdes * th * p.N0, 1.0);
        // to steady state
        for (int k = 0; k < 60000; ++k) CHECK_OK(aldkin_step(h, &Tf, dt, nullptr, nullptr));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &th));
        CHECK_CLOSE(th, a / (a + p.kdes), 1e-6);
        std::printf("  scheme %d: steady theta %.8f vs a/(a+kdes) %.8f\n", scheme, th, a / (a + p.kdes));
        aldkin_destroy(h);
    }
    // kdes on a patterned handle is refused
    aldkin_params q = tma_params(); q.kdes = 0.1;
    CHECK(aldkin_create(&q) == nullptr);
    CHECK(std::strstr(aldkin_last_error(nullptr), "kdes") != nullptr);
}

// An AR-100 hole at the constant face density T = 0.1 (3.79e19 1/m3, 0.3 Pa at 573 K), 1 ms steps, exponential
// scheme: the Python twin aldconf.solver.simulate_driven (same discretisation, same update) reaches thetaBottom 0.9 at
// t = 15.616 s with thetaBottom(15.616) = 0.900093212121 (computed with substeps_per_sample = 1 on a 1 ms grid).
// The Euler scheme (the boundary condition's) is checked to lie within 1 % of it.
void test_hole_constant_density_t90() {
    const double PY_T90 = 15.616, PY_THETA_AT_T90 = 0.900093212121;
    const double Tf = 0.1, dt = 1e-3;
    double t90[2] = {0, 0}, th_at[2] = {0, 0};
    for (int scheme = 1; scheme >= 0; --scheme) {
        aldkin_params p = tma_params(); p.scheme = scheme;
        aldkin_t* h = aldkin_create(&p);
        CHECK(h != nullptr);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        double bottom = 0.0, mid = 0.0, flat = 0.0;
        int k = 0;
        while (bottom < 0.9 && k < 40000) {
            CHECK_OK(aldkin_step(h, &Tf, dt, nullptr, nullptr));
            ++k;
            CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, &bottom));
            if (k == 15616) th_at[scheme] = bottom;
        }
        t90[scheme] = k * dt;
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_MID, &mid));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &flat));
        std::printf("  scheme %d: bottom t90 = %.3f s, thetaBottom(15.616 s) = %.12f, mid %.6f, flat %.6f\n",
                    scheme, t90[scheme], th_at[scheme], mid, flat);
        CHECK(mid > 0.9999 && flat > 0.9999);
        // the relative gas profile decays monotonically from the mouth
        std::vector<double> nrel(p.nCells);
        CHECK_OK(aldkin_hole_gas(h, 0, nrel.data()));
        CHECK(nrel[0] < 1.0 && nrel[0] > 0.0);
        bool mono = true;
        for (int i = 1; i < p.nCells; ++i) mono = mono && nrel[i] <= nrel[i - 1] + 1e-15;
        CHECK(mono);
        aldkin_destroy(h);
    }
    CHECK_CLOSE(t90[ALDKIN_SCHEME_EXPONENTIAL], PY_T90, 1e-3 + 1e-9);
    CHECK_CLOSE(th_at[ALDKIN_SCHEME_EXPONENTIAL], PY_THETA_AT_T90, 1e-6);
    CHECK(std::fabs(t90[ALDKIN_SCHEME_EULER] - PY_T90) / PY_T90 < 0.01);
}

// One Euler step of a patterned face equals the boundary condition's lines evaluated by hand (including the order:
// hole gas from the old coverages, K_eff from the new flat coverage and the old hole flux).
void test_hole_euler_literal_step() {
    aldkin_params p = tma_params(); p.nCells = 4;
    aldkin_t* h = aldkin_create(&p);
    CHECK_OK(aldkin_set_faces(h, 1, nullptr));
    const double thF0 = 0.3, thH0[4] = {0.4, 0.3, 0.2, 0.1};
    CHECK_OK(aldkin_set_theta(h, ALDKIN_THETA_FLAT, &thF0));
    CHECK_OK(aldkin_set_theta(h, ALDKIN_THETA_HOLE, thH0));
    const double Tf = 0.02, dt = 1e-3;
    // evaluate first (no advance): gas from the current coverages
    double K0;
    CHECK_OK(aldkin_evaluate(h, &Tf, nullptr, &K0));
    double nrel[4];
    CHECK_OK(aldkin_hole_gas(h, 0, nrel));
    double thH[4], thF;
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_HOLE, thH));
    CHECK(thH[0] == 0.4 && thH[3] == 0.1);
    // hand solve of the 4-cell tridiagonal system
    const int N = 4;
    const double Lh = p.AR * p.d, dx = Lh / N, DK = p.d * p.vbar / 3.0, a = DK / (dx * dx);
    double lower[4], diag[4], upper[4], rhs[4] = {0, 0, 0, 0}, cp[4], dp[4], nr[4];
    for (int i = 0; i < N; ++i) { lower[i] = -a; upper[i] = -a; diag[i] = 2.0 * a + (p.vbar / p.d) * (p.s * (1.0 - thH0[i]) + p.r); }
    diag[N - 1] = a + (diag[N - 1] - 2.0 * a);
    diag[0] += a;
    rhs[0] = 2.0 * a;
    cp[0] = upper[0] / diag[0]; dp[0] = rhs[0] / diag[0];
    for (int i = 1; i < N; ++i) { const double m = diag[i] - lower[i] * cp[i - 1]; cp[i] = (i < N - 1) ? upper[i] / m : 0.0; dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / m; }
    nr[N - 1] = dp[N - 1];
    for (int i = N - 2; i >= 0; --i) nr[i] = dp[i] - cp[i] * nr[i + 1];
    for (int i = 0; i < N; ++i) CHECK(nrel[i] == nr[i]);
    const double jhole = DK * (1.0 - nr[0]) / (0.5 * dx);
    const double phi = p.holeDensity * 3.14159265358979323846 * p.d * p.d / 4.0;
    CHECK(K0 == (1.0 - phi) * (p.vbar / 4.0) * (p.s * (1.0 - thF0) + p.r) + phi * jhole);
    // now one step
    double K1, flux1;
    CHECK_OK(aldkin_step(h, &Tf, dt, &flux1, &K1));
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &thF));
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_HOLE, thH));
    const double thF1 = std::fmin(1.0, thF0 + dt * (p.vbar / 4.0) * p.s * (1.0 - thF0) * p.c0 * Tf / p.N0);
    CHECK(thF == thF1);
    for (int i = 0; i < N; ++i) {
        const double e = std::fmin(1.0, thH0[i] + dt * (p.vbar / 4.0) * p.s * (1.0 - thH0[i]) * p.c0 * Tf * nr[i] / p.N0);
        CHECK(thH[i] == e);
    }
    CHECK(K1 == (1.0 - phi) * (p.vbar / 4.0) * (p.s * (1.0 - thF1) + p.r) + phi * jhole);
    CHECK(flux1 == K1 * p.c0 * Tf);
    // the mixed-BC mapping
    double Dd = 0.0084 * 4000.0, rv, rg, vf;
    CHECK_OK(aldkin_mixed_bc(h, &Dd, &rv, &rg, &vf));
    CHECK(rv == 0.0 && rg == 0.0 && vf == K1 / (K1 + Dd));
    // bottom and mid selectors
    double b, m;
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, &b));
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_MID, &m));
    CHECK(b == thH[3] && m == thH[2]);
    std::printf("  literal Euler step reproduced: K_eff %.10g -> %.10g m/s, thetaFlat %.10g\n", K0, K1, thF);
    aldkin_destroy(h);
}

// Radial tables: linear between points, clamped outside, a repeated radius is a step, NULL keeps the scalar,
// negatives clamp to 0, a density of 0 makes the face flat (no hole solved, hole coverages stay 0).
void test_profiles_interpolate() {
    aldkin_params p = tma_params();
    aldkin_t* h = aldkin_create(&p);
    const double radius[7] = {-1.0, 0.0, 0.05, 0.0875, 0.1, 0.12, 0.2};
    CHECK_OK(aldkin_set_faces(h, 7, radius));
    // ARProfile table ((0 100) (0.075 100) (0.1 60) (0.15 60)) as in paper/assumptions.md Night 3, part 1
    const double rAR[4] = {0.0, 0.075, 0.1, 0.15}, AR[4] = {100, 100, 60, 60};
    CHECK_OK(aldkin_set_profiles(h, 4, rAR, AR, nullptr));
    double ARf[7], dens[7];
    CHECK_OK(aldkin_face_pattern(h, ARf, dens));
    const double expAR[7] = {100, 100, 100, 80, 60, 60, 60};
    for (int i = 0; i < 7; ++i) CHECK_CLOSE(ARf[i], expAR[i], 1e-12);
    for (int i = 0; i < 7; ++i) CHECK(dens[i] == 1.78e13);              // NULL table: scalar kept
    std::printf("  AR at r = 0.0875 m: %.6g (expected 80)\n", ARf[3]);
    // density profile with a step at a repeated radius and a negative value: ((0 1.78e13) (0.1 1.78e13) (0.1 0) (0.15 -5))
    const double rD[4] = {0.0, 0.1, 0.1, 0.15}, D[4] = {1.78e13, 1.78e13, 0.0, -5.0};
    CHECK_OK(aldkin_set_profile(h, ALDKIN_PROFILE_DENSITY, 4, rD, D));    // independent radii, AR table kept
    CHECK_OK(aldkin_face_pattern(h, ARf, dens));
    // r = 0.1 lies in the segment (0.075, 0.1] of the AR table -> 60; density: r <= tab[1].first -> segment (0, 0.1] -> 1.78e13
    CHECK_CLOSE(ARf[4], 60.0, 1e-12);
    CHECK(dens[4] == 1.78e13);
    CHECK(dens[5] == 0.0);                                               // 0.12 in the (0.1, 0.15] segment: 0 -> -5 ramp, clamped at 0
    CHECK(dens[6] == 0.0);                                               // beyond the table: -5 clamped
    CHECK(dens[2] == 1.78e13);
    // faces with density 0 are flat: step a while, hole coverages stay 0, K_eff is the flat value
    const double T[7] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    double K[7], bottom[7], flat[7];
    for (int k = 0; k < 100; ++k) CHECK_OK(aldkin_step(h, T, 1e-3, nullptr, K));
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, bottom));
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, flat));
    CHECK(bottom[5] == 0.0 && bottom[6] == 0.0 && bottom[2] > 0.0);
    CHECK(K[5] == (p.vbar / 4.0) * (p.s * (1.0 - flat[5]) + p.r));
    CHECK(K[2] > K[5]);                                                  // the patterned face is the stronger sink
    // the AR-60 face (index 4) fills its bottom faster than the AR-100 faces
    CHECK(bottom[4] > bottom[2]);
    // clearing the profiles restores the scalars (coverages untouched)
    CHECK_OK(aldkin_set_profiles(h, 0, nullptr, nullptr, nullptr));
    CHECK_OK(aldkin_face_pattern(h, ARf, dens));
    CHECK(ARf[4] == 100.0 && dens[6] == 1.78e13);
    double flat2[7];
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, flat2));
    CHECK(flat2[3] == flat[3]);
    // descending radii are refused
    const double bad[2] = {0.1, 0.0};
    CHECK(aldkin_set_profiles(h, 2, bad, AR, nullptr) == ALDKIN_ERR_ARG);
    CHECK(std::strstr(aldkin_last_error(h), "ascending") != nullptr);
    CHECK(aldkin_set_profile(h, 5, 2, rAR, AR) == ALDKIN_ERR_ARG);
    // removing one table at a time
    CHECK_OK(aldkin_set_profile(h, ALDKIN_PROFILE_AR, 4, rAR, AR));
    CHECK_OK(aldkin_set_profile(h, ALDKIN_PROFILE_DENSITY, 4, rD, D));
    CHECK_OK(aldkin_set_profile(h, ALDKIN_PROFILE_AR, 0, nullptr, nullptr));
    CHECK_OK(aldkin_face_pattern(h, ARf, dens));
    CHECK(ARf[3] == 100.0 && dens[5] == 0.0);
    aldkin_destroy(h);
}

void test_zero_density_and_zero_holes_are_safe() {
    // zero and negative face density: nothing adsorbs, K_eff is finite and equals the bare-surface value
    {
        aldkin_params p = tma_params();
        aldkin_t* h = aldkin_create(&p);
        CHECK_OK(aldkin_set_faces(h, 2, nullptr));
        const double T[2] = {0.0, -0.5};
        double K[2], flux[2], flat[2], bottom[2];
        for (int k = 0; k < 10; ++k) CHECK_OK(aldkin_step(h, T, 1e-3, flux, K));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, flat));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, bottom));
        CHECK(flat[0] == 0.0 && flat[1] == 0.0 && bottom[0] == 0.0 && bottom[1] == 0.0);
        CHECK(std::isfinite(K[0]) && K[0] > 0.0 && K[0] == K[1]);
        CHECK(flux[0] == 0.0 && flux[1] == 0.0);
        double Dd[2] = {33.6, 33.6}, vf[2];
        CHECK_OK(aldkin_mixed_bc(h, Dd, nullptr, nullptr, vf));
        CHECK(vf[0] > 0.0 && vf[0] < 1.0);
        std::printf("  bare hole-patterned face at zero density: K_eff %.6g m/s, valueFraction %.6g\n", K[0], vf[0]);
        aldkin_destroy(h);
    }
    // holeDensity 0 and AR 0: flat Langmuir faces, phi = 0, no hole solved
    for (int variant = 0; variant < 2; ++variant) {
        aldkin_params p = tma_params();
        if (variant == 0) p.holeDensity = 0.0; else p.AR = 0.0;
        aldkin_t* h = aldkin_create(&p);
        CHECK(h != nullptr);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const double Tf = 0.2;
        double K, flat, bottom, mid;
        for (int k = 0; k < 50; ++k) CHECK_OK(aldkin_step(h, &Tf, 1e-3, nullptr, &K));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &flat));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, &bottom));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_MID, &mid));
        CHECK(bottom == 0.0 && mid == 0.0 && flat > 0.0);
        CHECK(K == (p.vbar / 4.0) * (p.s * (1.0 - flat) + p.r));
        std::vector<double> nrel(p.nCells, 1.0);
        CHECK_OK(aldkin_hole_gas(h, 0, nrel.data()));
        CHECK(nrel[0] == 0.0);
        aldkin_destroy(h);
    }
    // a single-cell hole and zero faces work
    {
        aldkin_params p = tma_params(); p.nCells = 1;
        aldkin_t* h = aldkin_create(&p);
        CHECK(h != nullptr);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const double Tf = 0.1; double K;
        CHECK_OK(aldkin_step(h, &Tf, 1e-3, nullptr, &K));
        CHECK(std::isfinite(K) && K > 0.0);
        CHECK_OK(aldkin_set_faces(h, 0, nullptr));
        CHECK_OK(aldkin_step(h, nullptr, 1e-3, nullptr, nullptr));
        aldkin_destroy(h);
    }
    // zero dt leaves the coverages alone
    {
        aldkin_params p = tma_params();
        aldkin_t* h = aldkin_create(&p);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const double Tf = 0.5; double b0, b1;
        for (int k = 0; k < 20; ++k) CHECK_OK(aldkin_step(h, &Tf, 1e-2, nullptr, nullptr));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, &b0));
        CHECK_OK(aldkin_step(h, &Tf, 0.0, nullptr, nullptr));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, &b1));
        CHECK(b0 == b1 && b0 > 0.0);
        aldkin_destroy(h);
    }
}

// Clone copies the state; set_theta/theta round-trip; a cloned handle continues identically; evaluate does not advance.
void test_clone_and_restart() {
    aldkin_params p = tma_params(); p.nCells = 50;
    aldkin_t* h = aldkin_create(&p);
    const double radius[2] = {0.0, 0.1};
    CHECK_OK(aldkin_set_faces(h, 2, radius));
    const double T[2] = {0.1, 0.05};
    for (int k = 0; k < 300; ++k) CHECK_OK(aldkin_step(h, T, 1e-3, nullptr, nullptr));
    std::vector<double> hole(2 * 50), flat(2);
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_HOLE, hole.data()));
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, flat.data()));
    aldkin_t* c = aldkin_clone(h);
    CHECK(c != nullptr);
    // a fresh handle restarted from the written state
    aldkin_t* r = aldkin_create(&p);
    CHECK_OK(aldkin_set_faces(r, 2, radius));
    CHECK_OK(aldkin_set_theta(r, ALDKIN_THETA_HOLE, hole.data()));
    CHECK_OK(aldkin_set_theta(r, ALDKIN_THETA_FLAT, flat.data()));
    double Kh[2], Kc[2], Kr[2];
    // evaluate is side-effect free
    CHECK_OK(aldkin_evaluate(h, T, nullptr, Kh));
    std::vector<double> hole2(2 * 50);
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_HOLE, hole2.data()));
    CHECK(hole2 == hole);
    for (int k = 0; k < 100; ++k) {
        CHECK_OK(aldkin_step(h, T, 1e-3, nullptr, Kh));
        CHECK_OK(aldkin_step(c, T, 1e-3, nullptr, Kc));
        CHECK_OK(aldkin_step(r, T, 1e-3, nullptr, Kr));
    }
    CHECK(Kh[0] == Kc[0] && Kh[1] == Kc[1] && Kh[0] == Kr[0] && Kh[1] == Kr[1]);
    std::vector<double> hh(100), hc(100), hr(100);
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_HOLE, hh.data()));
    CHECK_OK(aldkin_theta(c, ALDKIN_THETA_HOLE, hc.data()));
    CHECK_OK(aldkin_theta(r, ALDKIN_THETA_HOLE, hr.data()));
    CHECK(hh == hc && hh == hr);
    std::printf("  clone and restart continue bit-identically (K_eff %.10g m/s)\n", Kh[0]);
    aldkin_destroy(h); aldkin_destroy(c); aldkin_destroy(r);
}

void test_error_codes() {
    CHECK(aldkin_create(nullptr) == nullptr);
    CHECK(std::strlen(aldkin_last_error(nullptr)) > 0);
    aldkin_params p = tma_params(); p.vbar = -1;
    CHECK(aldkin_create(&p) == nullptr);
    CHECK(std::strstr(aldkin_last_error(nullptr), "vbar") != nullptr);
    p = tma_params(); p.nCells = 0;
    CHECK(aldkin_create(&p) == nullptr);
    p = tma_params(); p.scheme = 7;
    CHECK(aldkin_create(&p) == nullptr);
    p = tma_params(); p.s = 1.5;
    CHECK(aldkin_create(&p) == nullptr);
    p = tma_params();
    aldkin_t* h = aldkin_create(&p);
    const double Tf = 0.1;
    CHECK(aldkin_step(h, &Tf, 1e-3, nullptr, nullptr) == ALDKIN_ERR_STATE);     // faces not set
    CHECK(std::strstr(aldkin_last_error(h), "set_faces") != nullptr);
    CHECK_OK(aldkin_set_faces(h, 1, nullptr));
    CHECK(aldkin_step(h, nullptr, 1e-3, nullptr, nullptr) == ALDKIN_ERR_ARG);
    CHECK(aldkin_step(h, &Tf, -1.0, nullptr, nullptr) == ALDKIN_ERR_ARG);
    const double nan = std::nan("");
    CHECK(aldkin_step(h, &nan, 1e-3, nullptr, nullptr) == ALDKIN_ERR_ARG);
    double out;
    CHECK(aldkin_theta(h, 9, &out) == ALDKIN_ERR_ARG);
    CHECK(aldkin_theta(h, ALDKIN_THETA_FLAT, nullptr) == ALDKIN_ERR_ARG);
    CHECK(aldkin_hole_gas(h, 5, &out) == ALDKIN_ERR_ARG);
    CHECK(aldkin_set_theta(h, ALDKIN_THETA_MID, &out) == ALDKIN_ERR_ARG);
    CHECK(aldkin_set_faces(h, -1, nullptr) == ALDKIN_ERR_ARG);
    CHECK(aldkin_set_profiles(h, 2, nullptr, nullptr, nullptr) == ALDKIN_ERR_ARG);
    CHECK(aldkin_step(nullptr, &Tf, 1e-3, nullptr, nullptr) == ALDKIN_ERR_ARG);
    CHECK(aldkin_n_faces(nullptr) == -1);
    aldkin_destroy(nullptr);              // harmless
    // the handle still works after errors
    CHECK_OK(aldkin_step(h, &Tf, 1e-3, nullptr, nullptr));
    aldkin_destroy(h);
}


// ---- two-species option ---------------------------------------------------------------------------------------------

// TiCl4/HCl-like parameters on the AR-100 hole: s 0.1, N0 1.24e18, vbar 253 m/s (189.68 amu, 573 K), HCl vbarB 577 m/s
aldkin_params ticl4_params() {
    aldkin_params p = tma_params();
    p.s = 0.10; p.N0 = 1.24e18; p.vbar = 253.0; p.c0 = 3.79e20;
    p.nu = 0.5; p.vbarB = 577.0; p.sB = 0.0; p.kB = 0.0; p.kEtch = 0.0; p.cB0 = 0.0;
    return p;
}

// With sB = kEtch = 0 the two-species fields (any nu, vbarB, kB) leave the single-species path untouched: every
// coverage, K_eff and flux is bit-identical to the handle without them, for both models and both schemes.
void test_two_species_defaults_are_bit_identical() {
    for (int flat = 0; flat < 2; ++flat) {
        for (int scheme = 0; scheme < 2; ++scheme) {
            aldkin_params a = tma_params(); a.flatOnly = flat; a.scheme = scheme; if (flat) a.kdes = 0.3;
            aldkin_params b = a; b.nu = 4.0; b.vbarB = 577.0; b.kB = 10.0; b.cB0 = 1e18;   // inactive: sB = kEtch = 0
            aldkin_t* ha = aldkin_create(&a); aldkin_t* hb = aldkin_create(&b);
            CHECK(ha != nullptr && hb != nullptr);
            CHECK_OK(aldkin_set_faces(ha, 2, nullptr)); CHECK_OK(aldkin_set_faces(hb, 2, nullptr));
            const double T[2] = {0.1, 0.03};
            double Ka[2], Kb[2], Fa[2], Fb[2];
            for (int k = 0; k < 400; ++k) {
                CHECK_OK(aldkin_step(ha, T, 1e-3, Fa, Ka));
                CHECK_OK(aldkin_step(hb, T, 1e-3, Fb, Kb));
            }
            CHECK(Ka[0] == Kb[0] && Ka[1] == Kb[1] && Fa[0] == Fb[0] && Fa[1] == Fb[1]);
            std::vector<double> ta(2 * 200), tb(2 * 200), fa(2), fb(2), ba(2);
            CHECK_OK(aldkin_theta(ha, ALDKIN_THETA_HOLE, ta.data())); CHECK_OK(aldkin_theta(hb, ALDKIN_THETA_HOLE, tb.data()));
            CHECK_OK(aldkin_theta(ha, ALDKIN_THETA_FLAT, fa.data())); CHECK_OK(aldkin_theta(hb, ALDKIN_THETA_FLAT, fb.data()));
            CHECK(ta == tb && fa == fb);
            CHECK_OK(aldkin_theta(hb, ALDKIN_THETA_BLOCKED_FLAT, ba.data()));
            CHECK(ba[0] == 0.0 && ba[1] == 0.0);
            // byproduct reporting on the inactive path: everything produced escapes
            double J[2], b5[5];
            CHECK_OK(aldkin_byproduct_flux(hb, J)); CHECK(J[0] > 0.0);
            CHECK_OK(aldkin_byproduct_flux(ha, J)); CHECK(J[0] == 0.0);
            if (!flat) { CHECK_OK(aldkin_byproduct_balance(hb, 0, b5)); CHECK(b5[0] == b5[2] && b5[1] == 0.0 && b5[3] == 0.0); }
            aldkin_destroy(ha); aldkin_destroy(hb);
        }
    }
    std::printf("  inactive two-species fields: bit-identical to the single-species handle (4 model/scheme combinations)\n");
}

// kB -> infinity: blocked sites are released at once and the precursor coverages are the single-species ones to 1e-9.
void test_two_species_fast_desorption_limit() {
    aldkin_params a = ticl4_params(); a.scheme = ALDKIN_SCHEME_EXPONENTIAL;
    aldkin_params b = a; b.nu = 4.0; b.sB = 0.5; b.kB = 1e14;
    aldkin_t* ha = aldkin_create(&a); aldkin_t* hb = aldkin_create(&b);
    CHECK(ha != nullptr && hb != nullptr);
    CHECK_OK(aldkin_set_faces(ha, 1, nullptr)); CHECK_OK(aldkin_set_faces(hb, 1, nullptr));
    const double Tf = 0.1;
    double Ka, Kb;
    for (int k = 0; k < 2000; ++k) { CHECK_OK(aldkin_step(ha, &Tf, 1e-3, nullptr, &Ka)); CHECK_OK(aldkin_step(hb, &Tf, 1e-3, nullptr, &Kb)); }
    std::vector<double> ta(200), tb(200), bb(200);
    CHECK_OK(aldkin_theta(ha, ALDKIN_THETA_HOLE, ta.data())); CHECK_OK(aldkin_theta(hb, ALDKIN_THETA_HOLE, tb.data()));
    CHECK_OK(aldkin_theta(hb, ALDKIN_THETA_BLOCKED_HOLE, bb.data()));
    double maxdiff = 0.0, maxb = 0.0;
    for (int i = 0; i < 200; ++i) { maxdiff = std::fmax(maxdiff, std::fabs(ta[i] - tb[i])); maxb = std::fmax(maxb, bb[i]); }
    std::printf("  kB = 1e14: max |theta - single-species| = %.3e, max blocked coverage %.3e\n", maxdiff, maxb);
    CHECK(maxdiff <= 1e-9);
    CHECK(maxb <= 1e-9);
    CHECK_CLOSE(Ka, Kb, 1e-9 * std::fabs(Ka));
    aldkin_destroy(ha); aldkin_destroy(hb);
}

// Byproduct bookkeeping of the hole: produced + desorbed = escaped + re-adsorbed + etched after every evaluate, to
// 1e-9 relative, with and without etching, with and without a byproduct density at the mouth; the per-face flux
// equals the open-area blend.
void test_two_species_mass_balance() {
    const double sB[3] = {0.05, 0.05, 0.1}, kE[3] = {0.0, 1e-3, 0.0}, cB0[3] = {0.0, 0.0, 2e18};
    for (int v = 0; v < 3; ++v) {
        aldkin_params p = ticl4_params(); p.nu = 4.0; p.sB = sB[v]; p.kB = 10.0; p.kEtch = kE[v]; p.cB0 = cB0[v];
        p.scheme = ALDKIN_SCHEME_EXPONENTIAL;
        aldkin_t* h = aldkin_create(&p);
        CHECK(h != nullptr);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const double Tf = 0.1;
        double worst = 0.0, b5[5] = {0, 0, 0, 0, 0}, J;
        for (int k = 0; k < 500; ++k) {
            CHECK_OK(aldkin_step(h, &Tf, 1e-3, nullptr, nullptr));
            CHECK_OK(aldkin_byproduct_balance(h, 0, b5));
            const double in = b5[0] + b5[1], out = b5[2] + b5[3] + b5[4];
            worst = std::fmax(worst, std::fabs(in - out) / in);
        }
        CHECK(worst <= 1e-9);
        CHECK(b5[0] > 0.0 && b5[3] > 0.0);
        if (kE[v] == 0.0) CHECK(b5[4] == 0.0); else CHECK(b5[4] > 0.0);
        CHECK_OK(aldkin_byproduct_flux(h, &J));
        // flux = (1 - phi) flat part + phi escaped; the flat part at cB0 = 0 is nu x flat uptake + kB N0 thetaB_flat (no B at the face)
        double thF, thFB;
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_FLAT, &thF)); CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BLOCKED_FLAT, &thFB));
        const double phi = p.holeDensity * 3.14159265358979323846 * p.d * p.d / 4.0;
        const double freeF = 1.0 - thF - thFB;
        const double flat = p.nu * p.s * (p.vbar / 4.0) * freeF * p.c0 * Tf + p.kB * p.N0 * thFB
                            - (p.vbarB / 4.0) * (p.sB * freeF + p.kEtch * thF) * p.cB0;
        CHECK_CLOSE(J, (1.0 - phi) * flat + phi * b5[2], 1e-9 * std::fabs(J));
        if (cB0[v] > 0.0) CHECK(thFB > 0.0); else CHECK(thFB == 0.0);
        std::printf("  variant %d: worst closure %.2e; produced %.3e, desorbed %.3e, escaped %.3e, adsorbed %.3e, etched %.3e per m2 s\n",
                    v, worst, b5[0], b5[1], b5[2], b5[3], b5[4]);
        aldkin_destroy(h);
    }
    // the per-face mouth density override
    aldkin_params p = ticl4_params(); p.sB = 0.1; p.kB = 1.0; p.scheme = ALDKIN_SCHEME_EXPONENTIAL;
    aldkin_t* h = aldkin_create(&p);
    CHECK_OK(aldkin_set_faces(h, 2, nullptr));
    const double nBf[2] = {0.0, 3e18};
    CHECK_OK(aldkin_set_byproduct_mouth(h, nBf));
    const double T[2] = {0.1, 0.1};
    for (int k = 0; k < 50; ++k) CHECK_OK(aldkin_step(h, T, 1e-3, nullptr, nullptr));
    double bf[2];
    CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BLOCKED_FLAT, bf));
    CHECK(bf[0] == 0.0 && bf[1] > 0.0);
    CHECK_OK(aldkin_set_byproduct_mouth(h, nullptr));          // back to cB0 = 0
    aldkin_destroy(h);
}

// Re-adsorption slows the hole bottom; the Euler and exponential schemes agree; etching lowers the coverage further.
void test_two_species_blocking_slows_the_bottom() {
    double bottom[3], bottomE = 0.0;
    for (int v = 0; v < 3; ++v) {
        aldkin_params p = ticl4_params(); p.scheme = ALDKIN_SCHEME_EXPONENTIAL;
        if (v >= 1) { p.sB = 0.1; p.kB = 1.0; }
        if (v == 2) p.kEtch = 1e-3;
        aldkin_t* h = aldkin_create(&p);
        CHECK_OK(aldkin_set_faces(h, 1, nullptr));
        const double Tf = 0.1;
        for (int k = 0; k < 12000; ++k) CHECK_OK(aldkin_step(h, &Tf, 1e-3, nullptr, nullptr));
        CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BOTTOM, &bottom[v]));
        if (v == 1) {
            double bb;
            CHECK_OK(aldkin_theta(h, ALDKIN_THETA_BLOCKED_BOTTOM, &bb));
            CHECK(bb > 0.0);
            // Euler scheme at the same step
            aldkin_params q = p; q.scheme = ALDKIN_SCHEME_EULER;
            aldkin_t* he = aldkin_create(&q);
            CHECK_OK(aldkin_set_faces(he, 1, nullptr));
            for (int k = 0; k < 12000; ++k) CHECK_OK(aldkin_step(he, &Tf, 1e-3, nullptr, nullptr));
            CHECK_OK(aldkin_theta(he, ALDKIN_THETA_BOTTOM, &bottomE));
            aldkin_destroy(he);
        }
        aldkin_destroy(h);
    }
    std::printf("  bottom coverage at 12 s: no HCl action %.4f, blocking %.4f (Euler %.4f), blocking + etch %.4f\n",
                bottom[0], bottom[1], bottomE, bottom[2]);
    CHECK(bottom[0] > 0.5);
    CHECK(bottom[1] < bottom[0] - 1e-3);
    CHECK(bottom[2] < bottom[1]);
    CHECK(std::fabs(bottomE - bottom[1]) < 0.01);
}

void test_two_species_errors() {
    aldkin_params p = ticl4_params(); p.sB = 0.1;                 // vbarB set, fine
    aldkin_t* h = aldkin_create(&p); CHECK(h != nullptr); aldkin_destroy(h);
    p = ticl4_params(); p.sB = 0.1; p.vbarB = 0.0;
    CHECK(aldkin_create(&p) == nullptr);
    CHECK(std::strstr(aldkin_last_error(nullptr), "vbarB") != nullptr);
    p = ticl4_params(); p.sB = 1.5; CHECK(aldkin_create(&p) == nullptr);
    p = ticl4_params(); p.kEtch = -1.0; CHECK(aldkin_create(&p) == nullptr);
    p = ticl4_params(); p.nu = -1.0; CHECK(aldkin_create(&p) == nullptr);
    p = ticl4_params(); p.cB0 = -1.0; CHECK(aldkin_create(&p) == nullptr);
    p = ticl4_params(); p.sB = 0.1;
    h = aldkin_create(&p);
    const double nBf = 1.0;
    CHECK(aldkin_set_byproduct_mouth(h, &nBf) == ALDKIN_ERR_STATE);     // faces not set
    CHECK_OK(aldkin_set_faces(h, 1, nullptr));
    CHECK_OK(aldkin_set_byproduct_mouth(h, &nBf));
    double out[5];
    CHECK(aldkin_byproduct_balance(h, 3, out) == ALDKIN_ERR_ARG);
    CHECK(aldkin_hole_byproduct_gas(h, -1, out) == ALDKIN_ERR_ARG);
    CHECK(aldkin_byproduct_flux(h, nullptr) == ALDKIN_ERR_ARG);
    aldkin_destroy(h);
}

}  // namespace

int main(int argc, char** argv) {
    const std::string name = argc > 1 ? argv[1] : "";
    struct Entry { const char* n; void (*f)(); } tests[] = {
        {"version_and_defaults", test_version_and_defaults},
        {"flat_langmuir_analytic", test_flat_langmuir_analytic},
        {"flat_langmuir_euler_convergence", test_flat_langmuir_euler_convergence},
        {"flat_langmuir_desorption", test_flat_langmuir_desorption},
        {"hole_constant_density_t90", test_hole_constant_density_t90},
        {"hole_euler_literal_step", test_hole_euler_literal_step},
        {"profiles_interpolate", test_profiles_interpolate},
        {"zero_density_and_zero_holes_are_safe", test_zero_density_and_zero_holes_are_safe},
        {"clone_and_restart", test_clone_and_restart},
        {"error_codes", test_error_codes},
        {"two_species_defaults_are_bit_identical", test_two_species_defaults_are_bit_identical},
        {"two_species_fast_desorption_limit", test_two_species_fast_desorption_limit},
        {"two_species_mass_balance", test_two_species_mass_balance},
        {"two_species_blocking_slows_the_bottom", test_two_species_blocking_slows_the_bottom},
        {"two_species_errors", test_two_species_errors},
    };
    int ran = 0;
    for (const Entry& e : tests) {
        if (name.empty() || name == e.n) {
            std::printf("[%s]\n", e.n);
            e.f();
            ++ran;
        }
    }
    if (ran == 0) { std::printf("unknown test '%s'\n", name.c_str()); return 2; }
    std::printf("%s: %d failure(s)\n", name.empty() ? "all" : name.c_str(), failures);
    return failures == 0 ? 0 : 1;
}
