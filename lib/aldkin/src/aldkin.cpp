// aldkin: implementation. The arithmetic of the boundary conditions reactor/lib/aldWallKinetics (now in legacy/)
// is transcribed expression by expression, in the same evaluation order, so that a host calling aldkin_step once per
// time step reproduces them bit for bit. Comments marked [BC] quote the original line.
#include "aldkin.h"

#include <cmath>
#include <cstring>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;   // Foam::constant::mathematical::pi
constexpr double kSmall = 1.0e-15;                // Foam::SMALL for double

inline double fmin2(double a, double b) { return (a < b) ? a : b; }   // Foam::min(scalar, scalar)
inline double fmax2(double a, double b) { return (a > b) ? a : b; }   // Foam::max(scalar, scalar)

thread_local std::string g_create_error;

const char* version_string = "1.1.0";

}  // namespace

struct aldkin_t {
    aldkin_params p{};
    int nFaces = 0;
    bool facesSet = false;
    std::vector<double> radius;                 // face-centre radius (m)
    std::vector<std::pair<double, double>> ARProfile, holeDensityProfile;
    std::vector<double> ARField, holeDensityField;
    std::vector<double> thetaFlat;              // nFaces
    std::vector<double> thetaHole;              // nFaces * nCells, cell index fastest
    std::vector<double> nrel;                   // nFaces * nCells, last evaluated relative hole gas
    std::vector<double> Keff;                   // nFaces
    std::vector<double> refValue;               // nFaces (flat-only desorption)
    std::vector<double> Jdes;                   // nFaces
    // two-species option (version 1.1)
    bool twoSpecies = false;                    // p.sB > 0 || p.kEtch > 0
    std::vector<double> thetaFlatB;             // nFaces: flat sites blocked by the byproduct
    std::vector<double> thetaHoleB;             // nFaces * nCells
    std::vector<double> nB;                     // nFaces * nCells, last evaluated absolute byproduct density [1/m3]
    std::vector<double> nBface;                 // nFaces: byproduct density at the face / hole mouth [1/m3]
    std::vector<double> Jbyp;                   // nFaces: net byproduct flux per unit wafer area
    std::vector<double> bal;                    // nFaces * 5: produced, desorbed, escaped, adsorbed, etched (hole, per mouth area)
    // scratch for the tridiagonal solve
    std::vector<double> lower, diag, upper, rhs, cp, dp;
    std::string err;

    int fail(int code, const std::string& msg) { err = msg; return code; }

    bool hasHoles(int facei) const {            // [BC] holeDensityField_[facei] > 0 && ARField_[facei] > 0
        return !p.flatOnly && holeDensityField[facei] > 0 && ARField[facei] > 0;
    }

    static double interpolateProfile(const std::vector<std::pair<double, double>>& tab, double r) {
        const int n = static_cast<int>(tab.size());
        if (r <= tab[0].first) return tab[0].second;
        if (r >= tab[n - 1].first) return tab[n - 1].second;
        for (int i = 1; i < n; ++i) {
            if (r <= tab[i].first) {
                const double r0 = tab[i - 1].first, r1 = tab[i].first;
                const double v0 = tab[i - 1].second, v1 = tab[i].second;
                if (r1 <= r0) return v1;                   // repeated radius: step
                return v0 + (v1 - v0) * (r - r0) / (r1 - r0);
            }
        }
        return tab[n - 1].second;
    }

    // [BC] setPatternFields(): scalars everywhere, then the profiles at the face radii, clamped at 0
    void setPatternFields() {
        ARField.assign(nFaces, p.AR);
        holeDensityField.assign(nFaces, p.holeDensity);
        if (ARProfile.empty() && holeDensityProfile.empty()) return;
        for (int facei = 0; facei < nFaces; ++facei) {
            const double r = radius[facei];
            if (!ARProfile.empty()) ARField[facei] = interpolateProfile(ARProfile, r);
            if (!holeDensityProfile.empty()) holeDensityField[facei] = interpolateProfile(holeDensityProfile, r);
            if (ARField[facei] < 0) ARField[facei] = 0;
            if (holeDensityField[facei] < 0) holeDensityField[facei] = 0;
        }
    }

    // [BC] solveHole(): quasi-steady Knudsen diffusion with a Langmuir sink, mouth density 1, closed bottom.
    // Returns the mouth flux per unit mouth area per unit mouth density [m/s]; fills nrel for this face.
    // thB (two-species option only, else nullptr): the blocked coverage, which the free fraction in the sink excludes.
    double solveHole(int facei, double* nr, const double* thB = nullptr) {
        const int N = p.nCells;
        const double Lh = ARField[facei] * p.d;
        const double dx = Lh / N;
        const double DK = p.d * p.vbar / 3.0;              // Knudsen diffusivity, cylinder
        const double a = DK / (dx * dx);
        lower.assign(N, -a); diag.assign(N, 0.0); upper.assign(N, -a); rhs.assign(N, 0.0);
        const double* th = &thetaHole[static_cast<size_t>(facei) * N];
        if (thB) {
            for (int i = 0; i < N; ++i) diag[i] = 2.0 * a + (p.vbar / p.d) * (p.s * (1.0 - th[i] - thB[i]) + p.r);
        } else {
            for (int i = 0; i < N; ++i) {
                const double k = (p.vbar / p.d) * (p.s * (1.0 - th[i]) + p.r);   // (4/d)(vbar/4)[...]
                diag[i] = 2.0 * a + k;
            }
        }
        diag[N - 1] = a + (diag[N - 1] - 2.0 * a);         // Neumann at the bottom: only the west face conducts
        diag[0] += a;                                      // Dirichlet ghost at the mouth: n_ghost = 2 - n_0
        rhs[0] = 2.0 * a;
        // Thomas algorithm
        cp.assign(N, 0.0); dp.assign(N, 0.0);
        cp[0] = upper[0] / diag[0]; dp[0] = rhs[0] / diag[0];
        for (int i = 1; i < N; ++i) {
            const double m = diag[i] - lower[i] * cp[i - 1];
            cp[i] = (i < N - 1) ? upper[i] / m : 0.0;
            dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / m;
        }
        nr[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i) nr[i] = dp[i] - cp[i] * nr[i + 1];
        // flux into the mouth: D (n_mouth - n_0)/(dx/2) with n_mouth = 1
        return DK * (1.0 - nr[0]) / (0.5 * dx);
    }

    // ---- two-species option -------------------------------------------------------------------------------------
    // (exp(z) - 1)/z, 1 at z = 0
    static double e1(double z) { return (std::fabs(z) < 1e-12) ? 1.0 : std::expm1(z) / z; }
    // d/dz (exp(z) - 1)/z = (exp(z)(z - 1) + 1)/z^2, 1/2 at z = 0
    static double e1p(double z) {
        if (std::fabs(z) < 1e-3) return 0.5 + z / 3.0 + z * z / 8.0 + z * z * z / 30.0;
        return (std::exp(z) * (z - 1.0) + 1.0) / (z * z);
    }
    // Exact advance over dt of  d/dt [thP, thB] = M [thP, thB] + [aP, aB],  M = [[-(aP+e), -aP], [-aB, -(aB+kB)]]
    // (frozen gas), the update of aldconf.multispecies.advance_coverages: theta(dt) = theta0 + phi (M theta0 + b),
    // phi = psi(l1) I + (M - l1 I) psi[l1, l2], psi(l) = (exp(l dt) - 1)/l, l1 the eigenvalue of smaller magnitude
    // taken as det/l2 to avoid cancellation; psi[l1, l2] the divided difference (its derivative when l1 ~ l2).
    static void advance2(double& thP, double& thB, double aP, double aB, double e, double kB, double dt) {
        const double m11 = -(aP + e), m12 = -aP, m21 = -aB, m22 = -(aB + kB);
        const double tr = m11 + m22;
        const double det = aP * kB + e * aB + e * kB;
        const double disc = std::sqrt(fmax2(tr * tr - 4.0 * det, 0.0));
        const double l2 = 0.5 * (tr - disc);
        const double l1 = (l2 != 0.0) ? det / l2 : 0.0;
        const double psi1 = dt * e1(l1 * dt), psi2 = dt * e1(l2 * dt);
        const double dl = (l1 - l2) * dt;
        const double dd = (std::fabs(dl) < 1e-6) ? dt * dt * e1p(0.5 * (l1 + l2) * dt) : (psi1 - psi2) / (l1 - l2);
        const double r1 = m11 * thP + m12 * thB + aP;
        const double r2 = m21 * thP + m22 * thB + aB;
        const double phi11 = psi1 + (m11 - l1) * dd, phi12 = m12 * dd, phi21 = m21 * dd, phi22 = psi1 + (m22 - l1) * dd;
        double P = thP + phi11 * r1 + phi12 * r2;
        double B = thB + phi21 * r1 + phi22 * r2;
        P = fmin2(1.0, fmax2(0.0, P));
        B = fmin2(1.0 - P, fmax2(0.0, B));
        thP = P; thB = B;
    }

    // Quasi-steady absolute byproduct density in the hole of one face given the relative precursor profile nr (mouth
    // density c0 Tf) and the current coverages: -DKB n'' + (vbarB/d)[sB f + kEtch th] n = (4/d)[nu uptake + kB N0 thB],
    // Dirichlet nB0 at the mouth, closed bottom. Fills nb; fills the five bookkeeping rates per unit mouth area.
    void solveHoleB(int facei, const double* nr, double Tf, double nB0, double* nb, double* b5) {
        const int N = p.nCells;
        const double Lh = ARField[facei] * p.d;
        const double dx = Lh / N;
        const double DKB = p.d * p.vbarB / 3.0;
        const double a = DKB / (dx * dx);
        const double* th = &thetaHole[static_cast<size_t>(facei) * N];
        const double* thB = &thetaHoleB[static_cast<size_t>(facei) * N];
        lower.assign(N, -a); diag.assign(N, 0.0); upper.assign(N, -a); rhs.assign(N, 0.0);
        double produced = 0.0, desorbed = 0.0;
        for (int i = 0; i < N; ++i) {
            const double f = 1.0 - th[i] - thB[i];
            const double uptake = p.s * (p.vbar / 4.0) * f * p.c0 * Tf * nr[i];       // precursor adsorption per unit wall area
            const double src = (4.0 / p.d) * (p.nu * uptake + p.kB * p.N0 * thB[i]);
            diag[i] = 2.0 * a + (p.vbarB / p.d) * (p.sB * f + p.kEtch * th[i]);
            rhs[i] = src;
            produced += (4.0 / p.d) * p.nu * uptake * dx;
            desorbed += (4.0 / p.d) * p.kB * p.N0 * thB[i] * dx;
        }
        diag[N - 1] -= a;
        diag[0] += a;
        rhs[0] += 2.0 * a * nB0;
        cp.assign(N, 0.0); dp.assign(N, 0.0);
        cp[0] = upper[0] / diag[0]; dp[0] = rhs[0] / diag[0];
        for (int i = 1; i < N; ++i) {
            const double m = diag[i] - lower[i] * cp[i - 1];
            cp[i] = (i < N - 1) ? upper[i] / m : 0.0;
            dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / m;
        }
        nb[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i) nb[i] = dp[i] - cp[i] * nb[i + 1];
        double adsorbed = 0.0, etched = 0.0;
        for (int i = 0; i < N; ++i) {
            const double f = 1.0 - th[i] - thB[i];
            adsorbed += (4.0 / p.d) * (p.vbarB / 4.0) * p.sB * f * nb[i] * dx;
            etched += (4.0 / p.d) * (p.vbarB / 4.0) * p.kEtch * th[i] * nb[i] * dx;
        }
        b5[0] = produced; b5[1] = desorbed; b5[2] = DKB * 2.0 * (nb[0] - nB0) / dx; b5[3] = adsorbed; b5[4] = etched;
    }

    // The two-species pass (p.sB > 0 or p.kEtch > 0). Same order of operations as the single-species pass.
    void updateTwo(const double* Tface, double dt, bool advance, double* flux_out, double* keff_out) {
        const int N = p.nCells;
        const double vbar = p.vbar, vbarB = p.vbarB, s = p.s, r = p.r, c0 = p.c0, N0 = p.N0, kdes = p.kdes;
        const double sB = p.sB, kB = p.kB, kE = p.kEtch, nu = p.nu;
        const bool expo = (p.scheme == ALDKIN_SCHEME_EXPONENTIAL);
        for (int facei = 0; facei < nFaces; ++facei) {
            const double Tf = fmax2(Tface[facei], 0.0);
            const double nBf = nBface[facei];
            const bool holes = hasHoles(facei);
            double* nr = &nrel[static_cast<size_t>(facei) * N];
            double* nb = &nB[static_cast<size_t>(facei) * N];
            double* b5 = &bal[static_cast<size_t>(facei) * 5];
            const double phi = holes ? holeDensityField[facei] * kPi * p.d * p.d / 4.0 : 0.0;
            double jhole = 0.0;
            if (holes) {
                jhole = solveHole(facei, nr, &thetaHoleB[static_cast<size_t>(facei) * N]);
                solveHoleB(facei, nr, Tf, nBf, nb, b5);
            } else {
                for (int i = 0; i < N; ++i) { nr[i] = 0.0; nb[i] = 0.0; }
                for (int k = 0; k < 5; ++k) b5[k] = 0.0;
            }
            // flat surface rates (byproduct density at the face = nBf)
            const double aPf = (vbar / 4.0) * s * c0 * Tf / N0;
            const double aBf = (vbarB / 4.0) * sB * nBf / N0;
            const double ef = (vbarB / 4.0) * kE * nBf / N0;
            double& thF = thetaFlat[facei];
            double& thFB = thetaFlatB[facei];
            if (advance) {
                if (!expo) {
                    const double f = 1.0 - thF - thFB;
                    const double P = fmin2(1.0, fmax2(0.0, thF + dt * (aPf * f - (ef + (p.flatOnly ? kdes : 0.0)) * thF)));
                    const double B = fmin2(1.0 - P, fmax2(0.0, thFB + dt * (aBf * f - kB * thFB)));
                    thF = P; thFB = B;
                    if (holes) {
                        double* th = &thetaHole[static_cast<size_t>(facei) * N];
                        double* thB = &thetaHoleB[static_cast<size_t>(facei) * N];
                        for (int i = 0; i < N; ++i) {
                            const double fi = 1.0 - th[i] - thB[i];
                            const double aPi = (vbar / 4.0) * s * c0 * Tf * nr[i] / N0;
                            const double aBi = (vbarB / 4.0) * sB * nb[i] / N0;
                            const double ei = (vbarB / 4.0) * kE * nb[i] / N0;
                            const double Pi = fmin2(1.0, fmax2(0.0, th[i] + dt * (aPi * fi - ei * th[i])));
                            const double Bi = fmin2(1.0 - Pi, fmax2(0.0, thB[i] + dt * (aBi * fi - kB * thB[i])));
                            th[i] = Pi; thB[i] = Bi;
                        }
                    }
                } else {
                    advance2(thF, thFB, aPf, aBf, ef + (p.flatOnly ? kdes : 0.0), kB, dt);
                    if (holes) {
                        double* th = &thetaHole[static_cast<size_t>(facei) * N];
                        double* thB = &thetaHoleB[static_cast<size_t>(facei) * N];
                        for (int i = 0; i < N; ++i) {
                            const double aPi = (vbar / 4.0) * s * c0 * Tf * nr[i] / N0;
                            const double aBi = (vbarB / 4.0) * sB * nb[i] / N0;
                            const double ei = (vbarB / 4.0) * kE * nb[i] / N0;
                            advance2(th[i], thB[i], aPi, aBi, ei, kB, dt);
                        }
                    }
                }
            }
            const double freeF = 1.0 - thF - thFB;
            const double Kflat = (vbar / 4.0) * (s * freeF + r);
            const double K = (1.0 - phi) * Kflat + phi * jhole;
            Keff[facei] = K;
            if (p.flatOnly) {
                const double JdesT = kdes * thF * N0 / c0;
                refValue[facei] = JdesT / fmax2(K, kSmall);
                Jdes[facei] = kdes * thF * N0;
            } else {
                refValue[facei] = 0.0;
                Jdes[facei] = 0.0;
            }
            // net byproduct flux into the chamber per unit wafer area
            const double flatByp = nu * s * (vbar / 4.0) * freeF * c0 * Tf + kB * N0 * thFB
                                   - (vbarB / 4.0) * (sB * freeF + kE * thF) * nBf;
            Jbyp[facei] = (1.0 - phi) * flatByp + phi * b5[2];
            if (keff_out) keff_out[facei] = K;
            if (flux_out) flux_out[facei] = K * c0 * Tf - Jdes[facei];
        }
    }

    // One pass over the faces: evaluate (and, with advance, update the coverages first as the BCs do).
    void update(const double* Tface, double dt, bool advance, double* flux_out, double* keff_out) {
        if (twoSpecies) { updateTwo(Tface, dt, advance, flux_out, keff_out); return; }
        const int N = p.nCells;
        const double vbar = p.vbar, s = p.s, r = p.r, c0 = p.c0, N0 = p.N0, kdes = p.kdes;
        const bool expo = (p.scheme == ALDKIN_SCHEME_EXPONENTIAL);

        if (p.flatOnly) {
            // [BC] langmuirStickingFvPatchScalarField::updateCoeffs
            for (int facei = 0; facei < nFaces; ++facei) {
                double theta = thetaFlat[facei];
                const double Tf = fmax2(Tface[facei], 0.0);
                if (advance) {
                    if (!expo) {
                        // [BC] theta_ = min(1, max(0, theta_ + dt*((vbar_/4.0)*s_*(1.0 - theta_)*c0_*max(Tface, 0)/N0_ - kdes_*theta_)))
                        theta = fmin2(1.0, fmax2(0.0,
                            theta + dt * ((vbar / 4.0) * s * (1.0 - theta) * c0 * Tf / N0 - kdes * theta)));
                    } else {
                        // exact relaxation of dtheta/dt = a (1 - theta) - kdes theta at frozen gas
                        const double a = (vbar / 4.0) * s * c0 * Tf / N0;
                        const double tot = a + kdes;
                        if (tot > 0.0) {
                            const double thinf = a / tot;
                            theta = thinf + (theta - thinf) * std::exp(-tot * dt);
                        }
                    }
                    thetaFlat[facei] = theta;
                }
                // [BC] K = (vbar_/4.0)*(s_*(1.0 - theta_) + r_);  Jdes = kdes_*theta_*N0_/c0_;  refValue = Jdes/max(K, SMALL)
                const double K = (vbar / 4.0) * (s * (1.0 - theta) + r);
                const double JdesT = kdes * theta * N0 / c0;          // in T units (per c0)
                Keff[facei] = K;
                refValue[facei] = JdesT / fmax2(K, kSmall);
                Jdes[facei] = kdes * theta * N0;
                Jbyp[facei] = (p.nu > 0.0) ? p.nu * (vbar / 4.0) * s * (1.0 - theta) * c0 * Tf : 0.0;
                if (keff_out) keff_out[facei] = K;
                if (flux_out) flux_out[facei] = K * c0 * Tf - Jdes[facei];
            }
            return;
        }

        // [BC] patternedSurfaceFvPatchScalarField::updateCoeffs
        for (int facei = 0; facei < nFaces; ++facei) {
            const double Tf = fmax2(Tface[facei], 0.0);
            const bool holes = hasHoles(facei);
            double* nr = &nrel[static_cast<size_t>(facei) * N];
            // open-area fraction of this face; 0 on a bare (flat Langmuir) face
            const double phi = holes ? holeDensityField[facei] * kPi * p.d * p.d / 4.0 : 0.0;
            double jhole = 0.0, hole_prod = 0.0;
            if (holes) {
                jhole = solveHole(facei, nr);                                  // per unit mouth density
                if (p.nu > 0.0) {
                    // byproduct reporting only (no re-adsorption on this path: everything produced escapes), from the
                    // coverages the gas was solved with; does not enter K_eff or the coverages
                    const double dx = ARField[facei] * p.d / N;
                    const double* th = &thetaHole[static_cast<size_t>(facei) * N];
                    for (int i = 0; i < N; ++i) hole_prod += (4.0 / p.d) * p.nu * s * (vbar / 4.0) * (1.0 - th[i]) * c0 * Tf * nr[i] * dx;
                }
            } else {
                for (int i = 0; i < N; ++i) nr[i] = 0.0;
            }
            if (advance) {
                double& thF = thetaFlat[facei];
                if (!expo) {
                    // [BC] thetaFlat_ = min(1, thetaFlat_ + dt*(vbar_/4.0)*s_*(1.0 - thetaFlat_)*c0_*Tf/N0_)
                    thF = fmin2(1.0, thF + dt * (vbar / 4.0) * s * (1.0 - thF) * c0 * Tf / N0);
                    if (holes) {
                        double* th = &thetaHole[static_cast<size_t>(facei) * N];
                        for (int i = 0; i < N; ++i) {
                            // [BC] th = min(1, th + dt*(vbar_/4.0)*s_*(1.0 - th)*c0_*Tf*nrel[i]/N0_)
                            th[i] = fmin2(1.0, th[i] + dt * (vbar / 4.0) * s * (1.0 - th[i]) * c0 * Tf * nr[i] / N0);
                        }
                    }
                } else {
                    // aldconf.solver: theta <- 1 - (1 - theta) exp(-s (vbar/4) n dt / N0)
                    thF = 1.0 - (1.0 - thF) * std::exp(-s * (vbar / 4.0) * (c0 * Tf) * dt / N0);
                    if (holes) {
                        double* th = &thetaHole[static_cast<size_t>(facei) * N];
                        for (int i = 0; i < N; ++i) {
                            th[i] = 1.0 - (1.0 - th[i]) * std::exp(-s * (vbar / 4.0) * (c0 * Tf * nr[i]) * dt / N0);
                        }
                    }
                }
            }
            // [BC] Keff_ = (1.0 - phi)*(vbar_/4.0)*(s_*(1.0 - thetaFlat_) + r_) + phi*jhole
            const double K = (1.0 - phi) * (vbar / 4.0) * (s * (1.0 - thetaFlat[facei]) + r) + phi * jhole;
            Keff[facei] = K;
            refValue[facei] = 0.0;
            Jdes[facei] = 0.0;
            if (p.nu > 0.0) {
                double* b5 = &bal[static_cast<size_t>(facei) * 5];
                b5[0] = hole_prod; b5[1] = 0.0; b5[2] = hole_prod; b5[3] = 0.0; b5[4] = 0.0;
                Jbyp[facei] = (1.0 - phi) * p.nu * s * (vbar / 4.0) * (1.0 - thetaFlat[facei]) * c0 * Tf + phi * hole_prod;
            } else {
                Jbyp[facei] = 0.0;
            }
            if (keff_out) keff_out[facei] = K;
            if (flux_out) flux_out[facei] = K * c0 * Tf;
        }
    }
};

namespace {

bool finite_nonneg(double x) { return std::isfinite(x) && x >= 0.0; }

const char* validate(const aldkin_params& p) {
    if (!finite_nonneg(p.s) || p.s > 1.0) return "s must be in [0, 1]";
    if (!finite_nonneg(p.r) || p.r > 1.0) return "r must be in [0, 1]";
    if (!finite_nonneg(p.kdes)) return "kdes must be >= 0";
    if (!std::isfinite(p.vbar) || p.vbar <= 0.0) return "vbar must be > 0";
    if (!std::isfinite(p.N0) || p.N0 <= 0.0) return "N0 must be > 0";
    if (!std::isfinite(p.c0) || p.c0 <= 0.0) return "c0 must be > 0";
    if (p.flatOnly != 0 && p.flatOnly != 1) return "flatOnly must be 0 or 1";
    if (p.scheme != ALDKIN_SCHEME_EULER && p.scheme != ALDKIN_SCHEME_EXPONENTIAL) return "unknown scheme";
    if (!finite_nonneg(p.nu)) return "nu must be >= 0";
    if (!finite_nonneg(p.vbarB)) return "vbarB must be >= 0";
    if (!finite_nonneg(p.sB) || p.sB > 1.0) return "sB must be in [0, 1]";
    if (!finite_nonneg(p.kB)) return "kB must be >= 0";
    if (!finite_nonneg(p.kEtch) || p.kEtch > 1.0) return "kEtch must be in [0, 1]";
    if (!finite_nonneg(p.cB0)) return "cB0 must be >= 0";
    if ((p.sB > 0.0 || p.kEtch > 0.0) && p.vbarB <= 0.0) return "vbarB must be > 0 when sB > 0 or kEtch > 0 (two-species option)";
    if (!p.flatOnly) {
        if (p.kdes != 0.0) return "kdes is only available with flatOnly = 1 (the patterned model has no desorption term)";
        if (!std::isfinite(p.d) || p.d <= 0.0) return "d must be > 0";
        if (!std::isfinite(p.AR)) return "AR must be finite";
        if (!std::isfinite(p.holeDensity)) return "holeDensity must be finite";
        if (p.nCells < 1) return "nCells must be >= 1";
    }
    return nullptr;
}

template <class F>
int guarded(aldkin_t* h, F&& f) {
    try {
        return f();
    } catch (const std::bad_alloc&) {
        return h ? h->fail(ALDKIN_ERR_ALLOC, "allocation failed") : ALDKIN_ERR_ALLOC;
    } catch (const std::exception& e) {
        return h ? h->fail(ALDKIN_ERR_INTERNAL, std::string("internal error: ") + e.what()) : ALDKIN_ERR_INTERNAL;
    } catch (...) {
        return h ? h->fail(ALDKIN_ERR_INTERNAL, "internal error") : ALDKIN_ERR_INTERNAL;
    }
}

}  // namespace

extern "C" {

const char* aldkin_version(void) { return version_string; }

void aldkin_params_default(aldkin_params* p) {
    if (!p) return;
    // [BC] the dictionary-free constructors: s 1, r 0, kdes 0, vbar 410, N0 4e18, c0 1, d 1e-7, AR 50, holeDensity 0, nCells 100
    p->s = 1.0; p->r = 0.0; p->kdes = 0.0; p->vbar = 410.0; p->N0 = 4e18; p->c0 = 1.0;
    p->d = 1e-7; p->AR = 50.0; p->holeDensity = 0.0; p->nCells = 100;
    p->flatOnly = 0; p->scheme = ALDKIN_SCHEME_EULER;
    p->nu = 0.0; p->vbarB = 0.0; p->sB = 0.0; p->kB = 0.0; p->kEtch = 0.0; p->cB0 = 0.0;
}

aldkin_t* aldkin_create(const aldkin_params* p) {
    try {
        if (!p) { g_create_error = "aldkin_create: params is NULL"; return nullptr; }
        if (const char* m = validate(*p)) { g_create_error = std::string("aldkin_create: ") + m; return nullptr; }
        aldkin_t* h = new aldkin_t();
        h->p = *p;
        if (h->p.flatOnly) h->p.nCells = (h->p.nCells < 1) ? 1 : h->p.nCells;
        h->twoSpecies = (h->p.sB > 0.0 || h->p.kEtch > 0.0);
        return h;
    } catch (...) {
        g_create_error = "aldkin_create: allocation failed";
        return nullptr;
    }
}

aldkin_t* aldkin_clone(const aldkin_t* h) {
    try {
        if (!h) { g_create_error = "aldkin_clone: handle is NULL"; return nullptr; }
        return new aldkin_t(*h);
    } catch (...) {
        g_create_error = "aldkin_clone: allocation failed";
        return nullptr;
    }
}

void aldkin_destroy(aldkin_t* h) { delete h; }

const char* aldkin_last_error(const aldkin_t* h) { return h ? h->err.c_str() : g_create_error.c_str(); }

int aldkin_get_params(const aldkin_t* h, aldkin_params* p) {
    if (!h) return ALDKIN_ERR_ARG;
    if (!p) return const_cast<aldkin_t*>(h)->fail(ALDKIN_ERR_ARG, "aldkin_get_params: p is NULL");
    *p = h->p;
    return ALDKIN_OK;
}

int aldkin_set_faces(aldkin_t* h, int n_faces, const double* radius_m) {
    if (!h) return ALDKIN_ERR_ARG;
    return guarded(h, [&]() -> int {
        if (n_faces < 0) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_faces: n_faces < 0");
        const size_t n = static_cast<size_t>(n_faces), nc = static_cast<size_t>(h->p.nCells);
        h->nFaces = n_faces;
        h->facesSet = true;
        h->radius.assign(n, 0.0);
        if (radius_m) {
            for (size_t i = 0; i < n; ++i) {
                if (!std::isfinite(radius_m[i])) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_faces: radius not finite");
                h->radius[i] = radius_m[i];
            }
        }
        h->thetaFlat.assign(n, 0.0);
        h->thetaHole.assign(n * nc, 0.0);
        h->nrel.assign(n * nc, 0.0);
        h->Keff.assign(n, 0.0);
        h->refValue.assign(n, 0.0);
        h->Jdes.assign(n, 0.0);
        h->thetaFlatB.assign(n, 0.0);
        h->thetaHoleB.assign(n * nc, 0.0);
        h->nB.assign(n * nc, 0.0);
        h->nBface.assign(n, h->p.cB0);
        h->Jbyp.assign(n, 0.0);
        h->bal.assign(n * 5, 0.0);
        h->setPatternFields();
        return ALDKIN_OK;
    });
}

int aldkin_set_profile(aldkin_t* h, int which, int n, const double* r, const double* values) {
    if (!h) return ALDKIN_ERR_ARG;
    return guarded(h, [&]() -> int {
        if (which != ALDKIN_PROFILE_AR && which != ALDKIN_PROFILE_DENSITY)
            return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profile: which must be ALDKIN_PROFILE_AR or ALDKIN_PROFILE_DENSITY");
        if (n < 0) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profile: n < 0");
        std::vector<std::pair<double, double>> tab;
        if (n > 0) {
            if (!r || !values) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profile: r or values is NULL");
            for (int i = 0; i < n; ++i) {
                if (!std::isfinite(r[i])) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profile: radius not finite");
                if (i > 0 && r[i] < r[i - 1]) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profile: radii must be in ascending order");
                if (!std::isfinite(values[i])) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profile: value not finite");
            }
            tab.reserve(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) tab.emplace_back(r[i], values[i]);
        }
        if (which == ALDKIN_PROFILE_AR) h->ARProfile = std::move(tab); else h->holeDensityProfile = std::move(tab);
        h->setPatternFields();
        return ALDKIN_OK;
    });
}

int aldkin_set_profiles(aldkin_t* h, int n, const double* r, const double* AR, const double* density) {
    if (!h) return ALDKIN_ERR_ARG;
    if (n < 0) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profiles: n < 0");
    if (n > 0 && !r) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_profiles: r is NULL");
    int rc = aldkin_set_profile(h, ALDKIN_PROFILE_AR, AR ? n : 0, r, AR);
    if (rc != ALDKIN_OK) return rc;
    return aldkin_set_profile(h, ALDKIN_PROFILE_DENSITY, density ? n : 0, r, density);
}

int aldkin_n_faces(const aldkin_t* h) { return h ? h->nFaces : -1; }
int aldkin_n_cells(const aldkin_t* h) { return h ? h->p.nCells : -1; }

int aldkin_face_pattern(const aldkin_t* h, double* AR_out, double* density_out) {
    if (!h) return ALDKIN_ERR_ARG;
    const size_t n = static_cast<size_t>(h->nFaces);
    if (AR_out && n) std::memcpy(AR_out, h->ARField.data(), n * sizeof(double));
    if (density_out && n) std::memcpy(density_out, h->holeDensityField.data(), n * sizeof(double));
    return ALDKIN_OK;
}

namespace {
int do_update(aldkin_t* h, const double* T_face, double dt, bool advance, double* flux, double* keff, const char* who) {
    if (!h) return ALDKIN_ERR_ARG;
    return guarded(h, [&]() -> int {
        if (!h->facesSet) return h->fail(ALDKIN_ERR_STATE, std::string(who) + ": aldkin_set_faces has not been called");
        if (h->nFaces > 0 && !T_face) return h->fail(ALDKIN_ERR_ARG, std::string(who) + ": T_face is NULL");
        if (advance && !(std::isfinite(dt) && dt >= 0.0)) return h->fail(ALDKIN_ERR_ARG, std::string(who) + ": dt must be finite and >= 0");
        for (int i = 0; i < h->nFaces; ++i) {
            if (!std::isfinite(T_face[i])) return h->fail(ALDKIN_ERR_ARG, std::string(who) + ": T_face not finite");
        }
        h->update(T_face, dt, advance, flux, keff);
        return ALDKIN_OK;
    });
}
}  // namespace

int aldkin_step(aldkin_t* h, const double* T_face, double dt, double* uptake_flux_out, double* keff_out) {
    return do_update(h, T_face, dt, true, uptake_flux_out, keff_out, "aldkin_step");
}

int aldkin_evaluate(aldkin_t* h, const double* T_face, double* uptake_flux_out, double* keff_out) {
    return do_update(h, T_face, 0.0, false, uptake_flux_out, keff_out, "aldkin_evaluate");
}

int aldkin_theta(const aldkin_t* h, int which, double* out) {
    if (!h) return ALDKIN_ERR_ARG;
    aldkin_t* hh = const_cast<aldkin_t*>(h);
    if (!out) return hh->fail(ALDKIN_ERR_ARG, "aldkin_theta: out is NULL");
    const int n = h->nFaces, N = h->p.nCells;
    switch (which) {
        case ALDKIN_THETA_FLAT:
            for (int i = 0; i < n; ++i) out[i] = h->thetaFlat[i];
            return ALDKIN_OK;
        case ALDKIN_THETA_MID:
            for (int i = 0; i < n; ++i) out[i] = h->p.flatOnly ? 0.0 : h->thetaHole[static_cast<size_t>(i) * N + N / 2];
            return ALDKIN_OK;
        case ALDKIN_THETA_BOTTOM:
            for (int i = 0; i < n; ++i) out[i] = h->p.flatOnly ? 0.0 : h->thetaHole[static_cast<size_t>(i) * N + N - 1];
            return ALDKIN_OK;
        case ALDKIN_THETA_HOLE:
            for (size_t i = 0; i < static_cast<size_t>(n) * N; ++i) out[i] = h->p.flatOnly ? 0.0 : h->thetaHole[i];
            return ALDKIN_OK;
        case ALDKIN_THETA_BLOCKED_FLAT:
            for (int i = 0; i < n; ++i) out[i] = h->thetaFlatB[i];
            return ALDKIN_OK;
        case ALDKIN_THETA_BLOCKED_BOTTOM:
            for (int i = 0; i < n; ++i) out[i] = h->p.flatOnly ? 0.0 : h->thetaHoleB[static_cast<size_t>(i) * N + N - 1];
            return ALDKIN_OK;
        case ALDKIN_THETA_BLOCKED_HOLE:
            for (size_t i = 0; i < static_cast<size_t>(n) * N; ++i) out[i] = h->p.flatOnly ? 0.0 : h->thetaHoleB[i];
            return ALDKIN_OK;
        default:
            return hh->fail(ALDKIN_ERR_ARG, "aldkin_theta: unknown selector");
    }
}

int aldkin_set_theta(aldkin_t* h, int which, const double* in) {
    if (!h) return ALDKIN_ERR_ARG;
    if (!in) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_theta: in is NULL");
    const int n = h->nFaces, N = h->p.nCells;
    switch (which) {
        case ALDKIN_THETA_FLAT:
            for (int i = 0; i < n; ++i) h->thetaFlat[i] = in[i];
            return ALDKIN_OK;
        case ALDKIN_THETA_HOLE:
            if (h->p.flatOnly) return h->fail(ALDKIN_ERR_STATE, "aldkin_set_theta: flatOnly handle has no hole coverages");
            for (size_t i = 0; i < static_cast<size_t>(n) * N; ++i) h->thetaHole[i] = in[i];
            return ALDKIN_OK;
        case ALDKIN_THETA_BLOCKED_FLAT:
            for (int i = 0; i < n; ++i) h->thetaFlatB[i] = in[i];
            return ALDKIN_OK;
        case ALDKIN_THETA_BLOCKED_HOLE:
            if (h->p.flatOnly) return h->fail(ALDKIN_ERR_STATE, "aldkin_set_theta: flatOnly handle has no hole coverages");
            for (size_t i = 0; i < static_cast<size_t>(n) * N; ++i) h->thetaHoleB[i] = in[i];
            return ALDKIN_OK;
        default:
            return h->fail(ALDKIN_ERR_ARG, "aldkin_set_theta: selector must be ALDKIN_THETA_FLAT, _HOLE, _BLOCKED_FLAT or _BLOCKED_HOLE");
    }
}

int aldkin_keff(const aldkin_t* h, double* out) {
    if (!h) return ALDKIN_ERR_ARG;
    if (!out) return const_cast<aldkin_t*>(h)->fail(ALDKIN_ERR_ARG, "aldkin_keff: out is NULL");
    for (int i = 0; i < h->nFaces; ++i) out[i] = h->Keff[i];
    return ALDKIN_OK;
}

int aldkin_mixed_bc(const aldkin_t* h, const double* D_over_delta, double* refValue, double* refGrad, double* valueFraction) {
    if (!h) return ALDKIN_ERR_ARG;
    if (h->nFaces > 0 && !D_over_delta) return const_cast<aldkin_t*>(h)->fail(ALDKIN_ERR_ARG, "aldkin_mixed_bc: D_over_delta is NULL");
    for (int i = 0; i < h->nFaces; ++i) {
        if (refValue) refValue[i] = h->refValue[i];
        if (refGrad) refGrad[i] = 0.0;
        // [BC] valueFraction() = Keff_/(Keff_ + Dd)   (langmuir: K/(K + Dd))
        if (valueFraction) valueFraction[i] = h->Keff[i] / (h->Keff[i] + D_over_delta[i]);
    }
    return ALDKIN_OK;
}

int aldkin_desorption_flux(const aldkin_t* h, double* out) {
    if (!h) return ALDKIN_ERR_ARG;
    if (!out) return const_cast<aldkin_t*>(h)->fail(ALDKIN_ERR_ARG, "aldkin_desorption_flux: out is NULL");
    for (int i = 0; i < h->nFaces; ++i) out[i] = h->Jdes[i];
    return ALDKIN_OK;
}

int aldkin_set_byproduct_mouth(aldkin_t* h, const double* nB_face) {
    if (!h) return ALDKIN_ERR_ARG;
    return guarded(h, [&]() -> int {
        if (!h->facesSet) return h->fail(ALDKIN_ERR_STATE, "aldkin_set_byproduct_mouth: aldkin_set_faces has not been called");
        if (!nB_face) { h->nBface.assign(static_cast<size_t>(h->nFaces), h->p.cB0); return ALDKIN_OK; }
        for (int i = 0; i < h->nFaces; ++i) {
            if (!std::isfinite(nB_face[i])) return h->fail(ALDKIN_ERR_ARG, "aldkin_set_byproduct_mouth: value not finite");
        }
        for (int i = 0; i < h->nFaces; ++i) h->nBface[i] = fmax2(nB_face[i], 0.0);
        return ALDKIN_OK;
    });
}

int aldkin_byproduct_flux(const aldkin_t* h, double* out) {
    if (!h) return ALDKIN_ERR_ARG;
    if (!out) return const_cast<aldkin_t*>(h)->fail(ALDKIN_ERR_ARG, "aldkin_byproduct_flux: out is NULL");
    for (int i = 0; i < h->nFaces; ++i) out[i] = h->Jbyp[i];
    return ALDKIN_OK;
}

int aldkin_byproduct_balance(const aldkin_t* h, int face, double* out5) {
    if (!h) return ALDKIN_ERR_ARG;
    aldkin_t* hh = const_cast<aldkin_t*>(h);
    if (!out5) return hh->fail(ALDKIN_ERR_ARG, "aldkin_byproduct_balance: out5 is NULL");
    if (face < 0 || face >= h->nFaces) return hh->fail(ALDKIN_ERR_ARG, "aldkin_byproduct_balance: face out of range");
    for (int k = 0; k < 5; ++k) out5[k] = h->bal[static_cast<size_t>(face) * 5 + k];
    return ALDKIN_OK;
}

int aldkin_hole_byproduct_gas(const aldkin_t* h, int face, double* nB_out) {
    if (!h) return ALDKIN_ERR_ARG;
    aldkin_t* hh = const_cast<aldkin_t*>(h);
    if (!nB_out) return hh->fail(ALDKIN_ERR_ARG, "aldkin_hole_byproduct_gas: nB_out is NULL");
    if (face < 0 || face >= h->nFaces) return hh->fail(ALDKIN_ERR_ARG, "aldkin_hole_byproduct_gas: face out of range");
    const int N = h->p.nCells;
    for (int i = 0; i < N; ++i) nB_out[i] = h->p.flatOnly ? 0.0 : h->nB[static_cast<size_t>(face) * N + i];
    return ALDKIN_OK;
}

int aldkin_hole_gas(const aldkin_t* h, int face, double* nrel_out) {
    if (!h) return ALDKIN_ERR_ARG;
    aldkin_t* hh = const_cast<aldkin_t*>(h);
    if (!nrel_out) return hh->fail(ALDKIN_ERR_ARG, "aldkin_hole_gas: nrel_out is NULL");
    if (face < 0 || face >= h->nFaces) return hh->fail(ALDKIN_ERR_ARG, "aldkin_hole_gas: face out of range");
    const int N = h->p.nCells;
    for (int i = 0; i < N; ++i) nrel_out[i] = h->p.flatOnly ? 0.0 : h->nrel[static_cast<size_t>(face) * N + i];
    return ALDKIN_OK;
}

}  // extern "C"
