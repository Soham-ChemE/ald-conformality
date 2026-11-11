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

const char* version_string = "1.0.0";

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
    double solveHole(int facei, double* nr) {
        const int N = p.nCells;
        const double Lh = ARField[facei] * p.d;
        const double dx = Lh / N;
        const double DK = p.d * p.vbar / 3.0;              // Knudsen diffusivity, cylinder
        const double a = DK / (dx * dx);
        lower.assign(N, -a); diag.assign(N, 0.0); upper.assign(N, -a); rhs.assign(N, 0.0);
        const double* th = &thetaHole[static_cast<size_t>(facei) * N];
        for (int i = 0; i < N; ++i) {
            const double k = (p.vbar / p.d) * (p.s * (1.0 - th[i]) + p.r);   // (4/d)(vbar/4)[...]
            diag[i] = 2.0 * a + k;
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

    // One pass over the faces: evaluate (and, with advance, update the coverages first as the BCs do).
    void update(const double* Tface, double dt, bool advance, double* flux_out, double* keff_out) {
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
            double jhole = 0.0;
            if (holes) {
                jhole = solveHole(facei, nr);                                  // per unit mouth density
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
}

aldkin_t* aldkin_create(const aldkin_params* p) {
    try {
        if (!p) { g_create_error = "aldkin_create: params is NULL"; return nullptr; }
        if (const char* m = validate(*p)) { g_create_error = std::string("aldkin_create: ") + m; return nullptr; }
        aldkin_t* h = new aldkin_t();
        h->p = *p;
        if (h->p.flatOnly) h->p.nCells = (h->p.nCells < 1) ? 1 : h->p.nCells;
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
        default:
            return h->fail(ALDKIN_ERR_ARG, "aldkin_set_theta: selector must be ALDKIN_THETA_FLAT or ALDKIN_THETA_HOLE");
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
