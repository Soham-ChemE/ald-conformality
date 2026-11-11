/*
 * aldkin: multiscale ALD wall kinetics as a standalone library with a C ABI.
 *
 * Every wall face carries
 *   (a) a flat Langmuir surface with coverage thetaFlat, and
 *   (b) optionally a quasi-steady 1D Knudsen-diffusion model of a cylindrical hole (diameter d, aspect ratio AR,
 *       nCells cells, Langmuir sink on the wall) whose mouth sees the face gas density c0 * T_face.
 * Net uptake per unit wafer area is  j = K_eff * c0 * T_face,
 *   K_eff = (1 - phi) (vbar/4) [ s (1 - thetaFlat) + r ]  +  phi * j_hole / (c0 T_face),
 * phi = holeDensity * pi d^2 / 4 the open-area fraction and j_hole the Knudsen flux into the hole mouth. A face with
 * holeDensity 0 or AR 0 is a flat Langmuir surface (phi = 0, no hole solved). With flatOnly = 1 the hole model is
 * absent and a first-order desorption rate kdes may be given (the langmuirSticking model).
 *
 * Units. Gas densities enter and leave as the dimensionless T = n / c0 (c0 in 1/m3 is a parameter; pass c0 = 1 to
 * work with absolute number densities in 1/m3). vbar in m/s, N0 in 1/m2, d in m, holeDensity in 1/m2, kdes in 1/s,
 * dt in s. K_eff is a velocity in m/s; the uptake flux is in molecules/(m2 s).
 *
 * Call sequence per host solver time step (see README.md):
 *   once:        h = aldkin_create(&p); aldkin_set_faces(h, n, radius); [aldkin_set_profiles(h, ...)]
 *   per step:    aldkin_step(h, T_face, dt, flux, keff)        advances coverages, then evaluates K_eff
 *   per inner iteration (optional): aldkin_evaluate(h, T_face, flux, keff)   re-evaluates with the new T_face,
 *                coverages unchanged
 *   mixed BC:    aldkin_mixed_bc(h, D_over_delta, refValue, refGrad, valueFraction)
 *                   gives  refValue = J_des / max(K, 1e-15)  (0 without desorption), refGrad = 0,
 *                   valueFraction = K_eff / (K_eff + D/delta),   the exact quantities the OpenFOAM mixed
 *                   boundary condition uses:  T_face = f * refValue + (1 - f) * (T_cell + refGrad/deltaCoeffs).
 *                   Equivalently, for any solver: the wall flux balance  K_eff c_face - J_des = D (c_cell - c_face)/delta.
 *
 * Numerics are those of the OpenFOAM boundary conditions this library was extracted from (same discretisation of the
 * hole, same forward-Euler coverage update once per step, same blending, same clamping), so a host that calls
 * aldkin_step once per time step with the previous step's face value reproduces them bit for bit
 * (scheme ALDKIN_SCHEME_EULER). ALDKIN_SCHEME_EXPONENTIAL instead uses the exact exponential coverage update of
 * the Python twin aldconf.solver (theta <- 1 - (1 - theta) exp(-a dt)), unconditionally stable for any dt.
 *
 * Errors are return codes (ALDKIN_OK = 0); aldkin_last_error(h) gives the message of the last failure on that handle
 * (aldkin_last_error(NULL) the message of the last failed aldkin_create in this thread). No exception crosses the ABI.
 * Handles are independent; one handle must not be used from two threads at once.
 */
#ifndef ALDKIN_H
#define ALDKIN_H

#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(ALDKIN_BUILD_DLL)
#    define ALDKIN_API __declspec(dllexport)
#  elif defined(ALDKIN_DLL)
#    define ALDKIN_API __declspec(dllimport)
#  else
#    define ALDKIN_API
#  endif
#else
#  if defined(__GNUC__) || defined(__clang__)
#    define ALDKIN_API __attribute__((visibility("default")))
#  else
#    define ALDKIN_API
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ALDKIN_VERSION_MAJOR 1
#define ALDKIN_VERSION_MINOR 0
#define ALDKIN_VERSION_PATCH 0

/* return codes */
enum {
    ALDKIN_OK = 0,
    ALDKIN_ERR_ARG = 1,      /* null pointer, bad size or index */
    ALDKIN_ERR_PARAM = 2,    /* invalid physical parameters */
    ALDKIN_ERR_STATE = 3,    /* call out of sequence (e.g. step before set_faces) */
    ALDKIN_ERR_ALLOC = 4,    /* memory allocation failed */
    ALDKIN_ERR_INTERNAL = 5
};

/* coverage update scheme */
enum {
    ALDKIN_SCHEME_EULER = 0,        /* forward Euler, clamped to [0, 1]: the OpenFOAM boundary conditions */
    ALDKIN_SCHEME_EXPONENTIAL = 1   /* exact exponential relaxation at frozen gas: the aldconf Python solver */
};

/* selectors for aldkin_theta / aldkin_set_theta */
enum {
    ALDKIN_THETA_FLAT = 0,     /* n_faces values: flat-surface coverage (the only coverage when flatOnly = 1) */
    ALDKIN_THETA_MID = 1,      /* n_faces values: hole coverage at half depth (cell nCells/2) */
    ALDKIN_THETA_BOTTOM = 2,   /* n_faces values: hole coverage at the bottom (cell nCells - 1) */
    ALDKIN_THETA_HOLE = 3      /* n_faces * nCells values, cell index fastest: the whole hole profile */
};

/* selectors for aldkin_set_profile */
enum {
    ALDKIN_PROFILE_AR = 0,        /* aspect ratio versus radius */
    ALDKIN_PROFILE_DENSITY = 1    /* holes per m2 versus radius */
};

typedef struct aldkin_t aldkin_t;

typedef struct aldkin_params {
    double s;            /* sticking coefficient on an open site [-] */
    double r;            /* loss probability on every collision regardless of coverage (radical recombination) [-] */
    double kdes;         /* first-order desorption rate of the adsorbed layer [1/s]; flatOnly = 1 only; 0 = irreversible */
    double vbar;         /* mean thermal speed of the precursor [m/s] */
    double N0;           /* adsorption sites filled by one saturated half-cycle [1/m2] */
    double c0;           /* precursor number density at T = 1 [1/m3]; 1 to pass absolute densities */
    double d;            /* hole diameter [m] */
    double AR;           /* hole aspect ratio depth/d, wafer-wide default (0 = no holes) */
    double holeDensity;  /* holes per m2 of wafer, wafer-wide default (0 = no holes) */
    int nCells;          /* cells along the hole depth (>= 1) */
    int flatOnly;        /* 1: flat Langmuir surface only (langmuirSticking model, kdes allowed); 0: patterned */
    int scheme;          /* ALDKIN_SCHEME_EULER or ALDKIN_SCHEME_EXPONENTIAL */
} aldkin_params;

/* "major.minor.patch" */
ALDKIN_API const char* aldkin_version(void);

/* Fill p with the defaults of the OpenFOAM boundary conditions: s 1, r 0, kdes 0, vbar 410, N0 4e18, c0 1, d 1e-7,
 * AR 50, holeDensity 0, nCells 100, flatOnly 0, scheme Euler. */
ALDKIN_API void aldkin_params_default(aldkin_params* p);

/* New handle with no faces. NULL on failure (see aldkin_last_error(NULL)). */
ALDKIN_API aldkin_t* aldkin_create(const aldkin_params* p);

/* Deep copy of a handle including faces, profiles and coverages. NULL on failure. */
ALDKIN_API aldkin_t* aldkin_clone(const aldkin_t* h);

ALDKIN_API void aldkin_destroy(aldkin_t* h);

/* Message of the last failed call on h; with h = NULL the message of the last failed aldkin_create in this thread. */
ALDKIN_API const char* aldkin_last_error(const aldkin_t* h);

/* Read back the parameters the handle was created with. */
ALDKIN_API int aldkin_get_params(const aldkin_t* h, aldkin_params* p);

/* Define n_faces faces. radius_m[i] is the face-centre radius (m) used to look up the radial profiles; it may be
 * NULL (all radii 0) when no profiles are used. Resets every coverage to 0 and K_eff to 0. Can be called again. */
ALDKIN_API int aldkin_set_faces(aldkin_t* h, int n_faces, const double* radius_m);

/* Optional radial table for one quantity (ALDKIN_PROFILE_AR or ALDKIN_PROFILE_DENSITY): n (radius, value) pairs,
 * radii ascending, overriding the scalar AR or holeDensity face by face, linearly interpolated at the face radius
 * and clamped to the end values outside the table (a repeated radius makes a step); negative interpolated values
 * are clamped to 0. n = 0 removes that table (back to the scalar). The two tables are independent, as the OpenFOAM
 * entries ARProfile and holeDensityProfile are. The per-face values are recomputed whenever faces or profiles
 * change; coverages are not touched. */
ALDKIN_API int aldkin_set_profile(aldkin_t* h, int which, int n, const double* r, const double* values);

/* Both tables on one radius vector: AR and/or density may be NULL (that table is removed, the scalar is used);
 * n = 0 removes both. Equivalent to two aldkin_set_profile calls. */
ALDKIN_API int aldkin_set_profiles(aldkin_t* h, int n, const double* r, const double* AR, const double* density);

ALDKIN_API int aldkin_n_faces(const aldkin_t* h);
ALDKIN_API int aldkin_n_cells(const aldkin_t* h);

/* Per-face aspect ratio and hole density actually used (n_faces values each; either pointer may be NULL). */
ALDKIN_API int aldkin_face_pattern(const aldkin_t* h, double* AR_out, double* density_out);

/* Advance every face by one time step dt given T_face[i] = n_face/c0 (values < 0 are treated as 0), then evaluate.
 * Order of operations, as in the boundary conditions: the hole gas is solved with the coverages before the update,
 * the coverages are updated with that gas, K_eff is formed with the updated flat coverage and the pre-update hole
 * flux. Outputs (n_faces each, either may be NULL):
 *   uptake_flux_out  net uptake per unit wafer area, K_eff c0 T_face - J_des      [molecules/(m2 s)]
 *   keff_out         effective rate constant K_eff                                 [m/s]            */
ALDKIN_API int aldkin_step(aldkin_t* h, const double* T_face, double dt, double* uptake_flux_out, double* keff_out);

/* Evaluate K_eff and the flux for new face values without advancing the coverages (for the inner iterations of a
 * host that updates its boundary coefficients several times per time step). */
ALDKIN_API int aldkin_evaluate(aldkin_t* h, const double* T_face, double* uptake_flux_out, double* keff_out);

/* Coverages (see the ALDKIN_THETA_* selectors for the sizes). For flatOnly handles MID/BOTTOM/HOLE read as 0. */
ALDKIN_API int aldkin_theta(const aldkin_t* h, int which, double* out);

/* Overwrite coverages (restart): ALDKIN_THETA_FLAT (n_faces) or ALDKIN_THETA_HOLE (n_faces * nCells). */
ALDKIN_API int aldkin_set_theta(aldkin_t* h, int which, const double* in);

/* K_eff of the last step/evaluate (n_faces). */
ALDKIN_API int aldkin_keff(const aldkin_t* h, double* out);

/* Mixed-boundary-condition coefficients for the host's transported scalar T. D_over_delta[i] = D / delta_i
 * (diffusivity over the wall-normal cell distance; OpenFOAM: D * deltaCoeffs), n_faces values.
 *   refValue      J_des / max(K, 1e-15) in T units (flatOnly with kdes > 0), otherwise 0
 *   refGrad       0
 *   valueFraction K_eff / (K_eff + D/delta)
 * Any of the three output pointers may be NULL. */
ALDKIN_API int aldkin_mixed_bc(const aldkin_t* h, const double* D_over_delta,
                               double* refValue, double* refGrad, double* valueFraction);

/* Desorption flux J_des = kdes * theta * N0 per face [molecules/(m2 s)] (0 unless flatOnly with kdes > 0). */
ALDKIN_API int aldkin_desorption_flux(const aldkin_t* h, double* out);

/* Relative gas density profile n/n_mouth in the hole of one face from the last step/evaluate (nCells values;
 * zeros for a face without holes). */
ALDKIN_API int aldkin_hole_gas(const aldkin_t* h, int face, double* nrel_out);

#ifdef __cplusplus
}
#endif

#endif /* ALDKIN_H */
