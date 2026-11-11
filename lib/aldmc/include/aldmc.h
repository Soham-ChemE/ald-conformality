/*
 * aldmc: 3D feature-scale ALD solver on a voxel grid, test-particle Monte Carlo in the molecular-flow regime.
 *
 * Physics (molecular flow, Knudsen regime). Particles move in straight lines between wall collisions and never
 * collide with each other. They enter through the feature mouth (the open plane z = 0) with a cosine (Lambertian)
 * angular distribution and a uniform position over the open mouth voxels. On hitting a solid voxel face they stick
 * with probability s (1 - theta_face) (irreversible Langmuir adsorption on the face's own site coverage), are lost
 * with probability r (radical recombination, default 0), and otherwise re-emit diffusely (cosine law about the face
 * normal). Particles leaving through the mouth are lost to the gas phase ("escaped"). Coverage lives on the exposed
 * voxel faces: a face carries N0 * A_face sites, A_face = h^2 * w_face (see area correction below).
 *
 * Time stepping. The mouth flux is n(t) vbar / 4 per unit mouth area. One call of aldmc_step launches M particles
 * representing n vbar/4 * A_mouth * dt molecules, each carrying w = n vbar/4 A_mouth dt / M molecules. The coverage
 * field is frozen while the M particles fly (sticking is sampled with the coverage at the start of the step) and the
 * wall collisions of the step are counted per face. The default update (ALDMC_UPDATE_EXPONENTIAL) then applies the
 * exact Langmuir relaxation of the impinged dose,
 *     theta_f <- 1 - (1 - theta_f) exp( - s * hits_f * w / (N0 h^2) ),
 * which is the 1D solver's update with the Monte Carlo impingement in place of n vbar/4 dt; it uses every wall
 * collision, so its variance is about 1/s times smaller than counting sticking events. ALDMC_UPDATE_DEPOSIT is the
 * literal alternative: each sampled sticking event deposits w molecules on its face and theta_f is capped at 1.
 * Both are first order in dt; they agree as M and 1/dt grow. Film growth does not narrow the feature.
 *
 * Area correction (areaCorrection = 1, default). A voxel surface is a staircase: for a cylinder the voxel faces have
 * 4/pi times the area of the smooth wall. A smooth surface normal n is estimated for each exposed face from the
 * first moment of the occupancy in a neighbourhood of the face centre, and the face is given the weight
 * w_face = |n . e_face| (its projected area fraction): it carries N0 h^2 w_face sites, and sticking and loss per
 * collision are s w_face (1 - theta) and r w_face. In an equilibrium gas every face is hit at n vbar/4 h^2 whatever
 * its orientation, so the uptake per unit smooth area and the rate per site are the smooth surface's. Re-emission is
 * diffuse about the same smooth normal (emission = ALDMC_EMISSION_SMOOTH, directions pointing into the facet are
 * resampled), which restores the smooth wall's Knudsen transport: with re-emission about the facet normals
 * (ALDMC_EMISSION_FACET) the staircase's grazing micro-flights lower the axial transport of a 16-voxel cylinder by
 * several per cent, and Knudsen's formula for a rough tube (4A/P with the staircase perimeter) would give -21 %.
 * With areaCorrection = 0 every face has w = 1 and emission is about the facet normal (plain voxel faces).
 *
 * Threads. The M particles of a step are processed in fixed blocks; block b of step k uses its own random stream
 * seeded from (seed, k, b), and each thread accumulates integer collision counts per face which are merged by
 * integer addition. The result is therefore bitwise identical for any thread count and reproducible for a given
 * seed; threads = 0 uses the hardware concurrency.
 *
 * Grid. nx * ny * nz voxels of size h (m); occupancy[(k ny + j) nx + i] = 1 solid, 0 gas; z = k h increases into the
 * feature, the mouth is the plane z = 0 (gas voxels with k = 0 are mouth voxels). Particles crossing z = 0 outward
 * are counted "escaped"; in a non-periodic direction a particle crossing the domain boundary, or crossing z = nz h,
 * is counted "left" (geometries are expected to be closed by solid voxels there; periodicX/Y wrap instead).
 * Face directions: 0 +x, 1 -x, 2 +y, 3 -y, 4 +z, 5 -z is the direction from the gas voxel to the solid neighbour
 * (a floor face has direction 4).
 *
 * Errors are return codes (ALDMC_OK = 0); aldmc_last_error(h) gives the message of the last failure on that handle
 * (aldmc_last_error(NULL) the message of the last failed aldmc_create in this thread). No exception crosses the ABI.
 */
#ifndef ALDMC_H
#define ALDMC_H

#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(ALDMC_BUILD_DLL)
#    define ALDMC_API __declspec(dllexport)
#  elif defined(ALDMC_DLL)
#    define ALDMC_API __declspec(dllimport)
#  else
#    define ALDMC_API
#  endif
#else
#  if defined(__GNUC__) || defined(__clang__)
#    define ALDMC_API __attribute__((visibility("default")))
#  else
#    define ALDMC_API
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ALDMC_VERSION_MAJOR 1
#define ALDMC_VERSION_MINOR 0
#define ALDMC_VERSION_PATCH 0

enum aldmc_error {
    ALDMC_OK = 0,
    ALDMC_ERR_ARG = 1,      /* null pointer, bad size or index */
    ALDMC_ERR_PARAM = 2,    /* invalid physical parameters or grid */
    ALDMC_ERR_STATE = 3,    /* call out of sequence */
    ALDMC_ERR_ALLOC = 4,    /* memory allocation failed */
    ALDMC_ERR_INTERNAL = 5
};

enum aldmc_update {
    ALDMC_UPDATE_EXPONENTIAL = 0,   /* theta <- 1 - (1 - theta) exp(-s hits w / (N0 h^2)) from the impingement */
    ALDMC_UPDATE_DEPOSIT = 1        /* each sampled sticking event deposits w molecules, theta capped at 1 */
};

enum aldmc_emission {
    ALDMC_EMISSION_SMOOTH = 0,      /* cosine about the estimated smooth normal, directions into the facet resampled */
    ALDMC_EMISSION_FACET = 1        /* cosine about the voxel face normal (plain staircase) */
};

enum aldmc_direction { ALDMC_DIR_PX = 0, ALDMC_DIR_MX = 1, ALDMC_DIR_PY = 2, ALDMC_DIR_MY = 3, ALDMC_DIR_PZ = 4, ALDMC_DIR_MZ = 5 };

typedef struct aldmc_t aldmc_t;

typedef struct aldmc_params {
    double s;                 /* sticking probability on a bare site */
    double r;                 /* loss (recombination) probability per collision */
    double N0;                /* site density, 1/m2 */
    double vbar;              /* mean thermal speed, m/s */
    unsigned long long seed;  /* random seed */
    int threads;              /* worker threads; 0 = hardware concurrency */
    int periodicX;            /* 1: wrap in x */
    int periodicY;            /* 1: wrap in y */
    int areaCorrection;       /* 1: projected-area face weights from the estimated smooth normal; 0: w = 1 */
    int update;               /* aldmc_update */
    int normalRadius;         /* neighbourhood radius (voxels) of the normal estimate; 0 = default 5 */
    int emission;             /* aldmc_emission: re-emission about the smooth normal (default) or the facet normal */
    long long maxSegments;    /* per-particle cap on traversed voxel segments; 0 = default 1e8 (then "left") */
} aldmc_params;

ALDMC_API const char* aldmc_version(void);
ALDMC_API void aldmc_params_default(aldmc_params* p);

/* Create a solver on an nx*ny*nz grid of voxel size h (m). occupancy: nx*ny*nz bytes, 1 solid, 0 gas, index
   (k*ny + j)*nx + i. The exposed faces, the mouth voxels and the face weights are computed here. */
ALDMC_API aldmc_t* aldmc_create(int nx, int ny, int nz, double h, const unsigned char* occupancy, const aldmc_params* p);
ALDMC_API void aldmc_destroy(aldmc_t* h);
ALDMC_API const char* aldmc_last_error(const aldmc_t* h);
ALDMC_API int aldmc_get_params(const aldmc_t* h, aldmc_params* p);

/* Launch M particles for a mouth density n_mouth (1/m3) over dt (s) and advance the coverage. */
ALDMC_API int aldmc_step(aldmc_t* h, double n_mouth, double dt, long long M);

/* Coverage to 0, counters to 0, time to 0, step index to 0 (the random streams restart). */
ALDMC_API int aldmc_reset(aldmc_t* h);

ALDMC_API int aldmc_n_faces(const aldmc_t* h);
ALDMC_API int aldmc_n_mouth_voxels(const aldmc_t* h);
ALDMC_API int aldmc_mouth_area(const aldmc_t* h, double* out);      /* open mouth area, m2 */
ALDMC_API int aldmc_site_area(const aldmc_t* h, double* out);       /* sum of weighted face areas, m2 (sites / N0) */
ALDMC_API int aldmc_time(const aldmc_t* h, double* out);            /* accumulated dt, s */

/* Per-face data, arrays of length n_faces (any pointer may be NULL): the gas voxel's linear index, the direction
   (aldmc_direction), the weight w_face, and the face centre in metres. */
ALDMC_API int aldmc_face_info(const aldmc_t* h, long long* voxel, int* direction, double* weight,
                              double* cx, double* cy, double* cz);
ALDMC_API int aldmc_face_coverage(const aldmc_t* h, double* out);
ALDMC_API int aldmc_set_face_coverage(aldmc_t* h, const double* in);
/* Collisions and sticking events on each face in the last step (unsigned long long, length n_faces). */
ALDMC_API int aldmc_face_hits(const aldmc_t* h, unsigned long long* hits, unsigned long long* sticks);

/* Weighted mean coverage of the faces whose gas voxel lies in slice k (length nz; weight_out may be NULL and
   receives the weighted face area per slice in m2; slices without faces give coverage 0 and area 0). */
ALDMC_API int aldmc_slice_coverage(const aldmc_t* h, double* out, double* weight_out);
/* Weighted mean coverage of the floor faces (direction +z) in the deepest slice that has any. */
ALDMC_API int aldmc_bottom_coverage(const aldmc_t* h, double* out);

/* Cumulative particle counts since the last reset: escaped through the mouth, stuck, lost by recombination, left
   through another domain boundary (or hit the segment cap). Any pointer may be NULL. */
ALDMC_API int aldmc_counts(const aldmc_t* h, unsigned long long* escaped, unsigned long long* stuck,
                           unsigned long long* lost, unsigned long long* left);
/* The same for the last step only. */
ALDMC_API int aldmc_last_counts(const aldmc_t* h, unsigned long long* escaped, unsigned long long* stuck,
                                unsigned long long* lost, unsigned long long* left);
/* Voxel segments traversed in the last step (a measure of work). */
ALDMC_API int aldmc_last_segments(const aldmc_t* h, unsigned long long* out);

#ifdef __cplusplus
}
#endif

#endif /* ALDMC_H */
