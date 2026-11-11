/*
 * ald_wall_udf.c: ANSYS Fluent user-defined functions that put the aldkin multiscale ALD wall kinetics on a wall
 * zone of a transient species-transport case.
 *
 * STATUS: compiled and linked here against udf_stub.h and libaldkin as a syntax and link check only. It has NOT been
 * run inside Fluent; the hook shapes follow the Fluent Customization Manual and the usual species-UDF conventions
 * but every Fluent-specific line (zone id, species index, molar mass, parallel partitioning) must be checked by the
 * person who loads it. Edit the ALD_* defines below for the case.
 *
 * Hooks (all in this file):
 *   DEFINE_INIT          ald_wall_init        once at initialisation: create the handle, register the wall faces
 *                                             with their centre radii, install optional radial profiles
 *   DEFINE_ADJUST        ald_wall_adjust      every iteration: on the first iteration of a new time step call
 *                                             aldkin_step with the precursor number density on each face (from the
 *                                             previous step, as the OpenFOAM BC does) and the time step; on the
 *                                             further iterations aldkin_evaluate with the current face densities
 *   DEFINE_SR_RATE       ald_wall_sr_rate     wall surface-reaction rate for the precursor: rr = K_eff * C_face
 *                                             [kmol/(m2 s)], C_face = rho Y / MW; Fluent then applies the sink to
 *                                             the species equation as a Robin-type wall condition
 *   DEFINE_PROFILE       ald_wall_flux        alternative to the surface reaction: a specified species mass flux
 *                                             -K_eff rho Y [kg/(m2 s)] (explicit; use SR_RATE where available)
 *   DEFINE_EXECUTE_AT_END ald_wall_record     end of every time step: write thetaFlat, thetaMid, thetaBottom and
 *                                             K_eff into user-defined face memory 0..3 for contour plots
 *   DEFINE_ON_DEMAND     ald_wall_free        release the handle
 *
 * Units: the library is created with c0 = 1, so the "T_face" it receives is the absolute number density in 1/m3
 * (n = rho Y / MW * 1000 N_A with MW in kg/kmol) and K_eff comes back in m/s. The surface-reaction rate in
 * kmol/(m2 s) is K_eff * rho Y / MW. Set the kinetic parameters in SI as in the OpenFOAM cases.
 *
 * Linking: libaldkin must be linked into libudf. Either add lib/aldkin/build/libaldkin.a to the link line of the
 * generated libudf makefile (or Fluent's "User Defined > Functions > Compiled" with the library listed in the
 * makefile's LIBS), or put libaldkin.so/.dll beside libudf and add -laldkin. The header include/aldkin.h must be
 * on the include path (-I<repo>/lib/aldkin/include).
 *
 * Parallel: each compute node keeps its own handle for the faces it owns (begin_f_loop over the local partition);
 * nothing is exchanged between nodes because faces are independent. The host process does nothing (RP_HOST).
 */
#if defined(__has_include)
#  if __has_include("udf.h")
#    include "udf.h"
#    define ALD_IN_FLUENT 1
#  endif
#endif
#ifndef ALD_IN_FLUENT
#  include "udf_stub.h"
#endif

#include <math.h>
#include <stdlib.h>

#include "aldkin.h"

/* ---- case-specific settings: edit for the model ------------------------------------------------------------- */
#define ALD_WALL_ZONE_ID        7          /* zone id of the wafer wall (Boundary Conditions panel) */
#define ALD_PRECURSOR_SPECIES   0          /* index of the precursor in the species list (0 = first) */
#define ALD_PRECURSOR_MW        72.09      /* kg/kmol (TMA) */
#define ALD_S                   5e-3       /* sticking coefficient on an open site */
#define ALD_R                   0.0        /* loss probability per collision (recombination) */
#define ALD_VBAR                410.0      /* m/s mean thermal speed at the wall temperature */
#define ALD_N0                  4e18       /* sites per m2 */
#define ALD_D                   100e-9     /* hole diameter, m */
#define ALD_AR                  100.0      /* hole aspect ratio (0 = flat wafer) */
#define ALD_HOLE_DENSITY        1.78e13    /* holes per m2 (0 = flat wafer) */
#define ALD_NCELLS              200        /* cells along the hole */
#define ALD_FLAT_ONLY           0          /* 1: flat Langmuir only (then ALD_KDES may be > 0) */
#define ALD_KDES                0.0        /* 1/s first-order desorption (flat-only) */
#define ALD_AXIS                0          /* 0/1/2: the axis of the wafer (radius = distance from it) */
/* optional radial profiles (set ALD_N_PROFILE to 0 for a uniform pattern) */
#define ALD_N_PROFILE 0
#if ALD_N_PROFILE > 0
static const double ald_profile_r[ALD_N_PROFILE] = {0.0};
static const double ald_profile_AR[ALD_N_PROFILE] = {ALD_AR};
static const double ald_profile_density[ALD_N_PROFILE] = {ALD_HOLE_DENSITY};
#endif
/* ---------------------------------------------------------------------------------------------------------------- */

#define AVOGADRO_PER_KMOL 6.02214076e26

static aldkin_t* ald_h = NULL;
static int ald_nfaces = 0;
static int ald_last_timestep = -1;
static double* ald_T = NULL;          /* face number densities, 1/m3 */
static double* ald_keff = NULL;       /* m/s */

static void ald_fail(const char* what)
{
    Message("ald_wall_udf: %s: %s\n", what, aldkin_last_error(ald_h));
}

/* precursor number density on face f (1/m3) from the face mass fraction and density */
static double ald_face_density(face_t f, Thread* t)
{
    const double rho = F_R(f, t);
    const double Y = F_YI(f, t, ALD_PRECURSOR_SPECIES);
    double n = rho * Y / ALD_PRECURSOR_MW * AVOGADRO_PER_KMOL;
    (void)f; (void)t;                               /* the stub's accessors ignore them */
    return n > 0.0 ? n : 0.0;
}

DEFINE_INIT(ald_wall_init, d)
{
#if !RP_HOST
    Thread* t = Lookup_Thread(d, ALD_WALL_ZONE_ID);
    face_t f;
    aldkin_params p;
    double* radius;
    int i = 0;

    if (ald_h) { aldkin_destroy(ald_h); ald_h = NULL; }
    free(ald_T); free(ald_keff); ald_T = NULL; ald_keff = NULL;

    aldkin_params_default(&p);
    p.s = ALD_S; p.r = ALD_R; p.kdes = ALD_KDES; p.vbar = ALD_VBAR; p.N0 = ALD_N0;
    p.c0 = 1.0;                                     /* absolute number densities in and out */
    p.d = ALD_D; p.AR = ALD_AR; p.holeDensity = ALD_HOLE_DENSITY; p.nCells = ALD_NCELLS;
    p.flatOnly = ALD_FLAT_ONLY; p.scheme = ALDKIN_SCHEME_EULER;
    ald_h = aldkin_create(&p);
    if (!ald_h) { Message("ald_wall_udf: aldkin_create: %s\n", aldkin_last_error(NULL)); return; }

    ald_nfaces = THREAD_N_ELEMENTS_INT(t);
    radius = (double*)malloc((size_t)(ald_nfaces > 0 ? ald_nfaces : 1) * sizeof(double));
    ald_T = (double*)calloc((size_t)(ald_nfaces > 0 ? ald_nfaces : 1), sizeof(double));
    ald_keff = (double*)calloc((size_t)(ald_nfaces > 0 ? ald_nfaces : 1), sizeof(double));
    begin_f_loop(f, t)
    {
        real x[ND_ND];
        double r2 = 0.0;
        int k;
        F_CENTROID(x, f, t);
        for (k = 0; k < ND_ND; ++k) if (k != ALD_AXIS) r2 += (double)x[k] * (double)x[k];
        radius[i++] = sqrt(r2);
    }
    end_f_loop(f, t)
    if (aldkin_set_faces(ald_h, ald_nfaces, radius) != ALDKIN_OK) ald_fail("aldkin_set_faces");
    free(radius);
#if ALD_N_PROFILE > 0
    if (aldkin_set_profiles(ald_h, ALD_N_PROFILE, ald_profile_r, ald_profile_AR, ald_profile_density) != ALDKIN_OK)
        ald_fail("aldkin_set_profiles");
#endif
    ald_last_timestep = -1;
    Message("ald_wall_udf: aldkin %s, %d wall faces on zone %d, %d hole cells\n",
            aldkin_version(), ald_nfaces, ALD_WALL_ZONE_ID, ALD_NCELLS);
#endif
}

DEFINE_ADJUST(ald_wall_adjust, d)
{
#if !RP_HOST
    Thread* t;
    face_t f;
    int i = 0, rc;
    if (!ald_h) return;
    t = Lookup_Thread(d, ALD_WALL_ZONE_ID);
    if (THREAD_N_ELEMENTS_INT(t) != ald_nfaces) { Message("ald_wall_udf: face count changed; re-run ald_wall_init\n"); return; }
    begin_f_loop(f, t) { ald_T[i++] = ald_face_density(f, t); } end_f_loop(f, t)
    if (N_TIME != ald_last_timestep)
    {
        /* first iteration of a new time step: advance the coverages with the face densities of the previous step */
        rc = aldkin_step(ald_h, ald_T, (double)CURRENT_TIMESTEP, NULL, ald_keff);
        if (rc != ALDKIN_OK) { ald_fail("aldkin_step"); return; }
        ald_last_timestep = N_TIME;
    }
    else
    {
        rc = aldkin_evaluate(ald_h, ald_T, NULL, ald_keff);
        if (rc != ALDKIN_OK) { ald_fail("aldkin_evaluate"); return; }
    }
#endif
}

/* surface reaction rate of the precursor on the wall: rr = K_eff * C_face, kmol/(m2 s) */
DEFINE_SR_RATE(ald_wall_sr_rate, f, t, r, mw, yi, rr)
{
#if !RP_HOST
    double C;
    (void)r; (void)mw; (void)yi; (void)t;
    if (!ald_h || f < 0 || f >= ald_nfaces) { *rr = 0.0; return; }
    C = F_R(f, t) * F_YI(f, t, ALD_PRECURSOR_SPECIES) / ALD_PRECURSOR_MW;   /* kmol/m3 */
    *rr = C > 0.0 ? ald_keff[f] * C : 0.0;
#endif
}

/* alternative: specified species mass flux into the wall, kg/(m2 s), negative = leaving the gas */
DEFINE_PROFILE(ald_wall_flux, t, i)
{
#if !RP_HOST
    face_t f;
    (void)i;
    if (!ald_h) return;
    begin_f_loop(f, t)
    {
        const double rhoY = F_R(f, t) * F_YI(f, t, ALD_PRECURSOR_SPECIES);
        F_PROFILE(f, t, i) = (f >= 0 && f < ald_nfaces && rhoY > 0.0) ? -ald_keff[f] * rhoY : 0.0;
    }
    end_f_loop(f, t)
#endif
}

/* coverages and K_eff into user-defined face memory 0..3 (allocate 4 UDMs in the case) */
DEFINE_EXECUTE_AT_END(ald_wall_record)
{
#if !RP_HOST
    Domain* d = NULL;
    Thread* t;
    face_t f;
    double *flat, *mid, *bottom;
    int i = 0;
    if (!ald_h || ald_nfaces == 0) return;
#ifdef ALD_IN_FLUENT
    d = Get_Domain(1);
#endif
    t = Lookup_Thread(d, ALD_WALL_ZONE_ID);
    flat = (double*)malloc((size_t)ald_nfaces * sizeof(double));
    mid = (double*)malloc((size_t)ald_nfaces * sizeof(double));
    bottom = (double*)malloc((size_t)ald_nfaces * sizeof(double));
    aldkin_theta(ald_h, ALDKIN_THETA_FLAT, flat);
    aldkin_theta(ald_h, ALDKIN_THETA_MID, mid);
    aldkin_theta(ald_h, ALDKIN_THETA_BOTTOM, bottom);
    begin_f_loop(f, t)
    {
        if (i < ald_nfaces)
        {
            F_UDMI(f, t, 0) = (real)flat[i];
            F_UDMI(f, t, 1) = (real)mid[i];
            F_UDMI(f, t, 2) = (real)bottom[i];
            F_UDMI(f, t, 3) = (real)ald_keff[i];
        }
        ++i;
    }
    end_f_loop(f, t)
    free(flat); free(mid); free(bottom);
#endif
}

DEFINE_ON_DEMAND(ald_wall_free)
{
    if (ald_h) { aldkin_destroy(ald_h); ald_h = NULL; }
    free(ald_T); free(ald_keff); ald_T = NULL; ald_keff = NULL;
    ald_nfaces = 0;
    Message("ald_wall_udf: released\n");
}
