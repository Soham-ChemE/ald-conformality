/*
 * plain_c_example.c: one patterned wafer face driven through a recorded gas-density history, from plain C.
 *
 *   plain_c_example [centre_face_drive.csv]
 *
 * The CSV (written by export_centre_face_drive.py from reactor/showerhead_pulse_patterned_AR100_long) holds the
 * innermost wafer face's value of T = c/c0 every 50 ms over a 20 s TMA pulse and 2 s purge, together with the
 * boundary condition's own hole-bottom and flat coverages at those times. This program creates one face with the
 * case's parameters, interpolates the drive linearly to the solver's 1 ms step, calls aldkin_step once per step
 * (the forward-Euler scheme of the boundary condition) and prints the coverages, K_eff and the uptake flux every
 * second, then the time to 90 % hole-bottom coverage against the recorded one. The recorded run saw the true face
 * value at every 1 ms step, this one a 50 ms linear interpolation of it, so the two agree to about a write interval
 * (the program exits 1 if they differ by more than 0.1 s).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aldkin.h"

#define MAX_ROWS 100000

static int read_csv(const char* path, double* t, double* T, double* thetaB, int* n)
{
    FILE* f = fopen(path, "r");
    char line[512];
    *n = 0;
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 1; }
    while (fgets(line, sizeof line, f)) {
        double a, b, c, d;
        if (line[0] == '#' || line[0] == 't') continue;                 /* parameter line and column names */
        if (sscanf(line, "%lf,%lf,%lf,%lf", &a, &b, &c, &d) != 4) continue;
        if (*n >= MAX_ROWS) break;
        t[*n] = a; T[*n] = b; thetaB[*n] = c; ++*n;
    }
    fclose(f);
    return *n < 2;
}

static double interp(const double* x, const double* y, int n, double xq)
{
    int i;
    if (xq <= x[0]) return y[0];
    if (xq >= x[n - 1]) return y[n - 1];
    for (i = 1; i < n; ++i) {
        if (xq <= x[i]) return y[i - 1] + (y[i] - y[i - 1]) * (xq - x[i - 1]) / (x[i] - x[i - 1]);
    }
    return y[n - 1];
}

int main(int argc, char** argv)
{
    const char* path = argc > 1 ? argv[1] : "centre_face_drive.csv";
    static double t[MAX_ROWS], Tdrive[MAX_ROWS], thetaB_ref[MAX_ROWS];
    int n, k, nsteps, rc;
    aldkin_params p;
    aldkin_t* h;
    const double dt = 1e-3;                                             /* the OpenFOAM case's deltaT */
    const double radius = 0.00067;                                      /* face-centre radius of the innermost face, m */
    double Tf, flux, keff, flat, mid, bottom;
    double t90_lib = -1.0, t90_ref = -1.0;
    double next_print = 0.0;

    if (read_csv(path, t, Tdrive, thetaB_ref, &n)) return 2;
    for (k = 0; k < n; ++k) if (thetaB_ref[k] >= 0.9) { t90_ref = t[k]; break; }

    /* the case's kinetic parameters (the CSV's first line lists them) */
    aldkin_params_default(&p);
    p.s = 5e-3; p.r = 0.0; p.vbar = 410.0; p.N0 = 4e18; p.c0 = 3.79e20;
    p.d = 100e-9; p.AR = 100.0; p.holeDensity = 1.78e13; p.nCells = 200;
    p.flatOnly = 0; p.scheme = ALDKIN_SCHEME_EULER;

    h = aldkin_create(&p);
    if (!h) { fprintf(stderr, "aldkin_create: %s\n", aldkin_last_error(NULL)); return 3; }
    rc = aldkin_set_faces(h, 1, &radius);
    if (rc != ALDKIN_OK) { fprintf(stderr, "aldkin_set_faces: %s\n", aldkin_last_error(h)); return 3; }

    printf("aldkin %s: one AR-%g face, %d hole cells, drive %s (%d samples, %g..%g s)\n",
           aldkin_version(), p.AR, p.nCells, path, n, t[0], t[n - 1]);
    printf("%8s %10s %12s %12s %12s %12s %14s\n", "t (s)", "T_face", "thetaFlat", "thetaMid", "thetaBottom", "K_eff m/s", "flux 1/(m2 s)");

    nsteps = (int)floor(t[n - 1] / dt + 0.5);
    for (k = 1; k <= nsteps; ++k) {
        const double tnow = k * dt;
        /* the boundary condition advances with the face value of the previous step */
        Tf = interp(t, Tdrive, n, tnow - dt);
        rc = aldkin_step(h, &Tf, dt, &flux, &keff);
        if (rc != ALDKIN_OK) { fprintf(stderr, "aldkin_step: %s\n", aldkin_last_error(h)); return 4; }
        aldkin_theta(h, ALDKIN_THETA_BOTTOM, &bottom);
        if (t90_lib < 0.0 && bottom >= 0.9) t90_lib = tnow;
        if (tnow + 0.5 * dt >= next_print) {
            aldkin_theta(h, ALDKIN_THETA_FLAT, &flat);
            aldkin_theta(h, ALDKIN_THETA_MID, &mid);
            printf("%8.3f %10.5f %12.6f %12.6f %12.8f %12.5g %14.5g\n", tnow, Tf, flat, mid, bottom, keff, flux);
            next_print += 1.0;
        }
    }
    aldkin_destroy(h);

    printf("hole-bottom t90 = %.3f s from this run; %.3f s recorded by the boundary condition (50 ms writes)\n",
           t90_lib, t90_ref);
    if (t90_lib < 0.0 || t90_ref < 0.0 || fabs(t90_lib - t90_ref) > 0.1) {
        printf("DIFFERENCE LARGER THAN 0.1 s\n");
        return 1;
    }
    return 0;
}
