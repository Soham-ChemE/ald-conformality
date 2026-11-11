/*
 * udf_stub.h: the handful of Fluent UDF macros and types that ald_wall_udf.c uses, as a stand-in for Fluent's udf.h
 * so that the adapter can be compiled and linked here as a syntax and link check:
 *
 *   cc -I. -I../../lib/aldkin/include -c ald_wall_udf.c
 *   cc -shared -fPIC -I. -I../../lib/aldkin/include ald_wall_udf.c -L../../lib/aldkin/build -laldkin -o libudf_check.so
 *
 * Nothing in here does anything: the face loops run over an empty thread, the field accessors return 0, Message
 * prints to stdout. The real build inside Fluent uses udf.h (ald_wall_udf.c picks it when it is on the include path).
 * The shapes of the DEFINE_* macros follow the Fluent Customization Manual (DEFINE_INIT, DEFINE_ADJUST,
 * DEFINE_SR_RATE, DEFINE_PROFILE, DEFINE_EXECUTE_AT_END, DEFINE_ON_DEMAND).
 */
#ifndef ALD_UDF_STUB_H
#define ALD_UDF_STUB_H

#include <stdio.h>

typedef double real;
typedef int face_t;
typedef int cell_t;
typedef struct Thread_s { int id; int n; } Thread;
typedef struct Domain_s { int id; } Domain;
typedef struct Reaction_s { int id; } Reaction;

#define ND_ND 3
#define RP_HOST 0
#define RP_NODE 0
#define PARALLEL 0

#define DEFINE_INIT(name, d) void name(Domain* d)
#define DEFINE_ADJUST(name, d) void name(Domain* d)
#define DEFINE_SR_RATE(name, f, t, r, mw, yi, rr) void name(face_t f, Thread* t, Reaction* r, real* mw, real* yi, real* rr)
#define DEFINE_PROFILE(name, t, i) void name(Thread* t, int i)
#define DEFINE_EXECUTE_AT_END(name) void name(void)
#define DEFINE_ON_DEMAND(name) void name(void)

static Thread ald_stub_thread = {0, 0};
static real ald_stub_zero = 0.0;
static real ald_stub_centroid[ND_ND] = {0.0, 0.0, 0.0};

static Thread* Lookup_Thread(Domain* d, int id) { (void)d; ald_stub_thread.id = id; return &ald_stub_thread; }
#define THREAD_N_ELEMENTS_INT(t) ((t)->n)
#define THREAD_ID(t) ((t)->id)

#define begin_f_loop(f, t) for ((f) = 0; (f) < THREAD_N_ELEMENTS_INT(t); ++(f)) {
#define end_f_loop(f, t) }

#define F_CENTROID(x, f, t) do { (x)[0] = ald_stub_centroid[0]; (x)[1] = ald_stub_centroid[1]; (x)[2] = ald_stub_centroid[2]; (void)(f); (void)(t); } while (0)
#define F_R(f, t) (ald_stub_zero)
#define F_YI(f, t, i) (ald_stub_zero)
#define F_PROFILE(f, t, i) (ald_stub_zero)
#define F_UDMI(f, t, i) (ald_stub_zero)
#define F_T(f, t) (ald_stub_zero)

#define CURRENT_TIME (0.0)
#define CURRENT_TIMESTEP (0.0)
#define N_TIME (0)
#define N_ITER (0)

#define Message printf
#define Error(...) do { printf(__VA_ARGS__); } while (0)

#endif /* ALD_UDF_STUB_H */
