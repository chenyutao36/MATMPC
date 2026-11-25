#ifndef MATMPC_RECOVER_H
#define MATMPC_RECOVER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MATMPC_RECOVER_OK = 0,
    MATMPC_RECOVER_INVALID_DIMENSIONS,
    MATMPC_RECOVER_INVALID_BUFFER
} matmpc_recover_status;

typedef struct {
    size_t nx;
    size_t nu;
    size_t horizon;
    size_t state_bounds;
    size_t path_constraints;
    size_t terminal_constraints;
} matmpc_recover_dimensions;

typedef struct {
    const double *data;
    size_t count;
} matmpc_recover_buffer;

/* Each matrix block is stored column-major. Stage blocks are contiguous. */
typedef struct {
    matmpc_recover_buffer q;       /* (N + 1) state Hessians, nx by nx */
    matmpc_recover_buffer s;       /* N mixed Hessians, nx by nu */
    matmpc_recover_buffer a;       /* N state transition matrices, nx by nx */
    matmpc_recover_buffer b;       /* N input transition matrices, nx by nu */
    matmpc_recover_buffer cx;      /* state-bound Jacobian, nbx by nx */
    matmpc_recover_buffer cgx;     /* N path Jacobians, nc by nx */
    matmpc_recover_buffer cgn;     /* terminal Jacobian, ncN by nx */
    matmpc_recover_buffer gx;      /* state gradients, nx by N + 1 */
    matmpc_recover_buffer defect;  /* stage defects, nx by N */
    matmpc_recover_buffer ds0;     /* initial state deviation, nx */
    matmpc_recover_buffer du;      /* control increments, nu by N */
    matmpc_recover_buffer mu;      /* stage then terminal multipliers, N*nc + ncN */
    matmpc_recover_buffer mu_x;    /* state-bound multipliers at nodes 1..N, N*nbx */
} matmpc_recover_inputs;

/* Required number of double elements. dx also gives the length of lambda and gx. */
typedef struct {
    size_t dx;
    size_t du;
    size_t q;
    size_t s;
    size_t a;
    size_t cx;
    size_t cgx;
    size_t cgn;
    size_t defect;
    size_t mu;
    size_t mu_x;
} matmpc_recover_sizes;

/* The horizon and state dimension must be positive. Other dimensions may be zero. */
matmpc_recover_status matmpc_recover_get_sizes(const matmpc_recover_dimensions *dimensions,
                                             matmpc_recover_sizes *sizes);

/*
 * Recover state increments and costates in caller-owned nx by N + 1 outputs.
 * Outputs must not overlap one another or any input buffers. A zero-length
 * input may be null. The function uses no dynamic storage or retained state.
 */
matmpc_recover_status matmpc_recover(const matmpc_recover_dimensions *dimensions,
                                   const matmpc_recover_inputs *inputs,
                                   double *dx, size_t dx_count,
                                   double *lambda, size_t lambda_count);

#ifdef __cplusplus
}
#endif

#endif
