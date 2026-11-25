#ifndef MATMPC_CONDENSE_MATRIX_H
#define MATMPC_CONDENSE_MATRIX_H

#include "condense.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Each stage is a contiguous column-major block. Q uses its lower triangle. */
typedef struct {
    matmpc_condense_buffer q;   /* (N + 1) nx-by-nx blocks */
    matmpc_condense_buffer s;   /* N nx-by-nu blocks */
    matmpc_condense_buffer r;   /* N nu-by-nu blocks */
    matmpc_condense_buffer a;   /* N nx-by-nx blocks */
    matmpc_condense_buffer b;   /* N nx-by-nu blocks */
    matmpc_condense_buffer cx;  /* nbx-by-nx */
    matmpc_condense_buffer cgx; /* N nc-by-nx blocks */
    matmpc_condense_buffer cgu; /* N nc-by-nu blocks */
    matmpc_condense_buffer cgn; /* ncN-by-nx */
} matmpc_condense_matrix_inputs;

typedef struct {
    matmpc_condense_output hc;  /* N*nu-by-N*nu */
    matmpc_condense_output ccx; /* N*nbx-by-N*nu, states 1 through N */
    matmpc_condense_output ccg; /* (N*nc + ncN)-by-N*nu, stages then terminal */
} matmpc_condense_matrix_outputs;

/* Required lengths in doubles. Workspace holds G, W and one nu-by-nu block. */
typedef struct {
    size_t q, s, r, a, b, cx, cgx, cgu, cgn;
    size_t hc, ccx, ccg;
    size_t workspace;
} matmpc_condense_matrix_sizes;

matmpc_condense_status matmpc_condense_matrix_get_sizes(const matmpc_condense_dimensions *dimensions,
                                                      matmpc_condense_matrix_sizes *sizes);

/*
 * Compute condensed quadratic and constraint matrices with caller-owned storage.
 * Inputs, outputs, and scratch must not overlap. Empty buffers may have null
 * data. The routine performs no allocations and accepts previously used scratch.
 */
matmpc_condense_status matmpc_condense_matrix(const matmpc_condense_dimensions *dimensions,
                                            const matmpc_condense_matrix_inputs *inputs,
                                            const matmpc_condense_matrix_outputs *outputs,
                                            matmpc_condense_output workspace);

#ifdef __cplusplus
}
#endif

#endif
