#ifndef MATMPC_CONDENSE_AFFINE_H
#define MATMPC_CONDENSE_AFFINE_H

#include "condense.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Column-major matrices are stored in consecutive shooting-stage blocks. */
typedef struct {
    matmpc_condense_buffer q;              /* N+1 state Hessians, lower triangle */
    matmpc_condense_buffer s;              /* N state/control cross Hessians */
    matmpc_condense_buffer a;              /* N state transition matrices */
    matmpc_condense_buffer b;              /* N input transition matrices */
    matmpc_condense_buffer cx;             /* State-bound Jacobian, nbx by nx */
    matmpc_condense_buffer cgx;            /* N path-constraint state Jacobians */
    matmpc_condense_buffer cgn;            /* Terminal state Jacobian */
    matmpc_condense_buffer gx;             /* State gradients at nodes 0..N */
    matmpc_condense_buffer gu;             /* Control gradients per stage or per block */
    matmpc_condense_buffer defect;         /* N state transition defects */
    matmpc_condense_buffer ds0;            /* Initial state deviation */
    matmpc_condense_buffer lower_general;  /* N*nc path rows followed by ncN terminal rows */
    matmpc_condense_buffer upper_general;
    matmpc_condense_buffer lower_state;    /* N*nbx rows at state nodes 1..N */
    matmpc_condense_buffer upper_state;
} matmpc_condense_affine_inputs;

typedef struct {
    matmpc_condense_output gc;
    matmpc_condense_output lcc;
    matmpc_condense_output ucc;
    matmpc_condense_output lxc;
    matmpc_condense_output uxc;
} matmpc_condense_affine_outputs;

/* Required double counts. s covers B as well as S. */
typedef struct {
    size_t states;
    size_t controls;
    size_t q;
    size_t s;
    size_t a;
    size_t cx;
    size_t cgx;
    size_t cgn;
    size_t defect;
    size_t general;
    size_t state;
    size_t workspace;
} matmpc_condense_affine_sizes;

matmpc_condense_status matmpc_condense_affine_get_sizes(const matmpc_condense_dimensions *dimensions,
                                                       matmpc_condense_affine_sizes *sizes);

/*
 * Condense the gradient and bounds using caller-owned outputs and workspace.
 * The state dimension and horizon must be positive. Empty buffers may be null.
 * Outputs and workspace must not overlap any input or one another.
 */
matmpc_condense_status matmpc_condense_affine(const matmpc_condense_dimensions *dimensions,
                                             const matmpc_condense_affine_inputs *inputs,
                                             const matmpc_condense_affine_outputs *outputs,
                                             matmpc_condense_output workspace);

/* Use the same shooting grid with control increments shared by each block.
 * gu contains the accumulated direct gradients, one nu-vector per block.
 * The boundary array starts at zero and ends at the shooting horizon.
 */
matmpc_condense_status matmpc_condense_affine_blocked(const matmpc_condense_dimensions *dimensions,
                                                     const matmpc_condense_affine_inputs *inputs,
                                                     const matmpc_condense_affine_outputs *outputs,
                                                     matmpc_condense_output workspace,
                                                     const size_t *boundaries, size_t boundary_count);

#ifdef __cplusplus
}
#endif

#endif
