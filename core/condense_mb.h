#ifndef MATMPC_CONDENSE_MB_H
#define MATMPC_CONDENSE_MB_H

#include "condense_affine.h"
#include "condense_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    matmpc_condense_affine_inputs affine; /* gu contains one already aggregated vector per block */
    matmpc_condense_buffer r_block;       /* One accumulated and regularized nu-by-nu block per control block */
    matmpc_condense_buffer cgu;           /* N path-constraint control Jacobians */
} matmpc_condense_mb_inputs;

typedef struct {
    matmpc_condense_matrix_outputs matrices; /* r*nu control columns instead of N*nu */
    matmpc_condense_affine_outputs affine;   /* r*nu gradient components, uncompressed constraint rows */
} matmpc_condense_mb_outputs;

/* All counts are in doubles. Workspace is shared by the two sequential phases. */
typedef struct {
    matmpc_condense_affine_sizes affine;
    size_t controls, r_block, cgu;
    size_t hc, ccx, ccg, workspace;
} matmpc_condense_mb_sizes;

matmpc_condense_status matmpc_condense_mb_get_sizes(const matmpc_condense_dimensions *dimensions,
                                                   size_t block_count,
                                                   matmpc_condense_mb_sizes *sizes);

/*
 * Eliminate the state increments with one shared control increment per block.
 * Boundaries are strictly increasing shooting indices with endpoints 0 and N.
 * Callers provide nonoverlapping inputs, outputs and scratch. Storage is
 * column-major and each shooting stage is contiguous.
 */
matmpc_condense_status matmpc_condense_mb(const matmpc_condense_dimensions *dimensions,
                                         const size_t *boundaries, size_t boundary_count,
                                         const matmpc_condense_mb_inputs *inputs,
                                         const matmpc_condense_mb_outputs *outputs,
                                         matmpc_condense_output workspace);

#ifdef __cplusplus
}
#endif

#endif
