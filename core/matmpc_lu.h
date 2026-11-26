#ifndef MATMPC_LU_H
#define MATMPC_LU_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MATMPC_LU_OK = 0,
    MATMPC_LU_INVALID_ARGUMENT,
    MATMPC_LU_SINGULAR,
    MATMPC_LU_NONFINITE
} matmpc_lu_status;

/*
 * Partially pivot a square, column-major matrix in place. The strictly lower
 * triangle receives L with its unit diagonal omitted and the rest receives U.
 * pivots[k] is the zero-based row swapped with row k. A pivot of exactly zero
 * is singular. The input matrix and pivot storage may be partly overwritten
 * if the elimination encounters a singular pivot or overflows.
 */
matmpc_lu_status matmpc_lu_factor(size_t n, double *a, size_t lda, size_t *pivots);

/*
 * Solve A X = B (transposed=0) or A^T X = B (transposed=1) using the factors.
 * The solution overwrites the n-by-nrhs column-major RHS. A nonfinite result
 * leaves the RHS unspecified. These functions allocate no memory. The caller
 * supplies nonoverlapping buffers with at least n pivots and sufficient space
 * for all columns at the given leading dimensions. Empty work is a no-op.
 */
matmpc_lu_status matmpc_lu_solve(int transposed, size_t n, size_t nrhs,
                                  const double *lu, size_t lda, const size_t *pivots,
                                  double *rhs, size_t ldb);

#ifdef __cplusplus
}
#endif

#endif
