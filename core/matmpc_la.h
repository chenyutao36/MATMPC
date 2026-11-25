#ifndef MATMPC_LA_H
#define MATMPC_LA_H

#include <stddef.h>

/* Internal double-precision, column-major operations with unit vector strides.
 * A nonzero transpose flag selects the transposed operand. Symmetric inputs
 * use the stored lower triangle. Callers validate dimensions and buffers. */
void matmpc_la_gemv(int t, size_t m, size_t n, double alpha, const double *a,
                    size_t lda, const double *x, double beta, double *y);
void matmpc_la_symv_lower(size_t n, double alpha, const double *a, size_t lda,
                         const double *x, double beta, double *y);
void matmpc_la_gemm(int ta, int tb, size_t m, size_t n, size_t k, double alpha,
                    const double *a, size_t lda, const double *b, size_t ldb,
                    double beta, double *c, size_t ldc);
void matmpc_la_symm_left_lower(size_t m, size_t n, double alpha, const double *a,
                              size_t lda, const double *b, size_t ldb,
                              double beta, double *c, size_t ldc);

#endif
