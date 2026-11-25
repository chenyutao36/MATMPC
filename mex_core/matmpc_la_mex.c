#include "matmpc_la.h"
#include "blas.h"

static void scale(size_t m, size_t n, double beta, double *out, size_t ld)
{
    if (beta == 1.0)
        return;
    for (size_t j = 0; j < n; ++j)
        for (size_t i = 0; i < m; ++i)
            out[i + ld * j] = beta == 0.0 ? 0.0 : beta * out[i + ld * j];
}

void matmpc_la_gemv(int t, size_t m, size_t n, double alpha, const double *a,
                    size_t lda, const double *x, double beta, double *y)
{
    if (!m || !n || alpha == 0.0) {
        const size_t rows = t ? n : m;
        scale(rows, 1, beta, y, rows);
        return;
    }
    char trans = t ? 'T' : 'N';
    ptrdiff_t rows = (ptrdiff_t)m, cols = (ptrdiff_t)n, leading = (ptrdiff_t)lda, stride = 1;
    dgemv(&trans, &rows, &cols, &alpha, (double *)a, &leading,
          (double *)x, &stride, &beta, y, &stride);
}

void matmpc_la_symv_lower(size_t n, double alpha, const double *a, size_t lda,
                         const double *x, double beta, double *y)
{
    if (!n || alpha == 0.0) {
        scale(n, 1, beta, y, n);
        return;
    }
    char lower = 'L';
    ptrdiff_t count = (ptrdiff_t)n, leading = (ptrdiff_t)lda, stride = 1;
    dsymv(&lower, &count, &alpha, (double *)a, &leading,
          (double *)x, &stride, &beta, y, &stride);
}

void matmpc_la_gemm(int ta, int tb, size_t m, size_t n, size_t k, double alpha,
                    const double *a, size_t lda, const double *b, size_t ldb,
                    double beta, double *c, size_t ldc)
{
    if (!m || !n || !k || alpha == 0.0) {
        scale(m, n, beta, c, ldc);
        return;
    }
    char trans_a = ta ? 'T' : 'N', trans_b = tb ? 'T' : 'N';
    ptrdiff_t rows = (ptrdiff_t)m, cols = (ptrdiff_t)n, inner = (ptrdiff_t)k;
    ptrdiff_t lead_a = (ptrdiff_t)lda, lead_b = (ptrdiff_t)ldb, lead_c = (ptrdiff_t)ldc;
    dgemm(&trans_a, &trans_b, &rows, &cols, &inner, &alpha,
          (double *)a, &lead_a, (double *)b, &lead_b, &beta, c, &lead_c);
}

void matmpc_la_symm_left_lower(size_t m, size_t n, double alpha, const double *a,
                              size_t lda, const double *b, size_t ldb,
                              double beta, double *c, size_t ldc)
{
    if (!m || !n || alpha == 0.0) {
        scale(m, n, beta, c, ldc);
        return;
    }
    char left = 'L', lower = 'L';
    ptrdiff_t rows = (ptrdiff_t)m, cols = (ptrdiff_t)n;
    ptrdiff_t lead_a = (ptrdiff_t)lda, lead_b = (ptrdiff_t)ldb, lead_c = (ptrdiff_t)ldc;
    dsymm(&left, &lower, &rows, &cols, &alpha,
          (double *)a, &lead_a, (double *)b, &lead_b, &beta, c, &lead_c);
}
