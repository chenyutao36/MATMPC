#include "matmpc_la.h"

static double symmetric_lower(const double *a, size_t ld, size_t row, size_t col)
{
    return row >= col ? a[row + ld * col] : a[col + ld * row];
}

static void update(double *out, double beta, double value)
{
    *out = beta == 0.0 ? value : value + beta * (*out);
}

void matmpc_la_gemv(int t, size_t m, size_t n, double alpha, const double *a,
                    size_t lda, const double *x, double beta, double *y)
{
    const size_t rows = t ? n : m;
    const size_t inner = t ? m : n;
    for (size_t i = 0; i < rows; ++i) {
        double sum = 0.0;
        if (alpha != 0.0)
            for (size_t j = 0; j < inner; ++j)
                sum += (t ? a[j + lda * i] : a[i + lda * j]) * x[j];
        update(y + i, beta, inner == 0 ? 0.0 : alpha * sum);
    }
}

void matmpc_la_symv_lower(size_t n, double alpha, const double *a, size_t lda,
                         const double *x, double beta, double *y)
{
    for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        if (alpha != 0.0)
            for (size_t j = 0; j < n; ++j)
                sum += symmetric_lower(a, lda, i, j) * x[j];
        update(y + i, beta, alpha * sum);
    }
}

void matmpc_la_gemm(int ta, int tb, size_t m, size_t n, size_t k, double alpha,
                    const double *a, size_t lda, const double *b, size_t ldb,
                    double beta, double *c, size_t ldc)
{
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < m; ++i) {
            double sum = 0.0;
            if (alpha != 0.0)
                for (size_t p = 0; p < k; ++p) {
                    const double av = ta ? a[p + lda * i] : a[i + lda * p];
                    const double bv = tb ? b[j + ldb * p] : b[p + ldb * j];
                    sum += av * bv;
                }
            update(c + i + ldc * j, beta, k == 0 ? 0.0 : alpha * sum);
        }
    }
}

void matmpc_la_symm_left_lower(size_t m, size_t n, double alpha, const double *a,
                              size_t lda, const double *b, size_t ldb,
                              double beta, double *c, size_t ldc)
{
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < m; ++i) {
            double sum = 0.0;
            if (alpha != 0.0)
                for (size_t p = 0; p < m; ++p)
                    sum += symmetric_lower(a, lda, i, p) * b[p + ldb * j];
            update(c + i + ldc * j, beta, alpha * sum);
        }
    }
}
