#include "matmpc_lu.h"

#include <math.h>
#include <stdint.h>

static int matrix_extent_valid(size_t rows, size_t columns, size_t ld)
{
    const size_t limit = SIZE_MAX / sizeof(double);

    return rows == 0 || columns == 0 ||
           (ld >= rows && rows <= limit && columns - 1 <= (limit - rows) / ld);
}

static int matrix_is_finite(const double *matrix, size_t rows, size_t columns, size_t ld)
{
    for (size_t col = 0; col < columns; ++col) {
        for (size_t row = 0; row < rows; ++row) {
            if (!isfinite(matrix[row + col * ld])) {
                return 0;
            }
        }
    }
    return 1;
}

static void swap(double *first, double *second)
{
    double value = *first;
    *first = *second;
    *second = value;
}

matmpc_lu_status matmpc_lu_factor(size_t n, double *a, size_t lda, size_t *pivots)
{
    if (n == 0) {
        return MATMPC_LU_OK;
    }
    if (a == NULL || pivots == NULL || n > SIZE_MAX / sizeof(*pivots) ||
        !matrix_extent_valid(n, n, lda)) {
        return MATMPC_LU_INVALID_ARGUMENT;
    }
    if (!matrix_is_finite(a, n, n, lda)) {
        return MATMPC_LU_NONFINITE;
    }

    for (size_t k = 0; k < n; ++k) {
        double *column = a + k * lda;
        size_t pivot = k;
        double largest = fabs(column[k]);

        for (size_t row = k + 1; row < n; ++row) {
            const double candidate = fabs(column[row]);
            if (candidate > largest) {
                largest = candidate;
                pivot = row;
            }
        }
        if (largest == 0.0) {
            return MATMPC_LU_SINGULAR;
        }

        pivots[k] = pivot;
        if (pivot != k) {
            for (size_t col = 0; col < n; ++col) {
                swap(a + k + col * lda, a + pivot + col * lda);
            }
        }
        for (size_t row = k + 1; row < n; ++row) {
            /* Direct division is safe even when the reciprocal would overflow. */
            column[row] /= column[k];
        }
        for (size_t col = k + 1; col < n; ++col) {
            double *target = a + col * lda;
            for (size_t row = k + 1; row < n; ++row) {
                target[row] -= column[row] * target[k];
                if (!isfinite(target[row])) {
                    return MATMPC_LU_NONFINITE;
                }
            }
        }
    }
    return MATMPC_LU_OK;
}

static void solve_normal(size_t n, const double *lu, size_t lda, const size_t *pivots, double *rhs)
{
    for (size_t k = 0; k < n; ++k) {
        if (pivots[k] != k) {
            swap(rhs + k, rhs + pivots[k]);
        }
    }
    for (size_t k = 0; k < n; ++k) {
        const double *column = lu + k * lda;
        for (size_t row = k + 1; row < n; ++row) {
            rhs[row] -= column[row] * rhs[k];
        }
    }
    for (size_t end = n; end > 0; --end) {
        const size_t k = end - 1;
        const double *column = lu + k * lda;
        rhs[k] /= column[k];
        for (size_t row = 0; row < k; ++row) {
            rhs[row] -= column[row] * rhs[k];
        }
    }
}

static void solve_transpose(size_t n, const double *lu, size_t lda, const size_t *pivots, double *rhs)
{
    for (size_t k = 0; k < n; ++k) {
        const double *column = lu + k * lda;
        double value = rhs[k];
        for (size_t row = 0; row < k; ++row) {
            value -= column[row] * rhs[row];
        }
        rhs[k] = value / column[k];
    }
    for (size_t end = n; end > 0; --end) {
        const size_t k = end - 1;
        const double *column = lu + k * lda;
        double value = rhs[k];
        for (size_t row = k + 1; row < n; ++row) {
            value -= column[row] * rhs[row];
        }
        rhs[k] = value;
    }
    for (size_t end = n; end > 0; --end) {
        const size_t k = end - 1;
        if (pivots[k] != k) {
            swap(rhs + k, rhs + pivots[k]);
        }
    }
}

matmpc_lu_status matmpc_lu_solve(int transposed, size_t n, size_t nrhs,
                                  const double *lu, size_t lda, const size_t *pivots,
                                  double *rhs, size_t ldb)
{
    if (transposed != 0 && transposed != 1) {
        return MATMPC_LU_INVALID_ARGUMENT;
    }
    if (n == 0 || nrhs == 0) {
        return MATMPC_LU_OK;
    }
    if (lu == NULL || pivots == NULL || rhs == NULL || n > SIZE_MAX / sizeof(*pivots) ||
        !matrix_extent_valid(n, n, lda) || !matrix_extent_valid(n, nrhs, ldb)) {
        return MATMPC_LU_INVALID_ARGUMENT;
    }
    for (size_t k = 0; k < n; ++k) {
        if (pivots[k] < k || pivots[k] >= n) {
            return MATMPC_LU_INVALID_ARGUMENT;
        }
    }
    if (!matrix_is_finite(lu, n, n, lda) || !matrix_is_finite(rhs, n, nrhs, ldb)) {
        return MATMPC_LU_NONFINITE;
    }
    for (size_t k = 0; k < n; ++k) {
        if (lu[k + k * lda] == 0.0) {
            return MATMPC_LU_SINGULAR;
        }
    }

    for (size_t col = 0; col < nrhs; ++col) {
        if (transposed) {
            solve_transpose(n, lu, lda, pivots, rhs + col * ldb);
        } else {
            solve_normal(n, lu, lda, pivots, rhs + col * ldb);
        }
    }
    return matrix_is_finite(rhs, n, nrhs, ldb) ? MATMPC_LU_OK : MATMPC_LU_NONFINITE;
}
