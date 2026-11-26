#include "matmpc_lu.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void test_pivoting_and_padded_multiple_rhs(void)
{
    enum { N = 3, LDA = 5, LDB = 4, RHS = 2 };
    const double matrix[LDA * N] = {0, 1, 4, NAN, 19, 5, 0, 1, NAN, 29, 2, 3, 6, NAN, 39};
    const double expected[N * RHS] = {1, 2, -1, 3, -2, 4};
    double lu[LDA * N], rhs[LDB * RHS];
    size_t pivots[N];

    memcpy(lu, matrix, sizeof(lu));
    require(matmpc_lu_factor(N, lu, LDA, pivots) == MATMPC_LU_OK, "Pivoting failed");
    require(pivots[0] == 2 && pivots[1] == 2 && pivots[2] == 2, "Incorrect sequential row swaps");
    for (size_t col = 0; col < N; ++col) {
        require(isnan(lu[3 + col * LDA]) && lu[4 + col * LDA] == matrix[4 + col * LDA],
                "Factorization changed matrix padding");
    }

    for (int transpose = 0; transpose <= 1; ++transpose) {
        for (size_t col = 0; col < RHS; ++col) {
            for (size_t row = 0; row < N; ++row) {
                double value = 0;
                for (size_t k = 0; k < N; ++k) {
                    const double coefficient = transpose ? matrix[k + row * LDA] : matrix[row + k * LDA];
                    value += coefficient * expected[k + col * N];
                }
                rhs[row + col * LDB] = value;
            }
            rhs[N + col * LDB] = NAN;
        }
        require(matmpc_lu_solve(transpose, N, RHS, lu, LDA, pivots, rhs, LDB) == MATMPC_LU_OK,
                "A multiple-RHS solve failed");
        for (size_t col = 0; col < RHS; ++col) {
            for (size_t row = 0; row < N; ++row) {
                require(fabs(rhs[row + col * LDB] - expected[row + col * N]) <= 1e-12,
                        "A multiple-RHS solution does not reproduce the expected vector");
            }
            require(isnan(rhs[N + col * LDB]), "Solve changed RHS padding");
        }
    }
}

static void test_tiny_nonzero_pivots(void)
{
    double scale = DBL_MIN / 8.0;
    if (scale == 0.0) {
        scale = DBL_MIN;
    }
    double matrix[] = {scale / 2, scale, 0, scale};
    double rhs[] = {scale, 5 * scale};
    size_t pivots[2];

    require(matmpc_lu_factor(2, matrix, 2, pivots) == MATMPC_LU_OK,
            "A tiny nonsingular system must be accepted");
    require(matmpc_lu_solve(0, 2, 1, matrix, 2, pivots, rhs, 2) == MATMPC_LU_OK &&
            rhs[0] == 2 && rhs[1] == 3, "A tiny nonsingular system gave the wrong solution");
}

static void test_errors_and_empty_work(void)
{
    double singular[] = {1, 2, 2, 4};
    double nonfinite[] = {1, 0, NAN, 1};
    double overflow[] = {DBL_MAX, DBL_MAX, DBL_MAX, -DBL_MAX};
    double identity[] = {1, 0, 0, 1};
    const double zero_diagonal[] = {1, 0, 0, 0};
    size_t pivots[] = {0, 1};
    const size_t invalid_pivots[] = {0, 0};
    double rhs[] = {1, 2};
    double invalid_rhs[] = {1, INFINITY};

    require(matmpc_lu_factor(0, NULL, 0, NULL) == MATMPC_LU_OK, "Empty factorization failed");
    require(matmpc_lu_solve(1, 2, 0, NULL, 0, NULL, NULL, 0) == MATMPC_LU_OK, "Empty solve failed");
    require(matmpc_lu_factor(2, identity, 1, pivots) == MATMPC_LU_INVALID_ARGUMENT,
            "Reject undersized leading dimension");
    require(matmpc_lu_factor(2, identity, SIZE_MAX, pivots) == MATMPC_LU_INVALID_ARGUMENT,
            "Reject overflowing matrix extent");
    require(matmpc_lu_factor(2, singular, 2, pivots) == MATMPC_LU_SINGULAR, "Reject a zero pivot");
    require(matmpc_lu_factor(2, nonfinite, 2, pivots) == MATMPC_LU_NONFINITE, "Reject nonfinite input");
    require(matmpc_lu_factor(2, overflow, 2, pivots) == MATMPC_LU_NONFINITE, "Reject elimination overflow");

    pivots[0] = 0;
    pivots[1] = 1;
    require(matmpc_lu_solve(2, 2, 1, identity, 2, pivots, rhs, 2) == MATMPC_LU_INVALID_ARGUMENT,
            "Reject unknown transpose flags");
    require(matmpc_lu_solve(0, 2, 2, identity, 2, pivots, rhs, SIZE_MAX) == MATMPC_LU_INVALID_ARGUMENT,
            "Reject overflowing RHS extent");
    require(matmpc_lu_solve(0, 2, 1, identity, 2, invalid_pivots, rhs, 2) == MATMPC_LU_INVALID_ARGUMENT,
            "Reject invalid swap rows");
    require(matmpc_lu_solve(1, 2, 1, zero_diagonal, 2, pivots, rhs, 2) == MATMPC_LU_SINGULAR,
            "Reject singular factor before modifying RHS");
    require(matmpc_lu_solve(0, 2, 1, nonfinite, 2, pivots, rhs, 2) == MATMPC_LU_NONFINITE,
            "Reject nonfinite factor before modifying RHS");
    require(rhs[0] == 1 && rhs[1] == 2, "A rejected factor or dimension changed the RHS");
    require(matmpc_lu_solve(0, 2, 1, identity, 2, pivots, invalid_rhs, 2) == MATMPC_LU_NONFINITE &&
            invalid_rhs[0] == 1 && isinf(invalid_rhs[1]), "Reject nonfinite RHS before modifying it");

    double large_rhs = DBL_MAX;
    const double tiny_diagonal = DBL_MIN;
    require(matmpc_lu_solve(0, 1, 1, &tiny_diagonal, 1, pivots, &large_rhs, 1) == MATMPC_LU_NONFINITE,
            "Reject solve overflow");
}

int main(void)
{
    test_pivoting_and_padded_multiple_rhs();
    test_tiny_nonzero_pivots();
    test_errors_and_empty_work();
    return EXIT_SUCCESS;
}
