#include "matmpc_la.h"

#include <math.h>
#include <stdio.h>

int main(void)
{
    const double a[] = {1, 4, 2, 5, 3, 6};
    const double b[] = {7, 9, 11, 8, 10, 12};
    const double b_transposed[] = {7, 8, 9, 10, 11, 12};
    double product[] = {NAN, NAN, 19, NAN, NAN, 21};
    const double x[] = {1, 2, 3};
    const double xt[] = {2, -1};
    double vector[] = {3, 2, NAN};
    const double symmetric[] = {2, 3, NAN, 5};
    const double rhs[] = {1, 2, 3, 4};
    double result[] = {NAN, NAN, NAN, NAN};

    matmpc_la_gemm(0, 0, 2, 2, 3, 1, a, 2, b, 3, 0, product, 3);
    if (product[0] != 58 || product[1] != 139 || product[2] != 19 ||
        product[3] != 64 || product[4] != 154 || product[5] != 21) {
        fputs("Column-major GEMM or leading dimension is incorrect.\n", stderr);
        return 1;
    }
    matmpc_la_gemm(0, 1, 2, 2, 3, 1, a, 2, b_transposed, 2, 0, product, 3);
    if (product[0] != 58 || product[1] != 139 || product[3] != 64 || product[4] != 154) {
        fputs("Transposed GEMM operand is incorrect.\n", stderr);
        return 1;
    }
    matmpc_la_gemv(0, 2, 3, 2, a, 2, x, 3, vector);
    if (vector[0] != 37 || vector[1] != 70) {
        fputs("Normal GEMV or beta accumulation is incorrect.\n", stderr);
        return 1;
    }
    matmpc_la_gemv(1, 2, 3, 1, a, 2, xt, 0, vector);
    if (vector[0] != -2 || vector[1] != -1 || vector[2] != 0) {
        fputs("Transpose GEMV or beta=0 is incorrect.\n", stderr);
        return 1;
    }
    matmpc_la_symv_lower(2, 1, symmetric, 2, rhs, 0, vector);
    matmpc_la_symm_left_lower(2, 2, 1, symmetric, 2, rhs, 2, 0, result, 2);
    if (vector[0] != 8 || vector[1] != 13 ||
        result[0] != 8 || result[1] != 13 || result[2] != 18 || result[3] != 29) {
        fputs("Symmetric operations read the wrong triangle.\n", stderr);
        return 1;
    }
    return 0;
}
