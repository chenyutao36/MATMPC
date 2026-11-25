#include "condense_affine.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void verify(const double *actual, const double *expected, size_t count)
{
    size_t i;
    for (i = 0; i < count; ++i) {
        if (!(fabs(actual[i] - expected[i]) <= 1e-12)) {
            fprintf(stderr, "Element %lu: expected %.17g, got %.17g\n",
                    (unsigned long)i, expected[i], actual[i]);
            exit(EXIT_FAILURE);
        }
    }
}

static void test_multivariable_affine_condensing(void)
{
    const matmpc_condense_dimensions dimensions = {2, 3, 2, 2, 1, 2};
    const double q[] = {2, 1, 99, 3, 4, -1, 101, 2, 6, 2, -91, 5};
    const double s[] = {1, -1, 2, 0, 0, 3, 0, 2, 1, -1, 2, 0};
    const double a[] = {1, 0, 2, -1, 0.5, 1, 0, 2};
    const double b[] = {1, -1, 0, 3, 2, 0, 0, 1, 2, 0, -2, 1};
    const double cx[] = {2, 0.5, -1, 3};
    const double cgx[] = {1, 2, -1, 1.5};
    const double cgn[] = {1, 3, -2, 0.5};
    const double gx[] = {1, -1, -2, 3, 4, 1};
    const double gu[] = {0.25, -0.5, 2, 1, -2, 3};
    const double defect[] = {0.5, -0.5, 1, -2};
    const double ds0[] = {1, -2};
    const double lc[] = {-10, -5, -2, -9}, uc[] = {10, 8, 6, 13};
    const double lbx[] = {-4, -7, -3, -11}, ubx[] = {5, 10, 6, 7};
    const matmpc_condense_affine_inputs inputs = {
        {q, 12}, {s, 12}, {a, 8}, {b, 12}, {cx, 4}, {cgx, 4}, {cgn, 4},
        {gx, 6}, {gu, 6}, {defect, 4}, {ds0, 2}, {lc, 4}, {uc, 4}, {lbx, 4}, {ubx, 4}
    };

    /* Direct substitution gives L=[1,-2,-2.5,1.5,-0.25,-1.5],
       T1=[B0,0] and T2=[A1*B0,B1]. Q uses its stored lower triangle. */
    const double expected_gc[] = {-12, -15, -45.5, -3, -7, -8};
    const double expected_lcc[] = {-7, -9.75, -4.75, -7.5};
    const double expected_ucc[] = {13, 3.25, 3.25, 14.5};
    const double expected_lxc[] = {2.5, -10.25, -4, -6.375};
    const double expected_uxc[] = {11.5, 6.75, 5, 11.625};
    double gc[6], lcc[4] = {NAN, NAN, NAN, NAN}, ucc[4], lxc[4], uxc[4], storage[10];
    matmpc_condense_affine_outputs outputs = {{gc, 6}, {lcc, 4}, {ucc, 4}, {lxc, 4}, {uxc, 4}};
    matmpc_condense_output workspace = {storage, 10};
    matmpc_condense_affine_sizes sizes;

    require(matmpc_condense_affine_get_sizes(&dimensions, &sizes) == MATMPC_CONDENSE_OK &&
            sizes.workspace == 10 && sizes.general == 4 && sizes.state == 4,
            "Multivariable condensing counts are incorrect");
    require(matmpc_condense_affine(&dimensions, &inputs, &outputs, workspace) == MATMPC_CONDENSE_OK,
            "Multivariable affine condensing failed");
    verify(gc, expected_gc, 6);
    verify(lcc, expected_lcc, 4);
    verify(ucc, expected_ucc, 4);
    verify(lxc, expected_lxc, 4);
    verify(uxc, expected_uxc, 4);

    workspace.count = 9;
    gc[0] = 100;
    require(matmpc_condense_affine(&dimensions, &inputs, &outputs, workspace) == MATMPC_CONDENSE_INVALID_BUFFER &&
            gc[0] == 100, "An undersized workspace must be rejected before writing outputs");
    workspace.count = 10;
    outputs.lxc.count = 3;
    require(matmpc_condense_affine(&dimensions, &inputs, &outputs, workspace) == MATMPC_CONDENSE_INVALID_BUFFER &&
            gc[0] == 100, "An undersized bound output must be rejected before writing outputs");
}

static void test_single_stage_without_constraints(void)
{
    matmpc_condense_dimensions dimensions = {1, 1, 1, 0, 0, 0};
    const double q[] = {1, 2}, gx[] = {5, 6};
    const double s = 0.5, a = 3, b = 4, gu = -3, defect = -2, ds0 = 2;
    const matmpc_condense_affine_inputs inputs = {
        {q, 2}, {&s, 1}, {&a, 1}, {&b, 1}, {NULL, 0}, {NULL, 0}, {NULL, 0},
        {gx, 2}, {&gu, 1}, {&defect, 1}, {&ds0, 1}, {NULL, 0}, {NULL, 0}, {NULL, 0}, {NULL, 0}
    };
    double gc, storage[3];
    const matmpc_condense_affine_outputs outputs = {{&gc, 1}, {NULL, 0}, {NULL, 0}, {NULL, 0}, {NULL, 0}};
    const matmpc_condense_output workspace = {storage, 3};
    matmpc_condense_affine_sizes sizes;

    require(matmpc_condense_affine(&dimensions, &inputs, &outputs, workspace) == MATMPC_CONDENSE_OK &&
            gc == 54.0, "Single-stage gradient must be -3 + 0.5*2 + 4*(6 + 2*4)");
    dimensions.horizon = SIZE_MAX;
    require(matmpc_condense_affine_get_sizes(&dimensions, &sizes) == MATMPC_CONDENSE_INVALID_DIMENSIONS,
            "Horizon overflow must be rejected");
    dimensions.horizon = 1;
    dimensions.nx = SIZE_MAX / sizeof(double);
    require(matmpc_condense_affine_get_sizes(&dimensions, &sizes) == MATMPC_CONDENSE_INVALID_DIMENSIONS,
            "State product overflow must be rejected");
}

int main(void)
{
    test_multivariable_affine_condensing();
    test_single_stage_without_constraints();
    return EXIT_SUCCESS;
}
