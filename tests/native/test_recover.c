#include "recover.h"

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

static matmpc_recover_buffer buffer(const double *data, size_t count)
{
    matmpc_recover_buffer result = {data, count};
    return result;
}

static void check_trajectory(const double *actual, const double *expected, size_t count)
{
    size_t i;
    for (i = 0; i < count; ++i) {
        if (actual[i] != expected[i]) {
            fprintf(stderr, "Element %lu: expected %.17g, got %.17g\n",
                    (unsigned long)i, expected[i], actual[i]);
            exit(EXIT_FAILURE);
        }
    }
}

static void test_constrained_trajectory(void)
{
    const matmpc_recover_dimensions dimensions = {2, 3, 2, 1, 2, 1};
    const double q[] = {2, 1, 1, 3, 4, -1, -1, 2, 6, 2, 2, 5};
    const double s[] = {1, -1, 2, 0, 0, 3, 0, 2, 1, -1, 2, 0};
    const double a[] = {1, 0, 2, -1, 0.5, 1, 0, 2};
    const double b[] = {1, -1, 0, 3, 2, 0, 0, 1, 2, 0, -2, 1};
    const double cx[] = {2, -1};
    const double cgx[] = {1, -1, 2, 1, 0, 2, 3, -2};
    const double cgn[] = {1, -2};
    const double gx[] = {1, -1, -2, 3, 4, 1};
    const double defect[] = {0.5, -0.5, 1, -2};
    const double ds0[] = {1, -2};
    const double du[] = {0.5, -1, 2, -2, 0.5, 1};
    const double mu[] = {0.5, 1.25, -2, 0.5, 1.75};
    const double mu_x[] = {1.5, -0.5};
    const double expected_dx[] = {1, -2, 2, -2, 1, -5};
    const double expected_lambda[] = {11.375, -47.5, -10.125, -66, 0.75, -25};
    const matmpc_recover_inputs inputs = {
        {q, 12}, {s, 12}, {a, 8}, {b, 12}, {cx, 2}, {cgx, 8}, {cgn, 2},
        {gx, 6}, {defect, 4}, {ds0, 2}, {du, 6}, {mu, 5}, {mu_x, 2}
    };
    matmpc_recover_sizes sizes;
    double recovered_dx[6], recovered_lambda[6];

    require(matmpc_recover_get_sizes(&dimensions, &sizes) == MATMPC_RECOVER_OK,
            "Constrained problem dimensions must be accepted");
    require(sizes.dx == 6 && sizes.du == 6 && sizes.q == 12 && sizes.s == 12 && sizes.a == 8 &&
            sizes.cx == 2 && sizes.cgx == 8 && sizes.cgn == 2 && sizes.defect == 4 &&
            sizes.mu == 5 && sizes.mu_x == 2, "Constrained problem buffer sizes are incorrect");
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 6, recovered_lambda, 6) == MATMPC_RECOVER_OK,
            "Constrained recovery failed");
    check_trajectory(recovered_dx, expected_dx, 6);
    check_trajectory(recovered_lambda, expected_lambda, 6);
}

static void test_single_stage_and_rejected_buffers(void)
{
    matmpc_recover_dimensions dimensions = {1, 0, 1, 0, 0, 0};
    const double q[] = {1, 2};
    const double a[] = {3};
    const double gx[] = {5, 6};
    const double defect[] = {-2};
    const double ds0[] = {2};
    const double expected_dx[] = {2, 4};
    const double expected_lambda[] = {-49, 14};
    matmpc_recover_inputs inputs = {0};
    double recovered_dx[2], recovered_lambda[2];

    inputs.q = buffer(q, 2);
    inputs.a = buffer(a, 1);
    inputs.gx = buffer(gx, 2);
    inputs.defect = buffer(defect, 1);
    inputs.ds0 = buffer(ds0, 1);
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 2, recovered_lambda, 2) == MATMPC_RECOVER_OK,
            "Single-stage autonomous recovery failed");
    check_trajectory(recovered_dx, expected_dx, 2);
    check_trajectory(recovered_lambda, expected_lambda, 2);

    dimensions.horizon = 0;
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 2, recovered_lambda, 2) ==
            MATMPC_RECOVER_INVALID_DIMENSIONS, "Zero horizon must be rejected");
    dimensions.horizon = SIZE_MAX;
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 2, recovered_lambda, 2) ==
            MATMPC_RECOVER_INVALID_DIMENSIONS, "An overflowing horizon must be rejected");
    dimensions.horizon = 1;
    dimensions.nx = SIZE_MAX / sizeof(double);
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 2, recovered_lambda, 2) ==
            MATMPC_RECOVER_INVALID_DIMENSIONS, "Overflowing state products must be rejected");
    dimensions.nx = 1;

    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 1, recovered_lambda, 2) ==
            MATMPC_RECOVER_INVALID_BUFFER, "An undersized output must be rejected");
    inputs.q.count = 1;
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 2, recovered_lambda, 2) ==
            MATMPC_RECOVER_INVALID_BUFFER, "An undersized input must be rejected");
    inputs.q = buffer(NULL, 2);
    require(matmpc_recover(&dimensions, &inputs, recovered_dx, 2, recovered_lambda, 2) ==
            MATMPC_RECOVER_INVALID_BUFFER, "A null required input must be rejected");
    check_trajectory(recovered_dx, expected_dx, 2);
    check_trajectory(recovered_lambda, expected_lambda, 2);
}

int main(void)
{
    test_constrained_trajectory();
    test_single_stage_and_rejected_buffers();
    return EXIT_SUCCESS;
}
