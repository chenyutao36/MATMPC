#include "erk_kernel.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static const double rk4_a[] = {0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0};
static const double rk4_b[] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
static const double midpoint_a[] = {0, 0.5, 0, 0};
static const double midpoint_b[] = {0, 1};

static void require(bool condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static int rhs(void *context, const double *x, const double *u, const double *p, double *rate)
{
    if (context != NULL && *(const bool *)context)
        return 1;
    rate[0] = x[1] + sin(u[0]) + p[0] * x[0];
    rate[1] = cos(x[0]) - 0.4 * x[1] + u[0] * x[0];
    return 0;
}

static int forward(void *context, const double *x, const double *u, const double *p,
                   const double *Sx, const double *Su, double *dSx, double *dSu)
{
    (void)context;
    const double fx_10 = -sin(x[0]) + u[0];
    for (size_t column = 0; column < 2; ++column) {
        dSx[2 * column] = p[0] * Sx[2 * column] + Sx[2 * column + 1];
        dSx[2 * column + 1] = fx_10 * Sx[2 * column] - 0.4 * Sx[2 * column + 1];
    }
    dSu[0] = p[0] * Su[0] + Su[1] + cos(u[0]);
    dSu[1] = fx_10 * Su[0] - 0.4 * Su[1] + x[0];
    return 0;
}

static int adjoint(void *context, const double *x, const double *u, const double *p,
                   const double *lambda, double *product)
{
    (void)context;
    product[0] = p[0] * lambda[0] + (-sin(x[0]) + u[0]) * lambda[1];
    product[1] = lambda[0] - 0.4 * lambda[1];
    product[2] = cos(u[0]) * lambda[0] + x[0] * lambda[1];
    return 0;
}

static void reference_rk4(const double *x, const double *u, const double *p, double h, size_t steps, double *result)
{
    double k1[2], k2[2], k3[2], k4[2], stage[2];
    result[0] = x[0];
    result[1] = x[1];
    for (size_t step = 0; step < steps; ++step) {
        rhs(NULL, result, u, p, k1);
        for (size_t i = 0; i < 2; ++i)
            stage[i] = result[i] + h * 0.5 * k1[i];
        rhs(NULL, stage, u, p, k2);
        for (size_t i = 0; i < 2; ++i)
            stage[i] = result[i] + h * 0.5 * k2[i];
        rhs(NULL, stage, u, p, k3);
        for (size_t i = 0; i < 2; ++i)
            stage[i] = result[i] + h * k3[i];
        rhs(NULL, stage, u, p, k4);
        for (size_t i = 0; i < 2; ++i)
            result[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }
}

static void integrate_state(matmpc_erk_config config, const matmpc_erk_model *model,
                            const double *x, const double *u, const double *p, double *result,
                            matmpc_erk_buffer workspace)
{
    config.shape.forward = false;
    config.shape.adjoint = false;
    const matmpc_erk_input input = {{x, 2}, {u, 1}, p, {NULL, 0}, {NULL, 0}, {NULL, 0}};
    const matmpc_erk_output output = {{result, 2}, {NULL, 0}, {NULL, 0}, {NULL, 0}};
    require(matmpc_erk_integrate(&config, model, &input, &output, workspace) == MATMPC_ERK_OK,
            "state-only integration failed");
}

static void verify_sensitivities(matmpc_erk_config config)
{
    const double x[] = {0.7, -0.2}, u[] = {0.4}, p[] = {0.3};
    const double seed_x[] = {1, 0, 0, 1}, seed_u[] = {0, 0}, lambda[] = {0.35, -0.8};
    double x_final[2], Sx[4], Su[2], gradient[3], scratch[120];
    const matmpc_erk_model model = {NULL, rhs, forward, adjoint};
    const matmpc_erk_input input = {{x, 2}, {u, 1}, p, {seed_x, 4}, {seed_u, 2}, {lambda, 2}};
    const matmpc_erk_output output = {{x_final, 2}, {Sx, 4}, {Su, 2}, {gradient, 3}};
    size_t count = 0;
    require(matmpc_erk_workspace_size(&config.shape, &count) == MATMPC_ERK_OK && count <= 120,
            "workspace query failed");
    const matmpc_erk_buffer workspace = {scratch, count};
    require(matmpc_erk_integrate(&config, &model, &input, &output, workspace) == MATMPC_ERK_OK,
            "forward and adjoint integration failed");

    if (config.shape.stages == 4) {
        double expected[2];
        reference_rk4(x, u, p, config.step, config.shape.steps, expected);
        for (size_t i = 0; i < 2; ++i)
            require(fabs(expected[i] - x_final[i]) < 1e-14, "RK4 result changed from reference");
    }
    const double eps = 1e-6;
    for (size_t j = 0; j < 3; ++j) {
        double xp[] = {x[0], x[1], u[0]}, xm[] = {x[0], x[1], u[0]};
        double plus[2], minus[2];
        xp[j] += eps;
        xm[j] -= eps;
        integrate_state(config, &model, xp, xp + 2, p, plus, workspace);
        integrate_state(config, &model, xm, xm + 2, p, minus, workspace);

        const double *column = j < 2 ? Sx + 2 * j : Su;
        for (size_t k = 0; k < 2; ++k)
            require(fabs((plus[k] - minus[k]) / (2 * eps) - column[k]) < 3e-9,
                    "forward derivative mismatch");
        require(fabs(gradient[j] - (column[0] * lambda[0] + column[1] * lambda[1])) < 2e-14,
                "adjoint does not match transposed forward sensitivities");
    }

    const matmpc_erk_buffer too_small = {scratch, count - 1};
    require(matmpc_erk_integrate(&config, &model, &input, &output, too_small) == MATMPC_ERK_WORKSPACE_TOO_SMALL,
            "did not reject an undersized workspace");
    matmpc_erk_input bad_input = input;
    bad_input.u.count = 0;
    require(matmpc_erk_integrate(&config, &model, &bad_input, &output, workspace) == MATMPC_ERK_INVALID_ARGUMENT,
            "did not reject a wrong-sized control buffer");
    matmpc_erk_config bad_config = config;
    bad_config.step = 0;
    require(matmpc_erk_integrate(&bad_config, &model, &input, &output, workspace) == MATMPC_ERK_INVALID_ARGUMENT,
            "did not reject an invalid step size");

    bool fail = true;
    matmpc_erk_model failing_model = model;
    failing_model.context = &fail;
    require(matmpc_erk_integrate(&config, &failing_model, &input, &output, workspace) == MATMPC_ERK_MODEL_ERROR,
            "did not propagate a callback failure");
}

int main(void)
{
    const matmpc_erk_config rk4 = {{2, 1, 4, 3, true, true}, 0.08, {rk4_a, 16}, {rk4_b, 4}};
    const matmpc_erk_config midpoint = {{2, 1, 2, 3, true, true}, 0.08, {midpoint_a, 4}, {midpoint_b, 2}};
    verify_sensitivities(rk4);
    verify_sensitivities(midpoint);

    size_t count;
    const matmpc_erk_shape oversized = {SIZE_MAX, 1, 4, 2, true, true};
    require(matmpc_erk_workspace_size(&oversized, &count) == MATMPC_ERK_SIZE_OVERFLOW && count == 0,
            "did not detect overflowing workspace size");
    return EXIT_SUCCESS;
}
