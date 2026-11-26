#include "irk_kernel.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static const double root_15 = 3.8729833462074168852;
static const double gauss_b[] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};

static void require(bool condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void near(double actual, double expected, double tolerance, const char *message)
{
    if (!(fabs(actual - expected) <= tolerance)) {
        fprintf(stderr, "%s (actual %.17g, expected %.17g)\n", message, actual, expected);
        exit(EXIT_FAILURE);
    }
}

typedef enum {
    COUPLED,
    SINGULAR,
    BAD_RESIDUAL,
    BAD_JACOBIAN,
    CALLBACK_ERROR
} fixture_mode;

typedef struct {
    fixture_mode mode;
    size_t residual_calls;
    size_t jacobian_calls;
    size_t u_jacobian_calls;
} fixture;

static int residual(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, double *result)
{
    fixture *model = context;
    ++model->residual_calls;
    if (model->mode == CALLBACK_ERROR)
        return 1;
    result[0] = 0.1 * xdot[0] + 1.3 * xdot[1] - p[0] * x[0] - sin(x[1]) - u[0];
    result[1] = 1.4 * xdot[0] - 0.25 * xdot[1] - 0.5 * x[0] - x[1] * x[1] - u[1] * x[0];
    if (model->mode == BAD_RESIDUAL)
        result[0] = NAN;
    return 0;
}

static int jacobians(void *context, const double *x, const double *u, const double *p,
                     const double *xdot, double *Fx, double *Fu, double *Fxdot)
{
    (void)xdot;
    fixture *model = context;
    ++model->jacobian_calls;
    Fx[0] = -p[0];
    Fx[1] = -0.5 - u[1];
    Fx[2] = -cos(x[1]);
    Fx[3] = -2.0 * x[1];
    Fxdot[0] = 0.1;
    Fxdot[1] = 1.4;
    Fxdot[2] = 1.3;
    Fxdot[3] = -0.25;
    if (Fu != NULL) {
        ++model->u_jacobian_calls;
        Fu[0] = -1.0;
        Fu[1] = 0.0;
        Fu[2] = 0.0;
        Fu[3] = -x[0];
    }
    if (model->mode == SINGULAR) {
        for (size_t i = 0; i < 4; ++i) {
            Fx[i] = 0.0;
            Fxdot[i] = 0.0;
        }
    }
    if (model->mode == BAD_JACOBIAN)
        Fxdot[0] = INFINITY;
    return 0;
}

static void state_only(matmpc_irk_config config, const matmpc_irk_ode_model *model, const double *x,
                       const double *u, const double *p, double *result, matmpc_irk_workspace workspace)
{
    config.shape.forward = false;
    config.shape.adjoint = false;
    const matmpc_erk_input input = {{x, 2}, {u, 2}, p, {NULL, 0}, {NULL, 0}, {NULL, 0}};
    const matmpc_erk_output output = {{result, 2}, {NULL, 0}, {NULL, 0}, {NULL, 0}};
    require(matmpc_irk_ode_integrate(&config, model, &input, &output, workspace) == MATMPC_IRK_OK,
            "state-only nonlinear integration failed");
}

static void verify_nonlinear(matmpc_irk_config config)
{
    const double x[] = {0.7, -0.2}, u[] = {0.4, 0.6}, p[] = {0.3};
    const double identity[] = {1, 0, 0, 1}, zero[] = {0, 0, 0, 0}, lambda[] = {0.35, -0.8};
    matmpc_erk_input input = {{x, 2}, {u, 2}, p, {identity, 4}, {zero, 4}, {lambda, 2}};
    double final[2], state_jac[4], input_jac[4], gradient[4];
    const matmpc_erk_output output = {{final, 2}, {state_jac, 4}, {input_jac, 4}, {gradient, 4}};
    fixture callbacks = {0};
    const matmpc_irk_ode_model model = {&callbacks, residual, jacobians};
    size_t doubles, pivots;
    double scratch[512];
    size_t rows[32];
    require(matmpc_irk_ode_workspace_size(&config.shape, &doubles, &pivots) == MATMPC_IRK_OK &&
            doubles <= sizeof(scratch) / sizeof(scratch[0]) && pivots <= sizeof(rows) / sizeof(rows[0]),
            "nonlinear IRK workspace query failed");
    const matmpc_irk_workspace workspace = {{scratch, doubles}, rows, pivots};

    require(matmpc_irk_ode_integrate(&config, &model, &input, &output, workspace) == MATMPC_IRK_OK,
            "nonlinear IRK propagation failed");
    const size_t total_stages = config.shape.stages * config.shape.steps;
    require(callbacks.residual_calls == total_stages * config.newton_iterations &&
            callbacks.jacobian_calls == callbacks.residual_calls + total_stages &&
            callbacks.u_jacobian_calls == total_stages,
            "IRK changed the fixed Newton iteration and final linearization convention");

    const double eps = 1e-6;
    for (size_t col = 0; col < 4; ++col) {
        double perturbed_plus[] = {x[0], x[1], u[0], u[1]};
        double perturbed_minus[] = {x[0], x[1], u[0], u[1]};
        double plus[2], minus[2];
        perturbed_plus[col] += eps;
        perturbed_minus[col] -= eps;
        state_only(config, &model, perturbed_plus, perturbed_plus + 2, p, plus, workspace);
        state_only(config, &model, perturbed_minus, perturbed_minus + 2, p, minus, workspace);

        const double *jacobian = col < 2 ? state_jac : input_jac;
        const double *column = jacobian + (col % 2) * 2;
        for (size_t row = 0; row < 2; ++row)
            near(column[row], (plus[row] - minus[row]) / (2.0 * eps), 2e-9, "nonlinear IRK forward derivative");
        near(gradient[col], lambda[0] * column[0] + lambda[1] * column[1], 2e-14,
             "nonlinear IRK adjoint transpose identity");
    }

    double other[2], forward_x[4], forward_u[4], adjoint_only[4];
    state_only(config, &model, x, u, p, other, workspace);
    for (size_t i = 0; i < 2; ++i)
        near(other[i], final[i], 0.0, "state-only repeat changed state");
    matmpc_irk_config one_mode = config;
    one_mode.shape.adjoint = false;
    const matmpc_erk_output forward_output = {{other, 2}, {forward_x, 4}, {forward_u, 4}, {NULL, 0}};
    require(matmpc_irk_ode_integrate(&one_mode, &model, &input, &forward_output, workspace) == MATMPC_IRK_OK,
            "forward-only IRK propagation failed");
    one_mode = config;
    one_mode.shape.forward = false;
    const matmpc_erk_output adjoint_output = {{other, 2}, {NULL, 0}, {NULL, 0}, {adjoint_only, 4}};
    require(matmpc_irk_ode_integrate(&one_mode, &model, &input, &adjoint_output, workspace) == MATMPC_IRK_OK,
            "adjoint-only IRK propagation failed");
    for (size_t i = 0; i < 4; ++i) {
        near(forward_x[i], state_jac[i], 0.0, "forward-only state sensitivities");
        near(forward_u[i], input_jac[i], 0.0, "forward-only control sensitivities");
        near(adjoint_only[i], gradient[i], 0.0, "adjoint-only sensitivities");
    }

    const double seed_x[] = {0.8, -0.2, 0.4, 1.3};
    const double seed_u[] = {0.3, -0.1, -0.2, 0.2};
    input.Sx.data = seed_x;
    input.Su.data = seed_u;
    one_mode = config;
    one_mode.shape.adjoint = false;
    require(matmpc_irk_ode_integrate(&one_mode, &model, &input, &forward_output, workspace) == MATMPC_IRK_OK,
            "IRK with explicit sensitivity seeds failed");
    for (size_t col = 0; col < 2; ++col) {
        for (size_t row = 0; row < 2; ++row) {
            const size_t i = col * 2 + row;
            near(forward_x[i], state_jac[row] * seed_x[col * 2] + state_jac[row + 2] * seed_x[col * 2 + 1],
                 2e-14, "nonidentity IRK state seed");
            near(forward_u[i], input_jac[i] + state_jac[row] * seed_u[col * 2] +
                 state_jac[row + 2] * seed_u[col * 2 + 1], 2e-14, "nonzero IRK control seed");
        }
    }

    matmpc_irk_workspace undersized = workspace;
    --undersized.doubles.count;
    require(matmpc_irk_ode_integrate(&config, &model, &input, &output, undersized) == MATMPC_IRK_WORKSPACE_TOO_SMALL,
            "IRK double workspace size unchecked");
    undersized = workspace;
    --undersized.pivot_count;
    require(matmpc_irk_ode_integrate(&config, &model, &input, &output, undersized) == MATMPC_IRK_WORKSPACE_TOO_SMALL,
            "IRK pivot workspace size unchecked");

    fixture_mode failure_cases[] = {SINGULAR, BAD_RESIDUAL, BAD_JACOBIAN, CALLBACK_ERROR, COUPLED};
    const matmpc_irk_status expected_failures[] = {
        MATMPC_IRK_SINGULAR, MATMPC_IRK_NONFINITE, MATMPC_IRK_NONFINITE, MATMPC_IRK_MODEL_ERROR, MATMPC_IRK_NONCONVERGENCE
    };
    for (size_t i = 0; i < sizeof(failure_cases) / sizeof(failure_cases[0]); ++i) {
        matmpc_irk_config altered = config;
        altered.newton_iterations = 1;
        callbacks.mode = failure_cases[i];
        require(matmpc_irk_ode_integrate(&altered, &model, &input, &output, workspace) == expected_failures[i],
                "unexpected IRK failure status");
    }
    state_only(config, &model, x, u, p, other, workspace);
    require(isfinite(other[0]) && isfinite(other[1]), "failed IRK call contaminated subsequent workspace reuse");
    --input.u.count;
    require(matmpc_irk_ode_integrate(&config, &model, &input, &output, workspace) == MATMPC_IRK_INVALID_ARGUMENT,
            "IRK control shape unchecked");
}

static int linear_residual(void *context, const double *x, const double *u, const double *p,
                           const double *xdot, double *result)
{
    (void)context;
    (void)u;
    result[0] = xdot[0] - p[0] * x[0] - p[1];
    return 0;
}

static int linear_jacobians(void *context, const double *x, const double *u, const double *p,
                            const double *xdot, double *Fx, double *Fu, double *Fxdot)
{
    (void)context;
    (void)x;
    (void)u;
    (void)xdot;
    require(Fu == NULL, "a control Jacobian was requested for an autonomous ODE");
    Fx[0] = -p[0];
    Fxdot[0] = 1.0;
    return 0;
}

static void verify_scalar_gauss(matmpc_irk_config config)
{
    const double x = 0.24, seed = 1.0, lambda = -0.28;
    const double p[] = {-0.75, 0.7};
    const matmpc_erk_input input = {{&x, 1}, {NULL, 0}, p, {&seed, 1}, {NULL, 0}, {&lambda, 1}};
    double xf, jac, gradient;
    const matmpc_erk_output output = {{&xf, 1}, {&jac, 1}, {NULL, 0}, {&gradient, 1}};
    const matmpc_irk_ode_model model = {NULL, linear_residual, linear_jacobians};
    size_t doubles, pivots;
    double scratch[128];
    size_t rows[16];
    require(matmpc_irk_ode_workspace_size(&config.shape, &doubles, &pivots) == MATMPC_IRK_OK &&
            doubles <= sizeof(scratch) / sizeof(scratch[0]) && pivots <= sizeof(rows) / sizeof(rows[0]),
            "scalar IRK workspace query failed");
    const matmpc_irk_workspace workspace = {{scratch, doubles}, rows, pivots};

    const double base_step = 0.052;
    double jfk[3] = {0};
    const double steps[] = {0.031, 0.16};
    for (size_t i = 0; i < sizeof(steps) / sizeof(steps[0]); ++i) {
        config.step = steps[i];
        const double z = p[0] * config.step;
        /* Stability function of the three-stage Gauss method on a linear ODE. */
        const double r = (120.0 + 60.0 * z + 12.0 * z * z + z * z * z) /
                         (120.0 - 60.0 * z + 12.0 * z * z - z * z * z);
        const double final = pow(r, (double)config.shape.steps) * (x + p[1] / p[0]) - p[1] / p[0];
        for (size_t gain = 0; gain <= 2; ++gain) {
            config.sensitivity_weights = (matmpc_erk_view){NULL, 0};
            if (gain != 0) {
                for (size_t stage = 0; stage < 3; ++stage)
                    jfk[stage] = (double)gain * base_step * gauss_b[stage];
                config.sensitivity_weights = (matmpc_erk_view){jfk, 3};
            }
            config.base_step = base_step;
            require(matmpc_irk_ode_integrate(&config, &model, &input, &output, workspace) == MATMPC_IRK_OK,
                    "scalar IRK propagation failed");
            const double tangent_factor = gain == 2 ? 1.0 + 2.0 * (r - 1.0) : r;
            const double tangent = pow(tangent_factor, (double)config.shape.steps);
            near(xf, final, 2e-15, "Gauss IRK state versus Padé formula");
            near(jac, tangent, 2e-15, "Gauss IRK forward scaling");
            near(gradient, lambda * tangent, 2e-15, "Gauss IRK adjoint scaling");
        }
    }
    config.base_step = 0.0;
    require(matmpc_irk_ode_integrate(&config, &model, &input, &output, workspace) == MATMPC_IRK_INVALID_ARGUMENT,
            "IRK accepted custom sensitivity weights with a zero base step");
}

int main(void)
{
    const double gauss_a[] = {
        5.0 / 36.0, 5.0 / 36.0 + root_15 / 24.0, 5.0 / 36.0 + root_15 / 30.0,
        2.0 / 9.0 - root_15 / 15.0, 2.0 / 9.0, 2.0 / 9.0 + root_15 / 15.0,
        5.0 / 36.0 - root_15 / 30.0, 5.0 / 36.0 - root_15 / 24.0, 5.0 / 36.0
    };
    matmpc_irk_config config = {{2, 2, 3, 3, true, true}, 0.08, {gauss_a, 9}, {gauss_b, 3},
                               5, 1e-2, {NULL, 0}, 0.0};
    verify_nonlinear(config);
    config.shape.nx = 1;
    config.shape.nu = 0;
    verify_scalar_gauss(config);

    size_t doubles = 1, pivots = 1;
    const matmpc_erk_shape overflow = {SIZE_MAX, 1, 3, 2, false, true};
    require(matmpc_irk_ode_workspace_size(&overflow, &doubles, &pivots) == MATMPC_IRK_SIZE_OVERFLOW &&
            doubles == 0 && pivots == 0, "IRK did not reject workspace size overflow");
    return EXIT_SUCCESS;
}
