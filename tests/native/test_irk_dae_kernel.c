#include "irk_dae_kernel.h"

#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum { VALID, MODEL_FAILURE, SINGULAR_SYSTEM, NONFINITE_RESIDUAL } model_mode;

/* Reduced dynamics x' = [-2,1;0,-0.5] x + u. */
static int coupled_residual(void *context, const double *x, const double *u, const double *p,
                            const double *rate, const double *z, double *f, double *g)
{
    (void)p;
    model_mode mode = *(const model_mode *)context;
    if (mode == MODEL_FAILURE) {
        return 1;
    }
    f[0] = rate[0] - z[0];
    f[1] = rate[1] + 0.5 * x[1] - u[1];
    g[0] = mode == NONFINITE_RESIDUAL ? NAN : z[0] + 2.0 * x[0] - x[1] - u[0];
    return 0;
}

static int coupled_jacobian(void *context, const double *x, const double *u, const double *p,
                            const double *rate, const double *z, double *fx, double *fu,
                            double *fdot, double *fz, double *gx, double *gu, double *gz)
{
    (void)x;
    (void)u;
    (void)p;
    (void)rate;
    (void)z;
    const double slope = *(const model_mode *)context == SINGULAR_SYSTEM ? 0.0 : 1.0;
    fx[0] = fx[1] = fx[2] = 0.0;
    fx[3] = 0.5 * slope;
    fdot[0] = fdot[3] = slope;
    fdot[1] = fdot[2] = 0.0;
    fz[0] = -slope;
    fz[1] = 0.0;
    gx[0] = 2.0 * slope;
    gx[1] = -slope;
    gz[0] = slope;
    if (fu != NULL) {
        fu[0] = fu[1] = fu[2] = 0.0;
        fu[3] = -slope;
        gu[0] = -slope;
        gu[1] = 0.0;
    }
    return 0;
}

static matmpc_irk_workspace allocate_workspace(const matmpc_erk_shape *shape, size_t nz)
{
    size_t count, pivots;
    assert(matmpc_irk_dae_workspace_size(shape, nz, &count, &pivots) == MATMPC_IRK_OK);
    matmpc_irk_workspace workspace = {{malloc(count * sizeof(double)), count},
                                      malloc(pivots * sizeof(size_t)), pivots};
    assert(workspace.doubles.data != NULL && workspace.pivots != NULL);
    for (size_t i = 0; i < count; ++i) {
        workspace.doubles.data[i] = NAN;
    }
    for (size_t i = 0; i < pivots; ++i) {
        workspace.pivots[i] = SIZE_MAX;
    }
    return workspace;
}

static void free_workspace(matmpc_irk_workspace workspace)
{
    free(workspace.doubles.data);
    free(workspace.pivots);
}

static void near(double actual, double expected)
{
    assert(fabs(actual - expected) < 2e-12 * (1.0 + fabs(expected)));
}

/* Stability function of the three-stage Radau IIA method. */
static double radau_gain(double z)
{
    return (1.0 + (2.0 / 5.0) * z + z * z / 20.0) /
           (1.0 - (3.0 / 5.0) * z + 3.0 * z * z / 20.0 - z * z * z / 60.0);
}

static void test_coupled_radau_and_failure_reuse(void)
{
    const double root = sqrt(6.0);
    const double a[] = {(88.0 - 7.0 * root) / 360.0, (296.0 + 169.0 * root) / 1800.0,
                        (16.0 - root) / 36.0, (296.0 - 169.0 * root) / 1800.0,
                        (88.0 + 7.0 * root) / 360.0, (16.0 + root) / 36.0,
                        (-2.0 + 3.0 * root) / 225.0, (-2.0 - 3.0 * root) / 225.0, 1.0 / 9.0};
    const double b[] = {a[2], a[5], a[8]};
    const double x[] = {1.2, -0.7}, u[] = {0.9, 1.1}, z = 20.0;
    const double seed_x[] = {0.9, -0.4, 0.2, 1.1}, seed_u[] = {0.2, 0.4, 0.3, 0.5};
    const double lambda[] = {0.8, -1.3};
    matmpc_irk_config config = {{2, 2, 3, 4, true, true}, 0.13, {a, 9}, {b, 3}, 5, 1e-10, {NULL, 0}, 0.0};
    model_mode mode = VALID;
    const matmpc_irk_dae_model model = {&mode, coupled_residual, coupled_jacobian};
    const matmpc_irk_dae_input in = {{{x, 2}, {u, 2}, NULL, {seed_x, 4}, {seed_u, 4}, {lambda, 2}}, {&z, 1}};
    double xn[2], zn, sx[4], su[4], adj[4];
    const matmpc_irk_dae_output out = {{{xn, 2}, {sx, 4}, {su, 4}, {adj, 4}}, {&zn, 1}};
    matmpc_irk_workspace workspace = allocate_workspace(&config.shape, 1);

    const double p = pow(radau_gain(-2.0 * config.step), (double)config.shape.steps);
    const double q = pow(radau_gain(-0.5 * config.step), (double)config.shape.steps);
    const double g0 = (1.0 - p) / 2.0, g1 = (1.0 - q) / 0.5, c = 2.0 / 3.0;
    const double A[] = {p, 0.0, c * (q - p), q};
    const double B[] = {g0, 0.0, c * (g1 - g0), g1};

    for (size_t pass = 0; pass < 2; ++pass) {
        assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == MATMPC_IRK_OK);
        for (size_t i = 0; i < 2; ++i) {
            near(xn[i], A[i] * x[0] + A[i + 2] * x[1] + B[i] * u[0] + B[i + 2] * u[1]);
        }
        near(zn, -2.0 * xn[0] + xn[1] + u[0]);
        for (size_t j = 0; j < 2; ++j) {
            near(adj[j], A[2 * j] * lambda[0] + A[2 * j + 1] * lambda[1]);
            near(adj[2 + j], B[2 * j] * lambda[0] + B[2 * j + 1] * lambda[1]);
            for (size_t i = 0; i < 2; ++i) {
                near(sx[i + 2 * j], A[i] * seed_x[2 * j] + A[i + 2] * seed_x[2 * j + 1]);
                near(su[i + 2 * j], A[i] * seed_u[2 * j] + A[i + 2] * seed_u[2 * j + 1] + B[i + 2 * j]);
            }
        }
        const double saved = zn;
        const model_mode bad[] = {MODEL_FAILURE, SINGULAR_SYSTEM, NONFINITE_RESIDUAL};
        const matmpc_irk_status errors[] = {MATMPC_IRK_MODEL_ERROR, MATMPC_IRK_SINGULAR, MATMPC_IRK_NONFINITE};
        for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
            mode = bad[i];
            assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == errors[i]);
            assert(zn == saved);
        }
        mode = VALID;
        config.newton_iterations = 1;
        assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == MATMPC_IRK_NONCONVERGENCE);
        assert(zn == saved);
        config.newton_iterations = 5;
    }
    --workspace.pivot_count;
    assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == MATMPC_IRK_WORKSPACE_TOO_SMALL);
    ++workspace.pivot_count;
    --workspace.doubles.count;
    assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == MATMPC_IRK_WORKSPACE_TOO_SMALL);
    free_workspace(workspace);
}

static int scalar_residual(void *context, const double *x, const double *u, const double *p,
                           const double *rate, const double *z, double *f, double *g)
{
    (void)context;
    (void)p;
    *f = *rate - *z;
    *g = *z + 2.0 * *x - (u == NULL ? 0.0 : *u);
    return 0;
}

static int scalar_jacobian(void *context, const double *x, const double *u, const double *p,
                           const double *rate, const double *z, double *fx, double *fu,
                           double *fdot, double *fz, double *gx, double *gu, double *gz)
{
    (void)context;
    (void)x;
    (void)u;
    (void)p;
    (void)rate;
    (void)z;
    *fx = 0.0;
    *fdot = 1.0;
    *fz = -1.0;
    *gx = 2.0;
    *gz = 1.0;
    if (fu != NULL) {
        *fu = 0.0;
        *gu = -1.0;
    }
    return 0;
}

static void test_weights_and_algebraic_stage(void)
{
    const double x = 2.0, u = 1.0, z = -30.0, a = 1.0, b = 1.0, seed = 1.0, zero = 0.0, lambda = 2.0;
    const double weights[] = {0.125, 0.25};
    matmpc_irk_config config = {{1, 1, 1, 1, true, true}, 0.25, {&a, 1}, {&b, 1}, 5, 1e-12, {weights, 2}, 0.125};
    const matmpc_irk_dae_model model = {NULL, scalar_residual, scalar_jacobian};
    matmpc_irk_dae_input in = {{{&x, 1}, {&u, 1}, NULL, {&seed, 1}, {&zero, 1}, {&lambda, 1}}, {&z, 1}};
    double xn, zn, sx, su, adj[2];
    matmpc_irk_dae_output out = {{{&xn, 1}, {&sx, 1}, {&su, 1}, {adj, 2}}, {&zn, 1}};
    matmpc_irk_workspace workspace = allocate_workspace(&config.shape, 1);
    assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == MATMPC_IRK_OK);
    near(xn, 1.5);
    near(zn, -2.0);
    near(sx, 0.0);
    near(su, 0.5);
    near(adj[0], 0.0);
    near(adj[1], 1.0);
    free_workspace(workspace);

    /* A midpoint tableau is not stiffly accurate. Algebraic output is internal. */
    const double midpoint = 0.5;
    config.a.data = &midpoint;
    config.shape.nu = 0;
    config.shape.forward = false;
    config.shape.adjoint = false;
    in.differential.u = (matmpc_erk_view){NULL, 0};
    out.differential.Su = (matmpc_erk_buffer){NULL, 0};
    workspace = allocate_workspace(&config.shape, 1);
    assert(matmpc_irk_dae_integrate(&config, 1, &model, &in, &out, workspace) == MATMPC_IRK_OK);
    near(xn, 1.2);
    near(zn, -3.2);
    free_workspace(workspace);
}

int main(void)
{
    test_coupled_radau_and_failure_reuse();
    test_weights_and_algebraic_stage();
    matmpc_erk_shape overflow = {1, 1, SIZE_MAX, 1, true, true};
    size_t doubles, pivots;
    assert(matmpc_irk_dae_workspace_size(&overflow, 1, &doubles, &pivots) == MATMPC_IRK_SIZE_OVERFLOW);
    assert(matmpc_irk_dae_workspace_size(&overflow, 0, &doubles, &pivots) == MATMPC_IRK_INVALID_ARGUMENT);
    puts("Pure C IRK DAE tests passed.");
    return 0;
}
