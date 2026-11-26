#include "irk_dae_kernel.h"
#include "qp_gauss_newton.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { horizon = 2 };

typedef struct {
    matmpc_irk_config config;
    matmpc_irk_workspace workspace;
    bool dae;
    size_t interval;
    size_t residual_calls;
    double steps[horizon];
    double guesses[horizon];
} implicit_fixture;

static void require(bool condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void near(double actual, double expected, const char *message)
{
    if (!(fabs(actual - expected) <= 2e-14)) {
        fprintf(stderr, "%s (got %.17g, expected %.17g)\n", message, actual, expected);
        exit(EXIT_FAILURE);
    }
}

static int ode_residual(void *context, const double *x, const double *u, const double *p,
                        const double *xdot, double *r)
{
    (void)context;
    (void)x;
    (void)p;
    r[0] = xdot[0] - u[0];
    return 0;
}

static int ode_jacobian(void *context, const double *x, const double *u, const double *p,
                        const double *xdot, double *Fx, double *Fu, double *Fxdot)
{
    (void)context;
    (void)x;
    (void)u;
    (void)p;
    (void)xdot;
    Fx[0] = 0.0;
    Fxdot[0] = 1.0;
    if (Fu != NULL)
        Fu[0] = -1.0;
    return 0;
}

static int dae_residual(void *context, const double *x, const double *u, const double *p,
                        const double *xdot, const double *z, double *f, double *g)
{
    (void)x;
    (void)p;
    implicit_fixture *fixture = context;
    if (fixture->residual_calls++ == 0)
        fixture->guesses[fixture->interval] = z[0];
    f[0] = xdot[0] - z[0];
    g[0] = z[0] - u[0];
    return 0;
}

static int dae_jacobian(void *context, const double *x, const double *u, const double *p,
                        const double *xdot, const double *z, double *Fx, double *Fu,
                        double *Fxdot, double *Fz, double *Gx, double *Gu, double *Gz)
{
    (void)context;
    (void)x;
    (void)u;
    (void)p;
    (void)xdot;
    (void)z;
    Fx[0] = 0.0;
    Fxdot[0] = 1.0;
    Fz[0] = -1.0;
    Gx[0] = 0.0;
    Gz[0] = 1.0;
    if (Fu != NULL)
        Fu[0] = 0.0;
    if (Gu != NULL)
        Gu[0] = -1.0;
    return 0;
}

static int integrate_stage(void *context, double step, const double *x, const double *u, const double *p,
                            const double *z, double *xn, double *zn, double *A, double *B)
{
    static const double identity = 1.0, zero = 0.0;
    implicit_fixture *fixture = context;
    require(fixture->interval < horizon, "QP integrated an extra shooting interval");
    fixture->steps[fixture->interval] = step;
    fixture->residual_calls = 0;
    fixture->config.step = step;
    const matmpc_erk_input input = {{x, 1}, {u, 1}, p, {&identity, 1}, {&zero, 1}, {NULL, 0}};
    const matmpc_erk_output output = {{xn, 1}, {A, 1}, {B, 1}, {NULL, 0}};
    matmpc_irk_status status;
    if (fixture->dae) {
        require(z != NULL && zn != NULL, "QP omitted algebraic storage for DAE");
        const matmpc_irk_dae_model model = {fixture, dae_residual, dae_jacobian};
        const matmpc_irk_dae_input dae_input = {input, {z, 1}};
        const matmpc_irk_dae_output dae_output = {output, {zn, 1}};
        status = matmpc_irk_dae_integrate(&fixture->config, 1, &model, &dae_input, &dae_output,
                                           fixture->workspace);
    } else {
        require(z == NULL && zn == NULL, "QP passed algebraic storage to ODE");
        const matmpc_irk_ode_model model = {NULL, ode_residual, ode_jacobian};
        status = matmpc_irk_ode_integrate(&fixture->config, &model, &input, &output, fixture->workspace);
    }
    ++fixture->interval;
    return (int)status;
}

static int stage_cost(void *context, const double *x, const double *u, const double *p,
                      const double *reference, const double *weight, double *Jx, double *Ju,
                      double *gx, double *gu)
{
    (void)context;
    (void)p;
    const double error = x[0] + 2.0 * u[0] - reference[0];
    Jx[0] = weight[0];
    Ju[0] = 2.0 * weight[0];
    gx[0] = weight[0] * Jx[0] * error;
    gu[0] = weight[0] * Ju[0] * error;
    return 0;
}

static int terminal_cost(void *context, const double *x, const double *p, const double *reference,
                         const double *weight, double *J, double *gx)
{
    (void)context;
    (void)p;
    J[0] = weight[0];
    gx[0] = J[0] * J[0] * (x[0] - reference[0]);
    return 0;
}

static void verify_nonuniform_qp(bool dae)
{
    const double root_15 = sqrt(15.0);
    const double table[] = {
        5.0 / 36.0, 5.0 / 36.0 + root_15 / 24.0, 5.0 / 36.0 + root_15 / 30.0,
        2.0 / 9.0 - root_15 / 15.0, 2.0 / 9.0, 2.0 / 9.0 + root_15 / 15.0,
        5.0 / 36.0 - root_15 / 30.0, 5.0 / 36.0 - root_15 / 24.0, 5.0 / 36.0
    };
    const double weights[] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};
    const matmpc_irk_config irk = {{1, 1, 3, 2, true, false}, 0.125, {table, 9}, {weights, 3},
                                  5, 1e-2, {NULL, 0}, 0.0};
    double irk_scratch[128];
    size_t row_scratch[16];
    implicit_fixture fixture = {0};
    fixture.config = irk;
    fixture.dae = dae;
    size_t doubles, pivots;
    const matmpc_irk_status sizing = dae ? matmpc_irk_dae_workspace_size(&irk.shape, 1, &doubles, &pivots) :
                                          matmpc_irk_ode_workspace_size(&irk.shape, &doubles, &pivots);
    require(sizing == MATMPC_IRK_OK && doubles <= sizeof(irk_scratch) / sizeof(irk_scratch[0]) &&
            pivots <= sizeof(row_scratch) / sizeof(row_scratch[0]), "IRK workspace query for QP assembly failed");
    fixture.workspace = (matmpc_irk_workspace){{irk_scratch, doubles}, row_scratch, pivots};

    const double grid[] = {0.0, 1.0, 4.0};
    const size_t bounded_state = 0;
    const matmpc_qp_config config = {
        .integration = {{1, 1, 0, 0, false, false}, irk.step, {NULL, 0}, {NULL, 0}},
        .horizon = horizon, .ny = 1, .ny_terminal = 1,
        .bounded_states = &bounded_state, .nbx = 1,
        .stage_grid = {grid, horizon + 1},
        .integrate_stage = integrate_stage, .integration_context = &fixture, .nz = dae ? 1u : 0u
    };
    const matmpc_qp_model model = {.stage_cost = stage_cost, .terminal_cost = terminal_cost};
    const double x0 = 1.75, x[] = {1.5, -0.25, 0.5}, u[] = {2.0, -1.5}, z[] = {4.0, -4.0};
    const double y[] = {0.0, 0.0}, yN = 0.25, W[] = {1.0, 1.0}, WN = 1.0;
    const double lbu[] = {-5.0, -4.0}, ubu[] = {5.0, 4.0};
    const double lbx[] = {-2.0, -3.0}, ubx[] = {2.0, 3.0};
    const matmpc_qp_input input = {
        .x0 = {&x0, 1}, .x = {x, horizon + 1}, .u = {u, horizon}, .z = {dae ? z : NULL, dae ? horizon : 0},
        .reference = {y, horizon}, .reference_terminal = {&yN, 1}, .weight = {W, horizon}, .weight_terminal = {&WN, 1},
        .lower_u = {lbu, horizon}, .upper_u = {ubu, horizon}, .lower_x = {lbx, horizon}, .upper_x = {ubx, horizon}
    };
    double A[horizon], B[horizon], a[horizon], ds0;
    double Q[horizon + 1], S[horizon], R[horizon], gx[horizon + 1], gu[horizon];
    double lb_du[horizon], ub_du[horizon], lb_dx[horizon], ub_dx[horizon], z_out[horizon];
    const matmpc_qp_output output = {
        .A = {A, horizon}, .B = {B, horizon}, .a = {a, horizon}, .ds0 = {&ds0, 1},
        .Q = {Q, horizon + 1}, .S = {S, horizon}, .R = {R, horizon}, .gx = {gx, horizon + 1}, .gu = {gu, horizon},
        .lb_du = {lb_du, horizon}, .ub_du = {ub_du, horizon}, .lb_dx = {lb_dx, horizon}, .ub_dx = {ub_dx, horizon},
        .z_out = {dae ? z_out : NULL, dae ? horizon : 0}
    };
    size_t count;
    double qp_scratch[4];
    require(matmpc_qp_workspace_size(&config, &count) == MATMPC_QP_OK &&
            count <= sizeof(qp_scratch) / sizeof(qp_scratch[0]), "QP workspace query with IRK failed");
    const matmpc_erk_buffer workspace = {qp_scratch, count};
    require(matmpc_qp_gauss_newton(&config, &model, &input, &output, workspace) == MATMPC_QP_OK,
            dae ? "QP assembly with IRK DAE failed" : "QP assembly with IRK ODE failed");
    require(fixture.interval == horizon, "QP omitted a shooting interval");

    const double interval[] = {0.25, 0.75};
    for (size_t i = 0; i < horizon; ++i) {
        near(fixture.steps[i], interval[i] / (double)irk.shape.steps, "incorrect nonuniform IRK microstep");
        near(A[i], 1.0, "incorrect state Jacobian");
        near(B[i], interval[i], "incorrect control Jacobian");
        near(a[i], x[i] + interval[i] * u[i] - x[i + 1], "incorrect shifted shooting defect");
        near(lb_dx[i], lbx[i] - x[i + 1], "incorrect stage state lower bound");
        near(ub_dx[i], ubx[i] - x[i + 1], "incorrect stage state upper bound");
        near(Q[i], 1.0, "incorrect stage Q");
        near(S[i], 2.0, "incorrect stage S");
        near(R[i], 4.0, "incorrect stage R");
        if (dae) {
            near(fixture.guesses[i], z[i], "incorrect DAE initial guess");
            near(z_out[i], u[i], "incorrect terminal algebraic state in QP");
        }
    }
    near(ds0, x0 - x[0], "incorrect initial state shift");
    near(Q[horizon], 1.0, "incorrect terminal Q");
    near(gx[horizon], x[horizon] - yN, "incorrect terminal gradient");
}

int main(void)
{
    verify_nonuniform_qp(false);
    verify_nonuniform_qp(true);
    return EXIT_SUCCESS;
}
