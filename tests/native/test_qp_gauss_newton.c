#include "qp_gauss_newton.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void require(bool condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void check_array(const double *actual, const double *expected, size_t count, const char *name)
{
    for (size_t i = 0; i < count; ++i) {
        if (!(fabs(actual[i] - expected[i]) < 1e-12)) {
            fprintf(stderr, "%s[%zu]: expected %.16g but got %.16g\n", name, i, expected[i], actual[i]);
            exit(EXIT_FAILURE);
        }
    }
}

static int rhs(void *context, const double *x, const double *u, const double *p, double *rate)
{
    (void)context;
    (void)x;
    rate[0] = u[0] + p[0];
    rate[1] = 2 * u[1] - p[0];
    return 0;
}

static int forward(void *context, const double *x, const double *u, const double *p,
                   const double *Sx, const double *Su, double *dSx, double *dSu)
{
    (void)context;
    (void)x;
    (void)u;
    (void)p;
    (void)Sx;
    (void)Su;
    for (size_t i = 0; i < 4; ++i) {
        dSx[i] = 0.0;
        dSu[i] = 0.0;
    }
    dSu[0] = 1.0;
    dSu[3] = 2.0;
    return 0;
}

static int stage_cost(void *context, const double *x, const double *u, const double *p,
                      const double *y, const double *W, double *Jx, double *Ju, double *gx, double *gu)
{
    (void)context;
    const double r0 = x[0] + 2 * u[0] + p[0] - y[0];
    const double r1 = 3 * x[1] - u[1] - y[1];
    Jx[0] = sqrt(W[0]);
    Jx[1] = Jx[2] = 0.0;
    Jx[3] = 3 * sqrt(W[1]);
    Ju[0] = 2 * sqrt(W[0]);
    Ju[1] = Ju[2] = 0.0;
    Ju[3] = -sqrt(W[1]);
    gx[0] = W[0] * r0;
    gx[1] = 3 * W[1] * r1;
    gu[0] = 2 * W[0] * r0;
    gu[1] = -W[1] * r1;
    return 0;
}

static int terminal_cost(void *context, const double *x, const double *p, const double *y, const double *W,
                         double *J, double *gx)
{
    (void)context;
    J[0] = J[1] = sqrt(W[0]);
    gx[0] = gx[1] = W[0] * (x[0] + x[1] - p[0] - y[0]);
    return 0;
}

static int path_constraint(void *context, const double *x, const double *u, const double *p,
                           double *value, double *Cx, double *Cu)
{
    (void)context;
    value[0] = x[0] + u[1] + p[0];
    value[1] = x[1] - u[0];
    Cx[0] = Cx[3] = 1.0;
    Cx[1] = Cx[2] = 0.0;
    Cu[0] = Cu[3] = 0.0;
    Cu[1] = -1.0;
    Cu[2] = 1.0;
    return 0;
}

static int terminal_constraint(void *context, const double *x, const double *p, double *value, double *C)
{
    (void)context;
    value[0] = 2 * x[0] - x[1] + p[0];
    C[0] = 2.0;
    C[1] = -1.0;
    return 0;
}

int main(void)
{
    const double tableau[] = {0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0};
    const double weights[] = {1.0 / 6, 1.0 / 3, 1.0 / 3, 1.0 / 6};
    const double state_seed[] = {1, 0, 0, 1}, control_seed[] = {0, 0, 0, 0};
    const size_t indices[] = {1, 0};
    matmpc_qp_config config = {
        {{2, 2, 4, 2, true, false}, 0.125, {tableau, 16}, {weights, 4}},
        2, 1, 2, 1, 2, 1, indices, 2, {state_seed, 4}, {control_seed, 4}, 2.0
    };
    const matmpc_qp_model model = {
        {NULL, rhs, forward, NULL}, NULL, stage_cost, terminal_cost, path_constraint, terminal_constraint
    };
    const double x0[] = {1.2, 1.9}, x[] = {1, 2, 1.5, 1.5, 2, 1}, u[] = {2, -1, 0.5, 0.75}, p[] = {4, -2, 3};
    const double y[] = {1, 0, 2, -1}, W[] = {4, 9, 1, 4}, yN[] = {-2}, WN[] = {0.5};
    const double lbu[] = {-5, -10, -5, -10}, ubu[] = {5, 10, 5, 10};
    const double lbx[] = {-3, -7, -2, -6}, ubx[] = {3, 7, 2, 6};
    const double lb[] = {-1, -2, -3, -4, -5}, ub[] = {1, 2, 3, 4, 5};
    matmpc_qp_input input = {
        {x0, 2}, {x, 6}, {u, 4}, {p, 3}, {y, 4}, {yN, 1}, {W, 4}, {WN, 1},
        {lbu, 4}, {ubu, 4}, {lbx, 4}, {ubx, 4}, {lb, 5}, {ub, 5}
    };
    double A[8], B[8], a[4], ds0[2], Q[12], S[8], R[8], gx[6], gu[4];
    double Cgx[8], Cgu[8], CgN[2], lc[5], uc[5], lb_du[4], ub_du[4], lb_dx[4], ub_dx[4];
    matmpc_qp_output output = {
        {A, 8}, {B, 8}, {a, 4}, {ds0, 2}, {Q, 12}, {S, 8}, {R, 8}, {gx, 6}, {gu, 4},
        {Cgx, 8}, {Cgu, 8}, {CgN, 2}, {lc, 5}, {uc, 5}, {lb_du, 4}, {ub_du, 4}, {lb_dx, 4}, {ub_dx, 4}
    };
    size_t count;
    double scratch[128];
    require(matmpc_qp_workspace_size(&config, &count) == MATMPC_QP_OK && count <= 128, "workspace query failed");
    const matmpc_erk_buffer workspace = {scratch, count};
    require(matmpc_qp_gauss_newton(&config, &model, &input, &output, workspace) == MATMPC_QP_OK,
            "QP assembly failed");

    const double expected_A[] = {1, 0, 0, 1, 1, 0, 0, 1};
    const double expected_B[] = {0.25, 0, 0, 0.5, 0.25, 0, 0, 0.5};
    const double expected_a[] = {1, -1, -0.875, 1.375}, expected_ds0[] = {0.2, -0.1};
    const double expected_Q[] = {4, 0, 0, 81, 2, 0, 0, 36, 2, 0.5, 0.5, 2};
    const double expected_S[] = {8, 0, 0, -27, 2, 0, 0, -12};
    const double expected_R[] = {16, 0, 0, 9, 4, 0, 0, 4};
    const double expected_gx[] = {32, 189, -1.5, 57, 1, 1}, expected_gu[] = {64, -63, -3, -19};
    const double expected_Cgu[] = {0, -1, 1, 0, 0, -1, 1, 0}, expected_CgN[] = {2, -1};
    const double expected_lc[] = {-5, -2, -3.25, -5, -11}, expected_uc[] = {-3, 2, 2.75, 3, -1};
    const double expected_lb_du[] = {-7, -9, -5.5, -10.75}, expected_ub_du[] = {3, 11, 4.5, 9.25};
    const double expected_lb_dx[] = {-4.5, -8.5, -3, -8}, expected_ub_dx[] = {1.5, 5.5, 1, 4};
    check_array(A, expected_A, 8, "A");
    check_array(B, expected_B, 8, "B");
    check_array(a, expected_a, 4, "a");
    check_array(ds0, expected_ds0, 2, "ds0");
    check_array(Q, expected_Q, 12, "Q");
    check_array(S, expected_S, 8, "S");
    check_array(R, expected_R, 8, "R");
    check_array(gx, expected_gx, 6, "gx");
    check_array(gu, expected_gu, 4, "gu");
    check_array(Cgx, expected_A, 8, "Cgx");
    check_array(Cgu, expected_Cgu, 8, "Cgu");
    check_array(CgN, expected_CgN, 2, "CgN");
    check_array(lc, expected_lc, 5, "lc");
    check_array(uc, expected_uc, 5, "uc");
    check_array(lb_du, expected_lb_du, 4, "lb_du");
    check_array(ub_du, expected_ub_du, 4, "ub_du");
    check_array(lb_dx, expected_lb_dx, 4, "lb_dx");
    check_array(ub_dx, expected_ub_dx, 4, "ub_dx");

    const size_t invalid_indices[] = {1, 2};
    config.bounded_states = invalid_indices;
    require(matmpc_qp_gauss_newton(&config, &model, &input, &output, workspace) == MATMPC_QP_INVALID_ARGUMENT,
            "unsafe state index was accepted");
    config.bounded_states = indices;
    const matmpc_erk_buffer short_workspace = {scratch, count - 1};
    require(matmpc_qp_gauss_newton(&config, &model, &input, &output, short_workspace) == MATMPC_QP_WORKSPACE_TOO_SMALL,
            "short workspace was accepted");

    config.nc = config.nc_terminal = config.nbx = 0;
    const matmpc_erk_view empty_input = {NULL, 0};
    const matmpc_erk_buffer empty_output = {NULL, 0};
    input.lower_path = input.upper_path = input.lower_x = input.upper_x = empty_input;
    output.Cgx = output.Cgu = output.CgN = output.lc = output.uc = empty_output;
    output.lb_dx = output.ub_dx = empty_output;
    require(matmpc_qp_gauss_newton(&config, &model, &input, &output, workspace) == MATMPC_QP_OK,
            "zero-row constraints were rejected");

    config.horizon = SIZE_MAX;
    require(matmpc_qp_workspace_size(&config, &count) == MATMPC_QP_SIZE_OVERFLOW && count == 0,
            "overflowing horizon was accepted");
    return EXIT_SUCCESS;
}
