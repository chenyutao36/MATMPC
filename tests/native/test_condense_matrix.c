#include "condense_matrix.h"

#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static double lower_symmetric_entry(const double *matrix, size_t n, size_t row, size_t col)
{
    return row >= col ? matrix[row + col * n] : matrix[col + row * n];
}

/* Independently eliminate the dynamics by propagating each unit control. */
static void reference(const matmpc_condense_dimensions *d, const matmpc_condense_matrix_inputs *in,
                      const matmpc_condense_matrix_outputs *out)
{
    const size_t nx = d->nx, nu = d->nu, N = d->horizon;
    const size_t nbx = d->state_bounds, nc = d->path_constraints, ncN = d->terminal_constraints;
    const size_t states = (N + 1) * nx, controls = N * nu;
    const size_t bound_rows = N * nbx, general_rows = N * nc + ncN;
    double *trajectory = calloc(states * controls, sizeof(double));
    assert(controls == 0 || trajectory != NULL);

    for (size_t col = 0; col < controls; ++col) {
        double *path = trajectory + col * states;
        for (size_t stage = 0; stage < N; ++stage) {
            for (size_t row = 0; row < nx; ++row) {
                double value = 0.0;
                for (size_t k = 0; k < nx; ++k) {
                    value += in->a.data[stage * nx * nx + row + k * nx] * path[stage * nx + k];
                }
                if (stage == col / nu) {
                    value += in->b.data[stage * nx * nu + row + (col % nu) * nx];
                }
                path[(stage + 1) * nx + row] = value;
            }
            for (size_t row = 0; row < nbx; ++row) {
                double value = 0.0;
                for (size_t k = 0; k < nx; ++k) {
                    value += in->cx.data[row + k * nbx] * path[(stage + 1) * nx + k];
                }
                out->ccx.data[stage * nbx + row + col * bound_rows] = value;
            }
            for (size_t row = 0; row < nc; ++row) {
                double value = 0.0;
                for (size_t k = 0; k < nx; ++k) {
                    value += in->cgx.data[stage * nc * nx + row + k * nc] * path[stage * nx + k];
                }
                if (stage == col / nu) {
                    value += in->cgu.data[stage * nc * nu + row + (col % nu) * nc];
                }
                out->ccg.data[stage * nc + row + col * general_rows] = value;
            }
        }
        for (size_t row = 0; row < ncN; ++row) {
            double value = 0.0;
            for (size_t k = 0; k < nx; ++k) {
                value += in->cgn.data[row + k * ncN] * path[N * nx + k];
            }
            out->ccg.data[N * nc + row + col * general_rows] = value;
        }
    }

    /* Evaluate Q, S and R contributions independently of the G/W recurrence. */
    for (size_t col = 0; col < controls; ++col) {
        for (size_t row = 0; row < controls; ++row) {
            double value = 0.0;
            for (size_t stage = 0; stage <= N; ++stage) {
                const double *xL = trajectory + row * states + stage * nx;
                const double *xR = trajectory + col * states + stage * nx;
                const double *Q = in->q.data + stage * nx * nx;
                for (size_t i = 0; i < nx; ++i) {
                    for (size_t j = 0; j < nx; ++j) {
                        value += xL[i] * lower_symmetric_entry(Q, nx, i, j) * xR[j];
                    }
                    if (stage == col / nu) {
                        value += xL[i] * in->s.data[stage * nx * nu + i + (col % nu) * nx];
                    }
                    if (stage == row / nu) {
                        value += xR[i] * in->s.data[stage * nx * nu + i + (row % nu) * nx];
                    }
                }
                if (stage == row / nu && stage == col / nu) {
                    value += in->r.data[stage * nu * nu + (row % nu) + (col % nu) * nu];
                }
            }
            out->hc.data[row + col * controls] = value;
        }
    }
    free(trajectory);
}

static matmpc_condense_buffer make_input(size_t count, size_t seed)
{
    double *data = count > 0 ? malloc(count * sizeof(double)) : NULL;
    assert(count == 0 || data != NULL);
    for (size_t i = 0; i < count; ++i) {
        data[i] = ((double)((i * 3 + seed) % 19) - 9.0) / 16.0;
    }
    return (matmpc_condense_buffer){data, count};
}

static matmpc_condense_output make_output(size_t count)
{
    double *data = count > 0 ? malloc(count * sizeof(double)) : NULL;
    assert(count == 0 || data != NULL);
    for (size_t i = 0; i < count; ++i) {
        data[i] = NAN;
    }
    return (matmpc_condense_output){data, count};
}

static void check_output(matmpc_condense_output actual, matmpc_condense_output expected)
{
    assert(actual.count == expected.count);
    for (size_t i = 0; i < actual.count; ++i) {
        assert(isfinite(actual.data[i]));
        assert(fabs(actual.data[i] - expected.data[i]) <= 1e-12 * (1.0 + fabs(expected.data[i])));
    }
}

static void test_matrix(const matmpc_condense_dimensions *dimensions)
{
    matmpc_condense_matrix_sizes sizes;
    assert(matmpc_condense_matrix_get_sizes(dimensions, &sizes) == MATMPC_CONDENSE_OK);
    matmpc_condense_matrix_inputs in = {
        make_input(sizes.q, 1), make_input(sizes.s, 2), make_input(sizes.r, 3),
        make_input(sizes.a, 4), make_input(sizes.b, 5), make_input(sizes.cx, 6),
        make_input(sizes.cgx, 7), make_input(sizes.cgu, 8), make_input(sizes.cgn, 9)
    };
    matmpc_condense_matrix_outputs out = {
        make_output(sizes.hc), make_output(sizes.ccx), make_output(sizes.ccg)
    };
    matmpc_condense_matrix_outputs expected = {
        make_output(sizes.hc), make_output(sizes.ccx), make_output(sizes.ccg)
    };
    matmpc_condense_output work = make_output(sizes.workspace);

    if (work.count > 0) {
        matmpc_condense_output short_work = {work.data, work.count - 1};
        assert(matmpc_condense_matrix(dimensions, &in, &out, short_work) == MATMPC_CONDENSE_INVALID_BUFFER);
        assert(isnan(out.hc.data[0]));
    }
    for (size_t pass = 0; pass < 2; ++pass) {
        if (pass > 0 && in.b.count > 0) {
            free((void *)in.b.data);
            in.b = make_input(sizes.b, 15);
        }
        reference(dimensions, &in, &expected);
        assert(matmpc_condense_matrix(dimensions, &in, &out, work) == MATMPC_CONDENSE_OK);
        check_output(out.hc, expected.hc);
        check_output(out.ccx, expected.ccx);
        check_output(out.ccg, expected.ccg);
    }

    const matmpc_condense_buffer inputs[] = {in.q, in.s, in.r, in.a, in.b, in.cx, in.cgx, in.cgu, in.cgn};
    const matmpc_condense_output buffers[] = {out.hc, out.ccx, out.ccg, expected.hc, expected.ccx, expected.ccg, work};
    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); ++i) {
        free((void *)inputs[i].data);
    }
    for (size_t i = 0; i < sizeof(buffers) / sizeof(buffers[0]); ++i) {
        free(buffers[i].data);
    }
}

static void test_single_stage(void)
{
    const matmpc_condense_dimensions d = {1, 1, 1, 1, 1, 1};
    const double Q[] = {17, 7}, S = -4, R = 2, A = 8, B = 3, Cx = 5, Cgx = -6, Cgu = 11, CgN = 2;
    const matmpc_condense_matrix_inputs in = {
        {Q, 2}, {&S, 1}, {&R, 1}, {&A, 1}, {&B, 1}, {&Cx, 1}, {&Cgx, 1}, {&Cgu, 1}, {&CgN, 1}
    };
    double hc, ccx, ccg[2], scratch[] = {NAN, NAN, NAN};
    const matmpc_condense_matrix_outputs out = {{&hc, 1}, {&ccx, 1}, {ccg, 2}};
    assert(matmpc_condense_matrix(&d, &in, &out, (matmpc_condense_output){scratch, 3}) == MATMPC_CONDENSE_OK);
    assert(hc == 65.0 && ccx == 15.0 && ccg[0] == 11.0 && ccg[1] == 6.0);
}

int main(void)
{
    const matmpc_condense_dimensions cases[] = {
        {3, 2, 4, 2, 1, 3},
        {2, 3, 2, 1, 3, 1},
        {2, 1, 3, 0, 0, 0},
        {2, 0, 3, 1, 0, 1},
        {3, 2, 2, 0, 0, 2},
        {3, 1, 2, 2, 2, 0}
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        test_matrix(cases + i);
    }
    test_single_stage();

    matmpc_condense_matrix_sizes sizes;
    const matmpc_condense_dimensions overflow = {2, 1, SIZE_MAX - 1, 0, 0, 0};
    assert(matmpc_condense_matrix_get_sizes(&overflow, &sizes) == MATMPC_CONDENSE_INVALID_DIMENSIONS);
    assert(matmpc_condense_matrix_get_sizes(cases, NULL) == MATMPC_CONDENSE_INVALID_BUFFER);
    puts("Default-condensing matrix native tests passed.");
    return 0;
}
