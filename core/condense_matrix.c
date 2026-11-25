#include "condense_matrix.h"
#include "matmpc_la.h"

#include <stdint.h>
#include <string.h>

static int count_product(size_t left, size_t right, size_t *result)
{
    const size_t limit = SIZE_MAX / sizeof(double);

    if (left != 0 && right > limit / left) {
        return 0;
    }
    *result = left * right;
    return 1;
}

static int count_sum(size_t left, size_t right, size_t *result)
{
    const size_t limit = SIZE_MAX / sizeof(double);

    if (right > limit - left) {
        return 0;
    }
    *result = left + right;
    return 1;
}

matmpc_condense_status matmpc_condense_matrix_get_sizes(const matmpc_condense_dimensions *dimensions,
                                                      matmpc_condense_matrix_sizes *sizes)
{
    matmpc_condense_matrix_sizes required;
    size_t nx, nu, N, nbx, nc, ncN;
    size_t state_rows, input_rows, path_rows, bound_rows, general_rows;
    size_t input_block, transition_count, both_transitions;

    if (dimensions == NULL || dimensions->nx == 0 || dimensions->horizon == 0 ||
        dimensions->horizon == SIZE_MAX) {
        return MATMPC_CONDENSE_INVALID_DIMENSIONS;
    }
    if (sizes == NULL) {
        return MATMPC_CONDENSE_INVALID_BUFFER;
    }

    nx = dimensions->nx;
    nu = dimensions->nu;
    N = dimensions->horizon;
    nbx = dimensions->state_bounds;
    nc = dimensions->path_constraints;
    ncN = dimensions->terminal_constraints;
    if (!count_product(N, nx, &state_rows) ||
        !count_product(N, nu, &input_rows) ||
        !count_product(N, nc, &path_rows) ||
        !count_product(N, nbx, &bound_rows) ||
        !count_sum(path_rows, ncN, &general_rows) ||
        !count_product(N + 1, nx, &required.q) ||
        !count_product(required.q, nx, &required.q) ||
        !count_product(state_rows, nx, &required.a) ||
        !count_product(input_rows, nx, &required.s) ||
        !count_product(input_rows, nu, &required.r) ||
        !count_product(nbx, nx, &required.cx) ||
        !count_product(path_rows, nx, &required.cgx) ||
        !count_product(path_rows, nu, &required.cgu) ||
        !count_product(ncN, nx, &required.cgn) ||
        !count_product(input_rows, input_rows, &required.hc) ||
        !count_product(bound_rows, input_rows, &required.ccx) ||
        !count_product(general_rows, input_rows, &required.ccg) ||
        !count_product(nu, nu, &input_block) ||
        !count_product(N, required.s, &transition_count) ||
        !count_product(2, transition_count, &both_transitions) ||
        !count_sum(both_transitions, input_block, &required.workspace)) {
        return MATMPC_CONDENSE_INVALID_DIMENSIONS;
    }

    required.b = required.s;
    *sizes = required;
    return MATMPC_CONDENSE_OK;
}

static int readable(matmpc_condense_buffer buffer, size_t count)
{
    return count == 0 || (buffer.data != NULL && buffer.count >= count);
}

static int writable(matmpc_condense_output buffer, size_t count)
{
    return count == 0 || (buffer.data != NULL && buffer.count >= count);
}

static void store_block(double *target, size_t stride, const double *block, size_t rows, size_t columns)
{
    for (size_t col = 0; col < columns; ++col) {
        memcpy(target + col * stride, block + col * rows, rows * sizeof(double));
    }
}

static void store_transpose(double *target, size_t stride, const double *block, size_t width)
{
    for (size_t col = 0; col < width; ++col) {
        for (size_t row = 0; row < width; ++row) {
            target[row + col * stride] = block[col + row * width];
        }
    }
}

static void zero_row_prefix(double *matrix, size_t stride, size_t rows, size_t columns)
{
    for (size_t col = 0; col < columns; ++col) {
        memset(matrix + col * stride, 0, rows * sizeof(double));
    }
}

matmpc_condense_status matmpc_condense_matrix(const matmpc_condense_dimensions *dimensions,
                                            const matmpc_condense_matrix_inputs *inputs,
                                            const matmpc_condense_matrix_outputs *outputs,
                                            matmpc_condense_output workspace)
{
    matmpc_condense_matrix_sizes sizes;
    matmpc_condense_status status = matmpc_condense_matrix_get_sizes(dimensions, &sizes);
    size_t nx, nu, N, nbx, nc, ncN, xx, xu, uu, inputs_per_horizon, bound_rows, general_rows;
    double *G, *W, *Hi;

    if (status != MATMPC_CONDENSE_OK) {
        return status;
    }
    if (inputs == NULL || outputs == NULL ||
        !readable(inputs->q, sizes.q) || !readable(inputs->s, sizes.s) ||
        !readable(inputs->r, sizes.r) || !readable(inputs->a, sizes.a) ||
        !readable(inputs->b, sizes.b) || !readable(inputs->cx, sizes.cx) ||
        !readable(inputs->cgx, sizes.cgx) || !readable(inputs->cgu, sizes.cgu) ||
        !readable(inputs->cgn, sizes.cgn) || !writable(outputs->hc, sizes.hc) ||
        !writable(outputs->ccx, sizes.ccx) || !writable(outputs->ccg, sizes.ccg) ||
        !writable(workspace, sizes.workspace)) {
        return MATMPC_CONDENSE_INVALID_BUFFER;
    }
    if (sizes.hc == 0) {
        return MATMPC_CONDENSE_OK;
    }

    nx = dimensions->nx;
    nu = dimensions->nu;
    N = dimensions->horizon;
    nbx = dimensions->state_bounds;
    nc = dimensions->path_constraints;
    ncN = dimensions->terminal_constraints;
    xx = nx * nx;
    xu = nx * nu;
    uu = nu * nu;
    inputs_per_horizon = nu * N;
    bound_rows = nbx * N;
    general_rows = nc * N + ncN;
    G = workspace.data;
    W = G + N * sizes.b;
    Hi = W + N * sizes.b;

    for (size_t i = 0; i < N; ++i) {
        double *g = G + i * sizes.b;
        double *w = W + i * sizes.b;
        double *hc = outputs->hc.data + i * nu * inputs_per_horizon;

        /* g[j] maps control i to state j+1. */
        memcpy(g + i * xu, inputs->b.data + i * xu, xu * sizeof(double));
        for (size_t j = i + 1; j < N; ++j) {
            matmpc_la_gemm(0, 0, nx, nu, nx, 1.0, inputs->a.data + j * xx, nx,
                           g + (j - 1) * xu, nx, 0.0, g + j * xu, nx);
        }

        matmpc_la_symm_left_lower(nx, nu, 1.0, inputs->q.data + N * xx, nx,
                                 g + (N - 1) * xu, nx, 0.0, w + (N - 1) * xu, nx);
        for (size_t j = N - 1; j > i; --j) {
            matmpc_la_gemm(1, 0, nu, nu, nx, 1.0, inputs->s.data + j * xu, nx,
                           g + (j - 1) * xu, nx, 0.0, Hi, nu);
            matmpc_la_gemm(1, 0, nu, nu, nx, 1.0, inputs->b.data + j * xu, nx,
                           w + j * xu, nx, 1.0, Hi, nu);
            store_block(hc + j * nu, inputs_per_horizon, Hi, nu, nu);
            store_transpose(outputs->hc.data + j * nu * inputs_per_horizon + i * nu,
                            inputs_per_horizon, Hi, nu);

            matmpc_la_gemm(1, 0, nx, nu, nx, 1.0, inputs->a.data + j * xx, nx,
                           w + j * xu, nx, 0.0, w + (j - 1) * xu, nx);
            matmpc_la_symm_left_lower(nx, nu, 1.0, inputs->q.data + j * xx, nx,
                                     g + (j - 1) * xu, nx, 1.0, w + (j - 1) * xu, nx);
        }
        memcpy(Hi, inputs->r.data + i * uu, uu * sizeof(double));
        matmpc_la_gemm(1, 0, nu, nu, nx, 1.0, inputs->b.data + i * xu, nx,
                       w + i * xu, nx, 1.0, Hi, nu);
        store_block(hc + i * nu, inputs_per_horizon, Hi, nu, nu);

        if (nc > 0) {
            double *ccg = outputs->ccg.data + i * nu * general_rows;
            zero_row_prefix(ccg, general_rows, i * nc, nu);
            store_block(ccg + i * nc, general_rows, inputs->cgu.data + i * nc * nu, nc, nu);
            for (size_t j = i + 1; j < N; ++j) {
                matmpc_la_gemm(0, 0, nc, nu, nx, 1.0, inputs->cgx.data + j * nc * nx, nc,
                               g + (j - 1) * xu, nx, 0.0, ccg + j * nc, general_rows);
            }
        }
        if (ncN > 0) {
            double *terminal = outputs->ccg.data + i * nu * general_rows + N * nc;
            matmpc_la_gemm(0, 0, ncN, nu, nx, 1.0, inputs->cgn.data, ncN,
                           g + (N - 1) * xu, nx, 0.0, terminal, general_rows);
        }

        if (nbx > 0) {
            double *ccx = outputs->ccx.data + i * nu * bound_rows;
            zero_row_prefix(ccx, bound_rows, i * nbx, nu);
            for (size_t j = i; j < N; ++j) {
                matmpc_la_gemm(0, 0, nbx, nu, nx, 1.0, inputs->cx.data, nbx,
                               g + j * xu, nx, 0.0, ccx + j * nbx, bound_rows);
            }
        }
    }

    return MATMPC_CONDENSE_OK;
}
