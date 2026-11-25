#include "recover.h"

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

matmpc_recover_status matmpc_recover_get_sizes(const matmpc_recover_dimensions *dimensions,
                                             matmpc_recover_sizes *sizes)
{
    const size_t limit = SIZE_MAX / sizeof(double);
    matmpc_recover_sizes required;
    size_t path_multipliers;

    if (dimensions == NULL || dimensions->nx == 0 || dimensions->horizon == 0 ||
        dimensions->horizon == SIZE_MAX) {
        return MATMPC_RECOVER_INVALID_DIMENSIONS;
    }
    if (sizes == NULL) {
        return MATMPC_RECOVER_INVALID_BUFFER;
    }

    if (!count_product(dimensions->horizon + 1, dimensions->nx, &required.dx) ||
        !count_product(dimensions->horizon, dimensions->nu, &required.du) ||
        !count_product(required.dx, dimensions->nx, &required.q) ||
        !count_product(required.du, dimensions->nx, &required.s) ||
        !count_product(dimensions->horizon, dimensions->nx, &required.defect) ||
        !count_product(required.defect, dimensions->nx, &required.a) ||
        !count_product(dimensions->state_bounds, dimensions->nx, &required.cx) ||
        !count_product(required.defect, dimensions->path_constraints, &required.cgx) ||
        !count_product(dimensions->terminal_constraints, dimensions->nx, &required.cgn) ||
        !count_product(dimensions->horizon, dimensions->path_constraints, &path_multipliers) ||
        !count_product(dimensions->horizon, dimensions->state_bounds, &required.mu_x) ||
        dimensions->terminal_constraints > limit - path_multipliers) {
        return MATMPC_RECOVER_INVALID_DIMENSIONS;
    }

    required.mu = path_multipliers + dimensions->terminal_constraints;
    *sizes = required;
    return MATMPC_RECOVER_OK;
}

static int readable(matmpc_recover_buffer buffer, size_t count)
{
    return count == 0 || (buffer.data != NULL && buffer.count >= count);
}

static void add_product(double *result, const double *matrix, const double *vector,
                        size_t rows, size_t columns)
{
    size_t column, row;

    for (column = 0; column < columns; ++column) {
        const double *entries = matrix + column * rows;
        const double weight = vector[column];
        for (row = 0; row < rows; ++row) {
            result[row] += entries[row] * weight;
        }
    }
}

static void add_transpose_product(double *result, const double *matrix, const double *vector,
                                  size_t rows, size_t columns)
{
    size_t column, row;

    for (column = 0; column < columns; ++column) {
        const double *entries = matrix + column * rows;
        double value = 0.0;
        for (row = 0; row < rows; ++row) {
            value += entries[row] * vector[row];
        }
        result[column] += value;
    }
}

matmpc_recover_status matmpc_recover(const matmpc_recover_dimensions *dimensions,
                                   const matmpc_recover_inputs *inputs,
                                   double *dx, size_t dx_count,
                                   double *lambda, size_t lambda_count)
{
    matmpc_recover_sizes sizes;
    matmpc_recover_status status = matmpc_recover_get_sizes(dimensions, &sizes);
    size_t nx, nu, horizon, nbx, nc, ncN, state_bytes, stage, row;
    double *costate;

    if (status != MATMPC_RECOVER_OK) {
        return status;
    }

    nx = dimensions->nx;
    nu = dimensions->nu;
    horizon = dimensions->horizon;
    nbx = dimensions->state_bounds;
    nc = dimensions->path_constraints;
    ncN = dimensions->terminal_constraints;
    if (inputs == NULL || dx == NULL || lambda == NULL ||
        dx_count < sizes.dx || lambda_count < sizes.dx ||
        !readable(inputs->q, sizes.q) || !readable(inputs->s, sizes.s) ||
        !readable(inputs->a, sizes.a) || !readable(inputs->b, sizes.s) ||
        !readable(inputs->cx, sizes.cx) || !readable(inputs->cgx, sizes.cgx) ||
        !readable(inputs->cgn, sizes.cgn) || !readable(inputs->gx, sizes.dx) ||
        !readable(inputs->defect, sizes.defect) || !readable(inputs->ds0, nx) ||
        !readable(inputs->du, sizes.du) || !readable(inputs->mu, sizes.mu) ||
        !readable(inputs->mu_x, sizes.mu_x)) {
        return MATMPC_RECOVER_INVALID_BUFFER;
    }

    state_bytes = nx * sizeof(double);
    memcpy(dx, inputs->ds0.data, state_bytes);
    for (stage = 0; stage < horizon; ++stage) {
        const size_t state = stage * nx;
        double *next = dx + state + nx;

        memcpy(next, inputs->defect.data + state, state_bytes);
        add_product(next, inputs->a.data + state * nx, dx + state, nx, nx);
        if (nu > 0) {
            add_product(next, inputs->b.data + state * nu, inputs->du.data + stage * nu, nx, nu);
        }
    }

    costate = lambda + horizon * nx;
    memcpy(costate, inputs->gx.data + horizon * nx, state_bytes);
    add_product(costate, inputs->q.data + horizon * nx * nx, dx + horizon * nx, nx, nx);
    if (ncN > 0) {
        add_transpose_product(costate, inputs->cgn.data, inputs->mu.data + horizon * nc, ncN, nx);
    }
    if (nbx > 0) {
        add_transpose_product(costate, inputs->cx.data, inputs->mu_x.data + (horizon - 1) * nbx, nbx, nx);
    }

    for (stage = horizon; stage > 0;) {
        size_t state;

        --stage;
        state = stage * nx;
        costate = lambda + state;
        memcpy(costate, inputs->gx.data + state, state_bytes);
        add_product(costate, inputs->q.data + state * nx, dx + state, nx, nx);
        if (nu > 0) {
            add_product(costate, inputs->s.data + state * nu, inputs->du.data + stage * nu, nx, nu);
        }
        add_transpose_product(costate, inputs->a.data + state * nx, costate + nx, nx, nx);
        if (nc > 0) {
            add_transpose_product(costate, inputs->cgx.data + state * nc, inputs->mu.data + stage * nc, nc, nx);
        }
        if (stage > 0 && nbx > 0) {
            add_transpose_product(costate, inputs->cx.data, inputs->mu_x.data + (stage - 1) * nbx, nbx, nx);
        }
    }

    /* The initial-condition equality has the opposite sign of the dynamics. */
    for (row = 0; row < nx; ++row) {
        lambda[row] = -lambda[row];
    }
    return MATMPC_RECOVER_OK;
}
