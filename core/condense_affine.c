#include "condense_affine.h"
#include "matmpc_la.h"

#include <stdint.h>
#include <string.h>

static int count_product(size_t left, size_t right, size_t *result)
{
    if (left != 0 && right > (SIZE_MAX / sizeof(double)) / left) {
        return 0;
    }
    *result = left * right;
    return 1;
}

matmpc_condense_status matmpc_condense_affine_get_sizes(const matmpc_condense_dimensions *dimensions,
                                                       matmpc_condense_affine_sizes *sizes)
{
    const size_t limit = SIZE_MAX / sizeof(double);
    matmpc_condense_affine_sizes required;
    size_t path;

    if (dimensions == NULL || dimensions->nx == 0 || dimensions->horizon == 0 ||
        dimensions->horizon == SIZE_MAX) {
        return MATMPC_CONDENSE_INVALID_DIMENSIONS;
    }
    if (sizes == NULL) {
        return MATMPC_CONDENSE_INVALID_BUFFER;
    }

    if (!count_product(dimensions->horizon + 1, dimensions->nx, &required.states) ||
        !count_product(dimensions->horizon, dimensions->nu, &required.controls) ||
        !count_product(required.states, dimensions->nx, &required.q) ||
        !count_product(required.controls, dimensions->nx, &required.s) ||
        !count_product(dimensions->horizon, dimensions->nx, &required.defect) ||
        !count_product(required.defect, dimensions->nx, &required.a) ||
        !count_product(dimensions->state_bounds, dimensions->nx, &required.cx) ||
        !count_product(required.defect, dimensions->path_constraints, &required.cgx) ||
        !count_product(dimensions->terminal_constraints, dimensions->nx, &required.cgn) ||
        !count_product(dimensions->horizon, dimensions->path_constraints, &path) ||
        !count_product(dimensions->horizon, dimensions->state_bounds, &required.state) ||
        dimensions->terminal_constraints > limit - path ||
        required.defect > limit - required.states) {
        return MATMPC_CONDENSE_INVALID_DIMENSIONS;
    }

    required.general = path + dimensions->terminal_constraints;
    required.workspace = required.states + required.defect;
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

static void shift_bounds(const double *matrix, size_t rows, size_t nx, const double *state,
                         const double *lower, const double *upper,
                         double *shifted_lower, double *shifted_upper)
{
    size_t i;

    matmpc_la_gemv(0, rows, nx, -1.0, matrix, rows, state, 0.0, shifted_lower);
    for (i = 0; i < rows; ++i) {
        shifted_upper[i] = shifted_lower[i] + upper[i];
        shifted_lower[i] += lower[i];
    }
}

matmpc_condense_status matmpc_condense_affine(const matmpc_condense_dimensions *dimensions,
                                             const matmpc_condense_affine_inputs *inputs,
                                             const matmpc_condense_affine_outputs *outputs,
                                             matmpc_condense_output workspace)
{
    matmpc_condense_affine_sizes sizes;
    matmpc_condense_status status = matmpc_condense_affine_get_sizes(dimensions, &sizes);
    size_t nx, nu, horizon, nbx, nc, ncN, stage, state_bytes;
    double *trajectory, *costates;

    if (status != MATMPC_CONDENSE_OK) {
        return status;
    }

    nx = dimensions->nx;
    nu = dimensions->nu;
    horizon = dimensions->horizon;
    nbx = dimensions->state_bounds;
    nc = dimensions->path_constraints;
    ncN = dimensions->terminal_constraints;
    if (inputs == NULL || outputs == NULL ||
        !readable(inputs->q, sizes.q) || !readable(inputs->s, sizes.s) ||
        !readable(inputs->a, sizes.a) || !readable(inputs->b, sizes.s) ||
        !readable(inputs->cx, sizes.cx) || !readable(inputs->cgx, sizes.cgx) ||
        !readable(inputs->cgn, sizes.cgn) || !readable(inputs->gx, sizes.states) ||
        !readable(inputs->gu, sizes.controls) || !readable(inputs->defect, sizes.defect) ||
        !readable(inputs->ds0, nx) || !readable(inputs->lower_general, sizes.general) ||
        !readable(inputs->upper_general, sizes.general) || !readable(inputs->lower_state, sizes.state) ||
        !readable(inputs->upper_state, sizes.state) || !writable(outputs->gc, sizes.controls) ||
        !writable(outputs->lcc, sizes.general) || !writable(outputs->ucc, sizes.general) ||
        !writable(outputs->lxc, sizes.state) || !writable(outputs->uxc, sizes.state) ||
        !writable(workspace, sizes.workspace)) {
        return MATMPC_CONDENSE_INVALID_BUFFER;
    }

    trajectory = workspace.data;
    costates = trajectory + sizes.states;
    state_bytes = nx * sizeof(double);
    memcpy(trajectory, inputs->ds0.data, state_bytes);
    for (stage = 0; stage < horizon; ++stage) {
        const size_t offset = stage * nx;
        double *next = trajectory + offset + nx;
        memcpy(next, inputs->defect.data + offset, state_bytes);
        matmpc_la_gemv(0, nx, nx, 1.0, inputs->a.data + offset * nx, nx,
                      trajectory + offset, 1.0, next);
    }

    if (nu > 0) {
        double *terminal_costate = costates + (horizon - 1) * nx;
        memcpy(terminal_costate, inputs->gx.data + horizon * nx, state_bytes);
        matmpc_la_symv_lower(nx, 1.0, inputs->q.data + horizon * nx * nx, nx,
                            trajectory + horizon * nx, 1.0, terminal_costate);

        for (stage = horizon; stage > 0;) {
            size_t offset;
            const double *next_costate;
            double *gradient;

            --stage;
            offset = stage * nx;
            next_costate = costates + offset;
            gradient = outputs->gc.data + stage * nu;
            memcpy(gradient, inputs->gu.data + stage * nu, nu * sizeof(double));
            matmpc_la_gemv(1, nx, nu, 1.0, inputs->s.data + offset * nu, nx,
                          trajectory + offset, 1.0, gradient);
            matmpc_la_gemv(1, nx, nu, 1.0, inputs->b.data + offset * nu, nx,
                          next_costate, 1.0, gradient);

            if (stage > 0) {
                double *costate = costates + offset - nx;
                memcpy(costate, inputs->gx.data + offset, state_bytes);
                matmpc_la_symv_lower(nx, 1.0, inputs->q.data + offset * nx, nx,
                                    trajectory + offset, 1.0, costate);
                matmpc_la_gemv(1, nx, nx, 1.0, inputs->a.data + offset * nx, nx,
                              next_costate, 1.0, costate);
            }
        }
    }

    for (stage = 0; nc > 0 && stage < horizon; ++stage) {
        const size_t row = stage * nc;
        shift_bounds(inputs->cgx.data + row * nx, nc, nx, trajectory + stage * nx,
                     inputs->lower_general.data + row, inputs->upper_general.data + row,
                     outputs->lcc.data + row, outputs->ucc.data + row);
    }
    if (ncN > 0) {
        const size_t row = horizon * nc;
        shift_bounds(inputs->cgn.data, ncN, nx, trajectory + horizon * nx,
                     inputs->lower_general.data + row, inputs->upper_general.data + row,
                     outputs->lcc.data + row, outputs->ucc.data + row);
    }
    for (stage = 0; nbx > 0 && stage < horizon; ++stage) {
        const size_t row = stage * nbx;
        shift_bounds(inputs->cx.data, nbx, nx, trajectory + (stage + 1) * nx,
                     inputs->lower_state.data + row, inputs->upper_state.data + row,
                     outputs->lxc.data + row, outputs->uxc.data + row);
    }

    return MATMPC_CONDENSE_OK;
}
