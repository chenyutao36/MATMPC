#include "condense_mb.h"
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

matmpc_condense_status matmpc_condense_mb_get_sizes(const matmpc_condense_dimensions *dimensions,
                                                   size_t block_count,
                                                   matmpc_condense_mb_sizes *sizes)
{
    const size_t limit = SIZE_MAX / sizeof(double);
    matmpc_condense_mb_sizes required;
    size_t adjoints, matrix_workspace;
    matmpc_condense_status status = matmpc_condense_affine_get_sizes(dimensions, &required.affine);

    if (status != MATMPC_CONDENSE_OK) {
        return status;
    }
    if (dimensions->nu == 0 || block_count == 0 || block_count > dimensions->horizon) {
        return MATMPC_CONDENSE_INVALID_DIMENSIONS;
    }
    if (sizes == NULL) {
        return MATMPC_CONDENSE_INVALID_BUFFER;
    }

    if (!count_product(block_count, dimensions->nu, &required.controls) ||
        !count_product(required.controls, dimensions->nu, &required.r_block) ||
        !count_product(required.affine.general - dimensions->terminal_constraints,
                       dimensions->nu, &required.cgu) ||
        !count_product(required.controls, required.controls, &required.hc) ||
        !count_product(required.affine.state, required.controls, &required.ccx) ||
        !count_product(required.affine.general, required.controls, &required.ccg) ||
        !count_product(2, required.affine.s / dimensions->horizon, &adjoints) ||
        adjoints > limit - required.affine.s) {
        return MATMPC_CONDENSE_INVALID_DIMENSIONS;
    }

    matrix_workspace = required.affine.s + adjoints;
    required.workspace = matrix_workspace > required.affine.workspace ? matrix_workspace : required.affine.workspace;
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

static void copy_block(double *target, size_t stride, const double *source, size_t rows, size_t columns)
{
    for (size_t col = 0; col < columns; ++col) {
        memcpy(target + col * stride, source + col * rows, rows * sizeof(double));
    }
}

static void form_responses(const matmpc_condense_dimensions *dimensions, const matmpc_condense_mb_sizes *sizes,
                           const matmpc_condense_mb_inputs *inputs, const matmpc_condense_matrix_outputs *outputs,
                           size_t start, size_t end, size_t column, double *response)
{
    const size_t nx = dimensions->nx, nu = dimensions->nu, horizon = dimensions->horizon;
    const size_t nbx = dimensions->state_bounds, nc = dimensions->path_constraints;
    const size_t ncN = dimensions->terminal_constraints, stride = nx * nu;
    const matmpc_condense_affine_inputs *affine = &inputs->affine;

    memcpy(response + start * stride, affine->b.data + start * stride, stride * sizeof(double));
    for (size_t stage = start; stage < horizon; ++stage) {
        double *next = response + stage * stride;

        if (stage > start) {
            const double *previous = next - stride;
            matmpc_la_gemm(0, 0, nx, nu, nx, 1.0, affine->a.data + stage * nx * nx, nx,
                          previous, nx, 0.0, next, nx);
            if (stage < end) {
                const double *b = affine->b.data + stage * stride;
                for (size_t i = 0; i < stride; ++i) {
                    next[i] += b[i];
                }
            }
        }

        if (nbx > 0) {
            double *bounds = outputs->ccx.data + column * sizes->affine.state + stage * nbx;
            matmpc_la_gemm(0, 0, nbx, nu, nx, 1.0, affine->cx.data, nbx,
                          next, nx, 0.0, bounds, sizes->affine.state);
        }
        if (nc > 0) {
            double *path = outputs->ccg.data + column * sizes->affine.general + stage * nc;
            if (stage < end) {
                copy_block(path, sizes->affine.general, inputs->cgu.data + stage * nc * nu, nc, nu);
            }
            if (stage > start) {
                matmpc_la_gemm(0, 0, nc, nu, nx, 1.0, affine->cgx.data + stage * nc * nx, nc,
                              next - stride, nx, stage < end ? 1.0 : 0.0, path, sizes->affine.general);
            }
        }
    }
    if (ncN > 0) {
        double *terminal = outputs->ccg.data + column * sizes->affine.general + horizon * nc;
        matmpc_la_gemm(0, 0, ncN, nu, nx, 1.0, affine->cgn.data, ncN,
                      response + (horizon - 1) * stride, nx, 0.0, terminal, sizes->affine.general);
    }
}

static void form_hessian_column(const matmpc_condense_dimensions *dimensions, const matmpc_condense_mb_sizes *sizes,
                                const matmpc_condense_mb_inputs *inputs, const size_t *boundaries, size_t block_count,
                                size_t block, double *hessian, double *response)
{
    const size_t nx = dimensions->nx, nu = dimensions->nu, horizon = dimensions->horizon;
    const size_t stride = nx * nu, start = boundaries[block], end = boundaries[block + 1];
    const size_t column = block * nu, ld = sizes->controls;
    const matmpc_condense_affine_inputs *affine = &inputs->affine;
    double *costate = response + sizes->affine.s;
    double *previous_costate = costate + stride;
    size_t row_block = block_count - 1;

    copy_block(hessian + column * ld + column, ld, inputs->r_block.data + column * nu, nu, nu);
    matmpc_la_symm_left_lower(nx, nu, 1.0, affine->q.data + horizon * nx * nx, nx,
                              response + (horizon - 1) * stride, nx, 0.0, costate, nx);

    for (size_t stage = horizon; stage > start;) {
        --stage;
        if (stage < boundaries[row_block]) {
            --row_block;
        }
        double *entry = hessian + column * ld + row_block * nu;
        if (stage > start) {
            const double *state = response + (stage - 1) * stride;
            const double *mixed = affine->s.data + stage * stride;
            matmpc_la_gemm(1, 0, nu, nu, nx, 1.0, mixed, nx, state, nx, 1.0, entry, ld);
        }
        matmpc_la_gemm(1, 0, nu, nu, nx, 1.0, affine->b.data + stage * stride, nx,
                      costate, nx, 1.0, entry, ld);

        if (stage > start) {
            const double *state = response + (stage - 1) * stride;
            matmpc_la_gemm(1, 0, nx, nu, nx, 1.0, affine->a.data + stage * nx * nx, nx,
                          costate, nx, 0.0, previous_costate, nx);
            matmpc_la_symm_left_lower(nx, nu, 1.0, affine->q.data + stage * nx * nx, nx,
                                      state, nx, 1.0, previous_costate, nx);
            if (stage < end) {
                /* The shared input also enters the state gradient through S at this stage. */
                const double *mixed = affine->s.data + stage * stride;
                for (size_t i = 0; i < stride; ++i) {
                    previous_costate[i] += mixed[i];
                }
            }
            double *next = costate;
            costate = previous_costate;
            previous_costate = next;
        }
    }
}

matmpc_condense_status matmpc_condense_mb(const matmpc_condense_dimensions *dimensions,
                                         const size_t *boundaries, size_t boundary_count,
                                         const matmpc_condense_mb_inputs *inputs,
                                         const matmpc_condense_mb_outputs *outputs,
                                         matmpc_condense_output workspace)
{
    matmpc_condense_mb_sizes sizes;
    size_t block_count = boundary_count == 0 ? 0 : boundary_count - 1;
    matmpc_condense_status status = matmpc_condense_mb_get_sizes(dimensions, block_count, &sizes);

    if (status != MATMPC_CONDENSE_OK) {
        return status;
    }
    if (inputs == NULL || outputs == NULL || !writable(workspace, sizes.workspace) ||
        !readable(inputs->r_block, sizes.r_block) || !readable(inputs->cgu, sizes.cgu) ||
        !writable(outputs->matrices.hc, sizes.hc) || !writable(outputs->matrices.ccx, sizes.ccx) ||
        !writable(outputs->matrices.ccg, sizes.ccg)) {
        return MATMPC_CONDENSE_INVALID_BUFFER;
    }
    status = matmpc_condense_affine_blocked(dimensions, &inputs->affine, &outputs->affine,
                                            workspace, boundaries, boundary_count);
    if (status != MATMPC_CONDENSE_OK) {
        return status;
    }

    memset(outputs->matrices.hc.data, 0, sizes.hc * sizeof(double));
    if (sizes.ccx > 0) {
        memset(outputs->matrices.ccx.data, 0, sizes.ccx * sizeof(double));
    }
    if (sizes.ccg > 0) {
        memset(outputs->matrices.ccg.data, 0, sizes.ccg * sizeof(double));
    }
    for (size_t block = 0; block < block_count; ++block) {
        form_responses(dimensions, &sizes, inputs, &outputs->matrices,
                       boundaries[block], boundaries[block + 1], block * dimensions->nu, workspace.data);
        form_hessian_column(dimensions, &sizes, inputs, boundaries, block_count,
                            block, outputs->matrices.hc.data, workspace.data);
    }
    double *hessian = outputs->matrices.hc.data;
    for (size_t col = 0; col < sizes.controls; ++col) {
        for (size_t row = 0; row < col; ++row) {
            hessian[row + col * sizes.controls] = hessian[col + row * sizes.controls];
        }
    }

    return MATMPC_CONDENSE_OK;
}
