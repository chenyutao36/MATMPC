#include "irk_dae_kernel.h"
#include "matmpc_la.h"
#include "matmpc_lu.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    size_t nx, nu, nz, stages, slots;
    size_t width, nk, tableau, system, state_rhs, control_rhs;
    size_t xx, xu, xz, zu, zz, adjoint;
    size_t doubles, pivots;
} dae_layout;

typedef struct {
    double *fx, *fu, *fdot, *fz, *gx, *gu, *gz;
} dae_jacobians;

typedef struct {
    const matmpc_irk_config *config;
    const dae_layout *layout;
    const matmpc_irk_dae_model *model;
    const matmpc_irk_dae_input *input;
    double *stage_x;
    dae_jacobians jac;
} dae_system;

static int product(size_t a, size_t b, size_t *value)
{
    const size_t limit = SIZE_MAX / sizeof(double);
    if (a != 0 && b > limit / a) {
        return 0;
    }
    *value = a * b;
    return 1;
}

static int append(size_t *total, size_t count, size_t copies)
{
    const size_t limit = SIZE_MAX / sizeof(double);
    size_t length;
    if (!product(count, copies, &length) || length > limit - *total) {
        return 0;
    }
    *total += length;
    return 1;
}

static matmpc_irk_status layout_for(const matmpc_erk_shape *shape, size_t nz, dae_layout *d)
{
    *d = (dae_layout){0};
    if (shape == NULL || shape->nx == 0 || nz == 0 || shape->stages == 0 || shape->steps == 0) {
        return MATMPC_IRK_INVALID_ARGUMENT;
    }

    d->nx = shape->nx;
    d->nu = shape->nu;
    d->nz = nz;
    d->stages = shape->stages;
    d->slots = shape->adjoint ? shape->steps : 1;
    if (!append(&d->width, d->nx, 1) || !append(&d->width, nz, 1) ||
        !product(d->stages, d->width, &d->nk) ||
        !product(d->stages, d->stages, &d->tableau) ||
        !product(d->nk, d->nk, &d->system) ||
        !product(d->nx, d->nx, &d->xx) || !product(d->nx, d->nu, &d->xu) ||
        !product(d->nx, nz, &d->xz) || !product(nz, d->nu, &d->zu) ||
        !product(nz, nz, &d->zz) || !product(d->nk, d->slots, &d->pivots) ||
        d->pivots > SIZE_MAX / sizeof(size_t) ||
        !append(&d->doubles, d->nx, 2) || !append(&d->doubles, d->nk, 2) ||
        !append(&d->doubles, d->system, d->slots) ||
        !append(&d->doubles, d->xx, 2) || !append(&d->doubles, d->xz, 2) ||
        !append(&d->doubles, d->zz, 1)) {
        return MATMPC_IRK_SIZE_OVERFLOW;
    }
    if ((shape->forward || shape->adjoint) &&
        (!product(d->nk, d->nx, &d->state_rhs) || !product(d->nk, d->nu, &d->control_rhs) ||
         !append(&d->doubles, d->xu, 1) || !append(&d->doubles, d->zu, 1) ||
         !append(&d->doubles, d->state_rhs, d->slots) ||
         !append(&d->doubles, d->control_rhs, d->slots) ||
         !append(&d->doubles, d->state_rhs, 1))) {
        return MATMPC_IRK_SIZE_OVERFLOW;
    }
    if (shape->forward &&
        (!append(&d->doubles, d->xx, 1) || !append(&d->doubles, d->xu, 1) ||
         !append(&d->doubles, d->state_rhs, 1) || !append(&d->doubles, d->control_rhs, 1))) {
        return MATMPC_IRK_SIZE_OVERFLOW;
    }
    if (shape->adjoint &&
        (!append(&d->adjoint, d->nx, 1) || !append(&d->adjoint, d->nu, 1) ||
         !append(&d->doubles, d->adjoint, 1) || !append(&d->doubles, d->nk, 1))) {
        return MATMPC_IRK_SIZE_OVERFLOW;
    }
    return MATMPC_IRK_OK;
}

matmpc_irk_status matmpc_irk_dae_workspace_size(const matmpc_erk_shape *shape, size_t nz,
                                               size_t *double_count, size_t *pivot_count)
{
    dae_layout d;
    if (double_count == NULL || pivot_count == NULL) {
        return MATMPC_IRK_INVALID_ARGUMENT;
    }
    matmpc_irk_status status = layout_for(shape, nz, &d);
    if (status == MATMPC_IRK_OK) {
        *double_count = d.doubles;
        *pivot_count = d.pivots;
    }
    return status;
}

static bool input_view(matmpc_erk_view view, size_t count)
{
    return count == 0 || (view.data != NULL && view.count >= count);
}

static bool output_buffer(matmpc_erk_buffer buffer, size_t count)
{
    return count == 0 || (buffer.data != NULL && buffer.count >= count);
}

static bool all_finite(const double *values, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        if (!isfinite(values[i])) {
            return false;
        }
    }
    return true;
}

static double *take(double **cursor, size_t count)
{
    double *result = *cursor;
    *cursor += count;
    return result;
}

static matmpc_irk_status lu_status(matmpc_lu_status status)
{
    switch (status) {
    case MATMPC_LU_OK:
        return MATMPC_IRK_OK;
    case MATMPC_LU_SINGULAR:
        return MATMPC_IRK_SINGULAR;
    case MATMPC_LU_NONFINITE:
        return MATMPC_IRK_NONFINITE;
    default:
        return MATMPC_IRK_INVALID_ARGUMENT;
    }
}

static void store_system_stage(const dae_system *sys, size_t i, double *matrix,
                               double *state_jac, double *control_jac)
{
    const dae_layout *d = sys->layout;
    const dae_jacobians *j = &sys->jac;
    const size_t nx = d->nx, nz = d->nz;
    const size_t row = i * d->width, nk = d->nk;
    for (size_t k = 0; k < d->stages; ++k) {
        const double scale = sys->config->step * sys->config->a.data[i + k * d->stages];
        double *columns = matrix + k * d->width * nk + row;
        for (size_t c = 0; c < nx; ++c) {
            double *column = columns + c * nk;
            for (size_t r = 0; r < nx; ++r) {
                const size_t index = r + c * nx;
                column[r] = scale * j->fx[index];
                if (i == k) {
                    column[r] += j->fdot[index];
                }
            }
            for (size_t r = 0; r < nz; ++r) {
                column[nx + r] = scale * j->gx[r + c * nz];
            }
        }
        for (size_t c = 0; c < nz; ++c) {
            double *column = columns + (nx + c) * nk;
            if (i == k) {
                memcpy(column, j->fz + c * nx, nx * sizeof(double));
                memcpy(column + nx, j->gz + c * nz, nz * sizeof(double));
            } else {
                memset(column, 0, d->width * sizeof(double));
            }
        }
    }
    if (state_jac != NULL) {
        for (size_t c = 0; c < nx; ++c) {
            memcpy(state_jac + row + c * nk, j->fx + c * nx, nx * sizeof(double));
            memcpy(state_jac + row + c * nk + nx, j->gx + c * nz, nz * sizeof(double));
        }
    }
    if (control_jac != NULL) {
        for (size_t c = 0; c < d->nu; ++c) {
            memcpy(control_jac + row + c * nk, j->fu + c * nx, nx * sizeof(double));
            memcpy(control_jac + row + c * nk + nx, j->gu + c * nz, nz * sizeof(double));
        }
    }
}

static matmpc_irk_status form_system(const dae_system *sys, const double *state, const double *k,
                                     double *residual, double *matrix, double *state_jac,
                                     double *control_jac)
{
    const dae_layout *d = sys->layout;
    const matmpc_irk_config *config = sys->config;
    const matmpc_irk_dae_model *model = sys->model;
    const dae_jacobians *j = &sys->jac;
    const double *u = d->nu == 0 ? NULL : sys->input->differential.u.data;
    const double *p = sys->input->differential.p;
    for (size_t i = 0; i < d->stages; ++i) {
        const double *rate = k + i * d->width;
        const double *z = rate + d->nx;
        double *x = sys->stage_x;
        memcpy(x, state, d->nx * sizeof(double));
        for (size_t m = 0; m < d->stages; ++m) {
            const double scale = config->step * config->a.data[i + m * d->stages];
            if (scale != 0.0) {
                for (size_t q = 0; q < d->nx; ++q) {
                    x[q] += scale * k[m * d->width + q];
                }
            }
        }
        if (!all_finite(x, d->nx)) {
            return MATMPC_IRK_NONFINITE;
        }
        if (residual != NULL) {
            double *stage_residual = residual + i * d->width;
            if (model->residual(model->context, x, u, p, rate, z,
                                stage_residual, stage_residual + d->nx) != 0) {
                return MATMPC_IRK_MODEL_ERROR;
            }
            if (!all_finite(stage_residual, d->width)) {
                return MATMPC_IRK_NONFINITE;
            }
        }
        double *fu = control_jac != NULL && d->nu != 0 ? j->fu : NULL;
        double *gu = control_jac != NULL && d->nu != 0 ? j->gu : NULL;
        if (model->jacobian(model->context, x, u, p, rate, z,
                            j->fx, fu, j->fdot, j->fz, j->gx, gu, j->gz) != 0) {
            return MATMPC_IRK_MODEL_ERROR;
        }
        if (!all_finite(j->fx, d->xx) || !all_finite(j->fdot, d->xx) ||
            !all_finite(j->fz, d->xz) || !all_finite(j->gx, d->xz) ||
            !all_finite(j->gz, d->zz) ||
            (fu != NULL && (!all_finite(fu, d->xu) || !all_finite(gu, d->zu)))) {
            return MATMPC_IRK_NONFINITE;
        }
        store_system_stage(sys, i, matrix, state_jac, control_jac);
    }
    return MATMPC_IRK_OK;
}

static matmpc_irk_status newton_step(const dae_system *sys, const double *state, double *k,
                                     double *correction, double *matrix, size_t *pivots)
{
    const dae_layout *d = sys->layout;
    const matmpc_irk_config *config = sys->config;
    for (size_t i = 0; i < config->newton_iterations; ++i) {
        matmpc_irk_status status = form_system(sys, state, k, correction, matrix, NULL, NULL);
        if (status != MATMPC_IRK_OK) {
            return status;
        }
        status = lu_status(matmpc_lu_factor(d->nk, matrix, d->nk, pivots));
        if (status == MATMPC_IRK_OK) {
            status = lu_status(matmpc_lu_solve(0, d->nk, 1, matrix, d->nk, pivots, correction, d->nk));
        }
        if (status != MATMPC_IRK_OK) {
            return status;
        }
        for (size_t q = 0; q < d->nk; ++q) {
            k[q] -= correction[q];
        }
        if (!all_finite(k, d->nk)) {
            return MATMPC_IRK_NONFINITE;
        }
    }

    double norm = 0.0;
    for (size_t q = 0; q < d->nk; ++q) {
        norm = hypot(norm, correction[q]);
    }
    if (!isfinite(norm)) {
        return MATMPC_IRK_NONFINITE;
    }
    return norm > config->newton_tolerance ? MATMPC_IRK_NONCONVERGENCE : MATMPC_IRK_OK;
}

static matmpc_irk_status forward_step(const dae_layout *d, const double *system,
                                      const size_t *pivots, const double *jx, const double *ju,
                                      const double *weights, double scale,
                                      double *sx, double *su, double *rhs_x, double *rhs_u)
{
    const size_t nx = d->nx, nu = d->nu, nk = d->nk;
    matmpc_la_gemm(0, 0, nk, nx, nx, 1.0, jx, nk, sx, nx, 0.0, rhs_x, nk);
    matmpc_irk_status status = lu_status(matmpc_lu_solve(0, nk, nx, system, nk, pivots, rhs_x, nk));
    if (status != MATMPC_IRK_OK) {
        return status;
    }
    matmpc_la_gemm(0, 0, nx, nx, nk, -scale, weights, nx, rhs_x, nk, 1.0, sx, nx);
    if (nu != 0) {
        memcpy(rhs_u, ju, d->control_rhs * sizeof(double));
        matmpc_la_gemm(0, 0, nk, nu, nx, 1.0, jx, nk, su, nx, 1.0, rhs_u, nk);
        status = lu_status(matmpc_lu_solve(0, nk, nu, system, nk, pivots, rhs_u, nk));
        if (status != MATMPC_IRK_OK) {
            return status;
        }
        matmpc_la_gemm(0, 0, nx, nu, nk, -scale, weights, nx, rhs_u, nk, 1.0, su, nx);
    }
    return all_finite(sx, d->xx) && all_finite(su, d->xu) ? MATMPC_IRK_OK : MATMPC_IRK_NONFINITE;
}

static matmpc_irk_status adjoint_sweep(const dae_layout *d, const double *systems,
                                       const size_t *pivots, const double *state_jacs,
                                       const double *control_jacs, const double *weights,
                                       double scale, double *adjoint, double *rhs)
{
    const size_t nx = d->nx, nu = d->nu, nk = d->nk;
    for (size_t step = d->slots; step-- > 0;) {
        matmpc_la_gemv(1, nx, nk, -scale, weights, nx, adjoint, 0.0, rhs);
        matmpc_irk_status status = lu_status(matmpc_lu_solve(1, nk, 1, systems + step * d->system,
                                                          nk, pivots + step * nk, rhs, nk));
        if (status != MATMPC_IRK_OK) {
            return status;
        }
        matmpc_la_gemv(1, nk, nx, 1.0, state_jacs + step * d->state_rhs, nk, rhs, 1.0, adjoint);
        if (nu != 0) {
            matmpc_la_gemv(1, nk, nu, 1.0, control_jacs + step * d->control_rhs, nk, rhs, 1.0, adjoint + nx);
        }
        if (!all_finite(adjoint, d->adjoint)) {
            return MATMPC_IRK_NONFINITE;
        }
    }
    return MATMPC_IRK_OK;
}

matmpc_irk_status matmpc_irk_dae_integrate(const matmpc_irk_config *config, size_t nz,
                                          const matmpc_irk_dae_model *model,
                                          const matmpc_irk_dae_input *input,
                                          const matmpc_irk_dae_output *output,
                                          matmpc_irk_workspace workspace)
{
    if (config == NULL || model == NULL || input == NULL || output == NULL ||
        model->residual == NULL || model->jacobian == NULL) {
        return MATMPC_IRK_INVALID_ARGUMENT;
    }
    const matmpc_erk_shape *shape = &config->shape;
    const matmpc_erk_input *in = &input->differential;
    const matmpc_erk_output *out = &output->differential;
    dae_layout d;
    matmpc_irk_status status = layout_for(shape, nz, &d);
    if (status != MATMPC_IRK_OK) {
        return status;
    }
    const bool forward = shape->forward, adjoint = shape->adjoint, sensitivity = forward || adjoint;
    const bool custom_weights = sensitivity && config->sensitivity_weights.count != 0;
    if (!(config->step > 0.0) || !isfinite(config->step) || config->newton_iterations == 0 ||
        !(config->newton_tolerance >= 0.0) || !isfinite(config->newton_tolerance) ||
        !input_view(config->a, d.tableau) || !input_view(config->b, d.stages) ||
        !input_view(in->x, d.nx) || !input_view(in->u, d.nu) || !input_view(input->z, d.nz) ||
        !output_buffer(out->x, d.nx) || !output_buffer(output->z, d.nz) ||
        (forward && (!input_view(in->Sx, d.xx) || !input_view(in->Su, d.xu) ||
                     !output_buffer(out->Sx, d.xx) || !output_buffer(out->Su, d.xu))) ||
        (adjoint && (!input_view(in->lambda, d.nx) || !output_buffer(out->adjoint, d.adjoint))) ||
        (custom_weights && (!input_view(config->sensitivity_weights, d.state_rhs) ||
                            !(config->base_step > 0.0) || !isfinite(config->base_step)))) {
        return MATMPC_IRK_INVALID_ARGUMENT;
    }
    if (!output_buffer(workspace.doubles, d.doubles) || workspace.pivots == NULL ||
        workspace.pivot_count < d.pivots) {
        return MATMPC_IRK_WORKSPACE_TOO_SMALL;
    }
    const double scale = custom_weights ? config->step / config->base_step : 1.0;
    if (!isfinite(scale) || !all_finite(config->a.data, d.tableau) ||
        !all_finite(config->b.data, d.stages) || !all_finite(in->x.data, d.nx) ||
        !all_finite(in->u.data, d.nu) || !all_finite(input->z.data, d.nz) ||
        (forward && (!all_finite(in->Sx.data, d.xx) || !all_finite(in->Su.data, d.xu))) ||
        (adjoint && !all_finite(in->lambda.data, d.nx)) ||
        (custom_weights && !all_finite(config->sensitivity_weights.data, d.state_rhs))) {
        return MATMPC_IRK_NONFINITE;
    }

    double *cursor = workspace.doubles.data;
    double *state = take(&cursor, d.nx);
    double *stage_x = take(&cursor, d.nx);
    double *k = take(&cursor, d.nk);
    double *correction = take(&cursor, d.nk);
    double *systems = take(&cursor, d.system * d.slots);
    dae_jacobians jac;
    jac.fx = take(&cursor, d.xx);
    jac.fdot = take(&cursor, d.xx);
    jac.fz = take(&cursor, d.xz);
    jac.gx = take(&cursor, d.xz);
    jac.gz = take(&cursor, d.zz);
    jac.fu = sensitivity ? take(&cursor, d.xu) : NULL;
    jac.gu = sensitivity ? take(&cursor, d.zu) : NULL;
    const dae_system sys = {config, &d, model, input, stage_x, jac};
    double *state_jacs = sensitivity ? take(&cursor, d.state_rhs * d.slots) : NULL;
    double *control_jacs = sensitivity ? take(&cursor, d.control_rhs * d.slots) : NULL;
    double *default_weights = sensitivity ? take(&cursor, d.state_rhs) : NULL;
    const double *weights = custom_weights ? config->sensitivity_weights.data : default_weights;
    double *sx = forward ? take(&cursor, d.xx) : NULL;
    double *su = forward ? take(&cursor, d.xu) : NULL;
    double *rhs_x = forward ? take(&cursor, d.state_rhs) : NULL;
    double *rhs_u = forward ? take(&cursor, d.control_rhs) : NULL;
    double *adj = adjoint ? take(&cursor, d.adjoint) : NULL;
    double *adj_rhs = adjoint ? take(&cursor, d.nk) : NULL;

    memcpy(state, in->x.data, d.nx * sizeof(double));
    memset(k, 0, d.nk * sizeof(double));
    for (size_t stage = 0; stage < d.stages; ++stage) {
        memcpy(k + stage * d.width + d.nx, input->z.data, d.nz * sizeof(double));
    }
    if (sensitivity && !custom_weights) {
        memset(default_weights, 0, d.state_rhs * sizeof(double));
        for (size_t stage = 0; stage < d.stages; ++stage) {
            for (size_t q = 0; q < d.nx; ++q) {
                default_weights[q + (stage * d.width + q) * d.nx] = config->step * config->b.data[stage];
            }
        }
    }
    if (forward) {
        memcpy(sx, in->Sx.data, d.xx * sizeof(double));
        if (d.nu != 0) {
            memcpy(su, in->Su.data, d.xu * sizeof(double));
        }
    }
    if (adjoint) {
        memcpy(adj, in->lambda.data, d.nx * sizeof(double));
        memset(adj + d.nx, 0, d.nu * sizeof(double));
    }

    for (size_t step = 0; step < shape->steps; ++step) {
        const size_t slot = adjoint ? step : 0;
        double *system = systems + slot * d.system;
        size_t *pivots = workspace.pivots + slot * d.nk;
        status = newton_step(&sys, state, k, correction, system, pivots);
        if (status != MATMPC_IRK_OK) {
            return status;
        }
        if (sensitivity) {
            double *jx = state_jacs + slot * d.state_rhs;
            double *ju = control_jacs + slot * d.control_rhs;
            status = form_system(&sys, state, k, NULL, system, jx, ju);
            if (status == MATMPC_IRK_OK) {
                status = lu_status(matmpc_lu_factor(d.nk, system, d.nk, pivots));
            }
            if (status == MATMPC_IRK_OK && forward) {
                status = forward_step(&d, system, pivots, jx, ju, weights, scale, sx, su, rhs_x, rhs_u);
            }
            if (status != MATMPC_IRK_OK) {
                return status;
            }
        }
        for (size_t stage = 0; stage < d.stages; ++stage) {
            const double weight = config->step * config->b.data[stage];
            for (size_t q = 0; q < d.nx; ++q) {
                state[q] += weight * k[stage * d.width + q];
            }
        }
        if (!all_finite(state, d.nx)) {
            return MATMPC_IRK_NONFINITE;
        }
    }

    if (adjoint) {
        status = adjoint_sweep(&d, systems, workspace.pivots, state_jacs, control_jacs,
                               weights, scale, adj, adj_rhs);
        if (status != MATMPC_IRK_OK) {
            return status;
        }
    }
    memcpy(out->x.data, state, d.nx * sizeof(double));
    memcpy(output->z.data, k + (d.stages - 1) * d.width + d.nx, d.nz * sizeof(double));
    if (forward) {
        memcpy(out->Sx.data, sx, d.xx * sizeof(double));
        if (d.nu != 0) {
            memcpy(out->Su.data, su, d.xu * sizeof(double));
        }
    }
    if (adjoint) {
        memcpy(out->adjoint.data, adj, d.adjoint * sizeof(double));
    }
    return MATMPC_IRK_OK;
}
