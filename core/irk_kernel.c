#include "irk_kernel.h"
#include "matmpc_la.h"
#include "matmpc_lu.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    size_t tableau;
    size_t rates;
    size_t nx_nx;
    size_t nx_nu;
    size_t adjoint;
    size_t factor;
    size_t linearized_x;
    size_t linearized_u;
    size_t linearization;
    size_t tangent;
    size_t doubles;
    size_t pivots;
} ode_layout;

static bool multiply(size_t a, size_t b, size_t *product)
{
    if (b != 0 && a > SIZE_MAX / b)
        return false;
    *product = a * b;
    return true;
}

static bool add(size_t a, size_t b, size_t *sum)
{
    if (a > SIZE_MAX - b)
        return false;
    *sum = a + b;
    return true;
}

static matmpc_irk_status workspace_layout(const matmpc_erk_shape *shape, ode_layout *layout)
{
    if (shape == NULL || shape->nx == 0 || shape->stages == 0 || shape->steps == 0)
        return MATMPC_IRK_INVALID_ARGUMENT;

    const size_t nx = shape->nx;
    const size_t nu = shape->nu;
    const size_t stages = shape->stages;
    size_t scratch;
    if (!multiply(stages, stages, &layout->tableau) ||
        !multiply(stages, nx, &layout->rates) ||
        !multiply(nx, nx, &layout->nx_nx) || !multiply(nx, nu, &layout->nx_nu) ||
        !multiply(layout->rates, layout->rates, &layout->factor) ||
        !add(layout->rates, layout->nx_nx, &scratch) ||
        !multiply(2, scratch, &scratch) || !add(scratch, nx, &scratch))
        return MATMPC_IRK_SIZE_OVERFLOW;

    layout->linearization = layout->factor;
    if (shape->forward || shape->adjoint) {
        if (!add(nx, nu, &layout->adjoint) ||
            !multiply(layout->rates, nx, &layout->linearized_x) ||
            !multiply(layout->rates, nu, &layout->linearized_u) ||
            !add(layout->linearized_x, layout->linearized_u, &layout->tangent) ||
            !add(layout->factor, layout->tangent, &layout->linearization) ||
            !add(scratch, layout->nx_nu, &scratch))
            return MATMPC_IRK_SIZE_OVERFLOW;
    }
    if (shape->forward && !add(scratch, layout->tangent, &scratch))
        return MATMPC_IRK_SIZE_OVERFLOW;

    const size_t saved_steps = shape->adjoint ? shape->steps : 1;
    if (!multiply(layout->linearization, saved_steps, &layout->doubles) ||
        !add(layout->doubles, scratch, &layout->doubles) ||
        !multiply(layout->rates, saved_steps, &layout->pivots) ||
        layout->doubles > SIZE_MAX / sizeof(double) || layout->pivots > SIZE_MAX / sizeof(size_t))
        return MATMPC_IRK_SIZE_OVERFLOW;
    return MATMPC_IRK_OK;
}

matmpc_irk_status matmpc_irk_ode_workspace_size(const matmpc_erk_shape *shape,
                                               size_t *double_count, size_t *pivot_count)
{
    if (double_count == NULL || pivot_count == NULL || double_count == pivot_count)
        return MATMPC_IRK_INVALID_ARGUMENT;
    ode_layout layout = {0};
    const matmpc_irk_status status = workspace_layout(shape, &layout);
    *double_count = status == MATMPC_IRK_OK ? layout.doubles : 0;
    *pivot_count = status == MATMPC_IRK_OK ? layout.pivots : 0;
    return status;
}

static bool valid_view(matmpc_erk_view value, size_t count)
{
    return value.count == count && (count == 0 || value.data != NULL);
}

static bool valid_buffer(matmpc_erk_buffer value, size_t count)
{
    return value.count == count && (count == 0 || value.data != NULL);
}

static bool finite_values(const double *values, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        if (!isfinite(values[i]))
            return false;
    }
    return true;
}

static bool valid_config(const matmpc_irk_config *config, const ode_layout *layout)
{
    if (!isfinite(config->step) || config->step <= 0.0 ||
        config->newton_iterations == 0 || !isfinite(config->newton_tolerance) ||
        config->newton_tolerance < 0.0 ||
        !valid_view(config->a, layout->tableau) || !valid_view(config->b, config->shape.stages))
        return false;

    for (size_t i = 0; i < layout->tableau; ++i) {
        if (!isfinite(config->step * config->a.data[i]))
            return false;
    }
    for (size_t i = 0; i < config->shape.stages; ++i) {
        if (!isfinite(config->step * config->b.data[i]))
            return false;
    }

    if ((config->shape.forward || config->shape.adjoint) && config->sensitivity_weights.count != 0)
        return isfinite(config->base_step) && config->base_step > 0.0 &&
               isfinite(config->step / config->base_step) &&
               valid_view(config->sensitivity_weights, layout->linearized_x) &&
               finite_values(config->sensitivity_weights.data, layout->linearized_x);
    return true;
}

static void copy(double *destination, const double *source, size_t count)
{
    if (count != 0)
        memcpy(destination, source, count * sizeof(double));
}

static void axpy(double *destination, const double *source, size_t count, double scale)
{
    for (size_t i = 0; i < count; ++i)
        destination[i] += scale * source[i];
}

static double *partition(double **next, size_t count)
{
    double *start = *next;
    *next += count;
    return start;
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

typedef struct {
    double *stage_x;
    double *rates;
    double *correction;
    double *Fx;
    double *Fxdot;
    double *Fu;
    double *tangent;
} ode_scratch;

typedef struct {
    double *factor;
    double *x;
    double *u;
    size_t *pivots;
} ode_linearization;

static matmpc_irk_status assemble(const matmpc_irk_config *config, const ode_layout *layout,
                                 const matmpc_irk_ode_model *model, const matmpc_erk_input *input,
                                 const double *x, const ode_scratch *scratch,
                                 const ode_linearization *system, bool linearize)
{
    const size_t nx = config->shape.nx;
    const size_t nu = config->shape.nu;
    const size_t stages = config->shape.stages;
    const size_t order = layout->rates;
    double *Fu = linearize && nu != 0 ? scratch->Fu : NULL;

    for (size_t stage = 0; stage < stages; ++stage) {
        const size_t offset = stage * nx;
        const double *xdot = scratch->rates + offset;
        copy(scratch->stage_x, x, nx);
        for (size_t col = 0; col < stages; ++col) {
            const double a = config->step * config->a.data[col * stages + stage];
            if (a != 0.0)
                axpy(scratch->stage_x, scratch->rates + col * nx, nx, a);
        }
        if (!finite_values(scratch->stage_x, nx))
            return MATMPC_IRK_NONFINITE;
        if (!linearize) {
            double *residual = scratch->correction + offset;
            if (model->residual(model->context, scratch->stage_x, input->u.data, input->p, xdot, residual) != 0)
                return MATMPC_IRK_MODEL_ERROR;
            if (!finite_values(residual, nx))
                return MATMPC_IRK_NONFINITE;
        }
        if (model->jacobians(model->context, scratch->stage_x, input->u.data, input->p, xdot,
                             scratch->Fx, Fu, scratch->Fxdot) != 0)
            return MATMPC_IRK_MODEL_ERROR;
        if (!finite_values(scratch->Fx, layout->nx_nx) || !finite_values(scratch->Fxdot, layout->nx_nx) ||
            (Fu != NULL && !finite_values(Fu, layout->nx_nu)))
            return MATMPC_IRK_NONFINITE;

        for (size_t j = 0; j < stages; ++j) {
            const double a = config->step * config->a.data[j * stages + stage];
            for (size_t col = 0; col < nx; ++col) {
                double *block = system->factor + (j * nx + col) * order + offset;
                for (size_t row = 0; row < nx; ++row) {
                    const size_t element = col * nx + row;
                    block[row] = a * scratch->Fx[element];
                    if (j == stage)
                        block[row] += scratch->Fxdot[element];
                }
            }
        }

        if (linearize) {
            for (size_t col = 0; col < nx; ++col)
                copy(system->x + col * order + offset, scratch->Fx + col * nx, nx);
            for (size_t col = 0; col < nu; ++col)
                copy(system->u + col * order + offset, Fu + col * nx, nx);
        }
    }
    return lu_status(matmpc_lu_factor(order, system->factor, order, system->pivots));
}

static double correction_norm(const double *correction, size_t count)
{
    double norm = 0.0;
    for (size_t i = 0; i < count; ++i)
        norm = hypot(norm, correction[i]);
    return norm;
}

static void update_tangent(const matmpc_irk_config *config, const ode_layout *layout,
                           const double *stage_tangent, size_t columns, double *result)
{
    const size_t nx = config->shape.nx;
    if (config->sensitivity_weights.count != 0) {
        const double scale = -config->step / config->base_step;
        matmpc_la_gemm(0, 0, nx, columns, layout->rates, scale, config->sensitivity_weights.data,
                       nx, stage_tangent, layout->rates, 1.0, result, nx);
        return;
    }
    for (size_t col = 0; col < columns; ++col) {
        for (size_t stage = 0; stage < config->shape.stages; ++stage) {
            axpy(result + col * nx, stage_tangent + col * layout->rates + stage * nx,
                 nx, -config->step * config->b.data[stage]);
        }
    }
}

static matmpc_irk_status forward(const matmpc_irk_config *config, const ode_layout *layout,
                                 const ode_linearization *system, const matmpc_erk_output *output,
                                 double *tangent)
{
    const size_t nx = config->shape.nx;
    const size_t nu = config->shape.nu;
    const size_t order = layout->rates;
    double *tangent_u = tangent + layout->linearized_x;

    matmpc_la_gemm(0, 0, order, nx, nx, 1.0, system->x, order, output->Sx.data, nx,
                   0.0, tangent, order);
    if (nu != 0) {
        copy(tangent_u, system->u, layout->linearized_u);
        matmpc_la_gemm(0, 0, order, nu, nx, 1.0, system->x, order, output->Su.data, nx,
                       1.0, tangent_u, order);
    }
    const matmpc_irk_status status = lu_status(matmpc_lu_solve(0, order, layout->adjoint,
                                               system->factor, order, system->pivots, tangent, order));
    if (status != MATMPC_IRK_OK)
        return status;

    update_tangent(config, layout, tangent, nx, output->Sx.data);
    if (nu != 0)
        update_tangent(config, layout, tangent_u, nu, output->Su.data);
    if (!finite_values(output->Sx.data, layout->nx_nx) || !finite_values(output->Su.data, layout->nx_nu))
        return MATMPC_IRK_NONFINITE;
    return MATMPC_IRK_OK;
}

static matmpc_irk_status adjoint(const matmpc_irk_config *config, const ode_layout *layout,
                                 const ode_linearization *system, double *result, double *lambda)
{
    const size_t nx = config->shape.nx;
    const size_t nu = config->shape.nu;
    const size_t order = layout->rates;
    if (config->sensitivity_weights.count != 0) {
        const double scale = -config->step / config->base_step;
        matmpc_la_gemv(1, nx, order, scale, config->sensitivity_weights.data, nx, result, 0.0, lambda);
    } else {
        for (size_t stage = 0; stage < config->shape.stages; ++stage) {
            const double weight = -config->step * config->b.data[stage];
            for (size_t row = 0; row < nx; ++row)
                lambda[stage * nx + row] = weight * result[row];
        }
    }
    const matmpc_irk_status status = lu_status(matmpc_lu_solve(1, order, 1,
                                               system->factor, order, system->pivots, lambda, order));
    if (status != MATMPC_IRK_OK)
        return status;
    matmpc_la_gemv(1, order, nx, 1.0, system->x, order, lambda, 1.0, result);
    if (nu != 0)
        matmpc_la_gemv(1, order, nu, 1.0, system->u, order, lambda, 1.0, result + nx);
    return finite_values(result, layout->adjoint) ? MATMPC_IRK_OK : MATMPC_IRK_NONFINITE;
}

static ode_linearization load_system(double *data, size_t *pivots, const ode_layout *layout,
                                     size_t slot)
{
    double *factor = data + slot * layout->linearization;
    double *x = factor + layout->factor;
    ode_linearization system = {factor, x, x + layout->linearized_x, pivots + slot * layout->rates};
    return system;
}

matmpc_irk_status matmpc_irk_ode_integrate(const matmpc_irk_config *config,
                                         const matmpc_irk_ode_model *model,
                                         const matmpc_erk_input *input,
                                         const matmpc_erk_output *output,
                                         matmpc_irk_workspace workspace)
{
    if (config == NULL || model == NULL || input == NULL || output == NULL ||
        model->residual == NULL || model->jacobians == NULL)
        return MATMPC_IRK_INVALID_ARGUMENT;

    const matmpc_erk_shape *shape = &config->shape;
    ode_layout layout = {0};
    matmpc_irk_status status = workspace_layout(shape, &layout);
    if (status != MATMPC_IRK_OK)
        return status;

    const size_t nx = shape->nx;
    const size_t nu = shape->nu;
    if (!valid_config(config, &layout) || !valid_view(input->x, nx) ||
        !valid_view(input->u, nu) || !valid_buffer(output->x, nx) ||
        (shape->forward && (!valid_view(input->Sx, layout.nx_nx) || !valid_view(input->Su, layout.nx_nu) ||
                           !valid_buffer(output->Sx, layout.nx_nx) || !valid_buffer(output->Su, layout.nx_nu))) ||
        (shape->adjoint && (!valid_view(input->lambda, nx) || !valid_buffer(output->adjoint, layout.adjoint))))
        return MATMPC_IRK_INVALID_ARGUMENT;
    if (workspace.doubles.data == NULL || workspace.doubles.count < layout.doubles ||
        workspace.pivots == NULL || workspace.pivot_count < layout.pivots)
        return MATMPC_IRK_WORKSPACE_TOO_SMALL;
    if (!finite_values(input->x.data, nx) || !finite_values(input->u.data, nu) ||
        (shape->forward && (!finite_values(input->Sx.data, layout.nx_nx) ||
                           !finite_values(input->Su.data, layout.nx_nu))) ||
        (shape->adjoint && !finite_values(input->lambda.data, nx)))
        return MATMPC_IRK_NONFINITE;

    double *next = workspace.doubles.data;
    ode_scratch scratch;
    scratch.stage_x = partition(&next, nx);
    scratch.rates = partition(&next, layout.rates);
    scratch.correction = partition(&next, layout.rates);
    scratch.Fx = partition(&next, layout.nx_nx);
    scratch.Fxdot = partition(&next, layout.nx_nx);
    const bool linearize = shape->forward || shape->adjoint;
    scratch.Fu = linearize ? partition(&next, layout.nx_nu) : NULL;
    scratch.tangent = shape->forward ? partition(&next, layout.tangent) : NULL;
    double *systems = next;

    copy(output->x.data, input->x.data, nx);
    memset(scratch.rates, 0, layout.rates * sizeof(double));
    if (shape->forward) {
        copy(output->Sx.data, input->Sx.data, layout.nx_nx);
        copy(output->Su.data, input->Su.data, layout.nx_nu);
    }

    for (size_t step = 0; step < shape->steps; ++step) {
        const size_t slot = shape->adjoint ? step : 0;
        const ode_linearization system = load_system(systems, workspace.pivots, &layout, slot);

        /* The last iterate of each substep seeds the next substep. */
        for (size_t iteration = 0; iteration < config->newton_iterations; ++iteration) {
            status = assemble(config, &layout, model, input, output->x.data, &scratch, &system, false);
            if (status != MATMPC_IRK_OK)
                return status;
            status = lu_status(matmpc_lu_solve(0, layout.rates, 1, system.factor, layout.rates,
                                               system.pivots, scratch.correction, layout.rates));
            if (status != MATMPC_IRK_OK)
                return status;
            axpy(scratch.rates, scratch.correction, layout.rates, -1.0);
            if (!finite_values(scratch.rates, layout.rates))
                return MATMPC_IRK_NONFINITE;
        }
        if (correction_norm(scratch.correction, layout.rates) > config->newton_tolerance)
            return MATMPC_IRK_NONCONVERGENCE;

        if (linearize) {
            status = assemble(config, &layout, model, input, output->x.data, &scratch, &system, true);
            if (status != MATMPC_IRK_OK)
                return status;
        }
        if (shape->forward) {
            status = forward(config, &layout, &system, output, scratch.tangent);
            if (status != MATMPC_IRK_OK)
                return status;
        }

        for (size_t stage = 0; stage < shape->stages; ++stage) {
            axpy(output->x.data, scratch.rates + stage * nx, nx,
                 config->step * config->b.data[stage]);
        }
        if (!finite_values(output->x.data, nx))
            return MATMPC_IRK_NONFINITE;
    }

    if (shape->adjoint) {
        copy(output->adjoint.data, input->lambda.data, nx);
        for (size_t col = 0; col < nu; ++col)
            output->adjoint.data[nx + col] = 0.0;
        for (size_t step = shape->steps; step > 0; --step) {
            const ode_linearization system = load_system(systems, workspace.pivots, &layout, step - 1);
            status = adjoint(config, &layout, &system, output->adjoint.data, scratch.correction);
            if (status != MATMPC_IRK_OK)
                return status;
        }
    }
    return MATMPC_IRK_OK;
}
