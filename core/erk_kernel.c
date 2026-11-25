#include "erk_kernel.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    size_t tableau;
    size_t Sx;
    size_t Su;
    size_t rates;
    size_t trajectory;
    size_t adjoint_width;
    size_t count;
} workspace_layout;

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

static bool add_buffer(size_t count, size_t width, size_t *total)
{
    size_t length;
    return multiply(count, width, &length) && add(*total, length, total);
}

static matmpc_erk_status layout_workspace(const matmpc_erk_shape *shape, workspace_layout *layout)
{
    if (shape == NULL || shape->nx == 0 || shape->stages == 0 || shape->steps == 0)
        return MATMPC_ERK_INVALID_ARGUMENT;

    const size_t nx = shape->nx;
    const size_t nu = shape->nu;
    const size_t stages = shape->stages;
    if (!multiply(stages, stages, &layout->tableau) || !multiply(stages, nx, &layout->rates) ||
        !add(nx, layout->rates, &layout->count))
        return MATMPC_ERK_SIZE_OVERFLOW;

    if (shape->forward) {
        size_t buffers;
        if (!multiply(nx, nx, &layout->Sx) || !multiply(nx, nu, &layout->Su) ||
            !add(stages, 1, &buffers) || !add_buffer(buffers, layout->Sx, &layout->count) ||
            !add_buffer(buffers, layout->Su, &layout->count))
            return MATMPC_ERK_SIZE_OVERFLOW;
    }
    if (shape->adjoint) {
        if (!multiply(shape->steps, layout->rates, &layout->trajectory) ||
            !add(nx, nu, &layout->adjoint_width) || !add(layout->count, layout->trajectory, &layout->count) ||
            !add_buffer(stages, layout->adjoint_width, &layout->count) || !add(nx, layout->count, &layout->count))
            return MATMPC_ERK_SIZE_OVERFLOW;
    }
    if (layout->count > SIZE_MAX / sizeof(double) || layout->tableau > SIZE_MAX / sizeof(double))
        return MATMPC_ERK_SIZE_OVERFLOW;
    return MATMPC_ERK_OK;
}

matmpc_erk_status matmpc_erk_workspace_size(const matmpc_erk_shape *shape, size_t *count)
{
    if (count == NULL)
        return MATMPC_ERK_INVALID_ARGUMENT;
    workspace_layout layout = {0};
    const matmpc_erk_status status = layout_workspace(shape, &layout);
    *count = status == MATMPC_ERK_OK ? layout.count : 0;
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

static bool valid_tableau(const matmpc_erk_config *config, size_t tableau_count)
{
    const size_t stages = config->shape.stages;
    if (!isfinite(config->step) || config->step <= 0.0 || !valid_view(config->a, tableau_count) ||
        !valid_view(config->b, stages))
        return false;

    for (size_t col = 0; col < stages; ++col) {
        if (!isfinite(config->b.data[col]))
            return false;
        for (size_t row = 0; row < stages; ++row) {
            const double coefficient = config->a.data[col * stages + row];
            if (!isfinite(coefficient) || (row <= col && coefficient != 0.0))
                return false;
        }
    }
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

matmpc_erk_status matmpc_erk_integrate(const matmpc_erk_config *config, const matmpc_erk_model *model,
                                     const matmpc_erk_input *input, const matmpc_erk_output *output,
                                     matmpc_erk_buffer workspace)
{
    if (config == NULL || model == NULL || input == NULL || output == NULL)
        return MATMPC_ERK_INVALID_ARGUMENT;

    const matmpc_erk_shape *shape = &config->shape;
    workspace_layout layout = {0};
    const matmpc_erk_status status = layout_workspace(shape, &layout);
    if (status != MATMPC_ERK_OK)
        return status;
    if (workspace.data == NULL || workspace.count < layout.count)
        return MATMPC_ERK_WORKSPACE_TOO_SMALL;

    const size_t nx = shape->nx;
    const size_t stages = shape->stages;
    const size_t sx_size = layout.Sx;
    const size_t su_size = layout.Su;
    if (model->rhs == NULL || !valid_tableau(config, layout.tableau) ||
        !valid_view(input->x, nx) || !valid_view(input->u, shape->nu) || !valid_buffer(output->x, nx))
        return MATMPC_ERK_INVALID_ARGUMENT;
    if (shape->forward && (model->forward == NULL || !valid_view(input->Sx, sx_size) ||
        !valid_view(input->Su, su_size) || !valid_buffer(output->Sx, sx_size) || !valid_buffer(output->Su, su_size)))
        return MATMPC_ERK_INVALID_ARGUMENT;
    if (shape->adjoint && (model->adjoint == NULL || !valid_view(input->lambda, nx) ||
        !valid_buffer(output->adjoint, layout.adjoint_width)))
        return MATMPC_ERK_INVALID_ARGUMENT;

    double *stage_x = workspace.data;
    double *rates = stage_x + nx;
    double *stage_sx = rates + layout.rates;
    double *stage_su = stage_sx + sx_size;
    double *rates_sx = stage_su + su_size;
    double *rates_su = rates_sx + stages * sx_size;
    double *trajectory = rates_su + stages * su_size;
    double *adjoint_rates = trajectory + layout.trajectory;
    double *stage_lambda = adjoint_rates + stages * layout.adjoint_width;
    const double h = config->step;
    const double *a = config->a.data;
    const double *b = config->b.data;

    copy(output->x.data, input->x.data, nx);
    if (shape->forward) {
        copy(output->Sx.data, input->Sx.data, sx_size);
        copy(output->Su.data, input->Su.data, su_size);
    }
    for (size_t step = 0; step < shape->steps; ++step) {
        for (size_t stage = 0; stage < stages; ++stage) {
            copy(stage_x, output->x.data, nx);
            if (shape->forward) {
                copy(stage_sx, output->Sx.data, sx_size);
                copy(stage_su, output->Su.data, su_size);
            }
            for (size_t prev = 0; prev < stage; ++prev) {
                const double weight = h * a[prev * stages + stage];
                if (weight == 0.0)
                    continue;
                axpy(stage_x, rates + prev * nx, nx, weight);
                if (shape->forward) {
                    axpy(stage_sx, rates_sx + prev * sx_size, sx_size, weight);
                    axpy(stage_su, rates_su + prev * su_size, su_size, weight);
                }
            }
            if (model->rhs(model->context, stage_x, input->u.data, input->p, rates + stage * nx) != 0)
                return MATMPC_ERK_MODEL_ERROR;
            if (shape->forward && model->forward(model->context, stage_x, input->u.data, input->p, stage_sx,
                    stage_su, rates_sx + stage * sx_size, rates_su + stage * su_size) != 0)
                return MATMPC_ERK_MODEL_ERROR;
            if (shape->adjoint)
                copy(trajectory + (step * stages + stage) * nx, stage_x, nx);
        }
        for (size_t stage = 0; stage < stages; ++stage) {
            const double weight = h * b[stage];
            axpy(output->x.data, rates + stage * nx, nx, weight);
            if (shape->forward) {
                axpy(output->Sx.data, rates_sx + stage * sx_size, sx_size, weight);
                axpy(output->Su.data, rates_su + stage * su_size, su_size, weight);
            }
        }
    }
    if (shape->adjoint) {
        const size_t width = layout.adjoint_width;
        double *adjoint = output->adjoint.data;
        copy(adjoint, input->lambda.data, nx);
        for (size_t i = nx; i < width; ++i)
            adjoint[i] = 0.0;

        for (size_t step = shape->steps; step-- > 0;) {
            for (size_t stage = stages; stage-- > 0;) {
                for (size_t i = 0; i < nx; ++i)
                    stage_lambda[i] = h * b[stage] * adjoint[i];
                for (size_t next = stage + 1; next < stages; ++next) {
                    const double weight = h * a[stage * stages + next];
                    if (weight != 0.0)
                        axpy(stage_lambda, adjoint_rates + next * width, nx, weight);
                }
                if (model->adjoint(model->context, trajectory + (step * stages + stage) * nx, input->u.data,
                                   input->p, stage_lambda, adjoint_rates + stage * width) != 0)
                    return MATMPC_ERK_MODEL_ERROR;
            }
            for (size_t stage = 0; stage < stages; ++stage)
                axpy(adjoint, adjoint_rates + stage * width, width, 1.0);
        }
    }
    return MATMPC_ERK_OK;
}
