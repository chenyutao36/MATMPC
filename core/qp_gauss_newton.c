#include "qp_gauss_newton.h"

#include "matmpc_la.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    size_t x, u, parameters, references, state_bounds, path_bounds;
    size_t A, B, a, Q, R, Cgx, Cgu, CgN;
    size_t Jx, Ju, J_terminal, integration, count;
    size_t blocks, R_block, gu_block;
} qp_layout;

static bool multiply(size_t a, size_t b, size_t *count)
{
    if (b != 0 && a > (SIZE_MAX / sizeof(double)) / b)
        return false;
    *count = a * b;
    return true;
}

static bool add(size_t a, size_t b, size_t *count)
{
    const size_t maximum = SIZE_MAX / sizeof(double);
    if (a > maximum || b > maximum - a)
        return false;
    *count = a + b;
    return true;
}

static bool stages(size_t rows, size_t columns, size_t count, size_t *elements)
{
    size_t width;
    return multiply(rows, columns, &width) && multiply(width, count, elements);
}

static matmpc_qp_status layout_for(const matmpc_qp_config *config, qp_layout *layout)
{
    if (config == NULL || config->horizon == 0 || config->integration.shape.stages != 4 ||
        !config->integration.shape.forward || config->integration.shape.adjoint)
        return MATMPC_QP_INVALID_ARGUMENT;

    const size_t nx = config->integration.shape.nx;
    const size_t nu = config->integration.shape.nu;
    const size_t N = config->horizon;
    size_t states, path, jacobians;
    if (!add(N, 1, &states) || !multiply(nx, states, &layout->x) || !multiply(nu, N, &layout->u) ||
        !multiply(config->np, states, &layout->parameters) || !multiply(config->ny, N, &layout->references) ||
        !multiply(config->nbx, N, &layout->state_bounds) || !multiply(config->nc, N, &path) ||
        !add(path, config->nc_terminal, &layout->path_bounds) || !multiply(nx, N, &layout->a) ||
        !stages(nx, nx, N, &layout->A) || !stages(nx, nu, N, &layout->B) ||
        !stages(nx, nx, states, &layout->Q) || !stages(nu, nu, N, &layout->R) ||
        !stages(config->nc, nx, N, &layout->Cgx) || !stages(config->nc, nu, N, &layout->Cgu) ||
        !multiply(config->nc_terminal, nx, &layout->CgN) || !multiply(config->ny, nx, &layout->Jx) ||
        !multiply(config->ny, nu, &layout->Ju) || !multiply(config->ny_terminal, nx, &layout->J_terminal) ||
        !add(layout->Jx, layout->Ju, &jacobians))
        return MATMPC_QP_SIZE_OVERFLOW;

    const matmpc_erk_view blocks = config->control_blocks;
    if (blocks.count != 0) {
        if (config->stage_grid.count != 0 || nu == 0 || blocks.data == NULL || blocks.count < 2 || blocks.count > states ||
            blocks.data[0] != 0.0 || blocks.data[blocks.count - 1] != (double)N)
            return MATMPC_QP_INVALID_ARGUMENT;
        for (size_t k = 1; k < blocks.count; ++k) {
            const double boundary = blocks.data[k];
            if (!(boundary > blocks.data[k - 1] && boundary <= (double)N) || floor(boundary) != boundary)
                return MATMPC_QP_INVALID_ARGUMENT;
        }
        layout->blocks = blocks.count - 1;
        layout->R_block = layout->blocks * (nu * nu);
        layout->gu_block = layout->blocks * nu;
    }

    const matmpc_erk_status integration = matmpc_erk_workspace_size(&config->integration.shape, &layout->integration);
    if (integration == MATMPC_ERK_SIZE_OVERFLOW)
        return MATMPC_QP_SIZE_OVERFLOW;
    if (integration != MATMPC_ERK_OK)
        return MATMPC_QP_INVALID_ARGUMENT;

    const size_t cost_scratch = jacobians > layout->J_terminal ? jacobians : layout->J_terminal;
    if (!add(layout->integration, cost_scratch, &layout->count))
        return MATMPC_QP_SIZE_OVERFLOW;
    return MATMPC_QP_OK;
}

matmpc_qp_status matmpc_qp_workspace_size(const matmpc_qp_config *config, size_t *count)
{
    if (count == NULL)
        return MATMPC_QP_INVALID_ARGUMENT;
    qp_layout layout = {0};
    const matmpc_qp_status status = layout_for(config, &layout);
    *count = status == MATMPC_QP_OK ? layout.count : 0;
    return status;
}

static bool valid_input(matmpc_erk_view value, size_t count)
{
    return value.count == count && (count == 0 || value.data != NULL);
}

static bool valid_output(matmpc_erk_buffer value, size_t count)
{
    return value.count == count && (count == 0 || value.data != NULL);
}

static bool input_matches(const matmpc_qp_config *config, const qp_layout *layout, const matmpc_qp_input *input)
{
    return valid_input(input->x0, config->integration.shape.nx) && valid_input(input->x, layout->x) &&
           valid_input(input->u, layout->u) && valid_input(input->parameters, layout->parameters) &&
           valid_input(input->reference, layout->references) && valid_input(input->weight, layout->references) &&
           valid_input(input->reference_terminal, config->ny_terminal) && valid_input(input->weight_terminal, config->ny_terminal) &&
           valid_input(input->lower_u, layout->u) && valid_input(input->upper_u, layout->u) &&
           valid_input(input->lower_x, layout->state_bounds) && valid_input(input->upper_x, layout->state_bounds) &&
           valid_input(input->lower_path, layout->path_bounds) && valid_input(input->upper_path, layout->path_bounds);
}

static bool output_matches(const matmpc_qp_config *config, const qp_layout *layout, const matmpc_qp_output *output)
{
    return valid_output(output->A, layout->A) && valid_output(output->B, layout->B) && valid_output(output->a, layout->a) &&
           valid_output(output->ds0, config->integration.shape.nx) && valid_output(output->Q, layout->Q) &&
           valid_output(output->S, layout->B) && valid_output(output->R, layout->R) &&
           valid_output(output->gx, layout->x) && valid_output(output->gu, layout->u) &&
           valid_output(output->Cgx, layout->Cgx) && valid_output(output->Cgu, layout->Cgu) && valid_output(output->CgN, layout->CgN) &&
           valid_output(output->lc, layout->path_bounds) && valid_output(output->uc, layout->path_bounds) &&
           valid_output(output->lb_du, layout->u) && valid_output(output->ub_du, layout->u) &&
           valid_output(output->lb_dx, layout->state_bounds) && valid_output(output->ub_dx, layout->state_bounds) &&
           valid_output(output->R_block, layout->R_block) && valid_output(output->gu_block, layout->gu_block);
}

static bool valid_stage_grid(const matmpc_qp_config *config)
{
    const double base_step = config->integration.step;
    const matmpc_erk_view grid = config->stage_grid;
    if (!isfinite(base_step) || base_step <= 0.0)
        return false;
    if (grid.count == 0)
        return true;
    if (grid.data == NULL || grid.count != config->horizon + 1)
        return false;

    for (size_t i = 0; i < config->horizon; ++i) {
        const double step = base_step * (grid.data[i + 1] - grid.data[i]);
        if (!isfinite(step) || step <= 0.0)
            return false;
    }
    return true;
}

static void gram(size_t rows, size_t columns, size_t residuals, const double *left, const double *right, double *result)
{
    if (rows == 0 || columns == 0)
        return;
    if (residuals != 0) {
        matmpc_la_gemm(1, 0, rows, columns, residuals, 1.0, left, residuals,
                       right, residuals, 0.0, result, rows);
    } else {
        for (size_t i = 0; i < rows * columns; ++i)
            result[i] = 0.0;
    }
}

static void floor_diagonal(double *matrix, size_t n, double minimum)
{
    for (size_t i = 0; i < n; ++i) {
        if (matrix[i * (n + 1)] < minimum)
            matrix[i * (n + 1)] = minimum;
    }
}

static void aggregate_controls(const matmpc_qp_config *config, const qp_layout *layout, const matmpc_qp_output *output)
{
    const size_t nu = config->integration.shape.nu;
    const size_t width = nu * nu;
    const double *boundaries = config->control_blocks.data;
    for (size_t block = 0; block < layout->blocks; ++block) {
        const size_t first = (size_t)boundaries[block];
        const size_t end = (size_t)boundaries[block + 1];
        double *R = output->R_block.data + block * width;
        double *gu = output->gu_block.data + block * nu;
        memcpy(R, output->R.data + first * width, width * sizeof(*R));
        memcpy(gu, output->gu.data + first * nu, nu * sizeof(*gu));
        for (size_t stage = first + 1; stage < end; ++stage) {
            for (size_t j = 0; j < width; ++j)
                R[j] += output->R.data[stage * width + j];
            for (size_t j = 0; j < nu; ++j)
                gu[j] += output->gu.data[stage * nu + j];
        }
        floor_diagonal(R, nu, config->diagonal_floor);
    }
}

matmpc_qp_status matmpc_qp_gauss_newton(const matmpc_qp_config *config, const matmpc_qp_model *model,
                                        const matmpc_qp_input *input, const matmpc_qp_output *output,
                                        matmpc_erk_buffer workspace)
{
    if (config == NULL || model == NULL || input == NULL || output == NULL)
        return MATMPC_QP_INVALID_ARGUMENT;
    qp_layout layout = {0};
    const matmpc_qp_status status = layout_for(config, &layout);
    if (status != MATMPC_QP_OK)
        return status;
    if (workspace.data == NULL || workspace.count < layout.count)
        return MATMPC_QP_WORKSPACE_TOO_SMALL;

    const size_t nx = config->integration.shape.nx;
    const size_t nu = config->integration.shape.nu;
    const size_t nc = config->nc;
    const size_t N = config->horizon;
    if (!isfinite(config->diagonal_floor) || !valid_input(config->state_seed, nx * nx) ||
        !valid_input(config->control_seed, nx * nu) || (config->nbx != 0 && config->bounded_states == NULL) ||
        model->stage_cost == NULL || model->terminal_cost == NULL || (nc != 0 && model->path_constraint == NULL) ||
        (config->nc_terminal != 0 && model->terminal_constraint == NULL) ||
        !valid_stage_grid(config) || !input_matches(config, &layout, input) || !output_matches(config, &layout, output))
        return MATMPC_QP_INVALID_ARGUMENT;
    for (size_t j = 0; j < config->nbx; ++j) {
        if (config->bounded_states[j] >= nx)
            return MATMPC_QP_INVALID_ARGUMENT;
    }

    const matmpc_erk_buffer integration_workspace = {workspace.data, layout.integration};
    matmpc_erk_config integration = config->integration;
    const matmpc_erk_view grid = config->stage_grid;
    double *Jx = workspace.data + layout.integration;
    double *Ju = Jx + layout.Jx;
    for (size_t j = 0; j < nx; ++j)
        output->ds0.data[j] = input->x0.data[j] - input->x.data[j];

    for (size_t i = 0; i < N; ++i) {
        const double *x = input->x.data + i * nx;
        const double *u = nu != 0 ? input->u.data + i * nu : NULL;
        const double *p = config->np != 0 ? input->parameters.data + i * config->np : NULL;
        const double *y = config->ny != 0 ? input->reference.data + i * config->ny : NULL;
        const double *W = config->ny != 0 ? input->weight.data + i * config->ny : NULL;
        double *a = output->a.data + i * nx;
        const matmpc_erk_input step_input = {{x, nx}, {u, nu}, p, config->state_seed, config->control_seed, {NULL, 0}};
        const matmpc_erk_output step_output = {
            {a, nx}, {output->A.data + i * nx * nx, nx * nx},
            {nu != 0 ? output->B.data + i * nx * nu : NULL, nx * nu}, {NULL, 0}
        };
        if (grid.count != 0)
            integration.step = config->integration.step * (grid.data[i + 1] - grid.data[i]);
        const matmpc_erk_status step_status = matmpc_erk_integrate(&integration, &model->dynamics,
                                                                  &step_input, &step_output, integration_workspace);
        if (step_status != MATMPC_ERK_OK)
            return step_status == MATMPC_ERK_MODEL_ERROR ? MATMPC_QP_MODEL_ERROR : MATMPC_QP_INVALID_ARGUMENT;
        const double *next = x + nx;
        for (size_t j = 0; j < nx; ++j)
            a[j] -= next[j];

        double *Q = output->Q.data + i * nx * nx;
        double *S = nu != 0 ? output->S.data + i * nx * nu : NULL;
        double *R = nu != 0 ? output->R.data + i * nu * nu : NULL;
        if (model->stage_cost(model->context, x, u, p, y, W, Jx, Ju, output->gx.data + i * nx,
                              nu != 0 ? output->gu.data + i * nu : NULL) != 0)
            return MATMPC_QP_MODEL_ERROR;
        gram(nx, nx, config->ny, Jx, Jx, Q);
        gram(nx, nu, config->ny, Jx, Ju, S);
        gram(nu, nu, config->ny, Ju, Ju, R);
        floor_diagonal(Q, nx, config->diagonal_floor);
        if (layout.blocks == 0)
            floor_diagonal(R, nu, config->diagonal_floor);

        for (size_t j = 0; j < nu; ++j) {
            const size_t k = i * nu + j;
            output->lb_du.data[k] = input->lower_u.data[k] - u[j];
            output->ub_du.data[k] = input->upper_u.data[k] - u[j];
        }
        for (size_t j = 0; j < config->nbx; ++j) {
            const size_t k = i * config->nbx + j;
            const double value = next[config->bounded_states[j]];
            output->lb_dx.data[k] = input->lower_x.data[k] - value;
            output->ub_dx.data[k] = input->upper_x.data[k] - value;
        }
        if (nc != 0) {
            if (model->path_constraint(model->context, x, u, p, output->lc.data + i * nc,
                                       output->Cgx.data + i * nc * nx,
                                       nu != 0 ? output->Cgu.data + i * nc * nu : NULL) != 0)
                return MATMPC_QP_MODEL_ERROR;
            for (size_t j = 0; j < nc; ++j) {
                const size_t k = i * nc + j;
                output->uc.data[k] = input->upper_path.data[k] - output->lc.data[k];
                output->lc.data[k] = input->lower_path.data[k] - output->lc.data[k];
            }
        }
    }

    const double *x = input->x.data + N * nx;
    const double *p = config->np != 0 ? input->parameters.data + N * config->np : NULL;
    if (model->terminal_cost(model->context, x, p, input->reference_terminal.data, input->weight_terminal.data,
                             Jx, output->gx.data + N * nx) != 0)
        return MATMPC_QP_MODEL_ERROR;
    double *Q = output->Q.data + N * nx * nx;
    gram(nx, nx, config->ny_terminal, Jx, Jx, Q);
    floor_diagonal(Q, nx, config->diagonal_floor);
    if (config->nc_terminal != 0) {
        if (model->terminal_constraint(model->context, x, p, output->lc.data + N * nc, output->CgN.data) != 0)
            return MATMPC_QP_MODEL_ERROR;
        for (size_t j = N * nc; j < layout.path_bounds; ++j) {
            output->uc.data[j] = input->upper_path.data[j] - output->lc.data[j];
            output->lc.data[j] = input->lower_path.data[j] - output->lc.data[j];
        }
    }
    aggregate_controls(config, &layout, output);
    return MATMPC_QP_OK;
}
