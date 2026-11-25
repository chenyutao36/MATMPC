#include "erk.h"

#include "casadi_wrapper.h"
#include "erk_kernel.h"

#include <stdint.h>
#include <string.h>

static matmpc_erk_shape kernel_shape(const sim_opts *opts)
{
    const matmpc_erk_shape shape = {
        opts->nx, opts->nu, opts->num_stages, opts->num_steps, opts->forw_sens_flag, opts->adj_sens_flag
    };
    return shape;
}

static size_t matrix_count(size_t rows, size_t columns)
{
    if (columns != 0 && rows > (SIZE_MAX / sizeof(double)) / columns)
        mexErrMsgIdAndTxt("MATMPC:ERK:Workspace", "Runge-Kutta configuration is too large.");
    return rows * columns;
}

static size_t sum_count(size_t a, size_t b)
{
    if (b > SIZE_MAX / sizeof(double) - a)
        mexErrMsgIdAndTxt("MATMPC:ERK:Workspace", "Runge-Kutta configuration is too large.");
    return a + b;
}

sim_erk_workspace *sim_erk_workspace_create(sim_opts *opts)
{
    const matmpc_erk_shape shape = kernel_shape(opts);
    size_t scratch_count;
    if (matmpc_erk_workspace_size(&shape, &scratch_count) != MATMPC_ERK_OK)
        mexErrMsgIdAndTxt("MATMPC:ERK:Workspace", "Invalid Runge-Kutta workspace dimensions.");

    const size_t a = matrix_count(opts->num_stages, opts->num_stages);
    const size_t b = opts->num_stages;
    const size_t sx = matrix_count(opts->nx, opts->nx);
    const size_t su = matrix_count(opts->nx, opts->nu);
    const size_t coefficient_count = sum_count(sum_count(a, b), sum_count(sx, su));
    sim_erk_workspace *work = mxMalloc(sizeof(*work));
    work->A = mxMalloc(coefficient_count * sizeof(double));
    work->B = work->A + a;
    work->Sx = work->B + b;
    work->Su = work->Sx + sx;
    work->scratch = mxMalloc(scratch_count * sizeof(double));
    work->scratch_count = scratch_count;

    mexMakeMemoryPersistent(work->A);
    mexMakeMemoryPersistent(work->scratch);
    mexMakeMemoryPersistent(work);
    return work;
}

static void copy_coefficient(double *destination, const mxArray *mem, const char *name, size_t count)
{
    const mxArray *field = mxGetField(mem, 0, name);
    if (field == NULL || !mxIsDouble(field) || mxIsComplex(field) || mxIsSparse(field) ||
        mxGetNumberOfElements(field) != count)
        mexErrMsgIdAndTxt("MATMPC:ERK:Configuration", "Invalid Runge-Kutta configuration field '%s'.", name);
    if (count != 0)
        memcpy(destination, mxGetPr(field), count * sizeof(double));
}

void sim_erk_workspace_init(sim_opts *opts, const mxArray *mem, sim_erk_workspace *work)
{
    copy_coefficient(work->A, mem, "A_tab", matrix_count(opts->num_stages, opts->num_stages));
    copy_coefficient(work->B, mem, "B_tab", opts->num_stages);
    copy_coefficient(work->Sx, mem, "Sx", matrix_count(opts->nx, opts->nx));
    copy_coefficient(work->Su, mem, "Su", matrix_count(opts->nx, opts->nu));
}

/* The generated CasADi wrappers read their argument vectors. Their legacy API uses double **. */
static int model_rhs(void *context, const double *x, const double *u, const double *p, double *rate)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p, NULL};
    double *out[] = {rate};
    f_Fun(in, out);
    return 0;
}

static int model_forward(void *context, const double *x, const double *u, const double *p,
                         const double *sx, const double *su, double *dsx, double *dsu)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p, (double *)sx, (double *)su, NULL};
    double *out[] = {dsx, dsu};
    vde_Fun(in, out);
    return 0;
}

static int model_adjoint(void *context, const double *x, const double *u, const double *p,
                         const double *lambda, double *product)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p, (double *)lambda, NULL};
    double *out[] = {product};
    adj_ERK_Fun(in, out);
    return 0;
}

int sim_erk(sim_in *in, sim_out *out, sim_opts *opts, sim_erk_workspace *work)
{
    const matmpc_erk_shape shape = kernel_shape(opts);
    const size_t sx = opts->nx * opts->nx;
    const size_t su = opts->nx * opts->nu;
    const matmpc_erk_config config = {
        shape, opts->h, {work->A, opts->num_stages * opts->num_stages}, {work->B, opts->num_stages}
    };
    const matmpc_erk_model model = {NULL, model_rhs, model_forward, model_adjoint};
    const matmpc_erk_input input = {
        {in->x, opts->nx}, {in->u, opts->nu}, in->p, {work->Sx, sx}, {work->Su, su}, {in->lambda, opts->nx}
    };
    const matmpc_erk_output output = {
        {out->xn, opts->nx}, {out->Sx, sx}, {out->Su, su}, {out->adj_sens, opts->nx + opts->nu}
    };
    const matmpc_erk_buffer scratch = {work->scratch, work->scratch_count};
    const matmpc_erk_status status = matmpc_erk_integrate(&config, &model, &input, &output, scratch);
    if (status != MATMPC_ERK_OK)
        mexErrMsgIdAndTxt("MATMPC:ERK:Integration", "Explicit Runge-Kutta integration failed (status %d).", (int)status);
    return 0;
}

void sim_erk_workspace_free(sim_opts *opts, sim_erk_workspace *work)
{
    (void)opts;
    mxFree(work->A);
    mxFree(work->scratch);
    mxFree(work);
}
