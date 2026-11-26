#include "irk_mex.h"

#include <stdint.h>
#include <string.h>

static size_t multiply(size_t rows, size_t columns)
{
    if (columns != 0 && rows > (SIZE_MAX / sizeof(double)) / columns)
        mexErrMsgIdAndTxt("MATMPC:IRK:Workspace", "Implicit Runge-Kutta configuration is too large.");
    return rows * columns;
}

static size_t add(size_t a, size_t b)
{
    if (a > SIZE_MAX / sizeof(double) || b > SIZE_MAX / sizeof(double) - a)
        mexErrMsgIdAndTxt("MATMPC:IRK:Workspace", "Implicit Runge-Kutta configuration is too large.");
    return a + b;
}

static size_t weight_count(const sim_opts *opts)
{
    return multiply(opts->nx, multiply(opts->num_stages, add(opts->nx, opts->nz)));
}

sim_irk_workspace *sim_irk_mex_create(const sim_opts *opts, size_t doubles, size_t pivots)
{
    const size_t a = multiply(opts->num_stages, opts->num_stages);
    const size_t b = opts->num_stages;
    const size_t sx = multiply(opts->nx, opts->nx);
    const size_t su = multiply(opts->nx, opts->nu);
    const size_t count = add(add(add(a, b), add(sx, su)), weight_count(opts));
    sim_irk_workspace *work = mxMalloc(sizeof(*work));
    work->A = mxMalloc(count * sizeof(double));
    work->B = work->A + a;
    work->Sx = work->B + b;
    work->Su = work->Sx + sx;
    work->JFK = work->Su + su;
    work->scratch.doubles.data = mxMalloc(doubles * sizeof(double));
    work->scratch.doubles.count = doubles;
    work->scratch.pivots = mxMalloc(pivots * sizeof(size_t));
    work->scratch.pivot_count = pivots;
    mexMakeMemoryPersistent(work->A);
    mexMakeMemoryPersistent(work->scratch.doubles.data);
    mexMakeMemoryPersistent(work->scratch.pivots);
    mexMakeMemoryPersistent(work);
    return work;
}

static void copy_field(double *destination, const mxArray *mem, const char *name, size_t count)
{
    const mxArray *field = mxGetField(mem, 0, name);
    if (field == NULL || !mxIsDouble(field) || mxIsComplex(field) || mxIsSparse(field) ||
        mxGetNumberOfElements(field) != count)
        mexErrMsgIdAndTxt("MATMPC:IRK:Configuration", "Invalid implicit integrator field '%s'.", name);
    if (count != 0)
        memcpy(destination, mxGetPr(field), count * sizeof(double));
}

void sim_irk_mex_init(const sim_opts *opts, const mxArray *mem, sim_irk_workspace *work)
{
    copy_field(work->A, mem, "A_tab", multiply(opts->num_stages, opts->num_stages));
    copy_field(work->B, mem, "B_tab", opts->num_stages);
    copy_field(work->Sx, mem, "Sx", multiply(opts->nx, opts->nx));
    copy_field(work->Su, mem, "Su", multiply(opts->nx, opts->nu));
    copy_field(work->JFK, mem, "JFK", weight_count(opts));
    work->newton_iter = (int)mxGetScalar(mxGetField(mem, 0, "newton_iter"));
}

void sim_irk_mex_free(sim_irk_workspace *work)
{
    mxFree(work->A);
    mxFree(work->scratch.doubles.data);
    mxFree(work->scratch.pivots);
    mxFree(work);
}

matmpc_irk_config sim_irk_mex_config(const sim_opts *opts, const sim_irk_workspace *work)
{
    const matmpc_irk_config config = {
        {opts->nx, opts->nu, opts->num_stages, opts->num_steps, opts->forw_sens_flag, opts->adj_sens_flag},
        opts->h, {work->A, opts->num_stages * opts->num_stages}, {work->B, opts->num_stages},
        (size_t)work->newton_iter, 1e-2, {work->JFK, weight_count(opts)}, opts->base_h
    };
    return config;
}

void sim_irk_mex_check_status(matmpc_irk_status status, const char *integrator)
{
    switch (status) {
        case MATMPC_IRK_OK:
            break;
        case MATMPC_IRK_SINGULAR:
            mexErrMsgIdAndTxt("MATMPC:IRK:Singular", "%s collocation Jacobian is singular.", integrator);
            break;
        case MATMPC_IRK_NONFINITE:
            mexErrMsgIdAndTxt("MATMPC:IRK:Nonfinite", "%s produced nonfinite collocation data.", integrator);
            break;
        case MATMPC_IRK_NONCONVERGENCE:
            mexErrMsgIdAndTxt("MATMPC:IRK:Nonconvergence", "%s Newton iterations did not converge.", integrator);
            break;
        default:
            mexErrMsgIdAndTxt("MATMPC:IRK:Integration", "%s integration failed (status %d).", integrator, (int)status);
            break;
    }
}
