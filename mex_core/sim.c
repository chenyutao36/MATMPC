#include "mex.h"
#include "sim_workspace.h"
#include <limits.h>
#include <string.h>


sim_opts* sim_opts_create(const mxArray *mem)
{
    sim_opts* opts = (sim_opts*)mxMalloc(sizeof(sim_opts));
    mexMakeMemoryPersistent(opts);
       
    opts->nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    opts->nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    opts->nz = mxGetScalar( mxGetField(mem, 0, "nz") );
    opts->num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    opts->num_steps = mxGetScalar( mxGetField(mem, 0, "num_steps") );
    opts->h = mxGetScalar( mxGetField(mem, 0, "h") );
    opts->base_h = opts->h;
    
	return opts;
}

sim_in* sim_in_create(sim_opts *opts)
{
    sim_in* in = (sim_in*)mxCalloc(1, sizeof(sim_in));
    mexMakeMemoryPersistent(in);
    
    return in;
}

sim_out* sim_out_create(sim_opts *opts)
{
    sim_out* out = (sim_out*)mxCalloc(1, sizeof(sim_out));
    mexMakeMemoryPersistent(out);
    
    return out;
}

void sim_opts_free(sim_opts *opts)
{
    mxFree(opts);
}

void sim_in_free(sim_in *in)
{
    mxFree(in);
}

void sim_out_free(sim_out *out)
{
    mxFree(out);
}

static double require_scalar(const mxArray *mem, const char *name)
{
    const mxArray *value = mxGetField(mem, 0, name);
    if (value == NULL || !mxIsNumeric(value) || mxIsComplex(value) ||
        mxGetNumberOfElements(value) != 1 || !mxIsFinite(mxGetScalar(value))) {
        mexErrMsgIdAndTxt("MATMPC:Integrator:Configuration", "Invalid integrator configuration field '%s'.", name);
    }
    return mxGetScalar(value);
}

static const double *require_matrix(const mxArray *mem, const char *name, size_t count)
{
    const mxArray *value = mxGetField(mem, 0, name);
    if (value == NULL || !mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) ||
        mxGetNumberOfElements(value) != count) {
        mexErrMsgIdAndTxt("MATMPC:Integrator:Configuration", "Invalid integrator configuration field '%s'.", name);
    }
    return mxGetPr(value);
}

static bool same_matrix(const double *current, const double *configured, size_t count)
{
    return count == 0 || memcmp(current, configured, count * sizeof(double)) == 0;
}

void sim_workspace_free(sim_workspace *workspace)
{
    if (workspace->erk != NULL)
        sim_erk_workspace_free(workspace->opts, workspace->erk);
    if (workspace->irk_ode != NULL)
        sim_irk_ode_workspace_free(workspace->opts, workspace->irk_ode);
    if (workspace->irk_dae != NULL)
        sim_irk_dae_workspace_free(workspace->opts, workspace->irk_dae);
    if (workspace->opts != NULL)
        sim_opts_free(workspace->opts);
    if (workspace->in != NULL)
        sim_in_free(workspace->in);
    if (workspace->out != NULL)
        sim_out_free(workspace->out);
    memset(workspace, 0, sizeof(*workspace));
}

void sim_workspace_prepare(sim_workspace *workspace, const mxArray *mem, int method,
                           size_t nx, size_t nu, size_t nz, bool forward, bool adjoint)
{
    if (method < 1 || method > 3 || (nz != 0) != (method == 3)) {
        mexErrMsgIdAndTxt("MATMPC:Integrator:Configuration", "Choose an integrator compatible with the model's algebraic states.");
    }
    if (require_scalar(mem, "nx") != (double)nx || require_scalar(mem, "nu") != (double)nu ||
        require_scalar(mem, "nz") != (double)nz) {
        mexErrMsgIdAndTxt("MATMPC:Integrator:Configuration", "Integrator and model dimensions must agree.");
    }
    const double stages_value = require_scalar(mem, "num_stages");
    const double steps_value = require_scalar(mem, "num_steps");
    const double h = require_scalar(mem, "h");
    if (stages_value < 1 || stages_value > INT_MAX || stages_value != (int)stages_value ||
        steps_value < 1 || steps_value > INT_MAX || steps_value != (int)steps_value || h <= 0) {
        mexErrMsgIdAndTxt("MATMPC:Integrator:Configuration", "Stage count, step count and step size must be positive.");
    }

    const size_t stages = (size_t)stages_value;
    const size_t steps = (size_t)steps_value;
    const double *A = require_matrix(mem, "A_tab", stages * stages);
    const double *B = require_matrix(mem, "B_tab", stages);
    const double *Sx = require_matrix(mem, "Sx", nx * nx);
    const double *Su = require_matrix(mem, "Su", nx * nu);
    const double *JFK = NULL;
    int newton_iter = 0;
    if (method != 1) {
        const double iterations = require_scalar(mem, "newton_iter");
        if (iterations < 1 || iterations > INT_MAX || iterations != (int)iterations) {
            mexErrMsgIdAndTxt("MATMPC:Integrator:Configuration", "Newton iteration count must be a positive integer.");
        }
        newton_iter = (int)iterations;
        JFK = require_matrix(mem, "JFK", stages * nx * (nx + nz));
    }

    sim_opts *opts = workspace->opts;
    const bool recreate = opts == NULL || workspace->in == NULL || workspace->out == NULL ||
                          workspace->method != method || opts->nx != nx || opts->nu != nu ||
                          opts->nz != nz || opts->num_stages != stages || opts->num_steps != steps ||
                          opts->forw_sens_flag != forward || opts->adj_sens_flag != adjoint;
    if (recreate) {
        sim_workspace_free(workspace);
        opts = sim_opts_create(mem);
        workspace->opts = opts;
        workspace->method = method;
        opts->forw_sens_flag = forward;
        opts->adj_sens_flag = adjoint;
        workspace->in = sim_in_create(opts);
        workspace->out = sim_out_create(opts);
    }
    opts->base_h = h;
    opts->h = h;

    switch (method) {
        case 1: {
            const bool allocate = workspace->erk == NULL;
            if (allocate)
                workspace->erk = sim_erk_workspace_create(opts);
            sim_erk_workspace *work = workspace->erk;
            if (allocate || !same_matrix(A, work->A, stages * stages) || !same_matrix(B, work->B, stages) ||
                !same_matrix(Sx, work->Sx, nx * nx) || !same_matrix(Su, work->Su, nx * nu)) {
                sim_erk_workspace_init(opts, mem, work);
            }
            break;
        }
        case 2: {
            const bool allocate = workspace->irk_ode == NULL;
            if (allocate)
                workspace->irk_ode = sim_irk_ode_workspace_create(opts);
            sim_irk_ode_workspace *work = workspace->irk_ode;
            if (allocate || !same_matrix(A, work->A, stages * stages) || !same_matrix(B, work->B, stages) ||
                !same_matrix(Sx, work->Sx, nx * nx) || !same_matrix(Su, work->Su, nx * nu) ||
                !same_matrix(JFK, work->JFK, stages * nx * nx) || work->newton_iter != newton_iter) {
                sim_irk_ode_workspace_init(opts, mem, work);
            }
            break;
        }
        case 3: {
            const bool allocate = workspace->irk_dae == NULL;
            if (allocate)
                workspace->irk_dae = sim_irk_dae_workspace_create(opts);
            sim_irk_dae_workspace *work = workspace->irk_dae;
            if (allocate || !same_matrix(A, work->A, stages * stages) || !same_matrix(B, work->B, stages) ||
                !same_matrix(Sx, work->Sx, nx * nx) || !same_matrix(Su, work->Su, nx * nu) ||
                !same_matrix(JFK, work->JFK, stages * nx * (nx + nz)) || work->newton_iter != newton_iter) {
                sim_irk_dae_workspace_init(opts, mem, work);
            }
            break;
        }
    }
}

const double *sim_grid_indices(const mxArray *mem, size_t horizon, bool required)
{
    const mxArray *flag = mxGetField(mem, 0, "nonuniform_grid");
    if (!required && (flag == NULL || mxGetScalar(flag) == 0.0))
        return NULL;

    const double *grid = require_matrix(mem, "index_T", horizon + 1);
    for (size_t i = 0; i < horizon; ++i) {
        const double interval = grid[i + 1] - grid[i];
        if (!mxIsFinite(interval) || interval <= 0.0) {
            mexErrMsgIdAndTxt("MATMPC:Integrator:Grid", "Shooting grid intervals must be finite and positive.");
        }
    }
    return grid;
}

void sim_set_stage_step(sim_opts *opts, const double *grid, size_t stage)
{
    opts->h = opts->base_h;
    if (grid == NULL)
        return;

    opts->h *= grid[stage + 1] - grid[stage];
    if (!mxIsFinite(opts->h) || opts->h <= 0.0) {
        mexErrMsgIdAndTxt("MATMPC:Integrator:Grid", "Shooting interval produces an invalid integration step size.");
    }
}
