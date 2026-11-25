#include "mex.h"

#include "casadi_wrapper.h"
#include "sim_workspace.h"

#include <math.h>
#include <stdint.h>

static sim_workspace integrator;

static void exitFcn(void)
{
    sim_workspace_free(&integrator);
}

static size_t read_dimension(const mxArray *values, const char *name, size_t minimum, size_t maximum)
{
    const mxArray *field = mxGetField(values, 0, name);
    if (field != NULL && mxIsNumeric(field) && !mxIsComplex(field) && !mxIsSparse(field) &&
        mxGetNumberOfElements(field) == 1) {
        const double value = mxGetScalar(field);
        if (mxIsFinite(value) && value >= (double)minimum && value <= (double)maximum && floor(value) == value)
            return (size_t)value;
    }
    mexErrMsgIdAndTxt("MATMPC:SimulateSystem:Configuration", "Invalid simulation configuration field '%s'.", name);
    return 0;
}

static double *read_vector(const mxArray *value, size_t count, bool allow_placeholder, const char *name)
{
    const size_t length = mxGetNumberOfElements(value);
    const bool expected_length = length == count || (allow_placeholder && count == 0 && length == 1);
    const bool vector = length == 0 || (mxGetNumberOfDimensions(value) == 2 &&
                                       (mxGetM(value) == 1 || mxGetN(value) == 1));
    if (!mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) || !expected_length || !vector) {
        mexErrMsgIdAndTxt("MATMPC:SimulateSystem:Input", "Invalid real double vector for '%s'.", name);
    }
    return mxGetPr(value);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6 || nlhs > 2) {
        mexErrMsgIdAndTxt("MATMPC:SimulateSystem:Arguments",
                         "Call [xNext, zNext] = Simulate_System(x, u, z, parameters, mem, settings).");
    }
    const mxArray *mem = prhs[4];
    const mxArray *settings = prhs[5];
    if (!mxIsStruct(mem) || mxGetNumberOfElements(mem) != 1 ||
        !mxIsStruct(settings) || mxGetNumberOfElements(settings) != 1) {
        mexErrMsgIdAndTxt("MATMPC:SimulateSystem:Configuration", "mem and settings must be scalar structures.");
    }

    const size_t max_count = PTRDIFF_MAX / sizeof(double);
    const size_t nx = read_dimension(settings, "nx", 1, max_count);
    const size_t nu = read_dimension(settings, "nu", 0, max_count);
    const size_t nz = read_dimension(settings, "nz", 0, max_count);
    const int method = (int)read_dimension(mem, "sim_method", 1, 3);
    size_t model_dims[4];
    const int model_status = matmpc_model_dimensions(model_dims);
    const bool have_parameters = mxGetField(settings, 0, "np") != NULL;
    const size_t np = have_parameters ? read_dimension(settings, "np", 0, max_count) :
                     (model_status > 0 ? model_dims[2] : mxGetNumberOfElements(prhs[3]));
    if (model_status < 0 || (model_status > 0 &&
            (nx != model_dims[0] || nu != model_dims[1] || np != model_dims[2] || nz != model_dims[3]))) {
        mexErrMsgIdAndTxt("MATMPC:SimulateSystem:ModelDimensions",
                         "Simulation dimensions do not match the compiled model.");
    }

    sim_in in = {0};
    in.x = read_vector(prhs[0], nx, false, "x");
    in.u = read_vector(prhs[1], nu, false, "u");
    in.z = read_vector(prhs[2], nz, true, "z");
    in.p = read_vector(prhs[3], np, true, "parameters");

    if (integrator.opts == NULL)
        mexAtExit(exitFcn);
    sim_workspace_prepare(&integrator, mem, method, nx, nu, nz, false, false);

    mxArray *state = mxCreateDoubleMatrix(nx, 1, mxREAL);
    mxArray *algebraic = (nz != 0 || nlhs == 2) ? mxCreateDoubleMatrix(nz, 1, mxREAL) : NULL;
    sim_out out = {0};
    out.xn = mxGetPr(state);
    if (algebraic != NULL)
        out.zn = mxGetPr(algebraic);

    switch (method) {
        case 1:
            sim_erk(&in, &out, integrator.opts, integrator.erk);
            break;
        case 2:
            sim_irk_ode(&in, &out, integrator.opts, integrator.irk_ode);
            break;
        case 3:
            sim_irk_dae(&in, &out, integrator.opts, integrator.irk_dae);
            break;
    }

    /* MATLAB permits a first result to become ans when nlhs is zero. */
    plhs[0] = state;
    if (nlhs == 2)
        plhs[1] = algebraic;
    else if (algebraic != NULL)
        mxDestroyArray(algebraic);
}
