#include "mex.h"
#include "recover.h"

#include <math.h>
#include <stdint.h>

static int is_full_real_double(const mxArray *value)
{
    return value != NULL && mxIsDouble(value) && !mxIsComplex(value) && !mxIsSparse(value);
}

static size_t read_dimension(const mxArray *settings, const char *name)
{
    const mxArray *field = mxGetField(settings, 0, name);
    double value;

    if (!is_full_real_double(field) || mxGetNumberOfElements(field) != 1) {
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidDimensions",
                         "settings.%s must be a full real double scalar.", name);
    }
    value = mxGetScalar(field);
    if (!mxIsFinite(value) || value < 0.0 || value >= (double)SIZE_MAX || floor(value) != value) {
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidDimensions",
                         "settings.%s must be a representable nonnegative integer.", name);
    }
    return (size_t)value;
}

static matmpc_recover_buffer read_buffer(const mxArray *mem, const char *name,
                                         size_t count, size_t rows, int allow_vector)
{
    const mxArray *field = mxGetField(mem, 0, name);
    matmpc_recover_buffer buffer;

    if (!is_full_real_double(field) || mxGetNumberOfDimensions(field) != 2 ||
        mxGetNumberOfElements(field) != count ||
        (count != 0 && mxGetM(field) != rows &&
         !(allow_vector && (mxGetM(field) == 1 || mxGetN(field) == 1)))) {
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidBuffer",
                         "mem.%s must be a full real double matrix with the expected size.", name);
    }

    buffer.data = mxGetPr(field);
    buffer.count = count;
    return buffer;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const mxArray *mem;
    const mxArray *settings;
    matmpc_recover_dimensions dimensions;
    matmpc_recover_sizes sizes;
    matmpc_recover_inputs inputs;
    mxArray *dx;
    mxArray *lambda;
    matmpc_recover_status status;

    if (nrhs != 2 || nlhs != 2) {
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidArgumentCount",
                         "Use [dx, lambda_new] = Recover(mem, settings).");
    }
    mem = prhs[0];
    settings = prhs[1];
    if (!mxIsStruct(mem) || mxGetNumberOfElements(mem) != 1 ||
        !mxIsStruct(settings) || mxGetNumberOfElements(settings) != 1) {
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidStructure",
                         "mem and settings must be scalar structures.");
    }

    dimensions.nx = read_dimension(settings, "nx");
    dimensions.nu = read_dimension(settings, "nu");
    dimensions.horizon = read_dimension(settings, "N");
    dimensions.state_bounds = read_dimension(settings, "nbx");
    dimensions.path_constraints = read_dimension(settings, "nc");
    dimensions.terminal_constraints = read_dimension(settings, "ncN");
    if (matmpc_recover_get_sizes(&dimensions, &sizes) != MATMPC_RECOVER_OK) {
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidDimensions",
                         "The state dimension and horizon must be positive, and buffer sizes must be representable.");
    }

    inputs.q = read_buffer(mem, "Q", sizes.q, dimensions.nx, 0);
    inputs.s = read_buffer(mem, "S", sizes.s, dimensions.nx, 0);
    inputs.a = read_buffer(mem, "A", sizes.a, dimensions.nx, 0);
    inputs.b = read_buffer(mem, "B", sizes.s, dimensions.nx, 0);
    inputs.cx = read_buffer(mem, "Cx", sizes.cx, dimensions.state_bounds, 0);
    inputs.cgx = read_buffer(mem, "Cgx", sizes.cgx, dimensions.path_constraints, 0);
    inputs.cgn = read_buffer(mem, "CgN", sizes.cgn, dimensions.terminal_constraints, 0);
    inputs.gx = read_buffer(mem, "gx", sizes.dx, dimensions.nx, 0);
    inputs.defect = read_buffer(mem, "a", sizes.defect, dimensions.nx, 0);
    inputs.ds0 = read_buffer(mem, "ds0", dimensions.nx, dimensions.nx, 1);
    inputs.du = read_buffer(mem, "du", sizes.du, dimensions.nu, 1);
    inputs.mu = read_buffer(mem, "mu_new", sizes.mu, sizes.mu, 1);
    inputs.mu_x = read_buffer(mem, "mu_x_new", sizes.mu_x, sizes.mu_x, 1);

    dx = mxCreateDoubleMatrix(dimensions.nx, dimensions.horizon + 1, mxREAL);
    lambda = mxCreateDoubleMatrix(dimensions.nx, dimensions.horizon + 1, mxREAL);
    status = matmpc_recover(&dimensions, &inputs, mxGetPr(dx), sizes.dx, mxGetPr(lambda), sizes.dx);
    if (status != MATMPC_RECOVER_OK) {
        mxDestroyArray(dx);
        mxDestroyArray(lambda);
        mexErrMsgIdAndTxt("MATMPC:Recover:InvalidBuffer", "Unable to recover the state increments and costates.");
    }

    plhs[0] = dx;
    plhs[1] = lambda;
}
