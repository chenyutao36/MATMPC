#include "mex.h"
#include "condense_affine.h"
#include "condense_matrix.h"

#include <math.h>
#include <stdint.h>

enum output_index { HC, GC, CCX, CCG, LCC, UCC, LXC, UXC, OUTPUT_COUNT };
static const char *output_names[] = {"Hc", "gc", "Ccx", "Ccg", "lcc", "ucc", "lxc", "uxc"};
static double *scratch;
static size_t scratch_count;

static matmpc_condense_output output(double *data, size_t count)
{
    const matmpc_condense_output view = {data, count};
    return view;
}

static void free_scratch(void)
{
    mxFree(scratch);
    scratch = NULL;
    scratch_count = 0;
}

static matmpc_condense_output get_scratch(size_t count)
{
    if (scratch_count != count) {
        double *next = mxMalloc(count * sizeof(double));
        if (!scratch)
            mexAtExit(free_scratch);
        mexMakeMemoryPersistent(next);
        mxFree(scratch);
        scratch = next;
        scratch_count = count;
    }
    return output(scratch, scratch_count);
}

static size_t dimension(const mxArray *settings, const char *name)
{
    const mxArray *value = mxGetField(settings, 0, name);
    if (value && mxIsNumeric(value) && !mxIsComplex(value) && !mxIsSparse(value) &&
        mxGetNumberOfElements(value) == 1) {
        const double scalar = mxGetScalar(value);
        if (mxIsFinite(scalar) && scalar >= 0 && scalar < (double)PTRDIFF_MAX && floor(scalar) == scalar)
            return (size_t)scalar;
    }
    mexErrMsgIdAndTxt("MATMPC:Condensing:Dimensions", "Invalid settings.%s.", name);
    return 0;
}

static double *buffer(const mxArray *record, const char *name, size_t count)
{
    const mxArray *value = mxGetField(record, 0, name);
    if (!value || !mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) ||
        mxGetNumberOfElements(value) != count) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:Buffer", "mem.%s must be a real full double array with the expected size.", name);
    }
    return mxGetPr(value);
}

static matmpc_condense_buffer input(const mxArray *record, const char *name, size_t count)
{
    const matmpc_condense_buffer view = {buffer(record, name, count), count};
    return view;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2 || nlhs > 1)
        mexErrMsgIdAndTxt("MATMPC:Condensing:Arguments", "Call qp = Condensing(mem, settings).");
    const mxArray *mem = prhs[0];
    const mxArray *settings = prhs[1];
    if (!mxIsStruct(mem) || mxGetNumberOfElements(mem) != 1 ||
        !mxIsStruct(settings) || mxGetNumberOfElements(settings) != 1) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:Arguments", "mem and settings must be scalar structures.");
    }

    const matmpc_condense_dimensions dims = {
        dimension(settings, "nx"), dimension(settings, "nu"), dimension(settings, "N"),
        dimension(settings, "nbx"), dimension(settings, "nc"), dimension(settings, "ncN")
    };
    matmpc_condense_affine_sizes af;
    matmpc_condense_matrix_sizes mx;
    if (matmpc_condense_affine_get_sizes(&dims, &af) != MATMPC_CONDENSE_OK ||
        matmpc_condense_matrix_get_sizes(&dims, &mx) != MATMPC_CONDENSE_OK) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:Dimensions", "Invalid or overflowing QP dimensions.");
    }

    const matmpc_condense_matrix_inputs matrices = {
        input(mem, "Q", mx.q), input(mem, "S", mx.s), input(mem, "R", mx.r),
        input(mem, "A", mx.a), input(mem, "B", mx.b), input(mem, "Cx", mx.cx),
        input(mem, "Cgx", mx.cgx), input(mem, "Cgu", mx.cgu), input(mem, "CgN", mx.cgn)
    };
    const matmpc_condense_affine_inputs affine = {
        matrices.q, matrices.s, matrices.a, matrices.b, matrices.cx, matrices.cgx, matrices.cgn,
        input(mem, "gx", af.states), input(mem, "gu", af.controls), input(mem, "a", af.defect),
        input(mem, "ds0", dims.nx), input(mem, "lc", af.general), input(mem, "uc", af.general),
        input(mem, "lb_dx", af.state), input(mem, "ub_dx", af.state)
    };

    const size_t output_sizes[] = {mx.hc, af.controls, mx.ccx, mx.ccg,
                                   af.general, af.general, af.state, af.state};
    matmpc_condense_output outputs[OUTPUT_COUNT];
    mxArray *qp = NULL;
    if (nlhs == 1) {
        const size_t rows[] = {af.controls, af.controls, af.state, af.general,
                               af.general, af.general, af.state, af.state};
        const size_t cols[] = {af.controls, 1, af.controls, af.controls, 1, 1, 1, 1};
        qp = mxCreateStructMatrix(1, 1, OUTPUT_COUNT, output_names);
        for (int i = 0; i < OUTPUT_COUNT; ++i) {
            mxArray *value = mxCreateDoubleMatrix(rows[i], cols[i], mxREAL);
            mxSetFieldByNumber(qp, 0, i, value);
            outputs[i] = output(mxGetPr(value), output_sizes[i]);
        }
    } else {
        /* Preserve the former no-output MATLAB call for existing integrations. */
        for (int i = 0; i < OUTPUT_COUNT; ++i)
            outputs[i] = output(buffer(mem, output_names[i], output_sizes[i]), output_sizes[i]);
    }

    const matmpc_condense_matrix_outputs matrix_outputs = {outputs[HC], outputs[CCX], outputs[CCG]};
    const matmpc_condense_affine_outputs affine_outputs = {
        outputs[GC], outputs[LCC], outputs[UCC], outputs[LXC], outputs[UXC]
    };
    const matmpc_condense_output workspace = get_scratch(mx.workspace > af.workspace ? mx.workspace : af.workspace);
    if (matmpc_condense_matrix(&dims, &matrices, &matrix_outputs, workspace) != MATMPC_CONDENSE_OK ||
        matmpc_condense_affine(&dims, &affine, &affine_outputs, workspace) != MATMPC_CONDENSE_OK) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:Buffer", "Invalid buffers for QP condensation.");
    }
    if (nlhs == 1)
        plhs[0] = qp;
}
