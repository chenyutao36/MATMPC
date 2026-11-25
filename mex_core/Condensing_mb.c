#include "mex.h"
#include "condense_mb.h"

#include <math.h>
#include <stdint.h>

enum output_index { HC, GC, CCX, CCG, LCC, UCC, LXC, UXC, OUTPUT_COUNT };
static const char *output_names[] = {"Hc_r", "gc_r", "Ccx_r", "Ccg_r", "lcc", "ucc", "lxc", "uxc"};
static double *scratch;
static size_t scratch_count;

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
        if (scratch == NULL) {
            mexAtExit(free_scratch);
        }
        mexMakeMemoryPersistent(next);
        mxFree(scratch);
        scratch = next;
        scratch_count = count;
    }
    const matmpc_condense_output workspace = {scratch, scratch_count};
    return workspace;
}

static size_t dimension(const mxArray *record, const char *name)
{
    const mxArray *value = mxGetField(record, 0, name);
    if (value != NULL && mxIsNumeric(value) && !mxIsComplex(value) && !mxIsSparse(value) &&
        mxGetNumberOfElements(value) == 1) {
        const double scalar = mxGetScalar(value);
        if (mxIsFinite(scalar) && scalar >= 0 && scalar < (double)PTRDIFF_MAX && floor(scalar) == scalar) {
            return (size_t)scalar;
        }
    }
    mexErrMsgIdAndTxt("MATMPC:CondensingMB:Dimensions", "Invalid dimension field %s.", name);
    return 0;
}

static double *buffer(const mxArray *record, const char *name, size_t count)
{
    const mxArray *value = mxGetField(record, 0, name);
    if (value == NULL || !mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) ||
        mxGetNumberOfElements(value) != count) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Buffer", "mem.%s must be a real full double array with the expected size.", name);
    }
    return mxGetPr(value);
}

static matmpc_condense_buffer input(const mxArray *record, const char *name, size_t count)
{
    const matmpc_condense_buffer view = {buffer(record, name, count), count};
    return view;
}

static size_t *read_boundaries(const mxArray *mem, size_t count, size_t horizon)
{
    const double *raw = buffer(mem, "index_T", count);
    if (count > SIZE_MAX / sizeof(size_t)) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Dimensions", "The block partition is too large.");
    }
    size_t *boundaries = mxMalloc(count * sizeof(size_t));
    for (size_t i = 0; i < count; ++i) {
        if (!mxIsFinite(raw[i]) || raw[i] < 0 || raw[i] > (double)horizon || floor(raw[i]) != raw[i]) {
            mxFree(boundaries);
            mexErrMsgIdAndTxt("MATMPC:CondensingMB:Dimensions", "mem.index_T must contain integer shooting indices.");
        }
        boundaries[i] = (size_t)raw[i];
    }
    return boundaries;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2 || nlhs > 1) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Arguments", "Call qp = Condensing_mb(mem, settings).");
    }
    const mxArray *mem = prhs[0];
    const mxArray *settings = prhs[1];
    if (!mxIsStruct(mem) || mxGetNumberOfElements(mem) != 1 ||
        !mxIsStruct(settings) || mxGetNumberOfElements(settings) != 1) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Arguments", "mem and settings must be scalar structures.");
    }

    const matmpc_condense_dimensions dims = {
        dimension(settings, "nx"), dimension(settings, "nu"), dimension(settings, "N"),
        dimension(settings, "nbx"), dimension(settings, "nc"), dimension(settings, "ncN")
    };
    const size_t blocks = dimension(mem, "r");
    matmpc_condense_mb_sizes sizes;
    if (matmpc_condense_mb_get_sizes(&dims, blocks, &sizes) != MATMPC_CONDENSE_OK) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Dimensions", "Invalid or overflowing move-blocking dimensions.");
    }

    const matmpc_condense_affine_sizes *af = &sizes.affine;
    const matmpc_condense_mb_inputs inputs = {
        {
            input(mem, "Q", af->q), input(mem, "S", af->s),
            input(mem, "A", af->a), input(mem, "B", af->s), input(mem, "Cx", af->cx),
            input(mem, "Cgx", af->cgx), input(mem, "CgN", af->cgn), input(mem, "gx", af->states),
            input(mem, "gu_block", sizes.controls), input(mem, "a", af->defect), input(mem, "ds0", dims.nx),
            input(mem, "lc", af->general), input(mem, "uc", af->general),
            input(mem, "lb_dx", af->state), input(mem, "ub_dx", af->state)
        },
        input(mem, "R_block", sizes.r_block), input(mem, "Cgu", sizes.cgu)
    };

    const size_t counts[] = {sizes.hc, sizes.controls, sizes.ccx, sizes.ccg,
                             af->general, af->general, af->state, af->state};
    matmpc_condense_output fields[OUTPUT_COUNT];
    mxArray *qp = NULL;
    if (nlhs == 1) {
        const size_t rows[] = {sizes.controls, sizes.controls, af->state, af->general,
                               af->general, af->general, af->state, af->state};
        const size_t cols[] = {sizes.controls, 1, sizes.controls, sizes.controls, 1, 1, 1, 1};
        qp = mxCreateStructMatrix(1, 1, OUTPUT_COUNT, output_names);
        for (int i = 0; i < OUTPUT_COUNT; ++i) {
            mxArray *value = mxCreateDoubleMatrix(rows[i], cols[i], mxREAL);
            mxSetFieldByNumber(qp, 0, i, value);
            fields[i].data = mxGetPr(value);
            fields[i].count = counts[i];
        }
    } else {
        /* Retain the original no-output entry point for older MATLAB callers. */
        for (int i = 0; i < OUTPUT_COUNT; ++i) {
            fields[i].data = buffer(mem, output_names[i], counts[i]);
            fields[i].count = counts[i];
        }
    }

    const matmpc_condense_mb_outputs outputs = {
        {fields[HC], fields[CCX], fields[CCG]},
        {fields[GC], fields[LCC], fields[UCC], fields[LXC], fields[UXC]}
    };
    size_t *boundaries = read_boundaries(mem, blocks + 1, dims.horizon);
    const matmpc_condense_status status = matmpc_condense_mb(&dims, boundaries, blocks + 1,
                                                           &inputs, &outputs, get_scratch(sizes.workspace));
    mxFree(boundaries);
    if (status == MATMPC_CONDENSE_INVALID_DIMENSIONS) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Dimensions", "Invalid move-blocking partition.");
    }
    if (status != MATMPC_CONDENSE_OK) {
        mexErrMsgIdAndTxt("MATMPC:CondensingMB:Buffer", "Invalid buffers for move-blocking condensation.");
    }
    if (nlhs == 1) {
        plhs[0] = qp;
    }
}
