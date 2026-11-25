#include "mex.h"

/* Test-only stand-in for a legacy MEX gateway that overwrites array storage. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    (void)plhs;
    if (nlhs != 0 || nrhs != 3 || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        mxIsSparse(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) ||
        mxGetNumberOfElements(prhs[1]) != 1 || mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("MATMPCTest:OverwriteMatrix", "Expected a real full double matrix, index and value.");
    }
    const double index = mxGetScalar(prhs[1]);
    if (!(index >= 1 && index <= (double)mxGetNumberOfElements(prhs[0])) || index != (double)(size_t)index) {
        mexErrMsgIdAndTxt("MATMPCTest:OverwriteMatrix", "Index out of bounds.");
    }
    mxGetPr(prhs[0])[(size_t)index - 1] = mxGetScalar(prhs[2]);
}
