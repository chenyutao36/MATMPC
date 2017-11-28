#include "mex.h"

#include "blasfeo_target.h"
#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_aux_ext_dep.h"
#include "blasfeo_d_blas.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 
    double *A = mxGetPr(prhs[0]);
    double *B = mxGetPr(prhs[1]);
    
    int m = mxGetM(prhs[0]);
    int k = mxGetN(prhs[0]);
    int n = mxGetN(prhs[1]);
    
    if (k!=mxGetM(prhs[1]))
        mexErrMsgTxt("Matrix dimension mismatch");
    
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double *C = mxGetPr(plhs[0]);
    
    struct d_strmat sA;
    struct d_strmat sB;
    struct d_strmat sC;
    
    //
    d_allocate_strmat(m, k, &sA);
    d_cvt_mat2strmat(m, k, A, m, &sA, 0, 0);
    d_allocate_strmat(k, n, &sB);
    d_cvt_mat2strmat(k, n, B, k, &sB, 0, 0);
    d_allocate_strmat(m, n, &sC);
    d_cvt_mat2strmat(m, n, C, m, &sC, 0, 0);
    
    dgemm_nn_libstr(m, n, k, 1.0, &sA, 0, 0, &sB, 0, 0, 0.0, &sC, 0, 0, &sC, 0, 0);
    
    d_cvt_strmat2mat(m, n, &sC, 0, 0, C, m);
    
    d_free_strmat(&sA);
    d_free_strmat(&sB);
    d_free_strmat(&sC);
}