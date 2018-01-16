#include "mex.h"
#include "string.h"

// for builtin blas
#include "lapack.h"
#include "blas.h"

// for openblas
// #include "f77blas.h"
// #if !defined(_WIN32)
// #define dgemm dgemm_
// #define dgemv dgemv_
// #endif

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *z = mxGetPr( mxGetField(prhs[1], 0, "z") );
    double *xN = mxGetPr( mxGetField(prhs[1], 0, "xN") );
    double *lambda = mxGetPr( mxGetField(prhs[1], 0, "lambda") );
    double *mu = mxGetPr( mxGetField(prhs[1], 0, "mu") );
    double *muN = mxGetPr( mxGetField(prhs[1], 0, "muN") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[2], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[2], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[2], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[2], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[2], 0, "N") );
    
    mwSize nz = nx+nu;
    mwSize nw = N*nz;
    mwSize neq = (N+1)*nx;
    mwSize nineq = N*nc;
    
    double one_d = 1.0;
    mwSignedIndex one_i = 1;
    
    double *dz = mxGetPr( mxGetField(prhs[0], 0, "dz") );
    double *dxN = mxGetPr( mxGetField(prhs[0], 0, "dxN") );
    double *lambda_new = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu_new = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *muN_new = mxGetPr( mxGetField(prhs[0], 0, "muN_new") );
    
    double *q = mxGetPr( mxGetField(prhs[0], 0, "q") );
    daxpy(&nw, &one_d, dz, &one_i, q, &one_i);
    
    double alpha = 1.0;
    double inc = 1.0 - alpha;
     
    daxpy(&nw, &alpha, dz, &one_i, z, &one_i); 
    daxpy(&nx, &alpha, dxN, &one_i, xN, &one_i);
    
    dscal(&neq, &inc, lambda, &one_i);
    daxpy(&neq, &alpha, lambda_new, &one_i, lambda, &one_i);
    
    if (nc>0){
        dscal(&nineq, &inc, mu, &one_i);
        daxpy(&nineq, &alpha, mu_new, &one_i, mu, &one_i);
    }
    
    if (ncN>0){
        dscal(&ncN, &inc, muN, &one_i);
        daxpy(&ncN, &alpha, muN_new, &one_i, muN, &one_i);
    }
}
