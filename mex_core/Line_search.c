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
    double *z_old = mxGetPr( mxGetField(prhs[6], 0, "z") );
    double *xN_old = mxGetPr( mxGetField(prhs[6], 0, "xN") );
    double *lambda_old = mxGetPr( mxGetField(prhs[6], 0, "lambda") );
    double *mu_old = mxGetPr( mxGetField(prhs[6], 0, "mu") );
    double *muN_old = mxGetPr( mxGetField(prhs[6], 0, "muN") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[7], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[7], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[7], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[7], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[7], 0, "N") );
    
    mwSize nz = nx+nu;
    mwSize nw = N*nz;
    mwSize neq = (N+1)*nx;
    mwSize nineq = N*nc;
    
    mwSignedIndex one_i = 1;
       
    double *dz = mxGetPr(prhs[0]);
    double *dxN = mxGetPr(prhs[1]);
    double *lambda = mxGetPr(prhs[2]);
    double *mu = mxGetPr(prhs[3]);
    double *muN = mxGetPr(prhs[4]);
    
    double alpha = mxGetScalar(prhs[5]);
    double inc = 1.0 - alpha;
    
    plhs[0] = mxCreateDoubleMatrix(nz, N, mxREAL); // z
    plhs[1] = mxCreateDoubleMatrix(nx, 1, mxREAL); // xN
    plhs[2] = mxCreateDoubleMatrix(nx, N+1, mxREAL); // lambda
    plhs[3] = mxCreateDoubleMatrix(nc, N, mxREAL); // mu
    plhs[4] = mxCreateDoubleMatrix(ncN, 1, mxREAL); // muN
    
    double *z = mxGetPr(plhs[0]);
    double *xN = mxGetPr(plhs[1]);
    double *lambda_new = mxGetPr(plhs[2]);
    double *mu_new = mxGetPr(plhs[3]);
    double *muN_new = mxGetPr(plhs[4]);
    
    memcpy(&z[0], &z_old[0], nw*sizeof(double));
    daxpy(&nw, &alpha, dz, &one_i, z, &one_i);
    
    memcpy(&xN[0], &xN_old[0], nx*sizeof(double));
    daxpy(&nx, &alpha, dxN, &one_i, xN, &one_i);
    
    dscal(&neq, &alpha, lambda, &one_i);
    memcpy(&lambda_new[0], &lambda[0], neq*sizeof(double));
    daxpy(&neq, &inc, lambda_old, &one_i, lambda_new, &one_i);
    
    dscal(&nineq, &alpha, mu, &one_i);
    memcpy(&mu_new[0], &mu[0], nineq*sizeof(double));
    daxpy(&nineq, &inc, mu_old, &one_i, mu_new, &one_i);
    
    dscal(&ncN, &alpha, muN, &one_i);
    memcpy(&muN_new[0], &muN[0], ncN*sizeof(double));
    daxpy(&ncN, &inc, muN_old, &one_i, muN_new, &one_i);
    
}
