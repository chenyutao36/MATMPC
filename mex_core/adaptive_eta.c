#include "mex.h"
#include "string.h"
#include <stdio.h>
#include <math.h>

// for builtin blas
#include "blas.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
   
    double *q_pri = mxGetPr( mxGetField(prhs[0], 0, "dz") );
    double *q_dual = mxGetPr( mxGetField(prhs[0], 0, "q_dual") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A_sens") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B_sens") );
    
    double *threshold_pri = mxGetPr( mxGetField(prhs[0], 0, "threshold_pri") );
    double *threshold_dual = mxGetPr( mxGetField(prhs[0], 0, "threshold_dual") );
    double alpha = mxGetScalar( mxGetField(prhs[0], 0, "alpha") );
    double beta = mxGetScalar( mxGetField(prhs[0], 0, "beta") );
    double c1 = mxGetScalar( mxGetField(prhs[0], 0, "c1") );
    double rho_cmon = mxGetScalar( mxGetField(prhs[0], 0, "rho_cmon") );
    double tol = mxGetScalar( mxGetField(prhs[0], 0, "tol") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );    
    mwSize np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    
    mwSize m = N*nx;
    mwSize nz = nx+nu;
    mwSize n = N*nz;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
    
    double *V_pri = (double *)mxMalloc( m*sizeof(double) );
    double *V_dual = (double *)mxMalloc( n*sizeof(double) );
    
    mwIndex i,j;
    
    for (i=0;i<N;i++){
        dgemv(nTrans, &nx, &nx, &one_d, A+i*nx*nx, &nx, q_pri+i*nz, &one_i, &zero, V_pri+i*nx, &one_i);            
        dgemv(nTrans, &nx, &nu, &one_d, B+i*nx*nu, &nx, q_pri+i*nz+nx, &one_i, &one_d, V_pri+i*nx, &one_i);
        
        dgemv(Trans, &nx, &nx, &one_d, A+i*nx*nx, &nx, q_dual+(i+1)*nx, &one_i, &zero, V_dual+i*nz, &one_i);            
        dgemv(Trans, &nx, &nu, &one_d, B+i*nx*nu, &nx, q_dual+(i+1)*nx, &one_i, &zero, V_dual+i*nz+nx, &one_i);
    }
    
    double norm_V_pri = dnrm2(&m, V_pri, &one_i);
    double norm_V_dual = dnrm2(&n, V_dual, &one_i);
    
    *threshold_pri = (sqrt(c1)*tol)/(2*alpha*rho_cmon*norm_V_pri);
    *threshold_dual = (sqrt(1-c1)*tol)/(2*beta*rho_cmon*norm_V_dual);
    
    mxFree(V_pri);
    mxFree(V_dual);
}