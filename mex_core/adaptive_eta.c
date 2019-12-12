 #include "mex.h"	
#include "string.h"	
#include <stdio.h>	
#include <math.h>	
	
#include "blas.h"	
	
static double *dy=NULL;	
static bool mem_alloc = false;	
void exitFcn(){	
    if (mem_alloc){	
        mxFree(dy);	
    }	
}	
	
void	
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])	
{	
   	
    double *x_pri = mxGetPr( mxGetField(prhs[0], 0, "dx") );	
    double *u_pri = mxGetPr( mxGetField(prhs[0], 0, "du") );	
    double *q_dual = mxGetPr( mxGetField(prhs[0], 0, "q_dual") );	
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );	
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );	
    double *dmu = mxGetPr( mxGetField(prhs[0], 0, "dmu") );	
    	
    double *threshold_pri = mxGetPr( mxGetField(prhs[0], 0, "threshold_pri") );	
    double *threshold_dual = mxGetPr( mxGetField(prhs[0], 0, "threshold_dual") );	
    double *tol = mxGetPr( mxGetField(prhs[0], 0, "tol") );	
    	
    double alpha = mxGetScalar( mxGetField(prhs[0], 0, "alpha_cmon") );	
    double beta = mxGetScalar( mxGetField(prhs[0], 0, "beta_cmon") );	
    double c1 = mxGetScalar( mxGetField(prhs[0], 0, "c1") );	
    double rho_cmon = mxGetScalar( mxGetField(prhs[0], 0, "rho_cmon") );
    double gamma = mxGetScalar( mxGetField(prhs[0], 0, "gamma") );
    double tol_abs = mxGetScalar( mxGetField(prhs[0], 0, "tol_abs") );	
    double tol_ref = mxGetScalar( mxGetField(prhs[0], 0, "tol_ref") );
    
    double *shift_x = mxGetPr( mxGetField(prhs[0], 0, "shift_x") );
    double *shift_u = mxGetPr( mxGetField(prhs[0], 0, "shift_u") );
    double *V_pri = mxGetPr( mxGetField(prhs[0], 0, "V_pri") );
    double *V_dual = mxGetPr( mxGetField(prhs[0], 0, "V_dual") );
        	
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );	
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );	
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );	
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );	
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );	
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );    	
    size_t np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;	
    size_t ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );	
      	
    size_t nx_t = (N+1)*nx;
    size_t nu_t = N*nu;
    size_t m = N*nx;	
    size_t nv = nx+nu;	
    size_t n = N*nv;	
    size_t ns = (n+nx)+(m+nx)+(N*nu)+(N*nbx)+(N*nc+ncN);	
    char *nTrans = "N", *Trans="T";	
    double one_d = 1.0, zero = 0.0, minus_one_d = 1.0;	
    size_t one_i = 1;	
    	
    if (!mem_alloc){	
        dy = (double *)mxMalloc( ns*sizeof(double) );	
        mexMakeMemoryPersistent(dy);	
        	
        mem_alloc = true;     	
        mexAtExit(exitFcn);	
    }	
    	
    int i,j;	
    
    daxpy(&nx_t, &one_d, shift_x, &one_i, x_pri, &one_i);  
    daxpy(&nu_t, &one_d, shift_u, &one_i, u_pri, &one_i);  
    	
    for (i=0;i<N;i++){	
        dgemv(nTrans, &nx, &nx, &one_d, A+i*nx*nx, &nx, x_pri+i*nx, &one_i, &zero, V_pri+i*nx, &one_i);            	
        dgemv(nTrans, &nx, &nu, &one_d, B+i*nx*nu, &nx, u_pri+i*nu, &one_i, &one_d, V_pri+i*nx, &one_i);	
        	
        dgemv(Trans, &nx, &nx, &one_d, A+i*nx*nx, &nx, q_dual+(i+1)*nx, &one_i, &zero, V_dual+i*nv, &one_i);            	
        dgemv(Trans, &nx, &nu, &one_d, B+i*nx*nu, &nx, q_dual+(i+1)*nx, &one_i, &zero, V_dual+i*nv+nx, &one_i);	
        	
    }	
        	
    double norm_V_pri = dnrm2(&m, V_pri, &one_i);	
    double norm_V_dual = dnrm2(&n, V_dual, &one_i);	
        	
    memcpy(dy, x_pri, (N+1)*nx*sizeof(double));	
    memcpy(dy+(N+1)*nx, u_pri, N*nu*sizeof(double));	
    memcpy(dy+(N+1)*nx+N*nu, q_dual, (N+1)*nx*sizeof(double));	
    memcpy(dy+(N+1)*nx+N*nu+(N+1)*nx, dmu, (N*nu+N*nbx+N*nc+ncN)*sizeof(double));	
    *tol = tol_abs * sqrt(ns) + tol_ref * dnrm2(&ns, dy, &one_i); 	
    	
    *threshold_pri = (gamma*sqrt(c1)*(*tol))/(2*alpha*rho_cmon*norm_V_pri);	
    *threshold_dual = (gamma*sqrt(1-c1)*(*tol))/(beta*rho_cmon*norm_V_dual);	
    	
    	
} 