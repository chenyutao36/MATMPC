#include "mex.h"
#include "string.h"

// for builtin blas
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

    double *Qh = mxGetPr( mxGetField(prhs[0], 0, "Q_h") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A_sens") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B_sens") );
    double *Cx = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *CxN = mxGetPr( mxGetField(prhs[0], 0, "CxN") );
    double *Cu = mxGetPr( mxGetField(prhs[0], 0, "Cu") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    
    double *du = mxGetPr(prhs[2]);
    double *mu_vec = mxGetPr(prhs[3]);
    
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    
    mwSize nz = nx+nu;
    
    double *z = mxGetPr( mxGetField(prhs[0], 0, "dz") );
    double *xN = mxGetPr( mxGetField(prhs[0], 0, "dxN") );
    double *lambda = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *muN = mxGetPr( mxGetField(prhs[0], 0, "muN_new") );
   
    memcpy(&mu[0], &mu_vec[0], nc*N*sizeof(double));
    memcpy(&muN[0], &mu_vec[N*nc], ncN*sizeof(double));
 
    int i;
    double *tmp;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
       
    memcpy(&z[0],&ds0[0], nx*sizeof(double)); 
    
    for (i=0;i<N;i++){        
        if (i<N-1){
            memcpy(&z[(i+1)*nz],&a[i*nx], nx*sizeof(double));
            tmp = A + i*nx*nx;
            dgemv(nTrans,&nx,&nx,&one_d,tmp,&nx,z+i*nz,&one_i,&one_d,z+(i+1)*nz,&one_i);
            tmp = B + i*nx*nu;
            dgemv(nTrans,&nx,&nu,&one_d,tmp,&nx,du+i*nu,&one_i,&one_d,z+(i+1)*nz,&one_i);
        }else{
            memcpy(&xN[0],&a[i*nx], nx*sizeof(double));
            tmp = A + i*nx*nx;
            dgemv(nTrans,&nx,&nx,&one_d,tmp,&nx,z+i*nz,&one_i,&one_d,xN,&one_i);
            tmp = B + i*nx*nu;
            dgemv(nTrans,&nx,&nu,&one_d,tmp,&nx,du+i*nu,&one_i,&one_d,xN,&one_i);
        }
                
        memcpy(&z[i*nz+nx],&du[i*nu], nu*sizeof(double));
    }
    
    tmp = Qh + N*nx*nx;
    memcpy(&lambda[N*nx],&gx[N*nx], nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,tmp,&nx,xN,&one_i,&one_d,lambda+N*nx,&one_i);
    
    if (ncN>0){
        dgemv(Trans,&ncN,&nx,&one_d,CxN,&ncN,mu_vec+N*nc,&one_i,&one_d,lambda+N*nx,&one_i);
    }
    for (i=N-1;i>-1;i--){
        memcpy(&lambda[i*nx],&gx[i*nx], nx*sizeof(double));
        tmp = Qh + i*nx*nx;
        dgemv(nTrans,&nx,&nx,&one_d,tmp,&nx,z+i*nz,&one_i,&one_d,lambda+i*nx,&one_i);
        tmp = S + i*nx*nu;
        dgemv(nTrans,&nx,&nu,&one_d,tmp,&nx,du+i*nu,&one_i,&one_d,lambda+i*nx,&one_i);
        tmp = A + i*nx*nx;
        dgemv(Trans,&nx,&nx,&one_d,tmp,&nx,lambda+(i+1)*nx,&one_i,&one_d,lambda+i*nx,&one_i);
        
        if (nc>0){
            tmp = Cx + i*nc*nx;
            dgemv(Trans,&nc,&nx,&one_d,tmp,&nc,mu_vec+i*nc,&one_i,&one_d,lambda+i*nx,&one_i);
        }
    }
    
}