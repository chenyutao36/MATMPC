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
    double *Qh = mxGetPr(prhs[0]);
    double *S = mxGetPr(prhs[1]);    
    double *A = mxGetPr(prhs[2]);
    double *B = mxGetPr(prhs[3]);    
    double *Cx = mxGetPr(prhs[4]);
    double *CxN = mxGetPr(prhs[5]);
    double *a = mxGetPr(prhs[6]);
    double *gx = mxGetPr(prhs[7]);    
    double *du = mxGetPr(prhs[8]);
    double *ds0 = mxGetPr(prhs[9]);
    double *mu_vec = mxGetPr(prhs[10]);
    
    mwSize nx = mxGetScalar( mxGetField(prhs[11], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[11], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[11], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[11], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[11], 0, "N") );
    
    mwSize nz = nx+nu;
    
    plhs[0] = mxCreateDoubleMatrix(nz, N, mxREAL); //z
    plhs[1] = mxCreateDoubleMatrix(nx, 1, mxREAL); //xN
    plhs[2] = mxCreateDoubleMatrix(nx, N+1, mxREAL); //lambda
    plhs[3] = mxCreateDoubleMatrix(nc, N, mxREAL); //mu
    plhs[4] = mxCreateDoubleMatrix(ncN,1, mxREAL); //muN
    
    double *z = mxGetPr(plhs[0]);
    double *xN = mxGetPr(plhs[1]);
    double *lambda = mxGetPr(plhs[2]);
    double *mu = mxGetPr(plhs[3]);
    double *muN = mxGetPr(plhs[4]);
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