#include "mex.h"
#include "string.h"

#include "blas.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *Cgx = mxGetPr( mxGetField(prhs[0], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[0], 0, "Cgu") );
    double *CgN = mxGetPr( mxGetField(prhs[0], 0, "CgN") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
        
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
       
    double *dx = mxGetPr( mxGetField(prhs[0], 0, "dx") );
    double *du = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *lambda = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *mu_x = mxGetPr( mxGetField(prhs[0], 0, "mu_x_new") );
   
    int i,j;
        
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    size_t one_i = 1;
       
    memcpy(dx, ds0, nx*sizeof(double)); 
    
    for (i=0;i<N;i++){        
        memcpy(dx+(i+1)*nx, a+i*nx, nx*sizeof(double));          
        dgemv(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,dx+i*nx,&one_i,&one_d,dx+(i+1)*nx,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,B+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,dx+(i+1)*nx,&one_i);
    }
    
    memcpy(lambda+N*nx, gx+N*nx, nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,Q+N*nx*nx,&nx,dx+N*nx,&one_i,&one_d,lambda+N*nx,&one_i);
    
    if (ncN>0)
        dgemv(Trans,&ncN,&nx,&one_d,CgN,&ncN,mu+N*nc,&one_i,&one_d,lambda+N*nx,&one_i);
    if (nbx>0)
        dgemv(Trans,&nbx,&nx,&one_d,Cx,&nbx,mu_x+(N-1)*nbx,&one_i,&one_d,lambda+N*nx,&one_i);
    
    for (i=N-1;i>0;i--){
        memcpy(lambda+i*nx,gx+i*nx, nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,dx+i*nx,&one_i,&one_d,lambda+i*nx,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,S+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,lambda+i*nx,&one_i);
        dgemv(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,lambda+(i+1)*nx,&one_i,&one_d,lambda+i*nx,&one_i);
        
        if (nc>0)
            dgemv(Trans,&nc,&nx,&one_d,Cgx+i*nc*nx,&nc,mu+i*nc,&one_i,&one_d,lambda+i*nx,&one_i);
        if (nbx>0)
            dgemv(Trans,&nbx,&nx,&one_d,Cx,&nbx,mu_x+(i-1)*nbx,&one_i,&one_d,lambda+i*nx,&one_i);
    }
    
    // lambda_0 has a different sign
    for (j=0;j<nx;j++)
        lambda[j] = -1.0*gx[j];
    dgemv(nTrans,&nx,&nx,&minus_one_d,Q,&nx,dx,&one_i,&one_d,lambda,&one_i);
    dgemv(nTrans,&nx,&nu,&minus_one_d,S+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,lambda+i*nx,&one_i);
    dgemv(Trans,&nx,&nx,&minus_one_d,A+i*nx*nx,&nx,lambda+(i+1)*nx,&one_i,&one_d,lambda+i*nx,&one_i);
        
    if (nc>0)
        dgemv(Trans,&nc,&nx,&minus_one_d,Cgx,&nc,mu,&one_i,&one_d,lambda,&one_i);  
      
}