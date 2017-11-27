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

#define MAX(a,b) (((a)>(b))?(a):(b))

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    double *lambda = mxGetPr(prhs[0]);
    double *mu = mxGetPr(prhs[1]);
    double *muN = mxGetPr(prhs[2]);
    double *ds0 = mxGetPr(prhs[3]);
    
    double *z = mxGetPr( mxGetField(prhs[4], 0, "z") );
    double *xN = mxGetPr( mxGetField(prhs[4], 0, "xN") );
    double *y = mxGetPr( mxGetField(prhs[4], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[4], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[4], 0, "od") );
    double *Q = mxGetPr( mxGetField(prhs[4], 0, "W") );
    double *QN = mxGetPr( mxGetField(prhs[4], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[4], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[4], 0, "ub") );
    double *lbN = mxGetPr( mxGetField(prhs[4], 0, "lbN") );
    double *ubN = mxGetPr( mxGetField(prhs[4], 0, "ubN") );
    double *x0 = mxGetPr( mxGetField(prhs[4], 0, "x0") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[5], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[5], 0, "nu") );
    mwSize np = mxGetScalar( mxGetField(prhs[5], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[5], 0, "ny") );
    mwSize nyN = mxGetScalar( mxGetField(prhs[5], 0, "nyN") );
    mwSize nc = mxGetScalar( mxGetField(prhs[5], 0, "nc") ); nc *= 2;
    mwSize ncN = mxGetScalar( mxGetField(prhs[5], 0, "ncN") ); ncN *=2;
    mwSize N = mxGetScalar( mxGetField(prhs[5], 0, "N") );
      
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    mwSize nw = N*nz+nx;
    mwSize neq = (N+1)*nx;
    mwSize nineq = N*nc+ncN;
    char *nTrans = "N", *Trans="T", *Norm="O";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
    
    double **vec_out = (double **) mxMalloc(3 * sizeof(double*));
    vec_out[1] = (double *)mxCalloc(nz, sizeof(double));
    vec_out[2] = (double *)mxCalloc(nz, sizeof(double));
    
        
    double **vec_in = (double **) mxMalloc(6 * sizeof(double*));
    vec_in[3] = Q;
    
    double *L = (double *) mxCalloc( nw, sizeof(double));
    double *a = (double *) mxCalloc( neq, sizeof(double));
    memcpy(&a[0], &ds0[0], nx*sizeof(double));
    double *c = (double *) mxCalloc( nineq, sizeof(double));
    
    double *work;
    double KKT=0,eq_res=0,ineq_res=0;
    
    for (i=0;i<N;i++){
        vec_in[0] = z+i*nz;
        vec_in[1] = od+i*np;
        vec_in[2] = y+i*ny;
        vec_in[4] = lambda+(i+1)*nx;
        vec_in[5] = mu+i*nc;
           
        vec_out[0] = L+i*nz;
        adj_Fun(vec_in, vec_out);
        
        for (j=0;j<nx;j++){
            vec_out[1][j] -= lambda[i*nx+j];
        }
        
        daxpy(&nz, &one_d, vec_out[1], &one_i, vec_out[0], &one_i);
        daxpy(&nz, &one_d, vec_out[2], &one_i, vec_out[0], &one_i);
        
        vec_out[0] = a+(i+1)*nx;
        F_Fun(vec_in, vec_out);
        if (i < N-1){
            for (j=0;j<nx;j++)
                vec_out[0][j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                vec_out[0][j] -= xN[j];
        }
        
        vec_out[0] = c + i*nc;
        ineq_Fun(vec_in, vec_out);
        for (j=0;j<nc/2;j++){
            vec_out[0][nc/2+j] = lb[j] - vec_out[0][j];
            vec_out[0][j] -= ub[j];
        }
    }
    vec_in[0] = xN;
    vec_in[1] = od+N*np;
    vec_in[2] = yN;
    vec_in[3] = QN;
    vec_in[4] = muN;
    
    vec_out[0] = L+N*nz;
    adjN_Fun(vec_in, vec_out);
    daxpy(&nx, &one_d, vec_out[1], &one_i, vec_out[0], &one_i);
    
    vec_out[0] = c + N*nc;
    ineqN_Fun(vec_in, vec_out);
    for (j=0;j<ncN/2;j++){
        vec_out[0][ncN/2+j] = lbN[j] - vec_out[0][j];
        vec_out[0][j] -= ub[j];
    }
         
    eq_res = dlange(Norm, &neq, &one_i, a, &one_i, work);
    KKT = dlange(Norm, &nw, &one_i, L, &one_i, work);
    for (i=0;i<nineq;i++){
        ineq_res += MAX(c[i],0);
    }
    
    plhs[0] = mxCreateDoubleScalar(eq_res); // eq_res
    plhs[1] = mxCreateDoubleScalar(ineq_res); // ineq_res
    plhs[2] = mxCreateDoubleScalar(KKT); // KKT
    
    mxFree(vec_in);
    mxFree(vec_out);
    mxFree(L);
    mxFree(a);
    mxFree(c);
}