

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
    double *z = mxGetPr( mxGetField(prhs[0], 0, "z") );
    double *xN = mxGetPr( mxGetField(prhs[0], 0, "xN") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *QN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub") );
    double *lbN = mxGetPr( mxGetField(prhs[0], 0, "lbN") );
    double *ubN = mxGetPr( mxGetField(prhs[0], 0, "ubN") );
    double *x0 = mxGetPr( mxGetField(prhs[0], 0, "x0") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    mwSize nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); nc *= 2;
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") ); ncN *=2;
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    
    plhs[0] = mxCreateCellMatrix(1, N+1); // Qh
    plhs[1] = mxCreateCellMatrix(1, N); // S
    plhs[2] = mxCreateCellMatrix(1, N); // R
    plhs[3] = mxCreateCellMatrix(1, N); // A
    plhs[4] = mxCreateCellMatrix(1, N); // B
    plhs[5] = mxCreateCellMatrix(1, N+1); // Cx
    plhs[6] = mxCreateCellMatrix(1, N); // Cu
    plhs[7] = mxCreateDoubleMatrix(nx, N+1, mxREAL); //gx
    plhs[8] = mxCreateDoubleMatrix(nu, N, mxREAL); //gu
    plhs[9] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL); //c
    plhs[10] = mxCreateDoubleMatrix(nx,N, mxREAL); //a
    plhs[11] = mxCreateDoubleMatrix(nx,1,mxREAL); // ds0
    
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
      
    double *gx = mxGetPr(plhs[7]);
    double *gu = mxGetPr(plhs[8]);
    double *c = mxGetPr(plhs[9]);
    double *a = mxGetPr(plhs[10]);
    double *ds0 = mxGetPr(plhs[11]); 
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - z[i];
    
    // allocate memory
    mxArray ***Sens = (mxArray ***) mxMalloc(N * sizeof(mxArray**));
    mxArray ***Jac = (mxArray ***) mxMalloc((N+1) * sizeof(mxArray**));
    mxArray ***Cons = (mxArray ***) mxMalloc((N+1) * sizeof(mxArray**));
    mxArray **Qh = (mxArray **) mxMalloc((N+1) * sizeof(mxArray*));
    mxArray **S = (mxArray **) mxMalloc(N * sizeof(mxArray*));
    mxArray **R = (mxArray **) mxMalloc(N * sizeof(mxArray*));
    for(i=0;i<N;i++){
        Sens[i] = (mxArray **) mxMalloc(2 * sizeof(mxArray*));
        Sens[i][0] = mxCreateDoubleMatrix(nx, nx, mxREAL);
        Sens[i][1] = mxCreateDoubleMatrix(nx, nu, mxREAL);
        
        Jac[i] = (mxArray **) mxMalloc(2 * sizeof(mxArray*));
        Jac[i][0] = mxCreateDoubleMatrix(ny, nx, mxREAL);
        Jac[i][1] = mxCreateDoubleMatrix(ny, nu, mxREAL);
        
        Qh[i] = mxCreateDoubleMatrix(nx, nx, mxREAL);
        S[i] = mxCreateDoubleMatrix(nx, nu, mxREAL);
        R[i] = mxCreateDoubleMatrix(nu, nu, mxREAL);
        
        Cons[i] = (mxArray **) mxMalloc(2 * sizeof(mxArray*));
        Cons[i][0] = mxCreateDoubleMatrix(nc, nx, mxREAL);
        Cons[i][1] = mxCreateDoubleMatrix(nc, nu, mxREAL);
    }
    Jac[N] = (mxArray **) mxMalloc(sizeof(mxArray*));
    Jac[N][0] = mxCreateDoubleMatrix(nyN, nx, mxREAL);
    Qh[N] = mxCreateDoubleMatrix(nx, nx, mxREAL);
    Cons[N] = (mxArray **) mxMalloc(sizeof(mxArray*));
    Cons[N][0] = mxCreateDoubleMatrix(ncN, nx, mxREAL);
    
    double **vec_out = (double **) mxMalloc(2 * sizeof(double*));    
    double **vec_in = (double **) mxMalloc(4 * sizeof(double*));
    vec_in[3] = Q;
     
    // start loop
    for(i=0;i<N;i++){
        vec_in[0] = z+i*nz;
        vec_in[1] = od+i*np;
        vec_in[2] = y+i*ny;
        
        // integration
        vec_out[0] = a+i*nx;
        F_Fun(vec_in, vec_out);
        if (i < N-1){
            for (j=0;j<nx;j++)
                vec_out[0][j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                vec_out[0][j] -= xN[j];
        }
            
        // sensitiities
        D_Fun(vec_in, Sens[i]);  
        mxSetCell(plhs[3], i, Sens[i][0]);
        mxSetCell(plhs[4], i, Sens[i][1]);
        
        // Hessian
        Ji_Fun(vec_in, Jac[i]);
        dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, mxGetPr(Jac[i][0]), &ny, mxGetPr(Jac[i][0]), &ny, &zero, mxGetPr(Qh[i]), &nx);
        mxSetCell(plhs[0], i, Qh[i]);
        dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, mxGetPr(Jac[i][0]), &ny, mxGetPr(Jac[i][1]), &ny, &zero, mxGetPr(S[i]), &nx);
        mxSetCell(plhs[1], i, S[i]);
        dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, mxGetPr(Jac[i][1]), &ny, mxGetPr(Jac[i][1]), &ny, &zero, mxGetPr(R[i]), &nu);
        mxSetCell(plhs[2], i, R[i]);
        
        // gradient
        vec_out[0] = gx+i*nx;
        vec_out[1] = gu+i*nu;
        gi_Fun(vec_in, vec_out);
        
        // constraint residual
        vec_out[0] = c + i*nc;
        ineq_Fun(vec_in, vec_out);
        for (j=0;j<nc/2;j++){
            vec_out[0][nc/2+j] = lb[j] - vec_out[0][j];
            vec_out[0][j] -= ub[j];
        }
        
        // constraint Jacobian
        Ci_Fun(vec_in, Cons[i]);
        mxSetCell(plhs[5], i, Cons[i][0]);
        mxSetCell(plhs[6], i, Cons[i][1]);
    }
    
    // terminal data
    vec_in[0] = xN;
    vec_in[1] = od+N*np;
    vec_in[2] = yN;
    vec_in[3] = QN;
    
    JN_Fun(vec_in, Jac[N]);
    dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, mxGetPr(Jac[N][0]), &nyN, mxGetPr(Jac[N][0]), &nyN, &zero, mxGetPr(Qh[N]), &nx);
    mxSetCell(plhs[0], N, Qh[N]);
    
    vec_out[0] = gx+N*nx;
    gN_Fun(vec_in, vec_out);

    vec_out[0] = c + N*nc;
    ineqN_Fun(vec_in, vec_out);
    for (j=0;j<ncN/2;j++){
        vec_out[0][ncN/2+j] = lbN[j] - vec_out[0][j];
        vec_out[0][j] -= ub[j];
    }
    
    CN_Fun(vec_in, Cons[N]);
    mxSetCell(plhs[5], N, Cons[N][0]);
    
}