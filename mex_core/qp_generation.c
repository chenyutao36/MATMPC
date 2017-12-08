

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

static void *workspace = NULL;

void exitFcn_sim(){
    if (workspace!=NULL)
        mxFree(workspace);
}

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
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ubu") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    mwSize nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
    
    plhs[0] = mxCreateCellMatrix(1, N+1); // Qh
    plhs[1] = mxCreateCellMatrix(1, N); // S
    plhs[2] = mxCreateCellMatrix(1, N); // R
    plhs[3] = mxCreateCellMatrix(1, N); // A
    plhs[4] = mxCreateCellMatrix(1, N); // B
    plhs[5] = mxCreateCellMatrix(1, N+1); // Cx
    plhs[6] = mxCreateCellMatrix(1, N); // Cu
    plhs[7] = mxCreateDoubleMatrix(nx, N+1, mxREAL); //gx
    plhs[8] = mxCreateDoubleMatrix(nu, N, mxREAL); //gu 
    plhs[9] = mxCreateDoubleMatrix(nx,N, mxREAL); //a
    plhs[10] = mxCreateDoubleMatrix(nx,1,mxREAL); // ds0    
    plhs[11] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL); //lc
    plhs[12] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL); //uc
    plhs[13] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);   // lb_du
    plhs[14] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);   // ub_du   
    
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
      
    double *gx = mxGetPr(plhs[7]);
    double *gu = mxGetPr(plhs[8]);   
    double *a = mxGetPr(plhs[9]);
    double *ds0 = mxGetPr(plhs[10]);   
    double *lc = mxGetPr(plhs[11]);
    double *uc = mxGetPr(plhs[12]);
    double *lb_du = mxGetPr(plhs[13]);
    double *ub_du = mxGetPr(plhs[14]);
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - z[i];
    
    // allocate memory
    mxArray **Sens[N];
    mxArray **Jac[N+1];
    mxArray **Cons[N+1];
    mxArray *Qh[N+1];
    mxArray *S[N];
    mxArray *R[N];
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
    
    double *vec_in[4];
    double *vec_out[2];    
    vec_in[3] = Q;
   
    double *ode_in[3];

    if (workspace== NULL && sim_method!=0){
        int size = 0;
        if (sim_method == 1)
            size = sim_erk_calculate_workspace_size(prhs[2],true);
        if (sim_method ==2)
            size = sim_irk_calculate_workspace_size(prhs[2],true);    
        
        workspace = mxMalloc(size);
        mexMakeMemoryPersistent(workspace);        
        mexAtExit(exitFcn_sim);
    }
    
    // start loop
    for(i=0;i<N;i++){
        vec_in[0] = z+i*nz;
        vec_in[1] = od+i*np;
        vec_in[2] = y+i*ny;
        
        // control bounds
        for (j=0;j<nu;j++){
            lb_du[i*nu+j] = lbu[i*nu+j]-z[i*nz+nx+j];
            ub_du[i*nu+j] = ubu[i*nu+j]-z[i*nz+nx+j];
        }
        
        // integration
        vec_out[0] = a+i*nx;

        if (sim_method == 0){
            F_Fun(vec_in, vec_out);
            D_Fun(vec_in, Sens[i]);  
        }
        if (sim_method == 1){
            ode_in[0]=z+i*nz;
            ode_in[1]=z+i*nz+nx;
            ode_in[2]=od+i*np;          
            sim_erk(ode_in, vec_out, Sens[i], prhs[2], true, workspace);
        }
        if (sim_method == 2){
            ode_in[0]=z+i*nz;
            ode_in[1]=z+i*nz+nx;
            ode_in[2]=od+i*np;
            sim_irk(ode_in, vec_out, Sens[i], prhs[2], true, workspace);
        }
        
        if (i < N-1){
            for (j=0;j<nx;j++)
                vec_out[0][j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                vec_out[0][j] -= xN[j];
        }
            
        // sensitiities
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
        if (nc>0){        
            vec_out[0] = lc + i*nc;
            path_con_Fun(vec_in, vec_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[j] - vec_out[0][j];
                vec_out[0][j] = lb[j] - vec_out[0][j];            
            }
        
            // constraint Jacobian
            Ci_Fun(vec_in, Cons[i]);
            mxSetCell(plhs[5], i, Cons[i][0]);
            mxSetCell(plhs[6], i, Cons[i][1]);
        }
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

    if (ncN>0){
        vec_out[0] = lc + N*nc;
        path_con_N_Fun(vec_in, vec_out);
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ubN[j] - vec_out[0][j];
            vec_out[0][j] = lbN[j] - vec_out[0][j];            
        }

        CN_Fun(vec_in, Cons[N]);
        mxSetCell(plhs[5], N, Cons[N][0]);
    }
}