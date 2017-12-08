

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
static double *Jac[2];
static double *Jac_N;

static bool mem_alloc = false;

void exitFcn_sim(){
    if (mem_alloc){
        mxFree(workspace);
        mxFree(Jac[0]);
        mxFree(Jac[1]);
        mxFree(Jac_N);
    }
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
    
    plhs[0] = mxCreateDoubleMatrix(nx,nx*(N+1),mxREAL); // Qh
    plhs[1] = mxCreateDoubleMatrix(nx,nu*N,mxREAL); // S
    plhs[2] = mxCreateDoubleMatrix(nu,nu*N,mxREAL); // R
    plhs[3] = mxCreateDoubleMatrix(nx,nx*N,mxREAL); // A
    plhs[4] = mxCreateDoubleMatrix(nx,nu*N,mxREAL); // B
    plhs[5] = mxCreateDoubleMatrix(nc,nx*N,mxREAL); // Cx
    plhs[6] = mxCreateDoubleMatrix(nc,nu*N,mxREAL); // Cu
    plhs[7] = mxCreateDoubleMatrix(nx, N+1, mxREAL); //gx
    plhs[8] = mxCreateDoubleMatrix(nu, N, mxREAL); //gu 
    plhs[9] = mxCreateDoubleMatrix(nx,N, mxREAL); //a
    plhs[10] = mxCreateDoubleMatrix(nx,1,mxREAL); // ds0    
    plhs[11] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL); //lc
    plhs[12] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL); //uc
    plhs[13] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);   // lb_du
    plhs[14] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);   // ub_du     
    plhs[15] = mxCreateDoubleMatrix(ncN,nx,mxREAL); // CxN
    
    
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
      
    double *Qh = mxGetPr(plhs[0]);
    double *S = mxGetPr(plhs[1]);   
    double *R = mxGetPr(plhs[2]);
    double *A = mxGetPr(plhs[3]);   
    double *B = mxGetPr(plhs[4]);
    double *Cx = mxGetPr(plhs[5]);
    double *Cu = mxGetPr(plhs[6]);
    double *gx = mxGetPr(plhs[7]);
    double *gu = mxGetPr(plhs[8]);   
    double *a = mxGetPr(plhs[9]);
    double *ds0 = mxGetPr(plhs[10]);   
    double *lc = mxGetPr(plhs[11]);
    double *uc = mxGetPr(plhs[12]);
    double *lb_du = mxGetPr(plhs[13]);
    double *ub_du = mxGetPr(plhs[14]);
    double *CxN = mxGetPr(plhs[15]);
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - z[i];
    
    // allocate memory
    double *Sens[2];    
    double *Cons[2];
      
    double *vec_in[4];
    double *vec_out[2];    
    vec_in[3] = Q;
   
    double *ode_in[3];

    if (!mem_alloc){
        int size = 0;
        if (sim_method == 1)
            size = sim_erk_calculate_workspace_size(prhs[2],true);
        if (sim_method ==2)
            size = sim_irk_calculate_workspace_size(prhs[2],true);    
        
        workspace = mxMalloc(size);
        mexMakeMemoryPersistent(workspace); 
        
        Jac[0] = (double *) mxCalloc(ny*nx, sizeof(double));
        mexMakeMemoryPersistent(Jac[0]); 
        Jac[1] = (double *) mxCalloc(ny*nu, sizeof(double));
        mexMakeMemoryPersistent(Jac[1]); 
        Jac_N = (double *) mxCalloc(nyN*nx, sizeof(double));
        mexMakeMemoryPersistent(Jac_N);
        
        mem_alloc=true;
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
        Sens[0] = A + i*nx*nx;
        Sens[1] = B + i*nx*nu;

        if (sim_method == 0){
            F_Fun(vec_in, vec_out);
            D_Fun(vec_in, Sens);
        }
        if (sim_method == 1){
            ode_in[0]=z+i*nz;
            ode_in[1]=z+i*nz+nx;
            ode_in[2]=od+i*np;          
            sim_erk(ode_in, vec_out, Sens, prhs[2], true, workspace);
        }
        if (sim_method == 2){
            ode_in[0]=z+i*nz;
            ode_in[1]=z+i*nz+nx;
            ode_in[2]=od+i*np;
            sim_irk(ode_in, vec_out, Sens, prhs[2], true, workspace);
        }
        
        if (i < N-1){
            for (j=0;j<nx;j++)
                vec_out[0][j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                vec_out[0][j] -= xN[j];
        }

        
        // Hessian
        Ji_Fun(vec_in, Jac);
        dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, Jac[0], &ny, Jac[0], &ny, &zero, Qh+i*nx*nx, &nx);
        dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, Jac[0], &ny, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
        dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, Jac[1], &ny, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
        
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
            Cons[0] = Cx+i*nc*nx;
            Cons[1] = Cu+i*nc*nu;
            Ci_Fun(vec_in, Cons);
        }
    }
    
    // terminal data
    vec_in[0] = xN;
    vec_in[1] = od+N*np;
    vec_in[2] = yN;
    vec_in[3] = QN;
    
    JN_Fun(vec_in, Jac_N);
    dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, Jac_N, &nyN, Jac_N, &nyN, &zero, Qh+N*nx*nx, &nx);
    
    vec_out[0] = gx+N*nx;
    gN_Fun(vec_in, vec_out);

    if (ncN>0){
        vec_out[0] = lc + N*nc;
        path_con_N_Fun(vec_in, vec_out);
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ubN[j] - vec_out[0][j];
            vec_out[0][j] = lbN[j] - vec_out[0][j];            
        }

        CN_Fun(vec_in, CxN);
    }
    
}