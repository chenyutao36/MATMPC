
#include "mex.h"
#include "string.h"

#include "casadi_wrapper.h"
#include "sim.h"

#include "blas.h"

static void *workspace = NULL;
static double *Jac[2];
static double *Jac_N;

static bool mem_alloc = false;

void exitFcn_sim(){
    if (mem_alloc){
        if (workspace!=NULL)
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
    
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
      
    double *Qh = mxGetPr( mxGetField(prhs[2], 0, "Q_h") );
    double *S = mxGetPr( mxGetField(prhs[2], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[2], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[2], 0, "A_sens") );
    double *B = mxGetPr( mxGetField(prhs[2], 0, "B_sens") );
    double *Cx = mxGetPr( mxGetField(prhs[2], 0, "Cx") );
    double *Cu = mxGetPr( mxGetField(prhs[2], 0, "Cu") );
    double *gx = mxGetPr( mxGetField(prhs[2], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[2], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[2], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[2], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[2], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[2], 0, "uc") );
    double *lb_du = mxGetPr( mxGetField(prhs[2], 0, "lb_du") );
    double *ub_du = mxGetPr( mxGetField(prhs[2], 0, "ub_du") );
    double *CxN = mxGetPr( mxGetField(prhs[2], 0, "CxN") );
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - z[i];
    
    // allocate memory
    double *Sens[2];    
    double *Cons[2];
      
    double *vec_in[4];
    double *vec_out[2];    
    vec_in[3] = Q;
   
    double *ode_in[3];
    
    sim_opts opts;
    opts.forw_sens = true;
    
    int size;
    if (!mem_alloc){
        switch(sim_method){
            case 0:
                size = 0;
                break;
            case 1:
                size = sim_erk_calculate_workspace_size(prhs[2],&opts);
                break;
            case 2:
                size = sim_irk_calculate_workspace_size(prhs[2],&opts);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
        
        if (size > 0){
            workspace = mxMalloc(size);
            mexMakeMemoryPersistent(workspace); 
        }
        
        Jac[0] = (double *) mxMalloc(ny*nx * sizeof(double));
        mexMakeMemoryPersistent(Jac[0]); 
        Jac[1] = (double *) mxMalloc(ny*nu * sizeof(double));
        mexMakeMemoryPersistent(Jac[1]); 
        Jac_N = (double *) mxMalloc(nyN*nx * sizeof(double));
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
        
        switch(sim_method){
            case 0:
                F_Fun(vec_in, vec_out);
                D_Fun(vec_in, Sens);
                break;
            case 1:
                ode_in[0]=z+i*nz;
                ode_in[1]=z+i*nz+nx;
                ode_in[2]=od+i*np;          
                sim_erk(ode_in, vec_out, Sens, prhs[2], &opts, workspace);
                break;
            case 2:
                ode_in[0]=z+i*nz;
                ode_in[1]=z+i*nz+nx;
                ode_in[2]=od+i*np;
                sim_irk(ode_in, vec_out, Sens, prhs[2], &opts, workspace);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
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
            vec_in[0]=z+i*nz;
            vec_in[1]=z+i*nz+nx;
            vec_in[2]=od+i*np; 
            vec_out[0] = lc + i*nc;
            path_con_Fun(vec_in, vec_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - vec_out[0][j];
                vec_out[0][j] = lb[i*nc+j] - vec_out[0][j];            
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