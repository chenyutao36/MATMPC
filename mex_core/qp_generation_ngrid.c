
#include "mex.h"
#include "string.h"

#include "sim_workspace.h"
#include "casadi_wrapper.h"
#include "mpc_common.h"

#include "blas.h"

static sim_workspace integrator;
static bool mem_alloc = false;
static size_t workspace_shape[4];
static double *Hes[1];
static double *HesN[1];
static double *Jac[2];
static double *JacN[1];
static double *temp[3];

static void free_scratch(void){
    if (Hes[0] != NULL){
        mxFree(Hes[0]);
        Hes[0] = NULL;
    }
    if (HesN[0] != NULL){
        mxFree(HesN[0]);
        HesN[0] = NULL;
    }
    if (Jac[0] != NULL){
        mxFree(Jac[0]);
        Jac[0] = NULL;
    }
    if (Jac[1] != NULL){
        mxFree(Jac[1]);
        Jac[1] = NULL;
    }
    if (JacN[0] != NULL){
        mxFree(JacN[0]);
        JacN[0] = NULL;
    }
    if (temp[0] != NULL){
        mxFree(temp[0]);
        temp[0] = NULL;
    }
    if (temp[1] != NULL){
        mxFree(temp[1]);
        temp[1] = NULL;
    }
    if (temp[2] != NULL){
        mxFree(temp[2]);
        temp[2] = NULL;
    }
    mem_alloc = false;
}

static void exitFcn(void){
    free_scratch();
    sim_workspace_free(&integrator);
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *x = mxGetPr( mxGetField(prhs[0], 0, "x") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "u") );
    double *z = mxGetPr( mxGetField(prhs[0], 0, "z") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *W = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *WN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub") );
    double *x0 = mxGetPr( mxGetField(prhs[0], 0, "x0") );
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ubu") );
    double *lbx = mxGetPr( mxGetField(prhs[0], 0, "lbx") );
    double *ubx = mxGetPr( mxGetField(prhs[0], 0, "ubx") );
        
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nz = mxGetScalar( mxGetField(prhs[1], 0, "nz") );
    size_t np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    size_t ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    size_t nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    double *nbx_idx = mxGetPr( mxGetField(prhs[1], 0, "nbx_idx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
    
    int i=0,j=0;
    char *nTrans = "N", *Trans="T", *UPLO="L";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    size_t one_i = 1;
    int idx;
      
    double *Q = mxGetPr( mxGetField(prhs[2], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[2], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[2], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[2], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[2], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[2], 0, "Cx") );
    double *Cgx = mxGetPr( mxGetField(prhs[2], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[2], 0, "Cgu") );
    double *CgN = mxGetPr( mxGetField(prhs[2], 0, "CgN") );
    double *gx = mxGetPr( mxGetField(prhs[2], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[2], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[2], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[2], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[2], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[2], 0, "uc") );
    double *lb_du = mxGetPr( mxGetField(prhs[2], 0, "lb_du") );
    double *ub_du = mxGetPr( mxGetField(prhs[2], 0, "ub_du") );
    double *lb_dx = mxGetPr( mxGetField(prhs[2], 0, "lb_dx") );
    double *ub_dx = mxGetPr( mxGetField(prhs[2], 0, "ub_dx") );
    double *z_out = mxGetPr( mxGetField(prhs[2], 0, "z_out") );
    
    double reg = mxGetScalar( mxGetField(prhs[2], 0, "reg") );
    int hessian_type = mxGetScalar( mxGetField(prhs[2], 0, "hessian") );
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - x[i];
    
    // allocate memory
    double *Sens[2];    
    double *Cons[2];
      
    double *casadi_in[5];
    double *casadi_out[2];    
           
    if (integrator.opts == NULL)
        mexAtExit(exitFcn);
    sim_workspace_prepare(&integrator, prhs[2], sim_method, nx, nu, nz, true, false);
    const double *grid = sim_grid_indices(prhs[2], N, true);

    const size_t shape[] = { nx, nu, ny, nyN };
    if (!mem_alloc || memcmp(shape, workspace_shape, sizeof(shape)) != 0){
        free_scratch();
        
        Hes[0] = (double *) mxMalloc(ny*ny * sizeof(double));	
        mexMakeMemoryPersistent(Hes[0]); 	
        HesN[0] = (double *) mxMalloc(nyN*nyN * sizeof(double));	
        mexMakeMemoryPersistent(HesN[0]); 	
                       
        Jac[0] = (double *) mxMalloc(ny*nx * sizeof(double));	
        mexMakeMemoryPersistent(Jac[0]); 	
        Jac[1] = (double *) mxMalloc(ny*nu * sizeof(double));	
        mexMakeMemoryPersistent(Jac[1]); 	
        JacN[0] = (double *) mxMalloc(nyN*nx * sizeof(double));	
        mexMakeMemoryPersistent(JacN[0]);       
        
        temp[0] = (double *) mxMalloc(ny*nx * sizeof(double));	
        mexMakeMemoryPersistent(temp[0]); 	
        temp[1] = (double *) mxMalloc(ny*nu * sizeof(double));	
        mexMakeMemoryPersistent(temp[1]);        
        temp[2] = (double *) mxMalloc(nyN*nx * sizeof(double));	
        mexMakeMemoryPersistent(temp[2]); 
        
        memcpy(workspace_shape, shape, sizeof(shape));
        mem_alloc = true;
    }
        
    // start loop
    for(i=0;i<N;i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = od+i*np;
        casadi_in[3] = y+i*ny;
        casadi_in[4] = W+i*ny;
        
        // control bounds
        for (j=0;j<nu;j++){
            lb_du[i*nu+j] = lbu[i*nu+j]-u[i*nu+j];
            ub_du[i*nu+j] = ubu[i*nu+j]-u[i*nu+j];
        }
        
        // state bounds
        for (j=0;j<nbx;j++){
            idx = (int)nbx_idx[j]-1;
            lb_dx[i*nbx+j] = lbx[i*nbx+j]-x[(i+1)*nx+idx];
            ub_dx[i*nbx+j] = ubx[i*nbx+j]-x[(i+1)*nx+idx];
        }
        
        // integration    
        sim_set_stage_step(integrator.opts, grid, i);
        switch(sim_method){
            case 1:
                integrator.in->x = x+i*nx;
                integrator.in->u = u+i*nu;
                integrator.in->p = od+i*np;
                integrator.out->xn = a+i*nx;
                integrator.out->Sx = A + i*nx*nx;
                integrator.out->Su = B + i*nx*nu;
                sim_erk(integrator.in, integrator.out, integrator.opts, integrator.erk);
                break;
            case 2:
                integrator.in->x = x+i*nx;
                integrator.in->u = u+i*nu;
                integrator.in->p = od+i*np;
                integrator.in->z = z+i*nz;
                integrator.out->xn = a+i*nx;
                integrator.out->Sx = A + i*nx*nx;
                integrator.out->Su = B + i*nx*nu;
                sim_irk_ode(integrator.in, integrator.out, integrator.opts, integrator.irk_ode);
                break;
            case 3:
                integrator.in->x = x+i*nx;
                integrator.in->u = u+i*nu;
                integrator.in->p = od+i*np;
                integrator.in->z = z+i*nz;
                integrator.out->xn = a+i*nx;
                integrator.out->Sx = A + i*nx*nx;
                integrator.out->Su = B + i*nx*nu;
                integrator.out->zn = z_out + i*nz;
                sim_irk_dae(integrator.in, integrator.out, integrator.opts, integrator.irk_dae);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
       
        // equality residual        
        for (j=0;j<nx;j++)
            a[i*nx+j] -= x[(i+1)*nx+j];
       
        // Hessian
        Ji_Fun(casadi_in, Jac);
        switch(hessian_type){                
            case 0:                   
                dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, Jac[0], &ny, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
                dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, Jac[0], &ny, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
                dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, Jac[1], &ny, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                break;
            
            case 1:
                Hi_Fun(casadi_in, Hes);
                dgemm(Trans, nTrans, &nx, &ny, &ny, &one_d, Jac[0], &ny, Hes[0], &ny, &zero, temp[0], &nx);
                dgemm(nTrans, nTrans, &nx, &nx, &ny, &one_d, temp[0], &nx, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
                
                dgemm(nTrans, nTrans, &nx, &nu, &ny, &one_d, temp[0], &nx, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
                
                dgemm(Trans, nTrans, &nu, &ny, &ny, &one_d, Jac[1], &ny, Hes[0], &ny, &zero, temp[1], &nu);
                dgemm(nTrans, nTrans, &nu, &nu, &ny, &one_d, temp[1], &nu, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                
                break;
            default:
                mexErrMsgTxt("Please choose a supported Hessian type");
                break;
                            
        }
                
        regularization(nx, Q+i*nx*nx, reg);
        regularization(nu, R+i*nu*nu, reg);

        // gradient
        casadi_out[0] = gx+i*nx;
        casadi_out[1] = gu+i*nu;
        gi_Fun(casadi_in, casadi_out);
                        
        // constraint residual
        if (nc>0){  
            casadi_in[0]=x+i*nx;
            casadi_in[1]=u+i*nu;
            casadi_in[2]=od+i*np; 
            casadi_out[0] = lc + i*nc;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - casadi_out[0][j];
                casadi_out[0][j] = lb[i*nc+j] - casadi_out[0][j];            
            }
        
            // constraint Jacobian
            Cons[0] = Cgx+i*nc*nx;
            Cons[1] = Cgu+i*nc*nu;
            Ci_Fun(casadi_in, Cons);
        }
    }
    
    // terminal data
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    
    JN_Fun(casadi_in, JacN);
    switch(hessian_type){            
        case 0:
            dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, JacN[0], &nyN, JacN[0], &nyN, &zero, Q+N*nx*nx, &nx);
            break;
        case 1:
            HN_Fun(casadi_in, HesN);
            dgemm(Trans, nTrans, &nx, &nyN, &nyN, &one_d, JacN[0], &nyN, HesN[0], &nyN, &zero, temp[2], &nx);
            dgemm(nTrans, nTrans, &nx, &nx, &nyN, &one_d, temp[2], &nx, JacN[0], &nyN, &zero, Q+N*nx*nx, &nx);
            break;
        default:
            mexErrMsgTxt("Please choose a supported Hessian type");
            break;
            
    }
       
    regularization(nx, Q+N*nx*nx, reg);
    
        
    casadi_out[0] = gx+N*nx;
    gN_Fun(casadi_in, casadi_out);

    if (ncN>0){
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out);
        for (j=0;j<ncN;j++){
            uc[N*nc+j] = ub[N*nc+j] - casadi_out[0][j];
            casadi_out[0][j] = lb[N*nc+j] - casadi_out[0][j];            
        }

        CN_Fun(casadi_in, &CgN);
    }
    
}
