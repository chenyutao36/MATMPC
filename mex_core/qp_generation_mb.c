
#include "mex.h"
#include "string.h"

#include "sim.h"
#include "erk.h"
#include "irk.h"
#include "casadi_wrapper.h"
#include "mpc_common.h"

#include "blas.h"

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
static sim_irk_workspace *irk_workspace = NULL;
static bool mem_alloc = false;

static double *R_temp= NULL;
static double *gu_temp = NULL;

void exitFcn(){
    if (erk_workspace!=NULL)
        sim_erk_workspace_free(opts, erk_workspace);
    if (irk_workspace!=NULL)
        sim_irk_workspace_free(opts, irk_workspace);
    if (opts!=NULL)
        sim_opts_free(opts);
    if (in!=NULL)
        sim_in_free(in);
    if (out!=NULL)
        sim_out_free(out);    
    if (R_temp!=NULL)
        mxFree(R_temp);   
    if (gu_temp!=NULL)
        mxFree(gu_temp); 
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *x = mxGetPr( mxGetField(prhs[0], 0, "x") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "u") );
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
    size_t np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    size_t ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    size_t nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    double *nbx_idx = mxGetPr( mxGetField(prhs[1], 0, "nbx_idx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
    
    int i=0,j=0,k=0;
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
    double *index_T = mxGetPr( mxGetField(prhs[2], 0, "index_T") );
    size_t r = mxGetScalar( mxGetField(prhs[2], 0, "r") );
    
    int lin_obj = mxGetScalar( mxGetField(prhs[2], 0, "lin_obj") );
    double reg = mxGetScalar( mxGetField(prhs[2], 0, "reg") );
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - x[i];
    
    // allocate memory
    double *Sens[2];    
    double *Cons[2];
      
    double *casadi_in[5];
    double *casadi_out[2];    
    double *Hes[3];
       
    if (!mem_alloc){
        switch(sim_method){
            case 0:
                break;
            case 1:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, prhs[2], erk_workspace);                
                break;
            case 2:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_workspace = sim_irk_workspace_create(opts);               
                sim_irk_workspace_init(opts, prhs[2], irk_workspace);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }                      

        R_temp = (double *) mxMalloc(nu*nu * sizeof(double));
        mexMakeMemoryPersistent(R_temp);    
        gu_temp = (double *) mxMalloc(nu * sizeof(double));
        mexMakeMemoryPersistent(gu_temp);
        
        mem_alloc=true;
        mexAtExit(exitFcn);
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
        switch(sim_method){
            case 0:
                casadi_out[0] = a+i*nx;
                Sens[0] = A + i*nx*nx;
                Sens[1] = B + i*nx*nu;
                F_Fun(casadi_in, casadi_out);
                D_Fun(casadi_in, Sens);
                break;
            case 1:
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                out->xn = a+i*nx;
                out->Sx = A + i*nx*nx;
                out->Su = B + i*nx*nu;
                sim_erk(in, out, opts, erk_workspace);
                break;
            case 2:
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                out->xn = a+i*nx;
                out->Sx = A + i*nx*nx;
                out->Su = B + i*nx*nu;
                sim_irk(in, out, opts, irk_workspace);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
       
        // equality residual        
        for (j=0;j<nx;j++)
            a[i*nx+j] -= x[(i+1)*nx+j];
       
        // Hessian
        if (!lin_obj){
            Hes[0] = Q+i*nx*nx;
            Hes[1] = R+i*nu*nu;
            Hes[2] = S+i*nx*nu;
            Hi_Fun(casadi_in, Hes);
            
            if ( (int)(index_T[k])==i){
                memcpy(R_temp, R+i*nu*nu, nu*nu*sizeof(double));
            }  
            else{
                for (j=0;j<nu*nu;j++)
                      R_temp[j] += Hes[1][j];
            }  
            
            
            if ( (int)(index_T[k+1])==i+1){

                memcpy(R+k*nu*nu, R_temp, nu*nu*sizeof(double));
                regularization(nu, R+k*nu*nu, reg);
            }   
                
            regularization(nx, Q+i*nx*nx, reg);
            
        }
        
        // gradient
        casadi_out[0] = gx+i*nx;
        casadi_out[1] = gu_temp;
        gi_Fun(casadi_in, casadi_out);
        if ( (int)(index_T[k])==i){
                memcpy(gu+k*nu, gu_temp, nu*sizeof(double));
        }  
        else{
             for (j=0;j<nu;j++)
                *(gu+k*nu+j) += *(gu_temp+j);   
        }
          
                
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
            Cons[1] = Cgu+i*nu*nc;
            Ci_Fun(casadi_in, Cons);
            
            
        }
        if ( (int)(index_T[k+1])==i+1){
                k=k+1;
        }
    }
    
    // terminal data
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    
    if (!lin_obj){
        Hes[0]=Q+N*nx*nx;
        HN_Fun(casadi_in, Hes);
        regularization(nx, Q+N*nx*nx, reg);
    }
        
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