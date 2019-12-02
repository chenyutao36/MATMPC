#include "mex.h"
#include "string.h"

#include "casadi_wrapper.h"
#include "sim.h"
#include "erk.h"
#include "irk_ode.h"
#include "irk_dae.h"

#include "lapack.h"
#include "blas.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
static sim_irk_ode_workspace *irk_ode_workspace = NULL;
static sim_irk_dae_workspace *irk_dae_workspace = NULL;
static bool mem_alloc = false;

static double *casadi_out[3];
static double *L = NULL;
static double *eq_res_vec = NULL;
static double *lu, *uu, *lx, *ux, *tmp;

void exitFcn(){   
    if (erk_workspace!=NULL)
        sim_erk_workspace_free(opts, erk_workspace);
    if (irk_ode_workspace!=NULL)
        sim_irk_ode_workspace_free(opts, irk_ode_workspace);
    if (irk_dae_workspace!=NULL)
        sim_irk_dae_workspace_free(opts, irk_dae_workspace);
    if (opts!=NULL)
        sim_opts_free(opts);
    if (in!=NULL)
        sim_in_free(in);
    if (out!=NULL)
        sim_out_free(out);
    if (mem_alloc){
        mxFree(casadi_out[1]);
        mxFree(casadi_out[2]);
        mxFree(L);
        mxFree(eq_res_vec);
        mxFree(lu);
        mxFree(uu);
        mxFree(lx);
        mxFree(ux);
        mxFree(tmp);
    }
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{    
    double *x = mxGetPr( mxGetField(prhs[0], 0, "x") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "u") );
    double *z = mxGetPr( mxGetField(prhs[0], 0, "z") );
    double *lambda = mxGetPr( mxGetField(prhs[0], 0, "lambda") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu") );    
    double *mu_u = mxGetPr( mxGetField(prhs[0], 0, "mu_u") );
    double *mu_x = mxGetPr( mxGetField(prhs[0], 0, "mu_x") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *W = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *WN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub") );
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ubu") );
    double *lbx = mxGetPr( mxGetField(prhs[0], 0, "lbx") );
    double *ubx = mxGetPr( mxGetField(prhs[0], 0, "ubx") );
    
    double *ds0 = mxGetPr( mxGetField(prhs[2], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[2], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[2], 0, "uc") );
    double *obj = mxGetPr( mxGetField(prhs[2], 0, "obj") );
    double *z_out = mxGetPr( mxGetField(prhs[2], 0, "z_out") );
     
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
    double Ts_st = mxGetScalar(mxGetField(prhs[1], 0, "Ts_st") );
    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
          
    size_t i=0,j=0,l=0;
    size_t nv = nx+nu;
    size_t nw = N*nv+nx;
    size_t neq = (N+1)*nx;
    size_t nineq = N*nc+ncN;
    char *nTrans = "N", *Trans="T", *Norm="O";
    double one_d = 1.0, zero = 0.0;
    size_t one_i = 1;
    int idx;
    
    double *casadi_in[9];
        
    if (!mem_alloc){
        casadi_out[1] = (double *)mxMalloc(nv * sizeof(double));
        mexMakeMemoryPersistent(casadi_out[1]);
        casadi_out[2] = (double *)mxMalloc(nv * sizeof(double));
        mexMakeMemoryPersistent(casadi_out[2]);     
        L = (double *)mxMalloc( nw * sizeof(double));
        mexMakeMemoryPersistent(L);
        
        eq_res_vec = (double *)mxMalloc( neq * sizeof(double));
        mexMakeMemoryPersistent(eq_res_vec);
        
        lu = (double *)mxMalloc( N*nu * sizeof(double));
        mexMakeMemoryPersistent(lu);
        uu = (double *)mxMalloc( N*nu * sizeof(double));
        mexMakeMemoryPersistent(uu);
        lx = (double *)mxMalloc( N*nbx * sizeof(double)); 
        mexMakeMemoryPersistent(lx);
        ux = (double *)mxMalloc( N*nbx * sizeof(double));
        mexMakeMemoryPersistent(ux);
        
        tmp = (double *)mxCalloc(nbx,sizeof(double));
        mexMakeMemoryPersistent(tmp);
        
        switch(sim_method){
            case 1:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, prhs[2], erk_workspace);   
                break;
            case 2:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_ode_workspace = sim_irk_ode_workspace_create(opts);               
                sim_irk_ode_workspace_init(opts, prhs[2], irk_ode_workspace);
                break;
            case 3:
                opts = sim_opts_create(prhs[2]);
                opts->forw_sens_flag = false;
                opts->adj_sens_flag = true;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_dae_workspace = sim_irk_dae_workspace_create(opts);               
                sim_irk_dae_workspace_init(opts, prhs[2], irk_dae_workspace);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
         
        }   
                
        mem_alloc = true;     
        mexAtExit(exitFcn);
    }
    
    memcpy(eq_res_vec, ds0, nx*sizeof(double));
    
    double *work;
    double KKT=0, eq_res=0, ineq_res=0, OBJ=0;
    
    for (i=0;i<N;i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = od+i*np;
        casadi_in[3] = y+i*ny;
        casadi_in[4] = W+i*ny;
        casadi_in[5] = lambda+(i+1)*nx;
        if (i==0){
            casadi_in[6] = tmp;
        }
        else{
            casadi_in[6] = mu_x+(i-1)*nbx;
        }
        casadi_in[7] = mu_u+i*nu;
        casadi_in[8] = mu+i*nc;

        // objective value      
        casadi_out[0] = obj;
        obji_Fun(casadi_in, casadi_out);
        OBJ += obj[0];
              
        // call integrator to compute xend and adjoint sensitivities
        switch(sim_method){
            case 1:      
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->adj_sens = casadi_out[2];
                sim_erk(in, out, opts, erk_workspace);
                break;
            case 2:                         
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->adj_sens = casadi_out[2];
                sim_irk_ode(in, out, opts, irk_ode_workspace);
                break;
            case 3:                         
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                in->lambda = lambda+(i+1)*nx;
                out->xn = eq_res_vec+(i+1)*nx;
                out->zn = z_out+i*nz;
                out->adj_sens = casadi_out[2];
                sim_irk_dae(in, out, opts, irk_dae_workspace);
                break;
            default :
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }

        // compute gradient of Lagrangian
        if (i>0){
            for (j=0;j<nx;j++)
                casadi_out[2][j] -= lambda[i*nx+j];
        }else{
            for (j=0;j<nx;j++)
                casadi_out[2][j] += lambda[j];
        }
        casadi_out[0] = L+i*nv;
        adj_Fun(casadi_in, casadi_out);                      
        daxpy(&nv, &one_d, casadi_out[1], &one_i, casadi_out[0], &one_i); // dojb+dB'*mu
        daxpy(&nv, &one_d, casadi_out[2], &one_i, casadi_out[0], &one_i); // dobj+dB'*mu+dG'*lambda 
        
        // equality constraint residual
        for (j=0;j<nx;j++)
            eq_res_vec[(i+1)*nx+j] -= x[(i+1)*nx+j];
        
        // control bound residual
        for (j=0;j<nu;j++){
            lu[i*nu+j] = lbu[i*nu+j] - u[i*nu+j];
            uu[i*nu+j] = ubu[i*nu+j] - u[i*nu+j];
        }
        
        // state bound residual
        for (j=0;j<nbx;j++){
            idx = (int)nbx_idx[j]-1;
            lx[i*nbx+j] = lbx[i*nbx+j] - x[(i+1)*nx+idx];
            ux[i*nbx+j] = ubx[i*nbx+j] - x[(i+1)*nx+idx];
        }
        
        // general constraint residual
        if (nc>0){
            casadi_out[0] = lc + i*nc;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - casadi_out[0][j];
                casadi_out[0][j] = lb[i*nc+j] - casadi_out[0][j];            
            }
        }
    }
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    casadi_in[4] = mu_x+(N-1)*nbx;
    casadi_in[5] = mu+N*nc;
    
    casadi_out[0] = L+N*nv;
    adjN_Fun(casadi_in, casadi_out);
    for (j=0;j<nx;j++)
        casadi_out[0][j] -= lambda[N*nx+j];
    
    daxpy(&nx, &one_d, casadi_out[1], &one_i, casadi_out[0], &one_i);
               
    if (ncN>0){
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out); 
        for (j=0;j<ncN;j++){
            uc[N*nc+j] = ub[N*nc+j] - casadi_out[0][j];
            casadi_out[0][j] = lb[N*nc+j] - casadi_out[0][j];            
        }
    }
         
    eq_res = dlange(Norm, &neq, &one_i, eq_res_vec, &neq, work);
    KKT = dnrm2(&nw, L, &one_i);
    
    for (i=0;i<N*nu;i++)
        ineq_res += MAX(-1*uu[i],0) + MAX(lu[i],0);
    for (i=0;i<N*nbx;i++)
        ineq_res += MAX(-1*ux[i],0) + MAX(lx[i],0);
    for (i=0;i<nineq;i++)
        ineq_res += MAX(-1*uc[i],0) + MAX(lc[i],0);
    
    casadi_out[0] = obj;
    objN_Fun(casadi_in, casadi_out);
    OBJ += obj[0];
    
    OBJ *= Ts_st;

    memcpy(z, z_out, N*nz*sizeof(double));

    plhs[0] = mxCreateDoubleScalar(eq_res); // eq_res
    plhs[1] = mxCreateDoubleScalar(ineq_res); // ineq_res
    plhs[2] = mxCreateDoubleScalar(KKT); // KKT
    plhs[3] = mxCreateDoubleScalar(OBJ); // OBJ
}