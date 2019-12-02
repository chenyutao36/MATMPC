#include "mex.h"
#include <stdlib.h>
#include "string.h"

#include "sim.h"
#include "erk.h"
#include "mpc_common.h"
#include "casadi_wrapper.h"
#include "blas.h"

sim_erk_workspace* sim_erk_workspace_create(sim_opts *opts)
{
    sim_erk_workspace* work = (sim_erk_workspace*)mxMalloc(sizeof(sim_erk_workspace));
        
    work->A = mxMalloc(opts->num_stages*opts->num_stages*sizeof(double));
    work->B = mxMalloc(opts->num_stages*sizeof(double));
    work->Sx = mxMalloc(opts->nx*opts->nx*sizeof(double));
    work->Su = mxMalloc(opts->nx*opts->nu*sizeof(double));    
    work->xt = mxMalloc(opts->nx*sizeof(double));
    work->K = mxMalloc(opts->num_stages*opts->nx*sizeof(double));
    
    mexMakeMemoryPersistent(work->A);
    mexMakeMemoryPersistent(work->B);
    mexMakeMemoryPersistent(work->Sx);
    mexMakeMemoryPersistent(work->Su);
    mexMakeMemoryPersistent(work->xt);
    mexMakeMemoryPersistent(work->K);
    
    if (opts->forw_sens_flag){
        work->dKx = mxMalloc(opts->num_stages*opts->nx*opts->nx*sizeof(double));
        work->dKu = mxMalloc(opts->num_stages*opts->nx*opts->nu*sizeof(double));
        work->jacX_t = mxMalloc(opts->nx*opts->nx*sizeof(double));
        work->jacU_t = mxMalloc(opts->nx*opts->nu*sizeof(double));
        
        mexMakeMemoryPersistent(work->dKx);
        mexMakeMemoryPersistent(work->dKu);
        mexMakeMemoryPersistent(work->jacX_t);
        mexMakeMemoryPersistent(work->jacU_t);
    }

    if (opts->adj_sens_flag){
        work->X_traj = mxMalloc(opts->num_steps*opts->num_stages*opts->nx*sizeof(double));
        work->K_lambda = mxMalloc(opts->num_stages*(opts->nx+opts->nu)*sizeof(double));
        work->lambda_t = mxCalloc((opts->nx+opts->nu), sizeof(double));

        mexMakeMemoryPersistent(work->X_traj);
        mexMakeMemoryPersistent(work->K_lambda);
        mexMakeMemoryPersistent(work->lambda_t);
    }
    
    mexMakeMemoryPersistent(work);
    
    return work;
}

void sim_erk_workspace_init(sim_opts *opts, const mxArray *mem, sim_erk_workspace *work)
{
    size_t nx = opts->nx;
    size_t nu = opts->nu;
    size_t num_stages = opts->num_stages;
    
    double *A = mxGetPr( mxGetField(mem, 0, "A_tab"));
    double *B = mxGetPr( mxGetField(mem, 0, "B_tab"));
    double *Sx = mxGetPr( mxGetField(mem, 0, "Sx"));
    double *Su = mxGetPr( mxGetField(mem, 0, "Su"));
    
    memcpy(work->A, A, num_stages*num_stages*sizeof(double));
    memcpy(work->B, B, num_stages*sizeof(double));
    memcpy(work->Sx, Sx, nx*nx*sizeof(double));
    memcpy(work->Su, Su, nx*nu*sizeof(double));
}

int sim_erk(sim_in *in, sim_out *out, sim_opts *opts, sim_erk_workspace *workspace)
{    
    int istep, s, i, j;
    double a,b;

    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
    
    double h= opts->h;
    size_t nx = opts->nx;
    size_t nu = opts->nu;
    size_t num_stages = opts->num_stages;
    size_t num_steps = opts->num_steps;
    bool forw_sens_flag = opts->forw_sens_flag;
    bool adj_sens_flag = opts->adj_sens_flag;
        
    double *x = in->x;
    double *lambda=in->lambda;

    double *xn = out->xn;  
    double *Jac_x = out->Sx;
    double *Jac_u = out->Su;  
    double *adj_sens = out->adj_sens; 
               
    double *A = workspace->A;
    double *B = workspace->B;
    double *Sx = workspace->Sx;
    double *Su = workspace->Su;    
    double *xt = workspace->xt;
    double *K = workspace->K;
    double *dKx = workspace->dKx;
    double *dKu = workspace->dKu;
    double *jacX_t = workspace->jacX_t;
    double *jacU_t = workspace->jacU_t;
    double *X_traj = workspace->X_traj;
    double *K_lambda = workspace->K_lambda;
    double *lambda_t = workspace->lambda_t;

    size_t jx = nx*nx;
    size_t ju = nx*nu;
    size_t nw = nx+nu;
    
    double *vde_in[6];
    double *vde_out[2];
    double *ode_in[4];
    double *ode_out[1];
    double *adj_in[5];
    double *adj_out[2];
    
    // initialize
    memcpy(xn, x, nx*sizeof(double));    
    set_zeros(num_stages*nx, K);
    
    if (forw_sens_flag){
        memcpy(Jac_x, Sx, nx*nx*sizeof(double));
        memcpy(Jac_u, Su, nx*nu*sizeof(double));
    }
      
    ode_in[1] = in->u;
    ode_in[2] = in->p;
    ode_in[3] = 0;
    
    if (forw_sens_flag){
        vde_in[1] = in->u;
        vde_in[2] = in->p;
        vde_in[5] = 0;
    }

    if (adj_sens_flag){
        adj_in[1] = in->u;
        adj_in[2] = in->p;
        adj_in[4] = 0;

        for (j=0;j<nu;j++)
            adj_sens[nx+j] = 0.0;
    }
    
    // forward sweep
    for (istep = 0; istep < num_steps; istep++) {
        
        for (s = 0; s < num_stages; s++) {
            memcpy(xt, xn, nx*sizeof(double));
            
            if (forw_sens_flag){
                memcpy(jacX_t, Jac_x, nx*nx*sizeof(double));
                memcpy(jacU_t, Jac_u, nx*nu*sizeof(double));
            }
            for (j = 0; j < s; j++){
                a = A[j * num_stages + s];
                if (a!=0){
                    a *= h;                   
                    for(i=0;i<nx;i++)
                        xt[i]+=a*K[j*nx+i];
                    
                    if (forw_sens_flag){
                        for(i=0;i<jx;i++)
                            jacX_t[i]+=a*dKx[j*jx+i];
                        for(i=0;i<ju;i++)
                            jacU_t[i]+=a*dKu[j*ju+i];
                    }
                                 
                }
            }
            
            ode_in[0] = xt;
            ode_out[0] = K+s*nx;
            f_Fun(ode_in,ode_out);         
                 
            if(forw_sens_flag){
                vde_in[0] = xt;
                vde_in[3] = jacX_t;
                vde_in[4] = jacU_t;
                vde_out[0] = dKx+s*jx;
                vde_out[1] = dKu+s*ju;
                vde_Fun(vde_in, vde_out);
            }

            if (adj_sens_flag){
                memcpy(X_traj+istep*num_stages*nx+s*nx, xt, nx*sizeof(double));
            }
                      
        }
        
        for (s = 0; s < num_stages; s++){
            b = h * B[s];
            for(i=0;i<nx;i++)
                xn[i]+=b*K[s*nx+i];
            
            if(forw_sens_flag){
                for(i=0;i<jx;i++)
                    Jac_x[i]+=b*dKx[s*jx+i];
                for(i=0;i<ju;i++)
                    Jac_u[i]+=b*dKu[s*ju+i];
            }
        }
       
    }

    if (adj_sens_flag){
        memcpy(adj_sens, lambda, nx*sizeof(double));

        for (istep = num_steps-1; istep > -1; istep--) {
            for (s = num_stages-1; s > -1; s--){
                memcpy(lambda_t, adj_sens, nx*sizeof(double));
                if (s<num_stages-1){
                    a = A[(num_stages-s-2) * num_stages + (num_stages-s-1)]*h;
                    daxpy(&nx, &a, K_lambda+(s+1)*nw, &one_i, lambda_t, &one_i);
                }              
                adj_in[0] = X_traj+istep*num_stages*nx+s*nx;
                adj_in[3] = lambda_t;
                adj_out[0] = K_lambda+s*nw;
                adj_ERK_Fun(adj_in, adj_out);
            }
            
            for (s = 0;s < num_stages; s++){
                b = h * B[s];
                daxpy(&nw, &b, K_lambda+s*nw, &one_i, adj_sens, &one_i);
            }           
        }
    }

    return 0;      
}

void sim_erk_workspace_free(sim_opts *opts, sim_erk_workspace *work)
{
    mxFree(work->A);
    mxFree(work->B);
    mxFree(work->Sx);
    mxFree(work->Su);
    mxFree(work->xt);
    mxFree(work->K);
    
    if(opts->forw_sens_flag){
        mxFree(work->dKx);
        mxFree(work->dKu);
        mxFree(work->jacX_t);
        mxFree(work->jacU_t);
    }

    if(opts->adj_sens_flag){
        mxFree(work->X_traj);
        mxFree(work->K_lambda);
        mxFree(work->lambda_t);
    }
    
    mxFree(work);
}