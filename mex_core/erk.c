#include "mex.h"
#include <stdlib.h>
#include "string.h"

#include "sim.h"
#include "erk.h"
#include "mpc_common.h"
#include "casadi_wrapper.h"

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
    
    if (opts->forw_sens){
        work->dKx = mxMalloc(opts->num_stages*opts->nx*opts->nx*sizeof(double));
        work->dKu = mxMalloc(opts->num_stages*opts->nx*opts->nu*sizeof(double));
        work->jacX_t = mxMalloc(opts->nx*opts->nx*sizeof(double));
        work->jacU_t = mxMalloc(opts->nx*opts->nu*sizeof(double));
        
        mexMakeMemoryPersistent(work->dKx);
        mexMakeMemoryPersistent(work->dKu);
        mexMakeMemoryPersistent(work->jacX_t);
        mexMakeMemoryPersistent(work->jacU_t);
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
    
    double h= opts->h;
    size_t nx = opts->nx;
    size_t nu = opts->nu;
    size_t num_stages = opts->num_stages;
    size_t num_steps = opts->num_steps;
    bool forw_sens = opts->forw_sens;
        
    double *x = in->x;
    double *xn = out->xn;
    
    double *Jac_x = out->Sx;
    double *Jac_u = out->Su;   
               
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

    mwSize jx = nx*nx;
    mwSize ju = nx*nu;
    
    double *vde_in[5];
    double *vde_out[2];
    double *ode_in[3];
    double *ode_out[1];
    
    // initialize
    memcpy(xn, x, nx*sizeof(double));    
    set_zeros(num_stages*nx, K);
    
    if (forw_sens){
        memcpy(Jac_x, Sx, nx*nx*sizeof(double));
        memcpy(Jac_u, Su, nx*nu*sizeof(double));
    }
      
    ode_in[1] = in->u;
    ode_in[2] = in->p;
    
    if (forw_sens){
        vde_in[1] = in->u;
        vde_in[2] = in->p;
    }
    
    // forward sweep
    for (istep = 0; istep < num_steps; istep++) {
        
        for (s = 0; s < num_stages; s++) {
            memcpy(xt, xn, nx*sizeof(double));
            
            if (forw_sens){
                memcpy(jacX_t, Jac_x, nx*nx*sizeof(double));
                memcpy(jacU_t, Jac_u, nx*nu*sizeof(double));
            }
            for (j = 0; j < s; j++){
                a = A[j * num_stages + s];
                if (a!=0){
                    a *= h;                   
                    for(i=0;i<nx;i++)
                        xt[i]+=a*K[j*nx+i];
                    
                    if (forw_sens){
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
                 
            if(forw_sens){
                vde_in[0] = xt;
                vde_in[3] = jacX_t;
                vde_in[4] = jacU_t;
                vde_out[0] = dKx+s*jx;
                vde_out[1] = dKu+s*ju;
                vde_Fun(vde_in, vde_out);
            }
                      
        }
        
        for (s = 0; s < num_stages; s++){
            b = h * B[s];
            for(i=0;i<nx;i++)
                xn[i]+=b*K[s*nx+i];
            
            if(forw_sens){
                for(i=0;i<jx;i++)
                    Jac_x[i]+=b*dKx[s*jx+i];
                for(i=0;i<ju;i++)
                    Jac_u[i]+=b*dKu[s*ju+i];
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
    
    if(opts->forw_sens){
        mxFree(work->dKx);
        mxFree(work->dKu);
        mxFree(work->jacX_t);
        mxFree(work->jacU_t);
    }
    
    mxFree(work);
}