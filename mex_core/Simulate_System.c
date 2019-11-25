
#include "mex.h"
#include "string.h"

#include "sim.h"
#include "erk.h"
#include "irk_ode.h"
// #include "casadi_wrapper.h"
// #include "mpc_common.h"

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
static sim_irk_ode_workspace *irk_ode_workspace = NULL;
static bool mem_alloc = false;

void exitFcn(){
    if (erk_workspace!=NULL)
        sim_erk_workspace_free(opts, erk_workspace);
    if (irk_ode_workspace!=NULL)
        sim_irk_ode_workspace_free(opts, irk_ode_workspace);
    if (opts!=NULL)
        sim_opts_free(opts);
    if (in!=NULL)
        sim_in_free(in);
    if (out!=NULL)
        sim_out_free(out);
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *x = mxGetPr(prhs[0]);
    double *u = mxGetPr(prhs[1]);
    double *z = mxGetPr(prhs[2]);
    double *od = mxGetPr(prhs[3]);

    size_t nx = mxGetScalar( mxGetField(prhs[5], 0, "nx") );
    size_t nz = mxGetScalar( mxGetField(prhs[5], 0, "nz") );
    // size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    // size_t np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    // size_t ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    // size_t nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    // size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    // size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    // size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    // double *nbx_idx = mxGetPr( mxGetField(prhs[1], 0, "nbx_idx") );
    // size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );    
    int sim_method = mxGetScalar( mxGetField(prhs[4], 0, "sim_method") );

    plhs[0] = mxCreateDoubleMatrix(nx+nz, 1, mxREAL);
    double *x_out = mxGetPr(plhs[0]);
              
    if (!mem_alloc){
        opts = sim_opts_create(prhs[4]);
        opts->forw_sens_flag = false;
        opts->adj_sens_flag = false;
        in = sim_in_create(opts);              
        out = sim_out_create(opts);     
        switch(sim_method){
            case 1:                         
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, prhs[4], erk_workspace);                
                break;
            case 2:          
                irk_ode_workspace = sim_irk_ode_workspace_create(opts);               
                sim_irk_ode_workspace_init(opts, prhs[4], irk_ode_workspace);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }  
                
        mem_alloc=true;
        mexAtExit(exitFcn);
    }
    
        
     // integration                      
    switch(sim_method){
        case 1:
            in->x = x;
            in->u = u;
            in->p = od;
            out->xn = x_out;
            sim_erk(in, out, opts, erk_workspace);
            break;
        case 2:
            in->x = x;
            in->u = u;
            in->p = od;
            in->z = z;
            out->xn = x_out;
            sim_irk_ode(in, out, opts, irk_ode_workspace);
            break;
        default:
            mexErrMsgTxt("Please choose a supported integrator");
            break;
    }
        
}