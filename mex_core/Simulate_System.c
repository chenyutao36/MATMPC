
#include "mex.h"
#include "string.h"

#include "sim_workspace.h"

static sim_workspace integrator;

static void exitFcn(void){
    sim_workspace_free(&integrator);
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *x = mxGetPr(prhs[0]);
    double *u = mxGetPr(prhs[1]);
    double *z = mxGetPr(prhs[2]);
    double *od = mxGetPr(prhs[3]);

    size_t nx = mxGetScalar( mxGetField(prhs[5], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[5], 0, "nu") );
    size_t nz = mxGetScalar( mxGetField(prhs[5], 0, "nz") );
    int sim_method = mxGetScalar( mxGetField(prhs[4], 0, "sim_method") );

    plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL);
    double *x_out = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(nz, 1, mxREAL);
    double *z_out = mxGetPr(plhs[1]);
              
    if (integrator.opts == NULL)
        mexAtExit(exitFcn);
    sim_workspace_prepare(&integrator, prhs[4], sim_method, nx, nu, nz, false, false);

     // integration                      
    switch(sim_method){
        case 1:
            integrator.in->x = x;
            integrator.in->u = u;
            integrator.in->p = od;
            integrator.out->xn = x_out;
            sim_erk(integrator.in, integrator.out, integrator.opts, integrator.erk);
            break;
        case 2:
            integrator.in->x = x;
            integrator.in->u = u;
            integrator.in->p = od;
            integrator.in->z = z;
            integrator.out->xn = x_out;
            sim_irk_ode(integrator.in, integrator.out, integrator.opts, integrator.irk_ode);
            break;
        case 3:
            integrator.in->x = x;
            integrator.in->u = u;
            integrator.in->p = od;
            integrator.in->z = z;
            integrator.out->xn = x_out;
            integrator.out->zn = z_out;
            sim_irk_dae(integrator.in, integrator.out, integrator.opts, integrator.irk_dae);
            break;
        default:
            mexErrMsgTxt("Please choose a supported integrator");
            break;
    }
        
}