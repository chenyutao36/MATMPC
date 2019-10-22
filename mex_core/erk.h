#ifndef ERK_H_
#define ERK_H_

#include "mex.h"
#include "sim.h"

typedef struct{
    double *A;
    double *B;
    double *Sx;
    double *Su;
         
    double *xt;
    double *K;
    double *dKx;
    double *dKu;
    double *jacX_t;
    double *jacU_t;
    double *X_traj;
    double *K_lambda;
    double *lambda_t;
}sim_erk_workspace;

sim_erk_workspace* sim_erk_workspace_create(sim_opts *opts);

void sim_erk_workspace_init(sim_opts *opts, const mxArray *mem, sim_erk_workspace *work);

int sim_erk(sim_in *in, sim_out *out, sim_opts *opts, sim_erk_workspace *work);

void sim_erk_workspace_free(sim_opts *opts, sim_erk_workspace *work);

#endif