#ifndef IRK_ODE_H_
#define IRK_ODE_H_

#include "mex.h"
#include "sim.h"

typedef struct{
    double *A;
    double *B;
    double *Sx;
    double *Su;
    double *xt;
    double *K;
    double *rG;
    double *JGK;
    double *JGx;
    double *JGu;
    double *JKx;
    double *JKu;
    double *JFK;
    double *JGK_fact_traj;
    double *JGx_traj;
    double *JGu_traj;
    double *lambda_K;
    double **impl_ode_in;
    double **res_out;
    double **jac_x_out;
    double **jac_u_out;
    double **jac_xdot_out;
    size_t *IPIV;
    size_t *IPIV_traj;
    int newton_iter;
    double *tmp_nx_nx;
}sim_irk_ode_workspace;

sim_irk_ode_workspace* sim_irk_ode_workspace_create(sim_opts *opts);

void sim_irk_ode_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_ode_workspace *work);

int sim_irk_ode(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_ode_workspace *work);

void sim_irk_ode_workspace_free(sim_opts *opts, sim_irk_ode_workspace *work);

#endif