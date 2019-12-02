#ifndef IRK_DAE_H_
#define IRK_DAE_H_

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
    double **impl_dae_in;
    double **res_out;
    double **jac_f_x_out;
    double **jac_f_u_out;
    double **jac_f_xdot_out;
    double **jac_f_z_out;
    double **jac_g_x_out;
    double **jac_g_u_out;
    double **jac_g_z_out;
    size_t *IPIV;
    size_t *IPIV_traj;
    int newton_iter;
    double *tmp_nx_nx;
    double *tmp_nx_nz;
    double *tmp_block;
}sim_irk_dae_workspace;

sim_irk_dae_workspace* sim_irk_dae_workspace_create(sim_opts *opts);

void sim_irk_dae_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_dae_workspace *work);

int sim_irk_dae(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_dae_workspace *work);

void sim_irk_dae_workspace_free(sim_opts *opts, sim_irk_dae_workspace *work);

#endif