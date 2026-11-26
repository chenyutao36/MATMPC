#ifndef MATMPC_IRK_MEX_H
#define MATMPC_IRK_MEX_H

#include "irk_kernel.h"
#include "sim.h"

typedef struct {
    double *A, *B, *Sx, *Su, *JFK;
    int newton_iter;
    matmpc_irk_workspace scratch;
} sim_irk_workspace;

sim_irk_workspace *sim_irk_mex_create(const sim_opts *opts, size_t doubles, size_t pivots);
void sim_irk_mex_init(const sim_opts *opts, const mxArray *mem, sim_irk_workspace *work);
void sim_irk_mex_free(sim_irk_workspace *work);
matmpc_irk_config sim_irk_mex_config(const sim_opts *opts, const sim_irk_workspace *work);
void sim_irk_mex_check_status(matmpc_irk_status status, const char *integrator);

#endif
