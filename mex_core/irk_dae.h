#ifndef IRK_DAE_H_
#define IRK_DAE_H_

#include "irk_mex.h"

typedef sim_irk_workspace sim_irk_dae_workspace;

sim_irk_dae_workspace *sim_irk_dae_workspace_create(sim_opts *opts);
void sim_irk_dae_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_dae_workspace *work);
int sim_irk_dae(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_dae_workspace *work);
void sim_irk_dae_workspace_free(sim_opts *opts, sim_irk_dae_workspace *work);

#endif
