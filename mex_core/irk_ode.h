#ifndef IRK_ODE_H_
#define IRK_ODE_H_

#include "irk_mex.h"

typedef sim_irk_workspace sim_irk_ode_workspace;

sim_irk_ode_workspace *sim_irk_ode_workspace_create(sim_opts *opts);
void sim_irk_ode_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_ode_workspace *work);
int sim_irk_ode(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_ode_workspace *work);
void sim_irk_ode_workspace_free(sim_opts *opts, sim_irk_ode_workspace *work);

#endif
