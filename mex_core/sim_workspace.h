#ifndef SIM_WORKSPACE_H_
#define SIM_WORKSPACE_H_

#include "sim.h"
#include "erk.h"
#include "irk_ode.h"
#include "irk_dae.h"

typedef struct {
    sim_opts *opts;
    sim_in *in;
    sim_out *out;
    sim_erk_workspace *erk;
    sim_irk_ode_workspace *irk_ode;
    sim_irk_dae_workspace *irk_dae;
    int method;
} sim_workspace;

void sim_workspace_prepare(sim_workspace *workspace, const mxArray *mem, int method,
                           size_t nx, size_t nu, size_t nz, bool forward, bool adjoint);
void sim_workspace_free(sim_workspace *workspace);

#endif
