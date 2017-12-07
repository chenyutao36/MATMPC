#ifndef SIM_H_
#define SIM_H_

#include "mex.h"

void Block_Fill(mwSize m, mwSize n, double *Gi, double *G, mwSize idm, mwSize idn, mwSize ldG);

void set_zeros(int dim, double *A);

typedef struct{
    double *xt;
    double *K;
    double *dKx;
    double *dKu;
    double *jacX_t;
    double *jacU_t;
}sim_erk_workspace;

int sim_erk_calculate_workspace_size(mxArray *mem, bool forw_sens);

void *sim_erk_cast_workspace(mxArray *mem, bool forw_sens, void *raw_memory);

typedef struct{
    double *xt;
    double *K;
    double *rG;
    double *JGK;
    double *JGx;
    double *JGu;
    double *JKx;
    double *JKu;
    double **impl_ode_in;
    double **impl_ode_out;
    mwIndex *IPIV;
}sim_irk_workspace;

int sim_irk_calculate_workspace_size(mxArray *mem);

void *sim_irk_cast_workspace(mxArray *mem, void *raw_memory);

#endif