#ifndef SIM_H_
#define SIM_H_

#include "mex.h"

typedef struct{
    double *xt;
    double *K;
    double *dKx;
    double *dKu;
    double *jacX_t;
    double *jacU_t;
}sim_erk_workspace;

int sim_erk_calculate_workspace_size(const mxArray *mem, bool forw_sens);

void *sim_erk_cast_workspace(const mxArray *mem, bool forw_sens, void *raw_memory);

int sim_erk(double **in, double **out, double **Jac, const mxArray *mem, bool forw_sens, void *work_);

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

int sim_irk_calculate_workspace_size(const mxArray *mem, bool forw_sens);

void *sim_irk_cast_workspace(const mxArray *mem, bool forw_sens, void *raw_memory);

int sim_irk(double **in, double **out, double **Jac, const mxArray *mem, bool forw_sens, void *work_);

#endif