#ifndef RTI_STEP_COMMON_H_
#define RTI_STEP_COMMON_H_

// #include "mex.h"
// #include "stdlib.h"

typedef struct{
    int nx;
    int nu;
    int np;
    int ny;
    int nyN;
    int nc;
    int ncN;
    int N;    
}rti_step_dims;

typedef struct{
    double **Jac;
    double *Jac_N;
    
    double *L;
    double *w_vec;
    double *W_mat; 
    double *Hi;   
    double *Cci;
    double *CcN;
}rti_step_workspace;

int rti_step_calculate_workspace_size(rti_step_dims *dim);

void *rti_step_cast_workspace(rti_step_dims *dim, void *raw_memory);

#endif