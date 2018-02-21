#ifndef RTI_STEP_COMMON_H_
#define RTI_STEP_COMMON_H_

#include "mex.h"
#include "stdlib.h"

typedef struct{
    size_t nx;
    size_t nu;
    size_t np;
    size_t ny;
    size_t nyN;
    size_t nc;
    size_t ncN;
    size_t N;    
}rti_step_dims;

// void rti_step_init_dim(rti_step_dims *dim, size_t nx, size_t nu, size_t np,
//         size_t ny, size_t nyN, size_t nc, size_t ncN, size_t N);

typedef struct{
    double **Jac;
    double *Jac_N;
    
    double *L;
    double *w_vec;
    double *W_mat;
    double *w_vec_dup; 
    double *W_mat_dup;   
    double *Hi;   
    double *Cci;
    double *CcN;
    
//     mxArray **qpoases_in;
//     mxArray **qpoases_out;
}rti_step_workspace;

int rti_step_calculate_workspace_size(rti_step_dims *dim);

void *rti_step_cast_workspace(rti_step_dims *dim, void *raw_memory);

#endif