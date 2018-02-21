#include "stdlib.h"
#include "RTI_step_common.h"

// void rti_step_init_dim(rti_step_dims *dim, size_t nx, size_t nu, size_t np,
//         size_t ny, size_t nyN, size_t nc, size_t ncN, size_t N)
// {
//     dim->nx = nx;
//     dim->nu = nu;
//     dim->np = np;
//     dim->ny = ny;
//     dim->nyN = nyN;
//     dim->nc = nc;
//     dim->ncN = ncN;
//     dim->N = N;
// }

int rti_step_calculate_workspace_size(rti_step_dims *dim)
{
    size_t nx = dim->nx;
    size_t nu = dim->nu;
    size_t np = dim->np;
    size_t ny = dim->ny;
    size_t nyN = dim->nyN;
    size_t nc = dim->nc;
    size_t ncN = dim->ncN;   
    size_t N = dim->N;
    
    int size = sizeof(rti_step_workspace);
    
    size += 2*sizeof(double *); // Jac
    size += ny*nx*sizeof(double); // Jac[0]
    size += ny*nu*sizeof(double); // Jac[1]
    size += nyN*nx*sizeof(double); // Jac_N
    
    size += (N+1)*nx * sizeof(double); // L
    size += 2*nx * sizeof(double); // w_vec, w_vec_dup
    size += 2*nx*nu * sizeof(double); // W_mat, W_mat_dup
    size += nu*nu * sizeof(double); // Hi
    size += nc*nu * sizeof(double); // Cci
    size += ncN*nu * sizeof(double); // CcN
    
//     size += 7*sizeof(mxArray *); // qpoases_in
//     size += 9*sizeof(mxArray *); // qpoases_out
    
    return size;   
}

void *rti_step_cast_workspace(rti_step_dims *dim, void *raw_memory)
{
    size_t nx = dim->nx;
    size_t nu = dim->nu;
    size_t np = dim->np;
    size_t ny = dim->ny;
    size_t nyN = dim->nyN;
    size_t nc = dim->nc;
    size_t ncN = dim->ncN;   
    size_t N = dim->N;
    
    char *c_ptr = (char *)raw_memory;
    
    rti_step_workspace *workspace = (rti_step_workspace *) c_ptr;
    c_ptr += sizeof(rti_step_workspace);
    
    workspace->Jac = (double **)c_ptr;
    c_ptr += 2*sizeof(double *);
    
    workspace->Jac[0] = (double *)c_ptr;
    c_ptr += ny*nx*sizeof(double);
    
    workspace->Jac[1] = (double *)c_ptr;
    c_ptr += ny*nu*sizeof(double);
    
    workspace->Jac_N = (double *)c_ptr;
    c_ptr += nyN*nx*sizeof(double);
    
    workspace->L = (double *)c_ptr;
    c_ptr += (N+1)*nx*sizeof(double);
    
    workspace->w_vec = (double *)c_ptr;
    c_ptr += nx*sizeof(double);
    
    workspace->w_vec_dup = (double *)c_ptr;
    c_ptr += nx*sizeof(double);
    
    workspace->W_mat = (double *)c_ptr;
    c_ptr += nx*nu*sizeof(double);
    
    workspace->W_mat_dup = (double *)c_ptr;
    c_ptr += nx*nu*sizeof(double);
    
    workspace->Hi = (double *)c_ptr;
    c_ptr += nu*nu*sizeof(double);
    
    workspace->Cci = (double *)c_ptr;
    c_ptr += nc*nu*sizeof(double);
    
    workspace->CcN = (double *)c_ptr;
    c_ptr += ncN*nu*sizeof(double);
    
    return (void *)workspace;
}