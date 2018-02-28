#include "stdlib.h"
#include "RTI_step_common.h"

int rti_step_calculate_workspace_size(rti_step_dims *dim)
{
    int nx = dim->nx;
    int nu = dim->nu;
    int np = dim->np;
    int ny = dim->ny;
    int nyN = dim->nyN;
    int nc = dim->nc;
    int ncN = dim->ncN;   
    int N = dim->N;
    
    int size = sizeof(rti_step_workspace);
    
    size += 2*sizeof(double *); // Jac
    size += ny*nx*sizeof(double); // Jac[0]
    size += ny*nu*sizeof(double); // Jac[1]
    size += nyN*nx*sizeof(double); // Jac_N
    
    size += (N+1)*nx * sizeof(double); // L
    size += nx*nu*N*N * sizeof(double); // W_mat
    size += nx*N * sizeof(double); // w_vec
    size += nu*nu * sizeof(double); // Hi
    size += nc*nu * sizeof(double); // Cci
    size += ncN*nu * sizeof(double); // CcN
    
    return size;   
}

void *rti_step_cast_workspace(rti_step_dims *dim, void *raw_memory)
{
    int nx = dim->nx;
    int nu = dim->nu;
    int np = dim->np;
    int ny = dim->ny;
    int nyN = dim->nyN;
    int nc = dim->nc;
    int ncN = dim->ncN;   
    int N = dim->N;
    
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
    c_ptr += N*nx*sizeof(double);
       
    workspace->W_mat = (double *)c_ptr;
    c_ptr += nx*nu*N*N*sizeof(double);
        
    workspace->Hi = (double *)c_ptr;
    c_ptr += nu*nu*sizeof(double);
    
    workspace->Cci = (double *)c_ptr;
    c_ptr += nc*nu*sizeof(double);
    
    workspace->CcN = (double *)c_ptr;
    c_ptr += ncN*nu*sizeof(double);
    
    return (void *)workspace;
}