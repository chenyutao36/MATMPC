#include "mex.h"
#include "mpc_common.h"

void Block_Fill(size_t m, size_t n, double *Gi, double *G, size_t idm, size_t idn, size_t ldG){
       
    size_t i,j;
    size_t s;
    for (j=0;j<n;j++){
        s = idn*ldG + idm + j*ldG;
        for (i=0;i<m;i++){
            G[s+i] = Gi[j*m+i];
        }
    }
       
}

void Block_Fill_Trans(size_t m, size_t n, double *Gi, double *G, size_t idm, size_t idn, size_t ldG){
       
    size_t i,j;
    size_t s;
    for (j=0;j<m;j++){
        s = idn*ldG + idm + j*ldG;
        for (i=0;i<n;i++){
            G[s+i] = Gi[i*m+j];
        }
    }
       
}

void Block_Get(size_t m, size_t n, double *Gi, double *G, size_t idm, size_t idn, size_t ldG){
       
    size_t i,j;
    size_t s;
    for (j=0;j<n;j++){
        s = idn*ldG + idm + j*ldG;
        for (i=0;i<m;i++){
            Gi[j*m+i] = G[s+i];
        }
    }
       
}

void set_zeros(size_t dim, double *A){
    size_t i;
    for (i=0;i<dim;i++)
        A[i] = 0.0;
}

void print_matrix(double *A, size_t m, size_t n){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++)
            mexPrintf("%4.2f ", A[j*m+i]);
        mexPrintf("\n");
    }
}

void print_vector(double *x, size_t m){
    int i;
    for(i=0;i<m;i++){
        mexPrintf("%4.2f ", x[i]);
        mexPrintf("\n");
    }
}

void regularization(size_t n, double *A, double reg){
    int i;
    for (i=0;i<n;i++)
        if (A[i*n+i]<reg)
            A[i*n+i] = reg;
}