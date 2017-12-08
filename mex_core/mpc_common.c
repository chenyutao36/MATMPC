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

void set_zeros(size_t dim, double *A){
    size_t i;
    for (i=0;i<dim;i++)
        A[i] = 0.0;
}