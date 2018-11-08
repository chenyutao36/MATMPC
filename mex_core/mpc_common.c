#include "mex.h"
#include "mpc_common.h"
#include "string.h"

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

// void product_mb(double *P,double *M, size_t a, double *index, size_t r, size_t nu){    // M*T,   M is (a x (N*nu))
//     double *temp = mxCalloc(a*nu,sizeof(double));
//     double *ris  = mxCalloc(a*r*nu,sizeof(double));
//     int i,j,l;
//     
//             
//     for (i=0;i<r;i++){
//         memcpy(temp, M+(int)(index[i])*a*nu, a*nu*sizeof(double));
//         for (j=(int)(index[i])+1;j<(int)(index[i+1]);j++){
//             for (l=0;l<a*nu;l++){
//                 *(temp+l)=*(temp+l)+*(M+j*a*nu+l);
//             }
//         }
//         memcpy(ris+i*a*nu, temp, a*nu*sizeof(double));
//     }
//     memcpy(P,ris,a*r*nu*sizeof(double));
//     mxFree(temp);
//     mxFree(ris);
// }



// void product_mb_Trans(double *P, double *M, size_t b, double *index, size_t r,size_t nu){     // T'* M, M is (N*nu x b) 
//     double *temp = mxCalloc(b*nu,sizeof(double));
//     double *ris  = mxCalloc(b*r*nu,sizeof(double));
//     int i,j,l,k;
//     int N = (int)(index[r]);
//     
//     for (i=0;i<r;i++){
//         for (k=0;k<nu;k++){
//             for(l=0;l<b;l++){
//                 *(temp+l*nu+k)=*(M+l*N*nu+(int)(index[i])*nu+k);
//             }
//         }
//         for (j=(int)(index[i])+1;j<(int)(index[i+1]);j++){
//             for (k=0;k<nu;k++){
//                 for(l=0;l<b;l++){
//                     *(temp+l*nu+k)=*(temp+l*nu+k)+*(M+l*N*nu+j*nu+k);
//                 }
//             }
//         }
//         for (k=0;k<nu;k++){
//             for(l=0;l<b;l++){
//                 *(ris+l*r*nu+i*nu+k)=*(temp+l*nu+k);
//             }
//         }
//     }
//     memcpy(P,ris,b*r*nu*sizeof(double));
//     mxFree(ris);
//     mxFree(temp);    
// }