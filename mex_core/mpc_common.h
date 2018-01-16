#ifndef MPC_COMMON_H_
#define MPC_COMMON_H_

#include "stdlib.h"

// typedef struct{
//     double *z;
//     double xN;
//     double *y;
//     double *yN;
//     double *od;
//     double *Q;
//     double *QN;
//     double *lb;
//     double *ub;
//     double *lbN;
//     double *ubN;
//     double *x0;
//     double *lbu;
//     double *ubu;
// }mpc_input;
// 
// typedef struct{
//     size_t nx;
//     size_t nu;
//     size_t np;
//     size_t ny;
//     size_t nyN;
//     size_t nc; 
//     size_t ncN;
//     size_t N;
// }mpc_sizes;

void Block_Fill(size_t m, size_t n, double *Gi, double *G,
     size_t idm, size_t idn, size_t ldG);

void set_zeros(size_t dim, double *A);

void print_matrix(double *A, size_t m, size_t n);

void print_vector(double *x, size_t m);

#endif