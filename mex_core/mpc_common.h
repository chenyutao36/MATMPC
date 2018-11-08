#ifndef MPC_COMMON_H_
#define MPC_COMMON_H_

#include "stdlib.h"

void Block_Fill(size_t m, size_t n, double *Gi, double *G,
     size_t idm, size_t idn, size_t ldG);

void Block_Fill_Trans(size_t m, size_t n, double *Gi, double *G,
     size_t idm, size_t idn, size_t ldG);

void Block_Get(size_t m, size_t n, double *Gi, double *G, size_t idm, size_t idn, size_t ldG);

void set_zeros(size_t dim, double *A);

void print_matrix(double *A, size_t m, size_t n);

void print_vector(double *x, size_t m);

void regularization(size_t n, double *A, double reg);

// void product_mb(double *P,double *M, size_t a, double *index, size_t r,size_t nu);    // P = M*T,   M is (a x (N*nu))
// 
// void product_mb_Trans(double *P,double *M, size_t b, double *index, size_t r,size_t nu);     // P = T'* M, M is (N*nu x b) 

#endif