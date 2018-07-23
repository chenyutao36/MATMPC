#ifndef PARTIAL_CONDENSING_ROUTINES_H_
#define PARTIAL_CONDENSING_ROUTINES_H_

typedef struct
{
    double *G;
    double *L;
    double *I;
    
    double *block11; // nx,nx
    double *block12; // nx,nx
    double *block21; // nx,nu
    double *block22; // nx,nu
    double *block31; // nx,nu
    double *block32; // nu,nu
    double *block41; // nbx,nx
    double *block42; // nbx,nu
    double *block51; // nx,nx
    double *block52; // nx,nu
    double *block61; // nc,nx
    double *block62; // nc,nu
    double *vec1;
    double *vec2;
}partial_condensing_workspace;

partial_condensing_workspace* partial_condensing_workspace_allocate(size_t Npc, size_t nx, size_t nu, size_t nbx,
        size_t nc);

void partial_condensing_workspace_free(partial_condensing_workspace *work);

void compute_G(partial_condensing_workspace* work, double *Ap, double *Bp, double *A, double *B, size_t nx, size_t nu, size_t Npc);

void compute_L(partial_condensing_workspace* work, double *ap, double *A, double *a, size_t nx, size_t Npc);

void compute_H(double *H, double *Q, double *S, double *R, double *A, double *B, partial_condensing_workspace* work,
        size_t nx, size_t nu, size_t Npc);

void compute_g(double *g, double *Q, double *S, double *A, double *B, partial_condensing_workspace* work, double *gx, double *gu,
        size_t nx, size_t nu, size_t Npc);

void compute_Ccx(double *Ccx, double *Cx, partial_condensing_workspace* work,
        size_t nx, size_t nu, size_t nbx, size_t Npc);

void compute_ccx(double *lxc, double *uxc, double *lb_dx, double *ub_dx, double *Cx, partial_condensing_workspace* work, 
        size_t nx, size_t nu, size_t nbx, size_t Npc);

void compute_Ccg(double *Ccg, double *Cgx, double *Cgu, partial_condensing_workspace* work,
        size_t nx, size_t nu, size_t nc, size_t Npc);

void compute_ccg(double *lgc, double *ugc, double *lc, double *uc, double *Cgx, partial_condensing_workspace* work, 
        size_t nx, size_t nu, size_t nc, size_t Npc);

#endif