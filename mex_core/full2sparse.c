#include "mex.h"
#include "string.h"

#include "mpc_common.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{  
    /*Inputs*/
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *Cgx = mxGetPr( mxGetField(prhs[0], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[0], 0, "Cgu") );
    double *CgN = mxGetPr( mxGetField(prhs[0], 0, "CgN") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[0], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[0], 0, "uc") );
    double *lb_dx = mxGetPr( mxGetField(prhs[0], 0, "lb_dx") );
    double *ub_dx = mxGetPr( mxGetField(prhs[0], 0, "ub_dx") );
    double *lb_du = mxGetPr( mxGetField(prhs[0], 0, "lb_du") );
    double *ub_du = mxGetPr( mxGetField(prhs[0], 0, "ub_du") );

    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );

    double *H = mxGetPr( mxGetField(prhs[0], 0, "sparse_H") );
    double *g = mxGetPr( mxGetField(prhs[0], 0, "sparse_g") );  
    double *G = mxGetPr( mxGetField(prhs[0], 0, "sparse_G") );
    double *dG = mxGetPr( mxGetField(prhs[0], 0, "sparse_dG") );
    double *dB = mxGetPr( mxGetField(prhs[0], 0, "sparse_dB") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "sparse_ub") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "sparse_lb") );
    double *minus_eye = mxGetPr( mxGetField(prhs[0], 0, "sparse_minus_eye") );

    int i,j;
    size_t nz= nx+nu;
    size_t nw = N*nu + (N+1)*nx;
    size_t neq = (N+1)*nx;
    size_t nineq = N*nc+ncN;

    for (i=0;i<N;i++){
        Block_Fill(nx, nx, Q+i*nx*nx, H, i*nx, i*nx, nw);
        Block_Fill(nx, nu, S+i*nx*nu, H, i*nx, neq+i*nu, nw);
        Block_Fill_Trans(nx, nu, S+i*nx*nu, H, neq+i*nu, i*nx, nw);
        Block_Fill(nu, nu, R+i*nu*nu, H, neq+i*nu, neq+i*nu, nw);

        Block_Fill(nc, nx, Cgx+i*nc*nx, dB, i*nc, i*nx, nineq);
        Block_Fill(nc, nu, Cgu+i*nc*nu, dB, i*nc, neq+i*nu, nineq);

        Block_Fill(nx, nx, A+i*nx*nx, dG, (i+1)*nx, i*nx, neq);
        Block_Fill(nx, nx, minus_eye, dG, (i+1)*nx, (i+1)*nx, neq);
        Block_Fill(nx, nu, B+i*nx*nu, dG, (i+1)*nx, neq+i*nu, neq);

    }
    Block_Fill(nx, nx, Q+N*nx*nx, H, N*nx, N*nx, nw);
    Block_Fill(ncN, nx, CgN, dB, N*nc, N*nx, nineq);

    memcpy(g, gx, neq*sizeof(double));
    memcpy(g+neq, gu, N*nu*sizeof(double));

    for(j=0;j<nx;j++)
        G[j] = -ds0[j];
    memcpy(G+nx, a, N*nx*sizeof(double));

    memcpy(lb, lb_du, N*nu*sizeof(double));
    memcpy(ub, ub_du, N*nu*sizeof(double));
    memcpy(lb+N*nu, lb_dx, N*nbx*sizeof(double));
    memcpy(ub+N*nu, ub_dx, N*nbx*sizeof(double));
    memcpy(lb+N*nu+N*nbx, lc, nineq*sizeof(double));
    memcpy(ub+N*nu+N*nbx, uc, nineq*sizeof(double));
}