#include "mex.h"
#include "string.h"

#include "mpc_common.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{  
    /*mem2*/
    double *Q2 = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *Cgx = mxGetPr( mxGetField(prhs[0], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[0], 0, "Cgu") );
    double *gx2 = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[0], 0, "gu") );   
    double *lb_dx2 = mxGetPr( mxGetField(prhs[0], 0, "lb_dx") );
    double *ub_dx2 = mxGetPr( mxGetField(prhs[0], 0, "ub_dx") );

    /*mem */
    double *Q = mxGetPr( mxGetField(prhs[1], 0, "Q") );
    double *gx = mxGetPr( mxGetField(prhs[1], 0, "gx") );
    double *lb_dx = mxGetPr( mxGetField(prhs[1], 0, "lb_dx") );
    double *ub_dx = mxGetPr( mxGetField(prhs[1], 0, "ub_dx") );

    /*settings2*/
    size_t nx = mxGetScalar( mxGetField(prhs[2], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[2], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[2], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[2], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[2], 0, "nbx") );
    size_t N2 = mxGetScalar( mxGetField(prhs[2], 0, "N") );
    size_t Npc = mxGetScalar( mxGetField(prhs[2], 0, "Nc") );

    /*settings*/
    size_t N = mxGetScalar( mxGetField(prhs[3], 0, "N") );

    double *Hp = mxGetPr(prhs[4]);
    double *gp = mxGetPr(prhs[5]);
    double *Ccx = mxGetPr(prhs[6]);

    int i,j;
    size_t nz= nx+nu;
    size_t nw = N*nu + (N+1)*nx;
    size_t neq = (N+1)*nx;
    size_t nineq = N*nc+ncN;

    for (i=0;i<N2;i++){
        Block_Get(nx, nx, Q2+i*nx*nx, Hp, 0, i*nz, nz);
        Block_Get(nx, nu, S+i*nx*nu, Hp, 0, i*nz+nx, nz);
        Block_Get(nu, nu, R+i*nu*nu, Hp, nx, i*nz+nx, nz);

        memcpy(gx2+i*nx, gp+i*nz, nx*sizeof(double));
        memcpy(gu+i*nu, gp+i*nz+nx, nu*sizeof(double));

        memcpy(Cgx+i*nx*nc, Ccx+i*nc*nz, nx*nc*sizeof(double));
        memcpy(Cgu+i*nu*nc, Ccx+i*nc*nz+nc*nx, nu*nc*sizeof(double));

        memcpy(lb_dx2+i*nbx, lb_dx+i*Npc*nbx, nbx*sizeof(double));
        memcpy(ub_dx2+i*nbx, ub_dx+i*Npc*nbx, nbx*sizeof(double));

    }
    memcpy(Q2+N2*nx*nx, Q+N*nx*nx, nx*nx*sizeof(double));
    memcpy(gx2+N2*nx, gx+N*nx, nx*sizeof(double));

}