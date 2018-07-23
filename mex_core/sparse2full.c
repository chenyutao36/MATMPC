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
    double *lc2 = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *uc2 = mxGetPr( mxGetField(prhs[0], 0, "uc") );

    /*mem */
    double *Q = mxGetPr( mxGetField(prhs[1], 0, "Q") );
    double *gx = mxGetPr( mxGetField(prhs[1], 0, "gx") );
    double *lb_dx = mxGetPr( mxGetField(prhs[1], 0, "lb_dx") );
    double *ub_dx = mxGetPr( mxGetField(prhs[1], 0, "ub_dx") );
    double *lc = mxGetPr( mxGetField(prhs[1], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[1], 0, "uc") );

    /*settings2*/
    size_t nu2 = mxGetScalar( mxGetField(prhs[2], 0, "nu") );
    size_t nc2 = mxGetScalar( mxGetField(prhs[2], 0, "nc") );
    size_t N2 = mxGetScalar( mxGetField(prhs[2], 0, "N") );
    size_t Npc = mxGetScalar( mxGetField(prhs[2], 0, "Nc") );

    /*settings*/
    size_t nx = mxGetScalar( mxGetField(prhs[3], 0, "nx") );
    size_t nbx = mxGetScalar( mxGetField(prhs[3], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[3], 0, "N") );
    size_t nc = mxGetScalar( mxGetField(prhs[3], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[3], 0, "ncN") );

    double *Hp = mxGetPr(prhs[4]);
    double *gp = mxGetPr(prhs[5]);
    double *Ccx = mxGetPr(prhs[6]);
    double *Ccg = mxGetPr(prhs[7]);
    double *lxc = mxGetPr(prhs[8]);
    double *uxc = mxGetPr(prhs[9]);
    double *lgc = mxGetPr(prhs[10]);
    double *ugc = mxGetPr(prhs[11]);

    int i,j;
    size_t nz= nx+nu2;
    size_t nw = N*nu2 + (N+1)*nx;
    size_t neq = (N+1)*nx;
    
    for (i=0;i<N2;i++){
        Block_Get(nx, nx, Q2+i*nx*nx, Hp, 0, i*nz, nz);
        Block_Get(nx, nu2, S+i*nx*nu2, Hp, 0, i*nz+nx, nz);
        Block_Get(nu2, nu2, R+i*nu2*nu2, Hp, nx, i*nz+nx, nz);

        memcpy(gx2+i*nx, gp+i*nz, nx*sizeof(double));
        memcpy(gu+i*nu2, gp+i*nz+nx, nu2*sizeof(double));
        
        Block_Fill(nbx*(Npc-1), nx, Ccx+i*nbx*(Npc-1)*nz, Cgx, 0, i*nx, nc2);
        Block_Fill(nbx*(Npc-1), nu2, Ccx+i*nbx*(Npc-1)*nz+nbx*(Npc-1)*nx, Cgu, 0, i*nu2, nc2);           
        memcpy(lb_dx2+i*nbx, lb_dx+i*Npc*nbx, nbx*sizeof(double));
        memcpy(ub_dx2+i*nbx, ub_dx+i*Npc*nbx, nbx*sizeof(double));
      
        Block_Fill(nc*Npc, nx, Ccg+i*nc*Npc*nz, Cgx, nbx*(Npc-1), i*nx, nc2);
        Block_Fill(nc*Npc, nu2, Ccg+i*nc*Npc*nz+nc*Npc*nx, Cgu, nbx*(Npc-1), i*nu2, nc2);          
        memcpy(lc2+i*nc2, lxc+i*(Npc-1)*nbx, (Npc-1)*nbx*sizeof(double));
        memcpy(uc2+i*nc2, uxc+i*(Npc-1)*nbx, (Npc-1)*nbx*sizeof(double));
        memcpy(lc2+i*nc2+(Npc-1)*nbx, lgc+i*Npc*nc, Npc*nc*sizeof(double));
        memcpy(uc2+i*nc2+(Npc-1)*nbx, ugc+i*Npc*nc, Npc*nc*sizeof(double));

    }
    memcpy(Q2+N2*nx*nx, Q+N*nx*nx, nx*nx*sizeof(double));
    memcpy(gx2+N2*nx, gx+N*nx, nx*sizeof(double));
    
    memcpy(lc2+N2*nc2, lc+N*nc, ncN*sizeof(double));
    memcpy(uc2+N2*nc2, uc+N*nc, ncN*sizeof(double));
    
}