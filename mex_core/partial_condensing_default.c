#include "mex.h"
#include "string.h"
#include "partial_condensing_routines.h"
#include "mpc_common.h"

static partial_condensing_workspace *work = NULL;

void exitFcn(){
    if (work!=NULL)
        partial_condensing_workspace_free(work);
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *Cgx = mxGetPr( mxGetField(prhs[0], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[0], 0, "Cgu") );
//     double *CgN = mxGetPr( mxGetField(prhs[0], 0, "CgN") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[0], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[0], 0, "uc") );
    double *lb_dx = mxGetPr( mxGetField(prhs[0], 0, "lb_dx") );
    double *ub_dx = mxGetPr( mxGetField(prhs[0], 0, "ub_dx") );
    
    double *idxc = mxGetPr( mxGetField(prhs[0], 0, "idxc") );
        
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
//     size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    size_t N2 = mxGetScalar( mxGetField(prhs[1], 0, "N2") );
    
    size_t Npc = (size_t) N/N2;
            
    plhs[0] = mxCreateDoubleMatrix(nx+Npc*nu, (nx+Npc*nu)*N2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nx+Npc*nu, N2, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((Npc-1)*nbx, (nx+Npc*nu)*N2, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nx, nx*N2, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nx, Npc*nu*N2, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nx, N2, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(nbx*(Npc-1)*N2, 1, mxREAL);
    plhs[7] = mxCreateDoubleMatrix(nbx*(Npc-1)*N2, 1, mxREAL);
    
    plhs[8] = mxCreateDoubleMatrix(nc*Npc, (nx+Npc*nu)*N2, mxREAL);
    plhs[9] = mxCreateDoubleMatrix(nc*Npc*N2, 1, mxREAL);
    plhs[10] = mxCreateDoubleMatrix(nc*Npc*N2, 1, mxREAL);
    
    double *Hp = mxGetPr(plhs[0]);
    double *gp = mxGetPr(plhs[1]);
    double *Ccx = mxGetPr(plhs[2]);
    double *Ap = mxGetPr(plhs[3]);
    double *Bp = mxGetPr(plhs[4]);
    double *ap = mxGetPr(plhs[5]);
    double *lxc = mxGetPr(plhs[6]);
    double *uxc = mxGetPr(plhs[7]);
    
    double *Ccg = mxGetPr(plhs[8]);
    double *lgc = mxGetPr(plhs[9]);
    double *ugc = mxGetPr(plhs[10]);
            
    if (work==NULL){
        work = partial_condensing_workspace_allocate(Npc, nx, nu, nbx,nc);
        mexAtExit(exitFcn);
    }
    
    int i;
    
    for(i=0;i<N2;i++){        
        compute_G(work, Ap+i*nx*nx, Bp+i*nx*Npc*nu, A+i*Npc*nx*nx, B+i*Npc*nx*nu, nx, nu, Npc);              
        compute_H(Hp+i*(nx+Npc*nu)*(nx+Npc*nu), Q+i*Npc*nx*nx, S+i*Npc*nx*nu, R+i*Npc*nu*nu, A+i*Npc*nx*nx, B+i*Npc*nx*nu, work, nx, nu, Npc);       
        compute_Ccx(Ccx+i*(Npc-1)*nbx*(nx+Npc*nu), Cx, work, nx, nu, nbx, Npc); 
        compute_Ccg(Ccg+i*Npc*nc*(nx+Npc*nu), Cgx+i*Npc*nc*nx, Cgu+i*Npc*nc*nu, work, nx, nu, nc, Npc);       
        compute_L(work, ap+i*nx, A+i*Npc*nx*nx, a+i*Npc*nx, nx, Npc);    
        compute_g(gp+i*(nx+Npc*nu), Q+i*Npc*nx*nx, S+i*Npc*nx*nu, A+i*Npc*nx*nx, B+i*Npc*nx*nu, work, gx+i*Npc*nx, gu+i*Npc*nu, nx, nu, Npc);            
        compute_ccx(lxc+i*(Npc-1)*nbx, uxc+i*(Npc-1)*nbx, lb_dx+i*Npc*nbx, ub_dx+i*Npc*nbx, Cx, work, nx, nu, nbx, Npc);
        compute_ccg(lgc+i*Npc*nc, ugc+i*Npc*nc, lc+i*Npc*nc, uc+i*Npc*nc, Cgx+i*Npc*nc*nx, work, nx, nu, nc, Npc);
    }
    
}