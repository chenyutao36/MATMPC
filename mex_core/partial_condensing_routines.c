#include "mex.h"
#include "string.h"
#include "blas.h"
// #include "stdlib.h"
#include "partial_condensing_routines.h"
#include "mpc_common.h"

partial_condensing_workspace* partial_condensing_workspace_allocate(size_t Npc, size_t nx, size_t nu, size_t nbx,
        size_t nc)
{
    partial_condensing_workspace* work = mxMalloc(sizeof(partial_condensing_workspace));
    
    work->G = mxCalloc((Npc+1)*nx *(nx+Npc*nu), sizeof(double));
    work->L = mxCalloc((Npc+1)*nx, sizeof(double));
    work->I = mxCalloc( nx*nx, sizeof(double));
    
    work->block11 = mxCalloc(nx*nx, sizeof(double));
    work->block12 = mxCalloc(nx*nx, sizeof(double));
    work->block21 = mxCalloc(nx*nu, sizeof(double));
    work->block22 = mxCalloc(nx*nu, sizeof(double));
    work->block31 = mxCalloc(nx*nu, sizeof(double));
    work->block32 = mxCalloc(nu*nu, sizeof(double));
    work->block41 = mxCalloc(nbx*nx, sizeof(double));
    work->block42 = mxCalloc(nbx*nu, sizeof(double));
    work->block51 = mxCalloc(nx*nx, sizeof(double));
    work->block52 = mxCalloc(nx*nu, sizeof(double));
    work->block61 = mxCalloc(nc*nx, sizeof(double));
    work->block62 = mxCalloc(nc*nu, sizeof(double));
    work->vec1 = mxCalloc(nx, sizeof(double));
    work->vec2 = mxCalloc(nx, sizeof(double));
    
    mexMakeMemoryPersistent(work->G);
    mexMakeMemoryPersistent(work->L);
    mexMakeMemoryPersistent(work->I);
    mexMakeMemoryPersistent(work->block11);
    mexMakeMemoryPersistent(work->block12);
    mexMakeMemoryPersistent(work->block21);
    mexMakeMemoryPersistent(work->block22);
    mexMakeMemoryPersistent(work->block31);
    mexMakeMemoryPersistent(work->block32);
    mexMakeMemoryPersistent(work->block41);
    mexMakeMemoryPersistent(work->block42);
    mexMakeMemoryPersistent(work->block51);
    mexMakeMemoryPersistent(work->block52);
    mexMakeMemoryPersistent(work->block61);
    mexMakeMemoryPersistent(work->block62);
    mexMakeMemoryPersistent(work->vec1);
    mexMakeMemoryPersistent(work->vec2);
    
    mexMakeMemoryPersistent(work);
    
    return work;
}

void partial_condensing_workspace_free(partial_condensing_workspace *work)
{
    mxFree(work->G);
    mxFree(work->L);
    mxFree(work->block11);
    mxFree(work->block12);
    mxFree(work->block21);
    mxFree(work->block22);
    mxFree(work->block31);
    mxFree(work->block32);
    mxFree(work->block41);
    mxFree(work->block42);
    mxFree(work->block51);
    mxFree(work->block52);
    mxFree(work->block61);
    mxFree(work->block62);
    mxFree(work->vec1);
    mxFree(work->vec2);
    
    mxFree(work);
}


void compute_G(partial_condensing_workspace* work, double *Ap, double *Bp, double *A, double *B, size_t nx, size_t nu, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *G = work->G;
    double *I = work->I;      
    double *G_tmp_in1 = work->block11;
    double *G_tmp_out1 = work->block12;
    double *G_tmp_in2 = work->block21;
    double *G_tmp_out2 = work->block22;
    
    for (j=0;j<nx;j++)
        I[j*nx+j] = 1.0;
       
    for(i=0;i<=Npc;i++){
        if (i==0)
            Block_Fill(nx, nx, I, G, 0, 0, (Npc+1)*nx);
        else 
            Block_Fill(nx, nu, B+(i-1)*nx*nu, G, i*nx, nx+(i-1)*nu, (Npc+1)*nx);
        
        for (j=i+1;j<=Npc;j++){
            if (i==0){
                Block_Get(nx, nx, G_tmp_in1, G, (j-1)*nx, 0, (Npc+1)*nx);
                dgemm(nTrans, nTrans, &nx, &nx, &nx, &one_d, A+(j-1)*nx*nx, &nx, G_tmp_in1, &nx, &zero, G_tmp_out1, &nx);
                Block_Fill(nx, nx, G_tmp_out1, G, j*nx, 0, (Npc+1)*nx);
            }else{
                Block_Get(nx, nu, G_tmp_in2, G, (j-1)*nx, nx+(i-1)*nu, (Npc+1)*nx);
                dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, A+(j-1)*nx*nx, &nx, G_tmp_in2, &nx, &zero, G_tmp_out2, &nx);
                Block_Fill(nx, nu, G_tmp_out2, G, j*nx, nx+(i-1)*nu, (Npc+1)*nx);
            }
        }
    }
    
    Block_Get(nx, nx, G_tmp_in1, G, Npc*nx, 0, (Npc+1)*nx);
    memcpy(Ap, G_tmp_in1, nx*nx*sizeof(double));
    
    Block_Get(nx, Npc*nu, Bp, G, Npc*nx, nx, (Npc+1)*nx);
       

}

void compute_L(partial_condensing_workspace* work, double *ap, double *A, double *a, size_t nx, size_t Npc)
{
    int i;
    
    char *nTrans = "N";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *L = work->L;
        
    for(i=0;i<Npc;i++){
        memcpy(L+(i+1)*nx, a+i*nx, nx*sizeof(double)); 
        dgemv(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,L+i*nx,&one_i,&one_d,L+(i+1)*nx,&one_i);
    }
    
    memcpy(ap, L+Npc*nx, nx*sizeof(double));
}

void compute_H(double *H, double *Q, double *S, double *R, double *A, double *B, partial_condensing_workspace* work,
        size_t nx, size_t nu, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N", *Trans="T", *SIDE = "L", *UPLO = "L";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *G = work->G;
    double *W_tmp_in1 = work->block11;
    double *W_tmp_out1 = work->block12;
    double *W_tmp_in2 = work->block21;
    double *W_tmp_out2 = work->block22;
    double *H_tmp1 = work->block31;
    double *H_tmp2 = work->block32;
    double *G_tmp1 = work->block51;
    double *G_tmp2 = work->block52;
    
    Block_Fill(nu, nu, R+(Npc-1)*nu*nu, H, nx+(Npc-1)*nu, nx+(Npc-1)*nu, nx+Npc*nu);
    
    for(i=0;i<Npc;i++){
        if (i==0){
            Block_Get(nx, nx, G_tmp1, G, (Npc-1)*nx, 0, (Npc+1)*nx);
            dgemm(nTrans, nTrans, &nx, &nx, &nx, &one_d, Q+(Npc-1)*nx*nx, &nx, G_tmp1, &nx, &zero, W_tmp_in1, &nx);
            dgemm(Trans, nTrans, &nu, &nx, &nx, &one_d, S+(Npc-1)*nx*nu, &nx, G_tmp1, &nx, &zero, H_tmp1, &nu);  
            Block_Fill(nu, nx, H_tmp1, H, nx+(Npc-1)*nu, 0, nx+Npc*nu);
        }else{
            Block_Get(nx, nu, G_tmp2, G, (Npc-1)*nx, nx+(i-1)*nu, (Npc+1)*nx);
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, Q+(Npc-1)*nx*nx, &nx, G_tmp2, &nx, &zero, W_tmp_in2, &nx);
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, S+(Npc-1)*nx*nu, &nx, G_tmp2, &nx, &zero, H_tmp2, &nu);  
            Block_Fill(nu, nu, H_tmp2, H, nx+(Npc-1)*nu, nx+(i-1)*nu, nx+Npc*nu);
        }
        
        for(j=Npc-1;j>i;j--){
            if (i==0){
                Block_Get(nx, nx, G_tmp1, G, (j-1)*nx, 0, (Npc+1)*nx);
                dgemm(Trans, nTrans, &nu, &nx, &nx, &one_d, S+(j-1)*nx*nu, &nx, G_tmp1, &nx, &zero, H_tmp1, &nu);                     
                dgemm(Trans, nTrans, &nu, &nx, &nx, &one_d, B+(j-1)*nx*nu, &nx, W_tmp_in1, &nx, &one_d, H_tmp1, &nu);
                Block_Fill(nu, nx, H_tmp1, H, nx+(j-1)*nu, 0, nx+Npc*nu);                
                Block_Fill_Trans(nu, nx, H_tmp1, H, 0, nx+(j-1)*nu, nx+Npc*nu);
                
                dgemm(Trans, nTrans, &nx, &nx, &nx, &one_d, A+(j-1)*nx*nx, &nx, W_tmp_in1, &nx, &zero, W_tmp_out1, &nx); 
                dsymm(SIDE, UPLO, &nx, &nx, &one_d, Q+(j-1)*nx*nx, &nx, G_tmp1, &nx, &one_d, W_tmp_out1, &nx);
                memcpy(W_tmp_in1, W_tmp_out1, nx*nx*sizeof(double));
            }else{
                Block_Get(nx, nu, G_tmp2, G, (j-1)*nx, nx+(i-1)*nu, (Npc+1)*nx);
                dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, S+(j-1)*nx*nu, &nx, G_tmp2, &nx, &zero, H_tmp2, &nu);                     
                dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+(j-1)*nx*nu, &nx, W_tmp_in2, &nx, &one_d, H_tmp2, &nu);
                Block_Fill(nu, nu, H_tmp2, H, nx+(j-1)*nu, nx+(i-1)*nu, nx+Npc*nu);
                Block_Fill_Trans(nu, nu, H_tmp2, H, nx+(i-1)*nu, nx+(j-1)*nu, nx+Npc*nu);

                dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, A+(j-1)*nx*nx, &nx, W_tmp_in2, &nx, &zero, W_tmp_out2, &nx); 
                dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+(j-1)*nx*nx, &nx, G_tmp2, &nx, &one_d, W_tmp_out2, &nx);
                memcpy(W_tmp_in2, W_tmp_out2, nx*nu*sizeof(double));
            }         
       }
        
       if(i==0){
           Block_Fill(nx, nx, W_tmp_in1, H, 0, 0, nx+Npc*nu);          
       }else{
           memcpy(H_tmp2, R+(i-1)*nu*nu, nu*nu*sizeof(double));
           dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+(i-1)*nx*nu, &nx, W_tmp_in2, &nx, &one_d, H_tmp2, &nu);
           Block_Fill(nu, nu, H_tmp2, H, nx+(i-1)*nu, nx+(i-1)*nu, nx+Npc*nu);
       }
          
    }
    
}

void compute_g(double *g, double *Q, double *S, double *A, double *B, partial_condensing_workspace* work, double *gx, double *gu,
        size_t nx, size_t nu, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *L = work->L;
    double *w_in = work->vec1;
    double *w_out = work->vec2;
    
    memcpy(w_in, gx+(Npc-1)*nx, nx*sizeof(double));    
    dgemv(nTrans,&nx,&nx,&one_d,Q+(Npc-1)*nx*nx,&nx,L+(Npc-1)*nx,&one_i,&one_d,w_in,&one_i);
    
    memcpy(g+nx+(Npc-1)*nu, gu+(Npc-1)*nu, nu*sizeof(double));
    dgemv(Trans,&nx,&nu,&one_d,S+(Npc-1)*nx*nu,&nx,L+(Npc-1)*nx,&one_i,&one_d,g+nx+(Npc-1)*nu,&one_i);
    
    for(i=Npc-1;i>0;i--){
        memcpy(g+nx+(i-1)*nu, gu+(i-1)*nu, nu*sizeof(double));
        dgemv(Trans,&nx,&nu,&one_d,S+(i-1)*nx*nu,&nx,L+(i-1)*nx,&one_i,&one_d,g+nx+(i-1)*nu,&one_i);
        dgemv(Trans,&nx,&nu,&one_d,B+(i-1)*nx*nu,&nx,w_in,&one_i,&one_d,g+nx+(i-1)*nu,&one_i);
        
        memcpy(w_out, gx+(i-1)*nx, nx*sizeof(double));
        dgemv(Trans,&nx,&nx,&one_d,A+(i-1)*nx*nx,&nx,w_in,&one_i,&one_d,w_out,&one_i);
        dgemv(nTrans,&nx,&nx,&one_d,Q+(i-1)*nx*nx,&nx,L+(i-1)*nx,&one_i,&one_d,w_out,&one_i);        
        memcpy(w_in, w_out, nx*sizeof(double));
    }
    memcpy(g, w_in, nx*sizeof(double));
    
}

void compute_Ccx(double *Ccx, double *Cx, partial_condensing_workspace* work,
        size_t nx, size_t nu, size_t nbx, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
       
    double *G = work->G;
    double *C_tmp1 = work->block41;
    double *C_tmp2 = work->block42;
    double *G_tmp1 = work->block11;
    double *G_tmp2 = work->block12;
    
    for(i=0;i<Npc-1;i++){
        Block_Get(nx, nx, G_tmp1, G, (i+1)*nx, 0, (Npc+1)*nx);
        dgemm(nTrans, nTrans, &nbx, &nx, &nx, &one_d, Cx, &nbx, G_tmp1, &nx, &zero, C_tmp1, &nbx);
        Block_Fill(nbx, nx, C_tmp1, Ccx, i*nbx, 0, (Npc-1)*nbx);
        for(j=i+1;j<Npc;j++){
            Block_Get(nx, nu, G_tmp2, G, j*nx, nx+i*nu, (Npc+1)*nx);
            dgemm(nTrans, nTrans, &nbx, &nu, &nx, &one_d, Cx, &nbx, G_tmp2, &nx, &zero, C_tmp2, &nbx);
            Block_Fill(nbx, nu, C_tmp2, Ccx, (j-1)*nbx, nx+i*nu, (Npc-1)*nbx);
        }
    }
    
}

void compute_ccx(double *lxc, double *uxc, double *lb_dx, double *ub_dx, double *Cx, partial_condensing_workspace* work, 
        size_t nx, size_t nu, size_t nbx, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *L = work->L;
    
    for (i=0;i<Npc-1;i++){
        dgemv(nTrans,&nbx,&nx,&minus_one,Cx,&nbx,L+(i+1)*nx,&one_i,&zero,lxc+i*nbx,&one_i);
        for(j=0;j<nbx;j++){
            uxc[i*nbx+j] = lxc[i*nbx+j]+ ub_dx[i*nbx+j];
            lxc[i*nbx+j] += lb_dx[i*nbx+j];          
        }
    }
}

void compute_Ccg(double *Ccg, double *Cgx, double *Cgu, partial_condensing_workspace* work,
        size_t nx, size_t nu, size_t nc, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
       
    double *G = work->G;
    double *C_tmp1 = work->block61;
    double *C_tmp2 = work->block62;
    double *G_tmp1 = work->block11;
    double *G_tmp2 = work->block12;
    
    for(i=0;i<Npc;i++){
        Block_Get(nx, nx, G_tmp1, G, i*nx, 0, (Npc+1)*nx);
        dgemm(nTrans, nTrans, &nc, &nx, &nx, &one_d, Cgx+i*nc*nx, &nc, G_tmp1, &nx, &zero, C_tmp1, &nc);
        Block_Fill(nc, nx, C_tmp1, Ccg, i*nc, 0, Npc*nc);
        for(j=i+1;j<Npc;j++){
            Block_Get(nx, nu, G_tmp2, G, j*nx, nx+i*nu, (Npc+1)*nx);
            dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, Cgx+j*nc*nx, &nc, G_tmp2, &nx, &zero, C_tmp2, &nc);
            Block_Fill(nc, nu, C_tmp2, Ccg, j*nc, nx+i*nu, Npc*nc);
        }
        
        Block_Fill(nc, nu, Cgu+i*nc*nu, Ccg, i*nc, nx+i*nu, Npc*nc);
    }
    
}

void compute_ccg(double *lgc, double *ugc, double *lc, double *uc, double *Cgx, partial_condensing_workspace* work, 
        size_t nx, size_t nu, size_t nc, size_t Npc)
{
    int i,j;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *L = work->L;
    
    for (i=0;i<Npc;i++){
        dgemv(nTrans,&nc,&nx,&minus_one,Cgx+i*nc*nx,&nc,L+i*nx,&one_i,&zero,lgc+i*nc,&one_i);
        for(j=0;j<nc;j++){
            ugc[i*nc+j] = lgc[i*nc+j]+ uc[i*nc+j];
            lgc[i*nc+j] += lc[i*nc+j];          
        }
    }
}