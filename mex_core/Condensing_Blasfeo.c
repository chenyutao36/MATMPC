
#include "mex.h"
#include <stdbool.h>
#include "string.h"

#include "Condensing_Blasfeo.h"
#include "mpc_common.h"

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_kernel.h>
#include <blasfeo_d_blas.h>

static void *mem;
static Condensing_Blasfeo_workspace *workspace;
static bool mem_alloc = false;

int align_char_to(int num, char **c_ptr) {
    size_t s_ptr = (size_t)*c_ptr;
    s_ptr = (s_ptr + num - 1) / num * num;
    int offset = num - (int)(s_ptr - (size_t)(*c_ptr));
    *c_ptr = (char *)s_ptr;
    return offset;
}

void assign_blasfeo_dmat_mem(int m, int n, struct blasfeo_dmat *sA, char **ptr)
{
    blasfeo_create_dmat(m, n, sA, *ptr);
    *ptr += sA->memsize;
}

void assign_blasfeo_dvec_mem(int n, struct blasfeo_dvec *sv, char **ptr)
{
    blasfeo_create_dvec(n, sv, *ptr);
    *ptr += sv->memsize;
}

int Condensing_Blasfeo_workspace_calculate_size(int nx, int nu, int N, int nc, int ncN, int nbx)
{
    int size = sizeof(Condensing_Blasfeo_workspace);

    size += 11*sizeof(struct blasfeo_dmat);
    size += 6*sizeof(struct blasfeo_dvec); 
    
    size += blasfeo_memsize_dmat(nu+nx, nx*(N+1));       
    size += blasfeo_memsize_dmat(nu+nx, nx*N);        
    size += blasfeo_memsize_dmat(nu*(N+1), (nu+nx)*(N+1));        
    size += blasfeo_memsize_dmat(nu*N, nx*N);        
    size += blasfeo_memsize_dmat(nu, nu*N);         
    size += blasfeo_memsize_dmat(nc, nx*N);         
    size += blasfeo_memsize_dmat(ncN, nx);         
    size += blasfeo_memsize_dmat(nbx, nx);           
    size += blasfeo_memsize_dmat(N*nu, nu*N);       
    size += blasfeo_memsize_dmat(N*nc+ncN, nu*N);        
//     size += blasfeo_memsize_dmat((N+1)*nbx, nu*N);
    size += blasfeo_memsize_dmat(N*nbx, nu*N);
               
    size += blasfeo_memsize_dvec((nx+nu)*(N+1)); 
    size += blasfeo_memsize_dvec(nx*N);       
    size += blasfeo_memsize_dvec(nx*(N+1));       
    size += blasfeo_memsize_dvec((nu+nx)*(N+1));       
    size += blasfeo_memsize_dvec(N*nc+ncN);       
//     size += blasfeo_memsize_dvec((N+1)*nbx);
    size += blasfeo_memsize_dvec(N*nbx);

    size = (size + 64 - 1) / 64 * 64;
    size += 1 * 64;

    return size;
}

Condensing_Blasfeo_workspace* Condensing_Blasfeo_workspace_cast(int nx, int nu, int N, int nc, int ncN, int nbx, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    Condensing_Blasfeo_workspace *workspace = (Condensing_Blasfeo_workspace *) c_ptr;
    c_ptr += sizeof(Condensing_Blasfeo_workspace);

    workspace->StQ = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->BtAt = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->HtWt = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Gt = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Rt = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Cgx_dmat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->CgN_dmat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Cx_dmat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Hc_dmat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Ccg_dmat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    workspace->Ccx_dmat = (struct blasfeo_dmat *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dmat);
    
    workspace->rq = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    workspace->a_dvec = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    workspace->L = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    workspace->gcw = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    workspace->bg = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);
    workspace->bx = (struct blasfeo_dvec *)c_ptr;
    c_ptr += sizeof(struct blasfeo_dvec);

    align_char_to(64, &c_ptr);
    
    assign_blasfeo_dmat_mem(nu+nx, nx*(N+1), workspace->StQ, &c_ptr);       
    assign_blasfeo_dmat_mem(nu+nx, nx*N, workspace->BtAt, &c_ptr);        
    assign_blasfeo_dmat_mem(nu*(N+1), (nu+nx)*N, workspace->HtWt, &c_ptr);        
    assign_blasfeo_dmat_mem(nu*N, nx*N, workspace->Gt, &c_ptr);        
    assign_blasfeo_dmat_mem(nu, nu*N, workspace->Rt, &c_ptr);         
    assign_blasfeo_dmat_mem(nc, nx*N, workspace->Cgx_dmat, &c_ptr);         
    assign_blasfeo_dmat_mem(ncN, nx, workspace->CgN_dmat, &c_ptr);         
    assign_blasfeo_dmat_mem(nbx, nx, workspace->Cx_dmat, &c_ptr);           
    assign_blasfeo_dmat_mem(N*nu, nu*N, workspace->Hc_dmat, &c_ptr);       
    assign_blasfeo_dmat_mem(N*nc+ncN, nu*N, workspace->Ccg_dmat, &c_ptr);        
//     assign_blasfeo_dmat_mem((N+1)*nbx, nu*N, workspace->Ccx_dmat, &c_ptr);
    assign_blasfeo_dmat_mem(N*nbx, nu*N, workspace->Ccx_dmat, &c_ptr);
                
    assign_blasfeo_dvec_mem((nx+nu)*(N+1), workspace->rq, &c_ptr);
    assign_blasfeo_dvec_mem(nx*N, workspace->a_dvec, &c_ptr);         
    assign_blasfeo_dvec_mem(nx*(N+1), workspace->L, &c_ptr);       
    assign_blasfeo_dvec_mem((nu+nx)*(N+1), workspace->gcw, &c_ptr);       
    assign_blasfeo_dvec_mem(N*nc+ncN, workspace->bg, &c_ptr);       
//     assign_blasfeo_dvec_mem((N+1)*nbx, workspace->bx, &c_ptr);
    assign_blasfeo_dvec_mem(N*nbx, workspace->bx, &c_ptr);

    blasfeo_dgese(nu+nx, nx*(N+1), 0.0, workspace->StQ, 0, 0);       
    blasfeo_dgese(nu+nx, nx*N, 0.0, workspace->BtAt, 0, 0);        
    blasfeo_dgese(nu*(N+1), (nu+nx)*N, 0.0, workspace->HtWt, 0, 0);        
    blasfeo_dgese(nu*N, nx*N, 0.0, workspace->Gt, 0, 0);        
    blasfeo_dgese(nu, nu*N, 0.0, workspace->Rt, 0, 0);         
    blasfeo_dgese(nc, nx*N, 0.0, workspace->Cgx_dmat, 0, 0);         
    blasfeo_dgese(ncN, nx, 0.0, workspace->CgN_dmat, 0, 0);         
    blasfeo_dgese(nbx, nx, 0.0, workspace->Cx_dmat, 0, 0);           
    blasfeo_dgese(N*nu, nu*N, 0.0, workspace->Hc_dmat, 0, 0);       
    blasfeo_dgese(N*nc+ncN, nu*N, 0.0, workspace->Ccg_dmat, 0, 0);        
//     blasfeo_dgese((N+1)*nbx, nu*N, 0.0, workspace->Ccx_dmat, 0, 0);
    blasfeo_dgese(N*nbx, nu*N, 0.0, workspace->Ccx_dmat, 0, 0);
                
    blasfeo_dvecse((nx+nu)*(N+1), 0.0, workspace->rq, 0);
    blasfeo_dvecse(nx*N, 0.0, workspace->a_dvec, 0);         
    blasfeo_dvecse(nx*(N+1), 0.0, workspace->L, 0);       
    blasfeo_dvecse((nu+nx)*(N+1), 0.0, workspace->gcw, 0);       
    blasfeo_dvecse(N*nc+ncN, 0.0, workspace->bg, 0);       
//     blasfeo_dvecse((N+1)*nbx, 0.0, workspace->bx, 0);
    blasfeo_dvecse(N*nbx, 0.0, workspace->bx, 0);
    
//     mexPrintf("\npointer moved VS. size calculated = %d VS. %d \n", c_ptr- (char*)raw_memory, Condensing_Blasfeo_workspace_calculate_size(nx, nu, N, nc, ncN, nbx));

    return workspace;
}

void exitFcn(){
    if (mem_alloc){
        mxFree(mem);
    }
}

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
        
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    int ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
            
    /*Outputs*/
    double  *Hc, *gc, *Ccg, *Ccx, *lcc, *ucc, *lxc, *uxc; 
       
    Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc")  );
    gc = mxGetPr( mxGetField(prhs[0], 0, "gc")  );
    Ccg = mxGetPr( mxGetField(prhs[0], 0, "Ccg")  );
    Ccx = mxGetPr( mxGetField(prhs[0], 0, "Ccx")  );
    lcc = mxGetPr( mxGetField(prhs[0], 0, "lcc")  );
    ucc = mxGetPr( mxGetField(prhs[0], 0, "ucc")  );
    lxc = mxGetPr( mxGetField(prhs[0], 0, "lxc")  );
    uxc = mxGetPr( mxGetField(prhs[0], 0, "uxc")  );
    
    int iter = mxGetScalar( mxGetField(prhs[0], 0, "iter") );
    int hot_start = mxGetScalar( mxGetField(prhs[0], 0, "hot_start") );
    
    /*Allocate memory*/
    int i=0,j=0;
        
    if (!mem_alloc){   
        int size = Condensing_Blasfeo_workspace_calculate_size(nx, nu, N, nc, ncN, nbx);
        mem = mxCalloc(size,1);
        mexMakeMemoryPersistent(mem);
        workspace = Condensing_Blasfeo_workspace_cast(nx, nu, N, nc, ncN, nbx, mem);
        mem_alloc = true;       
        mexAtExit(exitFcn);
    }
    
    struct blasfeo_dmat *StQ = workspace->StQ;
    struct blasfeo_dmat *BtAt = workspace->BtAt;
    struct blasfeo_dmat *HtWt = workspace->HtWt;
    struct blasfeo_dmat *Gt = workspace->Gt;
    struct blasfeo_dmat *Rt = workspace->Rt;
    struct blasfeo_dmat *Cgx_dmat = workspace->Cgx_dmat;
    struct blasfeo_dmat *CgN_dmat = workspace->CgN_dmat;    
    struct blasfeo_dmat *Cx_dmat = workspace->Cx_dmat;    
    struct blasfeo_dmat *Hc_dmat = workspace->Hc_dmat;   
    struct blasfeo_dmat *Ccg_dmat = workspace->Ccg_dmat;   
    struct blasfeo_dmat *Ccx_dmat = workspace->Ccx_dmat;  
    struct blasfeo_dvec *rq = workspace->rq; 
    struct blasfeo_dvec *a_dvec = workspace->a_dvec;   
    struct blasfeo_dvec *L = workspace->L;   
    struct blasfeo_dvec *gcw = workspace->gcw;   
    struct blasfeo_dvec *bg = workspace->bg;  
    struct blasfeo_dvec *bx = workspace->bx;
          
    for(i=0;i<N;i++){
        blasfeo_pack_tran_dmat(nx, nu, S+i*nx*nu, nx, StQ, 0, i*nx);
        blasfeo_pack_dmat(nx, nx, Q+i*nx*nx, nx, StQ, nu, i*nx);            
        blasfeo_pack_tran_dmat(nu, nu, R+i*nu*nu, nu, Rt, 0, i*nu);
        blasfeo_pack_tran_dmat(nx, nu, B+i*nx*nu, nx, BtAt, 0, i*nx);
        blasfeo_pack_tran_dmat(nx, nx, A+i*nx*nx, nx, BtAt, nu, i*nx);
        blasfeo_pack_dmat(nc, nx, Cgx+i*nc*nx, nc, Cgx_dmat, 0, i*nx);        
        blasfeo_pack_dvec(nu, gu+i*nu, rq, i*(nx+nu));
        blasfeo_pack_dvec(nx, gx+i*nx, rq, i*(nx+nu)+nu);       
    }
    blasfeo_pack_dmat(nx, nx, Q+N*nx*nx, nx, StQ, nu, N*nx);
    blasfeo_pack_dmat(ncN, nx, CgN, ncN, CgN_dmat, 0, 0);
    blasfeo_pack_dmat(nbx, nx, Cx, nbx, Cx_dmat, 0, 0);    
    blasfeo_pack_dvec(nx, gx+N*nx, rq, N*(nx+nu)+nu);
    blasfeo_pack_dvec(nx*N, a, a_dvec, 0);
               
    /*Start the loop*/
        
    /* compute G */
    for(i=0;i<N;i++){
        blasfeo_dgecp(nu, nx, BtAt, 0, i*nx, Gt, i*nu, i*nx);
        for (j=i+1;j<N;j++){
            blasfeo_dgemm_nn(nu, nx, nx, 1.0, Gt, (j-1)*nu, i*nx, BtAt, nu, j*nx, 0.0, Gt, j*nu, i*nx, Gt, j*nu, i*nx);
        }
    }
               
    /* Compute Hc */
    for(i=0;i<N;i++){
        blasfeo_dgemm_nt(nu, nx, nx, 1.0, Gt, (N-1)*nu, i*nx, StQ, nu, N*nx, 0.0, HtWt, N*nu, i*(nu+nx)+nu, HtWt, N*nu, i*(nu+nx)+nu);
        for(j=N-1;j>i;j--){                       
            blasfeo_dgemm_nt(nu,nu+nx,nx,1.0,Gt,(j-1)*nu,i*nx,StQ,0,j*nx,0.0,HtWt,j*nu,i*(nu+nx),HtWt,j*nu,i*(nu+nx));
            blasfeo_dgemm_nt(nu,nu+nx,nx,1.0,HtWt,(j+1)*nu,i*(nx+nu)+nu,BtAt,0,j*nx,1.0,HtWt,j*nu,i*(nu+nx),HtWt,j*nu,i*(nu+nx));
            
            blasfeo_dgecp(nu, nu, HtWt, j*nu, i*(nx+nu), Hc_dmat, i*nu, j*nu);
        }          
        blasfeo_dgemm_nt(nu,nu,nx,1.0,HtWt,(i+1)*nu,i*(nx+nu)+nu,BtAt,0,i*nx,1.0,Rt,0,i*nu,Hc_dmat,i*nu,i*nu);
    }
    blasfeo_dtrtr_u(N*nu, Hc_dmat, 0, 0, Hc_dmat, 0, 0);
    blasfeo_unpack_dmat(N*nu,N*nu,Hc_dmat,0,0,Hc,N*nu);
            
    /* Compute Cc */
    if (nc>0){         
        for(i=0;i<N;i++){
            blasfeo_pack_dmat(nc, nu, Cgu+i*nc*nu, nc, Ccg_dmat, i*nc, i*nu);
            for(j=i+1;j<N;j++)   
                blasfeo_dgemm_nt(nc,nu,nx,1.0,Cgx_dmat,0,j*nx,Gt,(j-1)*nu,i*nx,0.0,Ccg_dmat,j*nc,i*nu,Ccg_dmat,j*nc,i*nu);                        
        }            
    }
    
    /* Compute CcN */
    if (ncN>0){          
        for(i=0;i<N;i++)                
            blasfeo_dgemm_nt(ncN,nu,nx,1.0,CgN_dmat,0,0,Gt,(N-1)*nu,i*nx,0.0,Ccg_dmat,N*nc,i*nu,Ccg_dmat,N*nc,i*nu);                                        
    }
        
    blasfeo_unpack_dmat(N*nc+ncN,N*nu,Ccg_dmat,0,0,Ccg,N*nc+ncN);
    
    /* Compute Ccx */
    if (nbx>0){         
        for(i=0;i<N;i++){
            for(j=i+1;j<=N;j++)             
//                    blasfeo_dgemm_nt(nbx,nu,nx,1.0,Cx_dmat,0,0,Gt,(j-1)*nu,i*nx,0.0,Ccx_dmat,j*nbx,i*nu,Ccx_dmat,j*nbx,i*nu); 
                blasfeo_dgemm_nt(nbx,nu,nx,1.0,Cx_dmat,0,0,Gt,(j-1)*nu,i*nx,0.0,Ccx_dmat,(j-1)*nbx,i*nu,Ccx_dmat,(j-1)*nbx,i*nu);
        }
//            blasfeo_unpack_dmat((N+1)*nbx,N*nu,Ccx_dmat,0,0,Ccx,(N+1)*nbx);
        blasfeo_unpack_dmat(N*nbx,N*nu,Ccx_dmat,0,0,Ccx,N*nbx);
    }  
             
    /* compute L */
    blasfeo_pack_dvec(nx, ds0, L, 0);
    for(i=0;i<N;i++){
        blasfeo_dgemv_t(nx, nx, 1.0, BtAt, nu, i*nx, L, i*nx, 1.0, a_dvec, i*nx, L, (i+1)*nx);
    }
     
    /* compute gc */
    blasfeo_dgemv_n(nx, nx, 1.0, StQ, nu, N*nx, L, N*nx, 1.0, rq, N*(nx+nu)+nu, gcw, N*(nu+nx)+nu);
    for(i=N-1;i>0;i--){
        blasfeo_dgemv_n(nu+nx, nx, 1.0, StQ, 0, i*nx, L, i*nx, 1.0, rq, i*(nx+nu), gcw, i*(nu+nx));
        blasfeo_dgemv_n(nu+nx, nx, 1.0, BtAt, 0, i*nx, gcw, (i+1)*(nx+nu)+nu, 1.0, gcw, i*(nx+nu), gcw, i*(nu+nx));
        blasfeo_unpack_dvec(nu, gcw, i*(nx+nu), gc+i*nu);
    }   
    blasfeo_dgemv_n(nu, nx, 1.0, StQ, 0, 0, L, 0, 1.0, rq, 0, gcw, 0);
    blasfeo_dgemv_n(nu, nx, 1.0, BtAt, 0, 0, gcw, (nx+nu)+nu, 1.0, gcw, 0, gcw, 0);
    blasfeo_unpack_dvec(nu, gcw, 0, gc);
   
    /* Compute cc */
    if (nc>0){                    
        for(i=0;i<N;i++){           
            blasfeo_dgemv_n(nc, nx, -1.0, Cgx_dmat, 0, i*nx, L, i*nx, 0.0, bg, i*nc, bg, i*nc);
            blasfeo_unpack_dvec(nc, bg, i*nc, lcc+i*nc);
            for(j=0;j<nc;j++){
                ucc[i*nc+j] = lcc[i*nc+j] + uc[i*nc+j];
                lcc[i*nc+j] += lc[i*nc+j];          
            }
        }        
    }   
    
    /* Compute ccN */
    if (ncN>0){    
        blasfeo_dgemv_n(ncN, nx, -1.0, CgN_dmat, 0, 0, L, N*nx, 0.0, bg, N*nc, bg, N*nc);
        blasfeo_unpack_dvec(ncN, bg, N*nc, lcc+N*nc);
        for(j=0;j<ncN;j++){
            ucc[N*nc+j] = lcc[N*nc+j] + uc[N*nc+j];
            lcc[N*nc+j] += lc[N*nc+j];          
        }
    }
    
    /* Compute ccx */
    if (nbx>0){                    
//         for(i=0;i<=N;i++){
        for(i=0;i<N;i++){
//             blasfeo_dgemv_n(nbx, nx, -1.0, Cx_dmat, 0, 0, L, i*nx, 0.0, bx, i*nbx, bx, i*nbx);
             blasfeo_dgemv_n(nbx, nx, -1.0, Cx_dmat, 0, 0, L, (i+1)*nx, 0.0, bx, i*nbx, bx, i*nbx);
            blasfeo_unpack_dvec(nbx, bx, i*nbx, lxc+i*nbx);
            for(j=0;j<nbx;j++){
                uxc[i*nbx+j] = lxc[i*nbx+j]+ ub_dx[i*nbx+j];
                lxc[i*nbx+j] += lb_dx[i*nbx+j];          
            }
        }        
    }   
    
}