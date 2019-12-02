
#include "mex.h"
#include <stdbool.h>
#include <stdio.h>
#include "string.h"
#include "mpc_common.h"

#include "blas.h"

static double *L = NULL,*Hc_temp=NULL,*G_temp=NULL,*G=NULL, *W_mat=NULL, *w_vec=NULL;
static double *Hi = NULL, *Cci = NULL, *Ccxi = NULL, *CcN = NULL;
static bool mem_alloc = false;

void exitFcn(){
    if (mem_alloc){
        mxFree(G);
        mxFree(G_temp);
        mxFree(Hc_temp);
        mxFree(W_mat);
        mxFree(w_vec);
        mxFree(L);
        mxFree(Hi);
        mxFree(Cci);   
        mxFree(CcN);
        mxFree(Ccxi);   
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
    double *index_T = mxGetPr( mxGetField(prhs[0], 0, "index_T") );
    size_t r = mxGetScalar( mxGetField(prhs[0], 0, "r") );
        
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
            
    /*Outputs*/
    double  *Hc, *gc, *Ccg, *Ccx, *lcc, *ucc, *lxc, *uxc; 
       
    Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc_r")  );
    gc = mxGetPr( mxGetField(prhs[0], 0, "gc_r")  );
    Ccg = mxGetPr( mxGetField(prhs[0], 0, "Ccg_r")  );
    Ccx = mxGetPr( mxGetField(prhs[0], 0, "Ccx_r")  );
    lcc = mxGetPr( mxGetField(prhs[0], 0, "lcc")  );
    ucc = mxGetPr( mxGetField(prhs[0], 0, "ucc")  );
    lxc = mxGetPr( mxGetField(prhs[0], 0, "lxc")  );
    uxc = mxGetPr( mxGetField(prhs[0], 0, "uxc")  );
    
    int iter = mxGetScalar( mxGetField(prhs[0], 0, "iter") );
    int hot_start = mxGetScalar( mxGetField(prhs[0], 0, "hot_start") );
 
    bool flag=false;
    
    /*Allocate memory*/
    int i=0,j=0,k=0,l=0;
     
    if (!mem_alloc){
        G = (double *)mxMalloc(nx*N*r*nu * sizeof(double));
        mexMakeMemoryPersistent(G);
        G_temp = (double *)mxMalloc(nx*N*nu * sizeof(double));
        mexMakeMemoryPersistent(G_temp);
        Hc_temp = (double *)mxMalloc(nu*nu*r*N* sizeof(double));
        mexMakeMemoryPersistent(Hc_temp);
        W_mat = (double *)mxMalloc(nx*N*nu * sizeof(double));
        mexMakeMemoryPersistent(W_mat);  
        w_vec = (double *)mxMalloc(nx*N * sizeof(double));
        mexMakeMemoryPersistent(w_vec);          
        L = (double *)mxMalloc((N+1)*nx * sizeof(double));
        mexMakeMemoryPersistent(L);        
        Hi = (double *)mxMalloc(nu*nu * sizeof(double));
        mexMakeMemoryPersistent(Hi);       
        Cci = (double *)mxMalloc(nc*nu * sizeof(double));
        mexMakeMemoryPersistent(Cci);       
        CcN = (double *)mxMalloc(ncN*nu * sizeof(double)); 
        mexMakeMemoryPersistent(CcN);
        Ccxi = (double *)mxMalloc(nbx*nu * sizeof(double));
        mexMakeMemoryPersistent(Ccxi);   
        
        mem_alloc = true;       
        mexAtExit(exitFcn);
    }
   
    char *nTrans = "N", *Trans="T", *SIDE = "L", *UPLO = "L";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
              
    /*Start the loop*/
        
    /* compute G */
    for(i=0;i<r;i++){
        k=(int)index_T[i];
        memcpy(G+(i*N+k)*nx*nu, B+k*nx*nu, nx*nu*sizeof(double));
        for(j=k+1;j<N;j++){
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, G+(i*N+j)*nx*nu, &nx);
            if (j<(int) index_T[i+1]){
                for(l=0;l<nx*nu;l++){
                    *(G+(i*N+j)*nx*nu+l) += *(B+j*nx*nu+l);
                }
            }
        }
    }
    
     /* Compute Hc */
    for(i=0;i<r;i++){
        dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+N*nx*nx, &nx, G+(i*N+N-1)*nx*nu, &nx, &zero, W_mat+(N-1)*nx*nu, &nx);
        for(j=N-1;j>(int)index_T[i];j--){        
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, S+j*nx*nu, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, Hi, &nu);                     
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+j*nx*nu, &nx, W_mat+(j)*nx*nu, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc_temp, j*nu, i*nu, N*nu);
       
            dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, W_mat+(j)*nx*nu, &nx, &zero, W_mat+(j-1)*nx*nu, &nx); 
            dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &one_d, W_mat+(j-1)*nx*nu, &nx);
        }
        memcpy(Hi,R+i*nu*nu,nu*nu*sizeof(double));
        dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+((int)index_T[i])*nx*nu, &nx, W_mat+((int)index_T[i])*nx*nu, &nx, &one_d, Hi, &nu);
        Block_Fill(nu, nu, Hi, Hc_temp, ((int)index_T[i])*nu, i*nu, N*nu);      
    }
    
    k=-1;
    
    for(j=0;j<N;j++){  
        if (j==(int)index_T[k+1]){
            k++;
            for (i=0; i< nu*(k+1);i++)
                for (l=0; l< nu;l++)    
                    *(Hc+i*r*nu+k*nu+l)=  *(Hc_temp+i*N*nu+j*nu+l); 
        }  
        else{
            for (i=0; i< nu*(k+1);i++)
                for (l=0; l< nu;l++)    
                    *(Hc+i*r*nu+k*nu+l)+=  *(Hc_temp+i*N*nu+j*nu+l); 
        }
    }
    
    for(i=0;i<r;i++)
        for(j=i+1;j<r;j++)
             for(l=0;l<nu;l++)
                 for(k=0;k<nu;k++)
                    *(Hc+r*nu*(j*nu+k)+i*nu+l) = *(Hc+r*nu*(i*nu+l)+j*nu+k);
                 
            
    /* Compute Cc */
    if (nc>0){  
        for(i=0;i<r;i++){   
            Block_Fill(nc, nu, Cgu+((int)(index_T[i]))*nc*nu, Ccg, ((int)(index_T[i]))*nc, i*nu, N*nc+ncN);              
            for(j=((int)(index_T[i]))+1;j<N;j++){
                if (j<(int)(index_T[i+1])){
                    memcpy(Cci,Cgu+(j)*nc*nu,nc*nu*sizeof(double));
                    dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, Cgx+j*nc*nx, &nc, G+(i*N+j-1)*nx*nu, &nx, &one_d, Cci, &nc);
                }    
                else{
                    dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, Cgx+j*nc*nx, &nc, G+(i*N+j-1)*nx*nu, &nx, &zero, Cci, &nc);
                }   
                Block_Fill(nc, nu, Cci, Ccg, j*nc, i*nu, N*nc+ncN); 
            }                
                          
        }
    }    
    if (ncN>0){
        for(i=0;i<r;i++){
            dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, CgN, &ncN, G+(i*N+N-1)*nx*nu, &nx, &zero, CcN, &ncN);
            Block_Fill(ncN, nu, CcN, Ccg, N*nc, i*nu, N*nc+ncN);
        }
    }
    
    
    
    /* Compute Ccx */
    if (nbx>0){         
        for(i=0;i<r;i++){
            for(j=index_T[i];j<N;j++){   
                dgemm(nTrans, nTrans, &nbx, &nu, &nx, &one_d, Cx, &nbx, G+(i*N+j)*nx*nu, &nx, &zero, Ccxi, &nbx);
                Block_Fill(nbx, nu, Ccxi, Ccx, j*nbx, i*nu, N*nbx);
            }  
        }  
    }
    
    
      /* compute L */
    memcpy(L,ds0, nx*sizeof(double)); 
    for(i=0;i<N;i++){
        memcpy(L+(i+1)*nx, a+i*nx, nx*sizeof(double)); 
        dgemv(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,L+i*nx,&one_i,&one_d,L+(i+1)*nx,&one_i);
    }
    
    /* compute gc */
    flag=false;    
    k=r-1;    
    memcpy(w_vec+(N-1)*nx,gx+N*nx,nx*sizeof(double));
    dsymv(UPLO, &nx, &one_d, Q+N*nx*nx, &nx, L+N*nx, &one_i, &one_d, w_vec+(N-1)*nx, &one_i);
    for(i=N-1;i>0;i--){
        if (k==0){
            flag=true;
        }
        if ( (int)(index_T[k+1])==i+1){
            memcpy(gc+k*nu,gu+k*nu,nu*sizeof(double));
            dgemv(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,L+i*nx,&one_i,&one_d,gc+k*nu,&one_i);
            dgemv(Trans,&nx,&nu,&one_d,B+i*nx*nu,&nx,w_vec+i*nx,&one_i,&one_d,gc+k*nu,&one_i);
        }   
        else{
            dgemv(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,L+i*nx,&one_i,&one_d,gc+k*nu,&one_i);
            dgemv(Trans,&nx,&nu,&one_d,B+i*nx*nu,&nx,w_vec+i*nx,&one_i,&one_d,gc+k*nu,&one_i);
        } 
        if ( (int)(index_T[k])==i){
            k=k-1;
        }    
        memcpy(w_vec+(i-1)*nx, gx+i*nx, nx*sizeof(double));
        dsymv(UPLO, &nx, &one_d, Q+i*nx*nx, &nx, L+i*nx, &one_i, &one_d, w_vec+(i-1)*nx, &one_i);
        dgemv(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,w_vec+i*nx,&one_i,&one_d,w_vec+(i-1)*nx,&one_i);
    }   
    if (!flag){
        memcpy(gc,gu,nu*sizeof(double));
    }
    dgemv(Trans,&nx,&nu,&one_d,S,&nx,L,&one_i,&one_d,gc,&one_i);
    dgemv(Trans,&nx,&nu,&one_d,B,&nx,w_vec,&one_i,&one_d,gc,&one_i);
   
 
    /* Compute cc */
    if (nc>0){                    
        for(i=0;i<N;i++){
            dgemv(nTrans,&nc,&nx,&minus_one,Cgx+i*nc*nx,&nc,L+i*nx,&one_i,&zero,lcc+i*nc,&one_i); 
            for(j=0;j<nc;j++){
                ucc[i*nc+j] = lcc[i*nc+j]+ uc[i*nc+j];
                lcc[i*nc+j] += lc[i*nc+j];          
            }
        }        
    }   
    
    /* Compute ccN */
    if (ncN>0){    
        dgemv(nTrans,&ncN,&nx,&minus_one,CgN,&ncN,L+N*nx,&one_i,&zero,lcc+N*nc,&one_i);
        for(j=0;j<ncN;j++){
            ucc[N*nc+j] = lcc[N*nc+j]+ uc[N*nc+j];
            lcc[N*nc+j] += lc[N*nc+j];          
        }
    }
    
    /* Compute ccx */
    if (nbx>0){                    
        for(i=0;i<N;i++){
            dgemv(nTrans,&nbx,&nx,&minus_one,Cx,&nbx,L+(i+1)*nx,&one_i,&zero,lxc+i*nbx,&one_i);
            for(j=0;j<nbx;j++){
                uxc[i*nbx+j] = lxc[i*nbx+j]+ ub_dx[i*nbx+j];
                lxc[i*nbx+j] += lb_dx[i*nbx+j];          
            }
        }        
    }   

    

}