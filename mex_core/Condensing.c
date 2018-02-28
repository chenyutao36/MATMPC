
#include "mex.h"
#include <stdbool.h>
#include "string.h"

#include "mpc_common.h"

#include "blas.h"

static double *L = NULL;
static double *Hi = NULL, *Cci = NULL, *CcN = NULL;
static bool mem_alloc = false;

void exitFcn(){
    if (mem_alloc){
        mxFree(L);
        mxFree(Hi);
        mxFree(Cci);   
        mxFree(CcN);
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
    double *Cu = mxGetPr( mxGetField(prhs[0], 0, "Cu") );
    double *CN = mxGetPr( mxGetField(prhs[0], 0, "CN") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[0], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[0], 0, "uc") );
        
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
            
    /*Outputs*/
    double  *Hc, *gc, *Cc, *lcc, *ucc, *G, *W_mat, *w_vec; 
    
    G = mxGetPr( mxGetField(prhs[0], 0, "G")  );
    W_mat = mxGetPr( mxGetField(prhs[0], 0, "W_mat")  );
    w_vec = mxGetPr( mxGetField(prhs[0], 0, "w_vec")  );
    
    Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc")  );
    gc = mxGetPr( mxGetField(prhs[0], 0, "gc")  );
    Cc = mxGetPr( mxGetField(prhs[0], 0, "Cc")  );
    lcc = mxGetPr( mxGetField(prhs[0], 0, "lcc")  );
    ucc = mxGetPr( mxGetField(prhs[0], 0, "ucc")  );
    
    int iter = mxGetScalar( mxGetField(prhs[0], 0, "iter") );
    int hot_start = mxGetScalar( mxGetField(prhs[0], 0, "hot_start") );
    int lin_obj = mxGetScalar( mxGetField(prhs[0], 0, "lin_obj") );  
    bool cond_save = (hot_start==1) && (lin_obj==1) ;
    
    /*Allocate memory*/
    int i=0,j=0;
     
    if (!mem_alloc){              
        L = (double *)mxMalloc((N+1)*nx * sizeof(double));
        mexMakeMemoryPersistent(L);        
        Hi = (double *)mxMalloc(nu*nu * sizeof(double));
        mexMakeMemoryPersistent(Hi);       
        Cci = (double *)mxMalloc(nc*nu * sizeof(double));
        mexMakeMemoryPersistent(Cci);       
        CcN = (double *)mxMalloc(ncN*nu * sizeof(double)); 
        mexMakeMemoryPersistent(CcN);
        
        mem_alloc = true;       
        mexAtExit(exitFcn);
    }
   
    char *nTrans = "N", *Trans="T", *SIDE = "L", *UPLO = "L";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
              
    /*Start the loop*/
        
    if (iter==1 || !cond_save){ // check if adjoint rti is used    
        /* compute G */
        for(i=0;i<N;i++){
            memcpy(G+(i*N+i)*nx*nu, B+i*nx*nu, nx*nu*sizeof(double));
            for (j=i+1;j<N;j++){
                dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, G+(i*N+j)*nx*nu, &nx);
            }
        }
        
        /* Compute Hc */
        for(i=0;i<N;i++){
            dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+N*nx*nx, &nx, G+(i*N+N-1)*nx*nu, &nx, &zero, W_mat+(i*N+N-1)*nx*nu, &nx);
            for(j=N-1;j>i;j--){        
                dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, S+j*nx*nu, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, Hi, &nu);                     
                dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+j*nx*nu, &nx, W_mat+(i*N+j)*nx*nu, &nx, &one_d, Hi, &nu);
                Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);
                Block_Fill_Trans(nu, nu, Hi, Hc, i*nu, j*nu, N*nu);

                dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, W_mat+(i*N+j)*nx*nu, &nx, &zero, W_mat+(i*N+j-1)*nx*nu, &nx); 
                dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &one_d, W_mat+(i*N+j-1)*nx*nu, &nx);
            }
            memcpy(Hi,R+i*nu*nu,nu*nu*sizeof(double));
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+i*nx*nu, &nx, W_mat+(i*N+i)*nx*nu, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc, i*nu, i*nu, N*nu);      
        }
        
        /* Compute Cc */
        if (nc>0){         
            for(i=0;i<N;i++){
                Block_Fill(nc, nu, Cu+i*nc*nu, Cc, i*nc, i*nu, N*nc+ncN);
                for(j=i+1;j<N;j++){   
                    dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, Cx+j*nc*nx, &nc, G+(i*N+j-1)*nx*nu, &nx, &zero, Cci, &nc);
                    Block_Fill(nc, nu, Cci, Cc, j*nc, i*nu, N*nc+ncN);
                }    
            }  
        }
        
        /* Compute CcN */
        if (ncN>0){          
            for(i=0;i<N;i++){                 
                dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, CN, &ncN, G+(i*N+N-1)*nx*nu, &nx, &zero, CcN, &ncN);
                Block_Fill(ncN, nu, CcN, Cc, N*nc, i*nu, N*nc+ncN);
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
    memcpy(w_vec+(N-1)*nx,gx+N*nx,nx*sizeof(double));
    dsymv(UPLO, &nx, &one_d, Q+N*nx*nx, &nx, L+N*nx, &one_i, &one_d, w_vec+(N-1)*nx, &one_i);
    for(i=N-1;i>0;i--){
        memcpy(gc+i*nu,gu+i*nu,nu*sizeof(double));
        dgemv(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,L+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
        dgemv(Trans,&nx,&nu,&one_d,B+i*nx*nu,&nx,w_vec+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
         
        memcpy(w_vec+(i-1)*nx, gx+i*nx, nx*sizeof(double));
        dsymv(UPLO, &nx, &one_d, Q+i*nx*nx, &nx, L+i*nx, &one_i, &one_d, w_vec+(i-1)*nx, &one_i);
        dgemv(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,w_vec+i*nx,&one_i,&one_d,w_vec+(i-1)*nx,&one_i);
    }   
    memcpy(gc,gu,nu*sizeof(double));
    dgemv(Trans,&nx,&nu,&one_d,S,&nx,L,&one_i,&one_d,gc,&one_i);
    dgemv(Trans,&nx,&nu,&one_d,B,&nx,w_vec,&one_i,&one_d,gc,&one_i);
   
    /* Compute cc */
    if (nc>0){                    
        for(i=0;i<N;i++){
            memcpy(lcc+i*nc,lc+i*nc,nc*sizeof(double));
            dgemv(nTrans,&nc,&nx,&minus_one,Cx+i*nc*nx,&nc,L+i*nx,&one_i,&one_d,lcc+i*nc,&one_i); 

            memcpy(ucc+i*nc,uc+i*nc,nc*sizeof(double));
            dgemv(nTrans,&nc,&nx,&minus_one,Cx+i*nc*nx,&nc,L+i*nx,&one_i,&one_d,ucc+i*nc,&one_i);
        }        
    }   
    
    /* Compute ccN */
    if (ncN>0){          
        memcpy(lcc+N*nc,lc+N*nc,ncN*sizeof(double));
        dgemv(nTrans,&ncN,&nx,&minus_one,CN,&ncN,L+N*nx,&one_i,&one_d,lcc+N*nc,&one_i);
        memcpy(ucc+N*nc,uc+N*nc,ncN*sizeof(double));
        dgemv(nTrans,&ncN,&nx,&minus_one,CN,&ncN,L+N*nx,&one_i,&one_d,ucc+N*nc,&one_i);
    }
    
}