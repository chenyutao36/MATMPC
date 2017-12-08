
#include "mex.h"
#include <stdbool.h>
#include "string.h"

#include "mpc_common.h"

// for builtin blas
#include "blas.h"

// for openblas
// #include "f77blas.h"
// #if !defined(_WIN32)
// #define dgemm dgemm_
// #define dgemv dgemv_
// #endif

static double *L = NULL, *w_vec = NULL, *W_mat = NULL, *w_vec_dup = NULL, *W_mat_dup = NULL;
static double *Hi = NULL, *Cci = NULL, *CcN = NULL;
static bool mem_alloc_cond = false;

void exitFcn(){
    if (mem_alloc_cond){
        mxFree(L);
        mxFree(w_vec);
        mxFree(W_mat);
        mxFree(w_vec_dup);
        mxFree(W_mat_dup);
        mxFree(Hi);
        mxFree(Cci);   
        mxFree(CcN);
    }
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{    
    /*Inputs*/
    double *A = mxGetPr(prhs[0]);
    double *B = mxGetPr(prhs[1]);    
    double *Qh = mxGetPr(prhs[2]);
    double *S = mxGetPr(prhs[3]);    
    double *R = mxGetPr(prhs[4]);
    double *Cx = mxGetPr(prhs[5]);
    double *Cu = mxGetPr(prhs[6]);
    double *CxN = mxGetPr(prhs[7]);
    double *ds0 = mxGetPr(prhs[8]);
    double *a = mxGetPr(prhs[9]);    
    double *gs = mxGetPr(prhs[10]);
    double *gu = mxGetPr(prhs[11]);    
    double *lc = mxGetPr(prhs[12]);
    double *uc = mxGetPr(prhs[13]);
    
    mwSize nx = mxGetScalar( mxGetField(prhs[14], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[14], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[14], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[14], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[14], 0, "N") );
        
    mwSize nz = nx+nu;
    
    /*Outputs*/
    double  *Hc, *gc, *Cc, *lcc, *ucc, *lcu, *ucu; 
    
    plhs[0] = mxCreateDoubleMatrix(N*nu, N*nu, mxREAL);
    Hc = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);
    gc = mxGetPr(plhs[1]);   
    plhs[2] = mxCreateDoubleMatrix(N*nc+ncN, N*nu, mxREAL);
    Cc = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL);
    lcc = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL);
    
    /*Allocate memory*/
    mwIndex i=0,j=0;
    double *cell; /* from cells */
     
    if (!mem_alloc_cond){       
        
        L = (double *)mxCalloc((N+1)*nx, sizeof(double));
        mexMakeMemoryPersistent(L);
        
        w_vec = (double *)mxCalloc(nx, sizeof(double));
        mexMakeMemoryPersistent(w_vec);
        
        W_mat = (double *)mxCalloc(nx*nu, sizeof(double));
        mexMakeMemoryPersistent(W_mat);
        
        w_vec_dup = (double *)mxCalloc(nx, sizeof(double));
        mexMakeMemoryPersistent(w_vec_dup);
        
        W_mat_dup = (double *)mxCalloc(nx*nu, sizeof(double));
        mexMakeMemoryPersistent(W_mat_dup);
        
        Hi = (double *)mxCalloc(nu*nu, sizeof(double));
        mexMakeMemoryPersistent(Hi);
        
        Cci = (double *)mxCalloc(nc*nu, sizeof(double));
        mexMakeMemoryPersistent(Cci);
        
        CcN = (double *)mxCalloc(ncN*nu,sizeof(double)); 
        mexMakeMemoryPersistent(CcN);

        mem_alloc_cond = true;       
        mexAtExit(exitFcn);
    }
   
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    mwSignedIndex one_i = 1; /* never use int for lapack and blas routines */
              
    /*Start the loop*/
    
    /* compute G */
    double *G[N*N];
    for(i=0;i<N;i++){
        G[i*N+i] = B+i*nx*nu; // Bi
        for (j=i+1;j<N;j++){
            cell = A+j*nx*nx;
            G[i*N+j] = (double *)mxCalloc(nx*nu, sizeof(double));
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, G[i*N+j-1], &nx, &zero, G[i*N+j], &nx);
        }
    }
         
    /* compute L */
    memcpy(&L[0],&ds0[0], nx*sizeof(double)); 
    for(i=0;i<N;i++){
        cell = A+i*nx*nx;
        memcpy(&L[(i+1)*nx],&a[i*nx], nx*sizeof(double)); 
        dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,L+i*nx,&one_i,&one_d,L+(i+1)*nx,&one_i);
    }
    
    /* compute gc */
    cell = Qh + N*nx*nx;
    memcpy(&w_vec[0],&gs[N*nx],nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,L+N*nx,&one_i,&one_d,w_vec,&one_i);
    for(i=N-1;i>0;i--){
        cell = S + i*nx*nu;
        memcpy(&gc[i*nu],&gu[i*nu],nu*sizeof(double));
        dgemv(Trans,&nx,&nu,&one_d,cell,&nx,L+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
        cell = B+i*nx*nu;
        dgemv(Trans,&nx,&nu,&one_d,cell,&nx,w_vec,&one_i,&one_d,gc+i*nu,&one_i);
         
        cell = Qh + i*nx*nx;
        memcpy(&w_vec_dup[0],&gs[i*nx],nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,L+i*nx,&one_i,&one_d,w_vec_dup,&one_i);
        cell = A + i*nx*nx;
        dgemv(Trans,&nx,&nx,&one_d,cell,&nx,w_vec,&one_i,&one_d,w_vec_dup,&one_i);
        memcpy(&w_vec[0],&w_vec_dup[0],nx*sizeof(double));
    }   
    cell = S;
    memcpy(&gc[0],&gu[0],nu*sizeof(double));
    dgemv(Trans,&nx,&nu,&one_d,cell,&nx,L,&one_i,&one_d,gc,&one_i);
    cell = B;
    dgemv(Trans,&nx,&nu,&one_d,cell,&nx,w_vec,&one_i,&one_d,gc,&one_i);
     
    /* Compute Hc (only the lower triangular part) */
    for(i=0;i<N;i++){
        cell = Qh + N*nx*nx;
        dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, G[i*N+N-1], &nx, &zero, W_mat, &nx);
        for(j=N-1;j>i;j--){
         
            cell = S + j*nx*nu;
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, G[i*N+j-1], &nx, &zero, Hi, &nu);         
            cell = B + j*nx*nu;
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);
             
            cell = A + j*nx*nx;
            dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &zero, W_mat_dup, &nx); 
            cell = Qh + j*nx*nx;
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, G[i*N+j-1], &nx, &one_d, W_mat_dup, &nx); 
            memcpy(&W_mat[0],&W_mat_dup[0],nx*nu*sizeof(double));
        }
        cell = R + i*nu*nu;
        memcpy(Hi,cell,nu*nu*sizeof(double));
        cell = B + i*nx*nu;
        dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &one_d, Hi, &nu);
        Block_Fill(nu, nu, Hi, Hc, i*nu, i*nu, N*nu);
    }    
    
    /* fill the upper triangular part of Hc (Hc is symmetric) */
    for(i=0;i<N*nu;i++){
        for(j=i+1;j<N*nu;j++)
            Hc[j*N*nu+i]=Hc[i*N*nu+j];
    }
    
    /* Compute Cc */
    if (nc>0){   
        for(i=0;i<N;i++){
            cell = Cu + i*nc*nu;
            Block_Fill(nc, nu, cell, Cc, i*nc, i*nu, N*nc+ncN);

            for(j=i+1;j<N;j++){   
                cell = Cx + j*nc*nx;
                dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, cell, &nc, G[i*N+j-1], &nx, &zero, Cci, &nc);
                Block_Fill(nc, nu, Cci, Cc, j*nc, i*nu, N*nc+ncN);
            }    
        }
    
        /* Compute cc */
        for(i=0;i<N;i++){
            cell = Cx + i*nc*nx;
            memcpy(&lcc[i*nc],&lc[i*nc],nc*sizeof(double));
            dgemv(nTrans,&nc,&nx,&minus_one,cell,&nc,L+i*nx,&one_i,&one_d,lcc+i*nc,&one_i); 

            memcpy(&ucc[i*nc],&uc[i*nc],nc*sizeof(double));
            dgemv(nTrans,&nc,&nx,&minus_one,cell,&nc,L+i*nx,&one_i,&one_d,ucc+i*nc,&one_i);
        }
    }
    
    /* Compute CcN and ccN */
    if (ncN>0){   
        cell = CxN;
        for(i=0;i<N;i++){                 
            dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, cell, &ncN, G[i*N+N-1], &nx, &zero, CcN, &ncN);
            Block_Fill(ncN, nu, CcN, Cc, N*nc, i*nu, N*nc+ncN);
        }

        memcpy(&lcc[N*nc],&lc[N*nc],ncN*sizeof(double));
        dgemv(nTrans,&ncN,&nx,&minus_one,cell,&ncN,L+N*nx,&one_i,&one_d,lcc+N*nc,&one_i);
        memcpy(&ucc[N*nc],&uc[N*nc],ncN*sizeof(double));
        dgemv(nTrans,&ncN,&nx,&minus_one,cell,&ncN,L+N*nx,&one_i,&one_d,ucc+N*nc,&one_i);
    }
    
        
    /* Free memory */
    for(i=0;i<N;i++){
        for (j=i+1;j<N;j++)
            mxFree(G[i*N+j]);
    }

}