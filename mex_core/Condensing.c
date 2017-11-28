
#include "mex.h"
#include <stdbool.h>
#include "string.h"

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

void Block_Fill(mwSize m, mwSize n, double *Gi, double *G, mwSize idm, mwSize idn, mwSize ldG){
       
    mwSize i,j;
    mwSize s;
    for (j=0;j<n;j++){
        s = idn*ldG + idm + j*ldG;
        for (i=0;i<m;i++){
            G[s+i] = Gi[j*m+i];
        }
    }
       
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{    
    /*Inputs*/
    if (!mxIsCell(prhs[0])) mexErrMsgTxt("Ak is not a cell array!"); /* A{} prhs[0] */
    if (!mxIsCell(prhs[1])) mexErrMsgTxt("Bk is not a cell array!"); /* B{} prhs[1] */ 
    if (!mxIsCell(prhs[2])) mexErrMsgTxt("Qk is not a cell array!"); /* Q{} prhs[2] */
    if (!mxIsCell(prhs[3])) mexErrMsgTxt("Sk is not a cell array!"); /* S{} prhs[3] */
    if (!mxIsCell(prhs[4])) mexErrMsgTxt("Rk is not a cell array!"); /* R{} prhs[4] */
    if (!mxIsCell(prhs[5])) mexErrMsgTxt("Cxk is not a cell array!"); /* Cx{} prhs[5] */
    if (!mxIsCell(prhs[6])) mexErrMsgTxt("Cuk is not a cell array!"); /* Cu{} prhs[6] */
    double *ds0 = mxGetPr(prhs[7]);
    double *a = mxGetPr(prhs[8]);    
    double *c = mxGetPr(prhs[9]);
    double *gs = mxGetPr(prhs[10]);
    double *gu = mxGetPr(prhs[11]);
    
    mwSize nx = mxGetScalar( mxGetField(prhs[12], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[12], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[12], 0, "nc") ); nc *= 2;
    mwSize ncN = mxGetScalar( mxGetField(prhs[12], 0, "ncN") ); ncN *=2;
    mwSize N = mxGetScalar( mxGetField(prhs[12], 0, "N") );
        
    mwSize nz = nx+nu;
    
    /*Outputs*/
    double  *Hc, *gc, *Cc_all, *cc; 
    
    plhs[0] = mxCreateDoubleMatrix(N*nu, N*nu, mxREAL);
    Hc = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);
    gc = mxGetPr(plhs[1]);   
    plhs[2] = mxCreateDoubleMatrix(N*nc+ncN, N*nu, mxREAL);
    Cc_all = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(N*nc+ncN, 1, mxREAL);
    cc = mxGetPr(plhs[3]);
    
    /*Allocate memory*/
    mwIndex i=0,j=0;
    const mxArray *cell_element;  
    double *cell; /* from cells */
    
    double **G = (double **) mxMalloc(N*N * sizeof(double*));
    
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
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1; /* never use int for lapack and blas routines */
              
    /*Start the loop*/
    
    /* compute G */
    for(i=0;i<N;i++){
        G[i*N+i] = mxGetPr(mxGetCell(prhs[1],i)); // Bi
        for (j=i+1;j<N;j++){
            cell_element = mxGetCell(prhs[0], j);           
            cell = mxGetPr(cell_element); /* Ai */
            G[i*N+j] = (double *)mxCalloc(nx*nu, sizeof(double));
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, G[i*N+j-1], &nx, &zero, G[i*N+j], &nx);
        }
    }
         
    /* compute L */
    memcpy(&L[0],&ds0[0], nx*sizeof(double)); 
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[0], i);       
        cell=mxGetPr(cell_element); /* Ai */
        memcpy(&L[(i+1)*nx],&a[i*nx], nx*sizeof(double)); 
        dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,L+i*nx,&one_i,&one_d,L+(i+1)*nx,&one_i);
    }
    
    /* compute gc */
    cell_element = mxGetCell(prhs[2], N);
    cell = mxGetPr(cell_element); /* Qi */
    memcpy(&w_vec[0],&gs[N*nx],nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,L+N*nx,&one_i,&one_d,w_vec,&one_i);
    for(i=N-1;i>0;i--){
        cell_element = mxGetCell(prhs[3], i);
        cell = mxGetPr(cell_element); /* Si */
        memcpy(&gc[i*nu],&gu[i*nu],nu*sizeof(double));
        dgemv(Trans,&nx,&nu,&one_d,cell,&nx,L+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
        cell_element = mxGetCell(prhs[1], i);
        cell = mxGetPr(cell_element); /* Bi */
        dgemv(Trans,&nx,&nu,&one_d,cell,&nx,w_vec,&one_i,&one_d,gc+i*nu,&one_i);
         
        cell_element = mxGetCell(prhs[2], i);
        cell = mxGetPr(cell_element); /* Qi */
        memcpy(&w_vec_dup[0],&gs[i*nx],nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,L+i*nx,&one_i,&one_d,w_vec_dup,&one_i);
        cell_element = mxGetCell(prhs[0], i);
        cell = mxGetPr(cell_element); /* Ai */
        dgemv(Trans,&nx,&nx,&one_d,cell,&nx,w_vec,&one_i,&one_d,w_vec_dup,&one_i);
        memcpy(&w_vec[0],&w_vec_dup[0],nx*sizeof(double));
    }   
    cell_element = mxGetCell(prhs[3], 0);
    cell = mxGetPr(cell_element); /* Si */
    memcpy(&gc[0],&gu[0],nu*sizeof(double));
    dgemv(Trans,&nx,&nu,&one_d,cell,&nx,L,&one_i,&one_d,gc,&one_i);
    cell_element = mxGetCell(prhs[1], 0);
    cell = mxGetPr(cell_element); /* Bi */
    dgemv(Trans,&nx,&nu,&one_d,cell,&nx,w_vec,&one_i,&one_d,gc,&one_i);
     
    /* Compute Hc (only the lower triangular part) */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[2], N);
        cell = mxGetPr(cell_element); /* Qi */
        dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, G[i*N+N-1], &nx, &zero, W_mat, &nx);
        for(j=N-1;j>i;j--){
            
            cell_element = mxGetCell(prhs[3], j);
            cell = mxGetPr(cell_element); /* Si */
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, G[i*N+j-1], &nx, &zero, Hi, &nu);
            cell_element = mxGetCell(prhs[1], j);
            cell = mxGetPr(cell_element); /* Bi */
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);
             
            cell_element = mxGetCell(prhs[0], j);
            cell = mxGetPr(cell_element); /* Ai */
            dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &zero, W_mat_dup, &nx); 
            cell_element = mxGetCell(prhs[2], j);
            cell = mxGetPr(cell_element); /* Qi */
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, G[i*N+j-1], &nx, &one_d, W_mat_dup, &nx); 
            memcpy(&W_mat[0],&W_mat_dup[0],nx*nu*sizeof(double));
        }
        cell_element = mxGetCell(prhs[4], i);
        cell = mxGetPr(cell_element); /* Ri */
        memcpy(Hi,cell,nu*nu*sizeof(double));
        cell_element = mxGetCell(prhs[1], i); 
        cell = mxGetPr(cell_element); /* Bi */
        dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &one_d, Hi, &nu);
        Block_Fill(nu, nu, Hi, Hc, i*nu, i*nu, N*nu);
    }    
    
    /* fill the upper triangular part of Hc (Hc is symmetric) */
    for(i=0;i<N*nu;i++){
        for(j=i+1;j<N*nu;j++)
            Hc[j*N*nu+i]=Hc[i*N*nu+j];
    }
    
    /* Compute Cc */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[6], i);
        cell = mxGetPr(cell_element); /* Cui */ 
        Block_Fill(nc, nu, cell, Cc_all, i*nc, i*nu, N*nc+ncN);
                
        for(j=i+1;j<N;j++){   
            cell_element = mxGetCell(prhs[5], j);
            cell = mxGetPr(cell_element);   /* Csi */
            dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, cell, &nc, G[i*N+j-1], &nx, &zero, Cci, &nc);
            Block_Fill(nc, nu, Cci, Cc_all, j*nc, i*nu, N*nc+ncN);
        }    
    }
    
    /* Compute cc */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[5], i);
        cell = mxGetPr(cell_element); /* Csi */       
        memcpy(&cc[i*nc],&c[i*nc],nc*sizeof(double));
        dgemv(nTrans,&nc,&nx,&one_d,cell,&nc,L+i*nx,&one_i,&one_d,cc+i*nc,&one_i);        
    }
    
    /* Compute CcN and ccN */
    cell_element = mxGetCell(prhs[5], N);
    cell = mxGetPr(cell_element);  /* CsN */      
    for(i=0;i<N;i++){                 
        dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, cell, &ncN, G[i*N+N-1], &nx, &zero, CcN, &ncN);
        Block_Fill(ncN, nu, CcN, Cc_all, N*nc, i*nu, N*nc+ncN);
    }
   
    memcpy(&cc[N*nc],&c[N*nc],ncN*sizeof(double));
    dgemv(nTrans,&ncN,&nx,&one_d,cell,&ncN,L+N*nx,&one_i,&one_d,cc+N*nc,&one_i);
    
        
    /* Free memory */
    for(i=0;i<N;i++){
        for (j=i+1;j<N;j++)
            mxFree(G[i*N+j]);
    }
    mxFree(G);

}