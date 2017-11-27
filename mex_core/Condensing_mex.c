
#include "mex.h"
// #include "casadi_wrapper.h"

// for builtin blas
#include "blas.h"

// for openblas
// #include "f77blas.h"
// #if !defined(_WIN32)
// #define dgemm dgemm_
// #define dgemv dgemv_
// #endif

#include "string.h"

/*
 *r: block row index
 *c: block column index
 *Gi: block data that is filled in
 *rows: number of rows of the block Gi
 *cols: number of columns of the block Gi
 *G: the original matrix
 *N: the number of row blocks in G
 *
 *These functions assume all blocks in G have the same dimension
*/
void Block_Fill(mwIndex r, mwIndex c, double *Gi, mwSize rows, mwSize cols, double *G, mwSize N){
        
    mwIndex m,n,spo,spi;
       
    spo = r*rows+c*cols*N*rows;
    for(n=0;n<cols;n++){
        spi=spo+ n*N*rows;
        for(m=0;m<rows;m++){
            G[spi+m]=Gi[m+n*rows];
        }
    }
}

void Block_Access(mwIndex r, mwIndex c, double *Gi, mwSize rows, mwSize cols, double *G, mwSize N){
       
    mwIndex m,n,spo,spi;
    
    spo = r*rows+c*cols*N*rows;
    for(n=0;n<cols;n++){
        spi=spo+ n*N*rows;
        for(m=0;m<rows;m++){
            Gi[m+n*rows]=G[spi+m];
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
    mwSize nx = mxGetScalar(prhs[12]);
    mwSize nu = mxGetScalar(prhs[13]);
    mwSize nc = mxGetScalar(prhs[14]);
    mwSize ncN = mxGetScalar(prhs[15]);
    mwSize N = mxGetScalar(prhs[16]);
        
    mwSize nz = nx+nu;
    
    /*Outputs*/
    double *G, *L, *gc, *Hc, *Cc, *cc, *CcN, *cN;
    plhs[0] = mxCreateDoubleMatrix(N*nx, N*nu, mxREAL);
    G = mxGetPr(plhs[0]);    
    plhs[1] = mxCreateDoubleMatrix((N+1)*nx, 1, mxREAL);
    L = mxGetPr(plhs[1]);   
    plhs[2] = mxCreateDoubleMatrix(N*nu, 1, mxREAL);
    gc = mxGetPr(plhs[2]);   
    plhs[3] = mxCreateDoubleMatrix(N*nu, N*nu, mxREAL);
    Hc = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(N*nc, N*nu, mxREAL);
    Cc = mxGetPr(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix(N*nc, 1, mxREAL);
    cc = mxGetPr(plhs[5]);
    plhs[6] = mxCreateDoubleMatrix(ncN, N*nu, mxREAL);
    CcN = mxGetPr(plhs[6]);
    plhs[7] = mxCreateDoubleMatrix(ncN, 1, mxREAL);
    cN = mxGetPr(plhs[7]);
    
    /*Allocate memory*/
    mwIndex i=0,j=0;
    const mxArray *cell_element;  
    double *cell; /* from cells */
    double *Gi, *Ci, *Li, *ai, *gsi, *gui, *w_vec, *W_mat, *Hi, *Cci, *ci, *CuN; /* Intermediate data */
    Gi = (double *)mxCalloc(nx*nu, sizeof(double));
    Ci = (double *)mxCalloc(nx*nu, sizeof(double));
    Li = (double *)mxCalloc(nx, sizeof(double));
    ai = (double *)mxCalloc(nx, sizeof(double));
    gsi = (double *)mxCalloc(nx, sizeof(double));    
    gui = (double *)mxCalloc(nu, sizeof(double)); 
    w_vec = (double *)mxCalloc(nx, sizeof(double));
    W_mat = (double *)mxCalloc(nx*nu, sizeof(double));   
    Hi = (double *)mxCalloc(nu*nu, sizeof(double));
    Cci = (double *)mxCalloc(nc*nu, sizeof(double));
    ci = (double *)mxCalloc(nc, sizeof(double));
    CuN = (double *)mxCalloc(ncN*nu, sizeof(double));   
    
    /*Initialization*/
    memcpy(&ai[0],&ds0[0], nx*sizeof(double));  

    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1; /* never use int for lapack and blas routines */
              
    /*Start the loop*/
    
    /* compute G */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[1], i);       
        cell = mxGetPr(cell_element); /* Bi */
        Block_Fill(i,i,cell,nx,nu,G,N); 
        
        for (j=i+1;j<N;j++){
            cell_element = mxGetCell(prhs[0], j);           
            cell = mxGetPr(cell_element); /* Ai */
            Block_Access(j-1,i,Gi,nx,nu,G,N);
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, Gi, &nx, &zero, Ci, &nx);
            Block_Fill(j,i,Ci,nx,nu,G,N);
        }
    }
         
    /* compute L */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[0], i);       
        cell=mxGetPr(cell_element); /* Ai */
        memcpy(&L[i*nx],&ai[0], nx*sizeof(double));
        memcpy(&ai[0],&a[i*nx], nx*sizeof(double));
        memcpy(&Li[0],&L[i*nx], nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,Li,&one_i,&one_d,ai,&one_i); 
    }
    memcpy(&L[N*nx],&ai[0], nx*sizeof(double));
    
    
    /* compute gc */
    cell_element = mxGetCell(prhs[2], N);
    cell = mxGetPr(cell_element); /* Qi */
    memcpy(&Li[0],&L[N*nx], nx*sizeof(double));
    memcpy(&gsi[0],&gs[N*nx],nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,Li,&one_i,&one_d,gsi,&one_i);
    memcpy(&w_vec[0],&gsi[0],nx*sizeof(double));
    for(i=N-1;i>0;i--){
        cell_element = mxGetCell(prhs[3], i);
        cell = mxGetPr(cell_element); /* Si */
        memcpy(&Li[0],&L[i*nx],nx*sizeof(double));
        memcpy(&gui[0],&gu[i*nu],nu*sizeof(double));
        dgemv(Trans,&nx,&nu,&one_d,cell,&nx,Li,&one_i,&one_d,gui,&one_i);
        cell_element = mxGetCell(prhs[1], i);
        cell = mxGetPr(cell_element); /* Bi */
        dgemv(Trans,&nx,&nu,&one_d,cell,&nx,w_vec,&one_i,&one_d,gui,&one_i);
        memcpy(&gc[i*nu],&gui[0],nu*sizeof(double));
         
        cell_element = mxGetCell(prhs[2], i);
        cell = mxGetPr(cell_element); /* Qi */
        memcpy(&gsi[0],&gs[i*nx],nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,cell,&nx,Li,&one_i,&one_d,gsi,&one_i);
        cell_element = mxGetCell(prhs[0], i);
        cell = mxGetPr(cell_element); /* Ai */
        dgemv(Trans,&nx,&nx,&one_d,cell,&nx,w_vec,&one_i,&one_d,gsi,&one_i);
        memcpy(&w_vec[0],&gsi[0],nx*sizeof(double));
    }   
    cell_element = mxGetCell(prhs[3], 0);
    cell = mxGetPr(cell_element); /* Si */
    memcpy(&Li[0],&L[0],nx*sizeof(double));
    memcpy(&gui[0],&gu[0],nu*sizeof(double));
    dgemv(Trans,&nx,&nu,&one_d,cell,&nx,Li,&one_i,&one_d,gui,&one_i);
    cell_element = mxGetCell(prhs[1], 0);
    cell = mxGetPr(cell_element); /* Bi */
    dgemv(Trans,&nx,&nu,&one_d,cell,&nx,w_vec,&one_i,&one_d,gui,&one_i);
    memcpy(&gc[0],&gui[0],nu*sizeof(double));
     
    /* Compute Hc (only the lower triangular part) */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[2], N);
        cell = mxGetPr(cell_element); /* Qi */
        Block_Access(N-1,i,Gi,nx,nu,G,N);
        dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, Gi, &nx, &zero, W_mat, &nx);
        for(j=N-1;j>i;j--){
            Block_Access(j-1,i,Gi,nx,nu,G,N);
            
            cell_element = mxGetCell(prhs[3], j);
            cell = mxGetPr(cell_element); /* Si */
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, Gi, &nx, &zero, Hi, &nu);
            cell_element = mxGetCell(prhs[1], j);
            cell = mxGetPr(cell_element); /* Bi */
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &one_d, Hi, &nu);
            Block_Fill(j,i,Hi,nu,nu,Hc,N);           
             
            cell_element = mxGetCell(prhs[0], j);
            cell = mxGetPr(cell_element); /* Ai */
            dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &zero, Ci, &nx); 
            cell_element = mxGetCell(prhs[2], j);
            cell = mxGetPr(cell_element); /* Qi */
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, cell, &nx, Gi, &nx, &one_d, Ci, &nx);                   
            memcpy(&W_mat[0],&Ci[0],nx*nu*sizeof(double));
        }
        cell_element = mxGetCell(prhs[4], i);
        cell = mxGetPr(cell_element); /* Ri */
        memcpy(Hi,cell,nu*nu*sizeof(double));
        cell_element = mxGetCell(prhs[1], i); 
        cell = mxGetPr(cell_element); /* Bi */
        dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, cell, &nx, W_mat, &nx, &one_d, Hi, &nu);
        Block_Fill(i,i,Hi,nu,nu,Hc,N);
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
        Block_Fill(i,i,cell,nc,nu,Cc,N);
                
        for(j=i+1;j<N;j++){   
            cell_element = mxGetCell(prhs[5], j);
            cell = mxGetPr(cell_element);   /* Csi */        
            Block_Access(j-1,i,Gi,nx,nu,G,N);
            dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, cell, &nc, Gi, &nx, &zero, Cci, &nc);
            Block_Fill(j,i,Cci,nc,nu,Cc,N);
        }    
    }
    
    /* Compute cc */
    for(i=0;i<N;i++){
        cell_element = mxGetCell(prhs[5], i);
        cell = mxGetPr(cell_element); /* Csi */ 
        memcpy(&Li[0],&L[i*nx],nx*sizeof(double));
        memcpy(&ci[0],&c[i*nc],nc*sizeof(double));
        dgemv(nTrans,&nc,&nx,&one_d,cell,&nc,Li,&one_i,&one_d,ci,&one_i);
        memcpy(&cc[i*nc],&ci[0],nc*sizeof(double));
    }
    
    /* Compute CcN and ccN */
    cell_element = mxGetCell(prhs[5], N);
    cell = mxGetPr(cell_element);  /* CsN */      
    for(i=0;i<N;i++){                 
        Block_Access(N-1,i,Gi,nx,nu,G,N);
        dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, cell, &ncN, Gi, &nx, &zero, CuN, &ncN);
        memcpy(&CcN[i*ncN*nu],&CuN[0],ncN*nu*sizeof(double));
    }
    
    memcpy(&Li[0],&L[N*nx],nx*sizeof(double));
    memcpy(&cN[0],&c[N*nc],ncN*sizeof(double));
    dgemv(nTrans,&ncN,&nx,&one_d,cell,&ncN,Li,&one_i,&one_d,cN,&one_i);
        
    /* Free memory */
    mxFree(Gi);
    mxFree(Ci);
    mxFree(Li);
    mxFree(ai);
    mxFree(gsi);
    mxFree(gui);
    mxFree(w_vec);
    mxFree(W_mat);
    mxFree(Hi);
    mxFree(Cci);
    mxFree(ci);
    mxFree(CuN);
       
}