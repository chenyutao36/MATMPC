#include "mex.h"
#include "string.h"
#include "blas.h"
#include "lapack.h"
#include <stdbool.h>

static double *K = NULL;
static double *rG = NULL;
static double *JGK = NULL, *JGx = NULL, *JGu = NULL, *JKx = NULL, *JKu = NULL;
static double *xt = NULL;
static double **impl_ode_out = NULL, **impl_ode_in = NULL;
static mwSize *IPIV = NULL;
static bool mem_alloc = false;

void exitFcn(){
    if (mem_alloc){
        mxFree(xt);
        mxFree(K);
        mxFree(rG);
        mxFree(JGK);
        mxFree(JGx);
        mxFree(JGu);
        mxFree(JKx);
        mxFree(JKu);
        mxFree(impl_ode_in);
        mxFree(impl_ode_out[1]);
        mxFree(impl_ode_out[2]);
        mxFree(impl_ode_out[3]);
        mxFree(impl_ode_out);
        mxFree(IPIV);
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

void sim_irk(double **in, double **out, mxArray **Jac, mxArray *mem){
    
    int istep, i, j, k, iter;
    double a,b;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
    
    
    double *x = in[0];
    double *xn = out[0];
   
    double *Jac_x = mxGetPr(Jac[0]);
    double *Jac_u = mxGetPr(Jac[1]);
    
    double *A = mxGetPr( mxGetField(mem, 0, "A"));
    double *B = mxGetPr( mxGetField(mem, 0, "B"));
    double h = mxGetScalar( mxGetField(mem, 0, "h") );
    mwSize num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    mwSize num_steps = mxGetScalar( mxGetField(mem, 0, "num_steps") );
    mwSize nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    double *Sx = mxGetPr( mxGetField(mem, 0, "Sx"));
    double *Su = mxGetPr( mxGetField(mem, 0, "Su"));
    mwSize newton_iter = mxGetScalar( mxGetField(mem, 0, "newton_iter") );
    double *JFK = mxGetPr( mxGetField(mem, 0, "JFK"));
        
    mwSize jx = nx*nx;
    mwSize ju = nx*nu;
    mwSize ldG = num_stages*nx;
    
    mwIndex INFO = 0;
    
    if (!mem_alloc){       
        K = (double *) mxCalloc(ldG, sizeof(double));
        mexMakeMemoryPersistent(K);
        
        rG = (double *) mxCalloc(ldG, sizeof(double));
        mexMakeMemoryPersistent(rG);
        
        JGK = (double *) mxCalloc(ldG*ldG, sizeof(double));
        mexMakeMemoryPersistent(JGK);
        
        JGx = (double *) mxCalloc(ldG*nx, sizeof(double));
        mexMakeMemoryPersistent(JGx);
        
        JKx = (double *) mxCalloc(ldG*nx, sizeof(double));
        mexMakeMemoryPersistent(JKx);
        
        JGu = (double *) mxCalloc(ldG*nu, sizeof(double));
        mexMakeMemoryPersistent(JGu);
        
        JKu = (double *) mxCalloc(ldG*nu, sizeof(double));
        mexMakeMemoryPersistent(JKu);
       
        xt = (double *)mxCalloc(nx, sizeof(double));
        mexMakeMemoryPersistent(xt);
        
        IPIV = (mwIndex *)mxMalloc(ldG * sizeof(mwIndex));
        for (i=0;i<ldG;i++)
            IPIV[i] = i+1;
        mexMakeMemoryPersistent(IPIV);
        
        impl_ode_in = (double **) mxMalloc(4*sizeof(double*));
        mexMakeMemoryPersistent(impl_ode_in);
        
        impl_ode_out = (double **) mxMalloc( 4 * sizeof(double*));
        impl_ode_out[1] = (double *) mxCalloc(jx, sizeof(double)); // df/d(x)
        impl_ode_out[2] = (double *) mxCalloc(ju, sizeof(double)); // df/d(u)
        impl_ode_out[3] = (double *) mxCalloc(jx, sizeof(double)); // df/d(x dot)
        mexMakeMemoryPersistent(impl_ode_out);
        mexMakeMemoryPersistent(impl_ode_out[1]);
        mexMakeMemoryPersistent(impl_ode_out[2]);
        mexMakeMemoryPersistent(impl_ode_out[3]);
        
        mem_alloc = true;
        
        mexAtExit(exitFcn);
    }
             
    memcpy(&xn[0], &x[0], nx*sizeof(double));
    memcpy(&Jac_x[0], &Sx[0], nx*nx*sizeof(double));
    memcpy(&Jac_u[0], &Su[0], nx*nu*sizeof(double));
    
    impl_ode_in[1] = in[1]; // u
    impl_ode_in[2] = in[2]; // p
        
    for (istep = 0; istep < num_steps; istep++) {
        
        for(iter=0; iter<newton_iter; iter++){
        
            for (i = 0; i < num_stages; i++) {
                
                memcpy(&xt[0], &xn[0], nx*sizeof(double));
                
                for (j = 0; j < num_stages; j++){
                    a = A[j * num_stages + i];
                    if (a!=0){
                        a *= h;                   
                        daxpy(&nx, &a, K+j*nx, &one_i, xt, &one_i); 
                    }
                }

                impl_ode_in[0] = xt;
                impl_ode_in[3] = K + i*nx; // xdot
                impl_ode_out[0] = rG + i*nx;

                impl_f_Fun(impl_ode_in, impl_ode_out);

                for (j=0; j<num_stages; j++){ //compute the block (ii,jj)th block = Jt
                    a = A[i + num_stages*j];
                    if (a!=0){
                        a *= h;
                        dscal(&jx, &a, impl_ode_out[1], &one_i);
                    }
                    if(j==i){
                        daxpy(&jx, &one_d, impl_ode_out[3], &one_i, impl_ode_out[1], &one_i);
                    }
                    // fill in the i-th, j-th block of JGK
                    Block_Fill(nx, nx, impl_ode_out[1], JGK, i*nx, j*nx, ldG);
                } // end j
            }// end i
            
            
            dgetrf(&ldG, &ldG, JGK, &ldG, IPIV, &INFO);	
                        
            dgetrs(nTrans, &ldG, &one_i, JGK, &ldG, IPIV, rG, &ldG, &INFO);
                        
            daxpy(&ldG, &minus_one_d, rG, &one_i, K, &one_i);
        }//end iter
        
        // evaluate forward sens
        for(i=0; i<num_stages; i++){  
            
            memcpy(&xt[0], &xn[0], nx*sizeof(double));
            
            for(j=0; j<num_stages; j++){ 
                a = A[i+num_stages*j];
                if(a!=0){
                   a *= h;
                   daxpy(&nx, &a, K+j*nx, &one_i, xt, &one_i); 
                }
            } 
           impl_ode_in[0] = xt;
           impl_ode_in[3] = K + i*nx;
           impl_ode_out[0] = rG + i*nx;
      
           impl_f_Fun(impl_ode_in, impl_ode_out);
           
           Block_Fill(nx, nx, impl_ode_out[1], JGx, i*nx, 0, ldG);
           Block_Fill(nx, nu, impl_ode_out[2], JKu, i*nx, 0, ldG);
                
           for (j=0; j<num_stages; j++){ //compute the block (ii,jj)th block = Jt
               a = A[i + num_stages*j];
               if (a!=0){
                   a *= h;
                   dscal(&jx, &a, impl_ode_out[1], &one_i);                       
               }
               if(j==i){
                   daxpy(&jx, &one_d, impl_ode_out[3], &one_i, impl_ode_out[1], &one_i);
               }
               // fill in the i-th, j-th block of JGK
               Block_Fill(nx, nx, impl_ode_out[1], JGK, i*nx, j*nx, ldG);
           } // end j
       }// end i
            
       dgetrf(&ldG, &ldG, JGK, &ldG, IPIV, &INFO);
       
       dgemm(nTrans, nTrans, &ldG, &nx, &nx, &one_d, JGx, &ldG, Jac_x, &nx, &zero, JKx, &ldG);
          
       dgemm(nTrans, nTrans, &ldG, &nu, &nx, &one_d, JGx, &ldG, Jac_u, &nx, &one_d, JKu, &ldG);
       
       dgetrs(nTrans, &ldG, &nx, JGK, &ldG, IPIV, JKx, &ldG, &INFO);
       dgetrs(nTrans, &ldG, &nu, JGK, &ldG, IPIV, JKu, &ldG, &INFO);
       
       dgemm(nTrans, nTrans, &nx, &nx, &ldG, &minus_one_d, JFK, &nx, JKx, &ldG, &one_d, Jac_x, &nx);
       dgemm(nTrans, nTrans, &nx, &nu, &ldG, &minus_one_d, JFK, &nx, JKu, &ldG, &one_d, Jac_u, &nx);
       
        // update xn
        for (i = 0; i < num_stages; i++){
            b = h * B[i];
            daxpy(&nx, &b, K+i*nx, &one_i, xn, &one_i);
        }
        
    }
    
}