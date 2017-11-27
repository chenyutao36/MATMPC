#include "mex.h"
#include "string.h"
#include "blas.h"
#include <stdbool.h>

static double **jac_out = NULL;
static double *K = NULL, *dKx = NULL, *dKu = NULL;
static double *jacX_t = NULL, *jacU_t = NULL, *xt = NULL;
static double **vde_in = NULL, **vde_out = NULL;

static bool mem_alloc_erk = false;

void exitFcn(){
    if (mem_alloc_erk){
        mxFree(jac_out[0]);
        mxFree(jac_out[1]);
        mxFree(jac_out);
        mxFree(K);
        mxFree(dKx);
        mxFree(dKu);
        mxFree(jacX_t);
        mxFree(jacU_t);
        mxFree(xt);
        mxFree(vde_in);
        mxFree(vde_out);
    }
}

void sim_erk(double **in, double **out, mxArray **Jac, mxArray *sim_mem){
    
    int istep, s, i, j;
    double a,b;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
    
    double *x = in[0];
    double *xn = out[0];
   
    double *Jac_x = mxGetPr(Jac[0]);
    double *Jac_u = mxGetPr(Jac[1]);
    
    double *A = mxGetPr( mxGetField(sim_mem, 0, "A_erk"));
    double *B = mxGetPr( mxGetField(sim_mem, 0, "B_erk"));
    double h = mxGetScalar( mxGetField(sim_mem, 0, "h") );
    mwSize num_stages = mxGetScalar( mxGetField(sim_mem, 0, "num_stages") );
    mwSize num_steps = mxGetScalar( mxGetField(sim_mem, 0, "num_steps") );
    mwSize nx = mxGetScalar( mxGetField(sim_mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(sim_mem, 0, "nu") );
    double *Sx = mxGetPr( mxGetField(sim_mem, 0, "Sx"));
    double *Su = mxGetPr( mxGetField(sim_mem, 0, "Su"));
        
    mwSize jx = nx*nx;
    mwSize ju = nx*nu;
    
    if (!mem_alloc_erk){
        K = (double *) mxCalloc(num_stages*nx, sizeof(double));
        mexMakeMemoryPersistent(K);
        
        dKx = (double *) mxCalloc(num_stages*nx*nx, sizeof(double));
        mexMakeMemoryPersistent(dKx);
        
        dKu = (double *) mxCalloc(num_stages*nx*nu, sizeof(double));
        mexMakeMemoryPersistent(dKu);
        
        xt = (double *)mxCalloc(nx, sizeof(double));
        mexMakeMemoryPersistent(xt);

        jac_out = (double **) mxMalloc( 2 * sizeof(double *));
        mexMakeMemoryPersistent(jac_out);      
        jac_out[0] = (double *)mxCalloc(nx*nx, sizeof(double));
        mexMakeMemoryPersistent(jac_out[0]);       
        jac_out[1] = (double *)mxCalloc(nx*nu, sizeof(double)); 
        mexMakeMemoryPersistent(jac_out[1]);
        
        jacX_t = (double *)mxCalloc(nx*nx, sizeof(double));
        mexMakeMemoryPersistent(jacX_t);
        jacU_t = (double *)mxCalloc(nx*nu, sizeof(double));
        mexMakeMemoryPersistent(jacU_t);
    
        vde_in = (double **) mxMalloc( 5 * sizeof(double*));
        mexMakeMemoryPersistent(vde_in);
        
        vde_out = (double **) mxMalloc( 3 * sizeof(double*));
        mexMakeMemoryPersistent(vde_out);
        
        mem_alloc_erk = true;
        
        mexAtExit(exitFcn);
    }
    
    memcpy(&xn[0], &x[0], nx*sizeof(double));
    memcpy(&Jac_x[0], &Sx[0], nx*nx*sizeof(double));
    memcpy(&Jac_u[0], &Su[0], nx*nu*sizeof(double));
      
    vde_in[1] = in[1];
    vde_in[2] = in[2];
    
    for (istep = 0; istep < num_steps; istep++) {
        
        for (s = 0; s < num_stages; s++) {
            memcpy(&xt[0], &xn[0], nx*sizeof(double));
            memcpy(&jacX_t[0], &Jac_x[0], nx*nx*sizeof(double));
            memcpy(&jacU_t[0], &Jac_u[0], nx*nu*sizeof(double));
            for (j = 0; j < s; j++){
                a = A[j * num_stages + s];
                if (a!=0){
                    a *= h;                   
                    daxpy(&nx, &a, K+j*nx, &one_i, xt, &one_i);

                    daxpy(&jx, &a, dKx+j*jx, &one_i, jacX_t, &one_i);

                    daxpy(&ju, &a, dKu+j*ju, &one_i, jacU_t, &one_i);
                                 
                }
            }
            
            vde_in[0] = xt;
            vde_in[3] = jacX_t;
            vde_in[4] = jacU_t;
            vde_out[0] = K+s*nx;
            vde_out[1] = dKx+s*jx;
            vde_out[2] = dKu+s*ju;
            vde_Fun(vde_in, vde_out);
          
        }
        
        for (s = 0; s < num_stages; s++){
            b = h * B[s];
            daxpy(&nx, &b, K+s*nx, &one_i, xn, &one_i);
            daxpy(&jx, &b, dKx+s*jx, &one_i, Jac_x, &one_i);
            daxpy(&ju, &b, dKu+s*ju, &one_i, Jac_u, &one_i);
        }
        
    }
    
}