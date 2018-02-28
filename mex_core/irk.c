
#include "mex.h"
#include <stdlib.h>
#include "string.h"

#include "sim.h"
#include "irk.h"
#include "mpc_common.h"
#include "casadi_wrapper.h"

#include "blas.h"
#include "lapack.h"

sim_irk_workspace* sim_irk_workspace_create(sim_opts *opts)
{
    sim_irk_workspace* work = (sim_irk_workspace*)mxMalloc(sizeof(sim_irk_workspace));
        
    work->impl_ode_in = mxMalloc(4*sizeof(double *));
    work->res_out = mxMalloc(sizeof(double *));
    work->jac_x_out = mxMalloc(sizeof(double *));
    work->jac_u_out = mxMalloc(sizeof(double *));
    work->jac_xdot_out = mxMalloc(sizeof(double *));
    
    mexMakeMemoryPersistent(work->impl_ode_in);
    mexMakeMemoryPersistent(work->res_out);
    mexMakeMemoryPersistent(work->jac_x_out);
    mexMakeMemoryPersistent(work->jac_u_out);
    mexMakeMemoryPersistent(work->jac_xdot_out);
    
    work->jac_x_out[0] = mxMalloc(opts->nx*opts->nx*sizeof(double));
    work->jac_u_out[0] = mxMalloc(opts->nx*opts->nu*sizeof(double));
    work->jac_xdot_out[0] = mxMalloc(opts->nx*opts->nx*sizeof(double));
    mexMakeMemoryPersistent(work->jac_x_out[0]);
    mexMakeMemoryPersistent(work->jac_u_out[0]);
    mexMakeMemoryPersistent(work->jac_xdot_out[0]);
       
    work->A = mxMalloc(opts->num_stages*opts->num_stages*sizeof(double));
    work->B = mxMalloc(opts->num_stages*sizeof(double));
    work->Sx = mxMalloc(opts->nx*opts->nx*sizeof(double));
    work->Su = mxMalloc(opts->nx*opts->nu*sizeof(double));    
    work->xt = mxMalloc(opts->nx*sizeof(double));
    work->K = mxMalloc(opts->num_stages*opts->nx*sizeof(double));
    work->rG = mxMalloc(opts->num_stages*opts->nx*sizeof(double));
    work->JGK = mxMalloc(opts->num_stages*opts->nx*opts->num_stages*opts->nx*sizeof(double));
    work->JFK = mxMalloc(opts->num_stages*opts->nx*opts->num_stages*opts->nx*sizeof(double));
    
    mexMakeMemoryPersistent(work->A);
    mexMakeMemoryPersistent(work->B);
    mexMakeMemoryPersistent(work->Sx);
    mexMakeMemoryPersistent(work->Su);
    mexMakeMemoryPersistent(work->xt);
    mexMakeMemoryPersistent(work->K);
    mexMakeMemoryPersistent(work->rG);
    mexMakeMemoryPersistent(work->JGK);
    mexMakeMemoryPersistent(work->JFK);
    
    if (opts->forw_sens){
        work->JGx = mxMalloc(opts->num_stages*opts->nx*opts->nx*sizeof(double));
        work->JGu = mxMalloc(opts->num_stages*opts->nx*opts->nu*sizeof(double));
        work->JKx = mxMalloc(opts->num_stages*opts->nx*opts->nx*sizeof(double));
        work->JKu = mxMalloc(opts->num_stages*opts->nx*opts->nu*sizeof(double));
        
        mexMakeMemoryPersistent(work->JGx);
        mexMakeMemoryPersistent(work->JGu);
        mexMakeMemoryPersistent(work->JKx);
        mexMakeMemoryPersistent(work->JKu);
    }
    
    work->IPIV = mxMalloc(opts->num_stages*opts->nx*sizeof(size_t));
    mexMakeMemoryPersistent(work->IPIV);
    
    mexMakeMemoryPersistent(work);
    
    return work;
}

void sim_irk_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_workspace *work)
{
    size_t nx = opts->nx;
    size_t nu = opts->nu;
    size_t num_stages = opts->num_stages;
    
    double *A = mxGetPr( mxGetField(mem, 0, "A_tab"));
    double *B = mxGetPr( mxGetField(mem, 0, "B_tab"));
    double *Sx = mxGetPr( mxGetField(mem, 0, "Sx"));
    double *Su = mxGetPr( mxGetField(mem, 0, "Su"));    
    double *JFK = mxGetPr( mxGetField(mem, 0, "JFK"));
    
    memcpy(work->A, A, num_stages*num_stages*sizeof(double));
    memcpy(work->B, B, num_stages*sizeof(double));
    memcpy(work->Sx, Sx, nx*nx*sizeof(double));
    memcpy(work->Su, Su, nx*nu*sizeof(double));    
    memcpy(work->JFK, JFK, num_stages*nx*nx*sizeof(double));
    
    work->newton_iter = mxGetScalar( mxGetField(mem, 0, "newton_iter") );
}

int sim_irk(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_workspace *workspace){    
    int istep, i, j, k, iter;
    double a,b;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
    
    double h= opts->h;
    size_t nx = opts->nx;
    size_t nu = opts->nu;
    size_t num_stages = opts->num_stages;
    size_t num_steps = opts->num_steps;
    bool forw_sens = opts->forw_sens;
      
    double *x = in->x;
    double *xn = out->xn;   
    double *Jac_x = out->Sx;
    double *Jac_u = out->Su;
    
    double *A = workspace->A;
    double *B = workspace->B;
    double *Sx = workspace->Sx;
    double *Su = workspace->Su;       
           
    double *xt = workspace->xt;
    double *K = workspace->K;
    double *rG = workspace->rG;
    double *JGK = workspace->JGK;
    double *JFK = workspace->JFK;
    double *JGx = workspace->JGx;
    double *JGu = workspace->JGu;
    double *JKx = workspace->JKx;
    double *JKu = workspace->JKu;
    double **impl_ode_in = workspace->impl_ode_in;
    double **res_out = workspace->res_out;
    double **jac_x_out = workspace->jac_x_out;
    double **jac_u_out = workspace->jac_u_out;
    double **jac_xdot_out = workspace->jac_xdot_out;
    size_t *IPIV = workspace->IPIV;
    int newton_iter = workspace->newton_iter;
    
    double measure;
       
    size_t jx = nx*nx;
    size_t ju = nx*nu;
    size_t ldG = num_stages*nx;  
    size_t INFO = 0;
            
    // initialize
    memcpy(xn, x, nx*sizeof(double));
    set_zeros(ldG, K);
    
    if (forw_sens){
        memcpy(Jac_x, Sx, nx*nx*sizeof(double));
        memcpy(Jac_u, Su, nx*nu*sizeof(double));
    }
    
    impl_ode_in[1] = in->u; // u
    impl_ode_in[2] = in->p; // p
      
    // forward sweep
    for (istep = 0; istep < num_steps; istep++) {
        
        for(iter=0; iter<newton_iter; iter++){
        
            for (i = 0; i < num_stages; i++) {
                
                memcpy(xt, xn, nx*sizeof(double));
                
                for (j = 0; j < num_stages; j++){
                    a = A[j * num_stages + i];
                    if (a!=0){
                        a *= h;                   
                        daxpy(&nx, &a, K+j*nx, &one_i, xt, &one_i); 
                    }
                }

                impl_ode_in[0] = xt;
                impl_ode_in[3] = K + i*nx; // xdot
                res_out[0] = rG+i*nx;

                impl_f_Fun(impl_ode_in, res_out);
                impl_jac_x_Fun(impl_ode_in, jac_x_out);
                impl_jac_xdot_Fun(impl_ode_in, jac_xdot_out);
                              
                for (j=0; j<num_stages; j++){ //compute the block (ii,jj)th block = Jt
                    a = A[i + num_stages*j];
                    if (a!=0){
                        a *= h;
                        dscal(&jx, &a, jac_x_out[0], &one_i);
                    }
                    if(j==i){
                        daxpy(&jx, &one_d, jac_xdot_out[0], &one_i, jac_x_out[0], &one_i);
                    }
                    // fill in the i-th, j-th block of JGK
                    Block_Fill(nx, nx, jac_x_out[0], JGK, i*nx, j*nx, ldG);
                } // end j
            }// end i             
            
            dgetrf(&ldG, &ldG, JGK, &ldG, IPIV, &INFO);	
                      
            dgetrs(nTrans, &ldG, &one_i, JGK, &ldG, IPIV, rG, &ldG, &INFO);
                                    
            daxpy(&ldG, &minus_one_d, rG, &one_i, K, &one_i);
        }//end iter
        
        measure = dnrm2(&ldG, rG, &one_i);
        if (measure>1E-2){
            mexPrintf("Implicit ODE residual is %5.3e\n",measure);
            mexErrMsgTxt("Implicit IRK has not converged");
        }
               
        // evaluate forward sens
        if (forw_sens){
            for(i=0; i<num_stages; i++){  

                memcpy(xt, xn, nx*sizeof(double));

                for(j=0; j<num_stages; j++){ 
                    a = A[i+num_stages*j];
                    if(a!=0){
                       a *= h;
                       daxpy(&nx, &a, K+j*nx, &one_i, xt, &one_i); 
                    }
                } 
               impl_ode_in[0] = xt;
               impl_ode_in[3] = K + i*nx;

               impl_jac_x_Fun(impl_ode_in, jac_x_out);
               impl_jac_u_Fun(impl_ode_in, jac_u_out);
               impl_jac_xdot_Fun(impl_ode_in, jac_xdot_out);

               Block_Fill(nx, nx, jac_x_out[0], JGx, i*nx, 0, ldG);
               Block_Fill(nx, nu, jac_u_out[0], JKu, i*nx, 0, ldG);

               for (j=0; j<num_stages; j++){ //compute the block (ii,jj)th block = Jt
                   a = A[i + num_stages*j];
                   if (a!=0){
                       a *= h;
                       dscal(&jx, &a, jac_x_out[0], &one_i);   
                   }
                   if(j==i){
                       daxpy(&jx, &one_d, jac_xdot_out[0], &one_i, jac_x_out[0], &one_i);
                   }
                   // fill in the i-th, j-th block of JGK
                   Block_Fill(nx, nx, jac_x_out[0], JGK, i*nx, j*nx, ldG);
               } // end j
           }// end i

           dgetrf(&ldG, &ldG, JGK, &ldG, IPIV, &INFO);

           dgemm(nTrans, nTrans, &ldG, &nx, &nx, &one_d, JGx, &ldG, Jac_x, &nx, &zero, JKx, &ldG);
           dgemm(nTrans, nTrans, &ldG, &nu, &nx, &one_d, JGx, &ldG, Jac_u, &nx, &one_d, JKu, &ldG);

           dgetrs(nTrans, &ldG, &nx, JGK, &ldG, IPIV, JKx, &ldG, &INFO);
           dgetrs(nTrans, &ldG, &nu, JGK, &ldG, IPIV, JKu, &ldG, &INFO);

           dgemm(nTrans, nTrans, &nx, &nx, &ldG, &minus_one_d, JFK, &nx, JKx, &ldG, &one_d, Jac_x, &nx);
           dgemm(nTrans, nTrans, &nx, &nu, &ldG, &minus_one_d, JFK, &nx, JKu, &ldG, &one_d, Jac_u, &nx);
        }
       
        // update xn
        for (i = 0; i < num_stages; i++){
            b = h * B[i];
            daxpy(&nx, &b, K+i*nx, &one_i, xn, &one_i);
        }
        
    }
    
    return 0;
    
}

void sim_irk_workspace_free(sim_opts *opts, sim_irk_workspace *work)
{
    mxFree(work->A);
    mxFree(work->B);
    mxFree(work->Sx);
    mxFree(work->Su);
    
    mxFree(work->xt);
    mxFree(work->K);
    mxFree(work->rG);
    mxFree(work->JGK);
    mxFree(work->JFK);
    
    if (opts->forw_sens){
        mxFree(work->JGx);
        mxFree(work->JGu);
        mxFree(work->JKx);
        mxFree(work->JKu);
    }
    
    mxFree(work->impl_ode_in);
    mxFree(work->res_out);
    mxFree(work->jac_x_out[0]);
    mxFree(work->jac_u_out[0]);
    mxFree(work->jac_xdot_out[0]);
    mxFree(work->jac_x_out);
    mxFree(work->jac_u_out);
    mxFree(work->jac_xdot_out);
    
    mxFree(work->IPIV);
    
    mxFree(work);
}