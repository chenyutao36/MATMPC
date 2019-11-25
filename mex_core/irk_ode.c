
#include "mex.h"
#include <stdlib.h>
#include "string.h"

#include "sim.h"
#include "irk_ode.h"
#include "mpc_common.h"
#include "casadi_wrapper.h"

#include "blas.h"
#include "lapack.h"

sim_irk_ode_workspace* sim_irk_ode_workspace_create(sim_opts *opts)
{
    sim_irk_ode_workspace* work = (sim_irk_ode_workspace*)mxMalloc(sizeof(sim_irk_ode_workspace));
        
    work->impl_ode_in = mxMalloc(5*sizeof(double *));
    work->res_out = mxMalloc(sizeof(double *));
    work->jac_x_out = mxMalloc(sizeof(double *));
    work->jac_u_out = mxMalloc(sizeof(double *));
    work->jac_xdot_out = mxMalloc(sizeof(double *));
    work->tmp_nx_nx = mxMalloc(opts->nx*opts->nx*sizeof(double));
    
    mexMakeMemoryPersistent(work->impl_ode_in);
    mexMakeMemoryPersistent(work->res_out);
    mexMakeMemoryPersistent(work->jac_x_out);
    mexMakeMemoryPersistent(work->jac_u_out);
    mexMakeMemoryPersistent(work->jac_xdot_out);
    mexMakeMemoryPersistent(work->tmp_nx_nx);
    
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

    work->IPIV = mxMalloc(opts->num_stages*opts->nx*sizeof(size_t));
    mexMakeMemoryPersistent(work->IPIV);
    
    if (opts->forw_sens_flag || opts->adj_sens_flag){
        work->JGx = mxMalloc(opts->num_stages*opts->nx*opts->nx*sizeof(double));
        work->JGu = mxMalloc(opts->num_stages*opts->nx*opts->nu*sizeof(double));
        work->JKx = mxMalloc(opts->num_stages*opts->nx*opts->nx*sizeof(double));
        work->JKu = mxMalloc(opts->num_stages*opts->nx*opts->nu*sizeof(double));
                
        mexMakeMemoryPersistent(work->JGx);
        mexMakeMemoryPersistent(work->JGu);
        mexMakeMemoryPersistent(work->JKx);
        mexMakeMemoryPersistent(work->JKu);
        
    }

    if (opts->adj_sens_flag){
        work->JGK_fact_traj = mxMalloc(opts->num_stages*opts->nx*opts->num_stages*opts->nx*opts->num_steps*sizeof(double));
        work->JGx_traj = mxMalloc(opts->num_stages*opts->nx*opts->nx*opts->num_steps*sizeof(double));
        work->JGu_traj = mxMalloc(opts->num_stages*opts->nx*opts->nu*opts->num_steps*sizeof(double));
        work->IPIV_traj = mxMalloc(opts->num_stages*opts->nx *opts->num_steps*sizeof(size_t));
        work->lambda_K = mxCalloc(opts->num_stages*opts->nx, sizeof(double));

        mexMakeMemoryPersistent(work->JGK_fact_traj);
        mexMakeMemoryPersistent(work->JGx_traj);
        mexMakeMemoryPersistent(work->JGu_traj);
        mexMakeMemoryPersistent(work->IPIV_traj);
        mexMakeMemoryPersistent(work->lambda_K);
    }
    
    mexMakeMemoryPersistent(work);
    
    return work;
}

void sim_irk_ode_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_ode_workspace *work)
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

int sim_irk_ode(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_ode_workspace *workspace){    
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
    bool forw_sens_flag = opts->forw_sens_flag;
    bool adj_sens_flag = opts->adj_sens_flag;
      
    double *x = in->x;
    double *z = in->z;
    double *lambda = in->lambda;

    double *xn = out->xn;   
    double *Jac_x = out->Sx;
    double *Jac_u = out->Su;
    double *adj_sens = out->adj_sens;
    
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
    double *JGK_fact_traj = workspace->JGK_fact_traj;
    double *JGx_traj = workspace->JGx_traj;
    double *JGu_traj = workspace->JGu_traj;
    double *lambda_K = workspace->lambda_K;
    double **impl_ode_in = workspace->impl_ode_in;
    double **res_out = workspace->res_out;
    double **jac_x_out = workspace->jac_x_out;
    double **jac_u_out = workspace->jac_u_out;
    double **jac_xdot_out = workspace->jac_xdot_out;
    size_t *IPIV = workspace->IPIV;
    size_t *IPIV_traj = workspace->IPIV_traj;
    int newton_iter = workspace->newton_iter;

    double *tmp_nx_nx = workspace->tmp_nx_nx;
    
    double measure;
       
    size_t jx = nx*nx;
    size_t ju = nx*nu;
    size_t ldG = num_stages*nx;  
    size_t INFO = 0;
    size_t nw = nx+nu;
            
    // initialize
    memcpy(xn, x, nx*sizeof(double));
    set_zeros(ldG, K);
    
    if (forw_sens_flag){
        memcpy(Jac_x, Sx, nx*nx*sizeof(double));
        memcpy(Jac_u, Su, nx*nu*sizeof(double));
    }

    if (adj_sens_flag){
        memcpy(adj_sens, lambda, nx*sizeof(double));

        // initialize the u-part to be zero
        for (j=0;j<nu;j++)
            adj_sens[nx+j] = 0.0;
    }
    
    impl_ode_in[1] = in->u; // u
    impl_ode_in[2] = in->p; // p
    impl_ode_in[4] = in->z; // z
      
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
                impl_jac_f_x_Fun(impl_ode_in, jac_x_out);
                impl_jac_f_xdot_Fun(impl_ode_in, jac_xdot_out);
                                                            
                for (j=0; j<num_stages; j++){ //compute the block (ii,jj)th block = Jt
                    a = A[j * num_stages + i];
                    memcpy(tmp_nx_nx, jac_x_out[0], nx*nx*sizeof(double));
                    a *= h;
                    dscal(&jx, &a, tmp_nx_nx, &one_i);                   
                    if(j==i){
                        daxpy(&jx, &one_d, jac_xdot_out[0], &one_i, tmp_nx_nx, &one_i);
                    }
                    // fill in the i-th, j-th block of JGK
                    Block_Fill(nx, nx, tmp_nx_nx, JGK, i*nx, j*nx, ldG);
                } // end j
            }// end i  
                  
            dgetrf(&ldG, &ldG, JGK, &ldG, IPIV, &INFO);	

            if (INFO < 0){
                mexPrintf("The %d th argument had an illegal value\n", -1.0*INFO);
                mexErrMsgTxt("The linear solver failed in the integrator!");
            }
            if (INFO > 0){
                mexPrintf("The factorized U(%d,%d) is zero\n", INFO, INFO);
                mexErrMsgTxt("The linear solver failed in the integrator!");
            }
                      
            dgetrs(nTrans, &ldG, &one_i, JGK, &ldG, IPIV, rG, &ldG, &INFO);
                                    
            daxpy(&ldG, &minus_one_d, rG, &one_i, K, &one_i);

        }//end Newton iter
      
        measure = dnrm2(&ldG, rG, &one_i);
        if (measure>1E-2){
            mexPrintf("Integrator residual is %5.3e\n",measure);
            mexErrMsgTxt("The integrator has not converged!");
        }
               
        // after L steps, evaluate derivatives dG/dK, dG/dx, dG/du at K^[L]
        // evaluate forward sens
        if (forw_sens_flag || adj_sens_flag){

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
                
                impl_jac_f_x_Fun(impl_ode_in, jac_x_out);
                impl_jac_f_xdot_Fun(impl_ode_in, jac_xdot_out);
                impl_jac_f_u_Fun(impl_ode_in, jac_u_out);
                Block_Fill(nx, nx, jac_x_out[0], JGx, i*nx, 0, ldG);                   
                if (adj_sens_flag){
                    Block_Fill(nx, nu, jac_u_out[0], JGu, i*nx, 0, ldG);
                }else{
                    Block_Fill(nx, nu, jac_u_out[0], JKu, i*nx, 0, ldG);
                }
                                              
                for (j=0; j<num_stages; j++){ //compute the block (ii,jj)th block = Jt
                    a = A[j * num_stages + i];
                    memcpy(tmp_nx_nx, jac_x_out[0], nx*nx*sizeof(double));
                    a *= h;
                    dscal(&jx, &a, tmp_nx_nx, &one_i);
                    if(j==i){
                        daxpy(&jx, &one_d, jac_xdot_out[0], &one_i, tmp_nx_nx, &one_i);
                    }
                    // fill in the i-th, j-th block of JGK
                    Block_Fill(nx, nx, tmp_nx_nx, JGK, i*nx, j*nx, ldG);
                } // end j
            }// end i             
            
            dgetrf(&ldG, &ldG, JGK, &ldG, IPIV, &INFO);	

            if (INFO < 0){
                mexPrintf("The %d th argument had an illegal value\n", -1.0*INFO);
                mexErrMsgTxt("The linear solver failed in the integrator!");
            }
            if (INFO > 0){
                mexPrintf("The U(%d,%d) is zero\n", INFO, INFO);
                mexErrMsgTxt("The linear solver failed in the integrator!");
            }

            if (forw_sens_flag){
                dgemm(nTrans, nTrans, &ldG, &nx, &nx, &one_d, JGx, &ldG, Jac_x, &nx, &zero, JKx, &ldG);
                dgemm(nTrans, nTrans, &ldG, &nu, &nx, &one_d, JGx, &ldG, Jac_u, &nx, &one_d, JKu, &ldG);

                dgetrs(nTrans, &ldG, &nx, JGK, &ldG, IPIV, JKx, &ldG, &INFO);
                dgetrs(nTrans, &ldG, &nu, JGK, &ldG, IPIV, JKu, &ldG, &INFO);

                dgemm(nTrans, nTrans, &nx, &nx, &ldG, &minus_one_d, JFK, &nx, JKx, &ldG, &one_d, Jac_x, &nx);
                dgemm(nTrans, nTrans, &nx, &nu, &ldG, &minus_one_d, JFK, &nx, JKu, &ldG, &one_d, Jac_u, &nx);
            }
        } 

        // store data at each integration step
        if (adj_sens_flag){
            memcpy(JGK_fact_traj+istep*num_stages*nx*num_stages*nx, JGK, num_stages * nx * num_stages * nx * sizeof(double));
            memcpy(JGx_traj+istep*num_stages*nx*nx, JGx, num_stages * nx * nx * sizeof(double));
            memcpy(JGu_traj+istep*num_stages*nx*nu, JGu, num_stages * nx * nu * sizeof(double));
            memcpy(IPIV_traj+istep*num_stages*nx, IPIV, num_stages * nx* sizeof(size_t));
        }

        // update xn
        for (i = 0; i < num_stages; i++){
            b = h * B[i];
            daxpy(&nx, &b, K+i*nx, &one_i, xn, &one_i);
        }     
        
    }// end num_steps

    // backward sweep
    if (adj_sens_flag){
        for (istep=num_steps-1;istep>-1;istep--){
            dgemv(Trans, &nx, &ldG, &minus_one_d, JFK, &nx, adj_sens, &one_i, &zero, lambda_K, &one_i);
            dgetrs(Trans, &ldG, &one_i, JGK_fact_traj+istep*num_stages*nx*num_stages*nx, &ldG, IPIV_traj+istep*num_stages*nx, lambda_K, &ldG, &INFO);
            dgemv(Trans, &ldG, &nx, &one_d, JGx_traj+istep*num_stages*nx*nx, &ldG, lambda_K, &one_i, &one_d, adj_sens, &one_i);
            dgemv(Trans, &ldG, &nu, &one_d, JGu_traj+istep*num_stages*nx*nu, &ldG, lambda_K, &one_i, &one_d, adj_sens+nx, &one_i);
        }
    }
    
    return 0;
    
}

void sim_irk_ode_workspace_free(sim_opts *opts, sim_irk_ode_workspace *work)
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

    mxFree(work->tmp_nx_nx);
    
    if (opts->forw_sens_flag || opts->adj_sens_flag){
        mxFree(work->JGx);
        mxFree(work->JGu);
        mxFree(work->JKx);
        mxFree(work->JKu);        
    }

    if (opts->adj_sens_flag){
        mxFree(work->JGK_fact_traj);
        mxFree(work->JGx_traj);
        mxFree(work->JGu_traj);
        mxFree(work->IPIV_traj);
        mxFree(work->lambda_K);
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