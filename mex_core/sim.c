#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "string.h"

#include "sim.h"
#include "mpc_common.h"
#include "blas.h"
#include "lapack.h"

int sim_erk_calculate_workspace_size(const mxArray *mem, bool forw_sens)
{
    mwSize nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    mwSize num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );

    int size = sizeof(sim_erk_workspace);
     
    size += nx * sizeof(double); // xt
    size += num_stages*nx * sizeof(double); // K
    
    if (forw_sens){
        size += num_stages*nx*nx * sizeof(double); // dKx
        size += num_stages*nx*nu * sizeof(double); // dKu
        size += nx*nx * sizeof(double); // jacX_t
        size += nx*nu * sizeof(double); // jacU_t
    }

    return size;
}

void *sim_erk_cast_workspace(const mxArray *mem, bool forw_sens, void *raw_memory){
    mwSize nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    mwSize num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    
    char *c_ptr = (char *)raw_memory;
    
    sim_erk_workspace *workspace = (sim_erk_workspace *) c_ptr;
    c_ptr += sizeof(sim_erk_workspace);
    
    workspace->xt = (double *)c_ptr;
    c_ptr += nx * sizeof(double);
    
    workspace->K = (double *)c_ptr;
    c_ptr += num_stages*nx * sizeof(double);
    
    if (forw_sens){
    
        workspace->dKx = (double *)c_ptr;
        c_ptr += num_stages*nx*nx * sizeof(double);

        workspace->dKu = (double *)c_ptr;
        c_ptr += num_stages*nx*nu * sizeof(double);

        workspace->jacX_t = (double *)c_ptr;
        c_ptr += nx*nx * sizeof(double);

        workspace->jacU_t = (double *)c_ptr;
        c_ptr += nx*nu * sizeof(double);   
    }
    
    assert((char*)raw_memory + sim_erk_calculate_workspace_size(mem, forw_sens) >= c_ptr);

    return (void *)workspace;

}

int sim_erk(double **in, double **out, double **Jac, const mxArray *mem, bool forw_sens, void *work_){    
    int istep, s, i, j;
    double a,b;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
    
    double *x = in[0];
    double *xn = out[0];
    
    double *A = mxGetPr( mxGetField(mem, 0, "A"));
    double *B = mxGetPr( mxGetField(mem, 0, "B"));
    double h = mxGetScalar( mxGetField(mem, 0, "h") );
    mwSize num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    mwSize num_steps = mxGetScalar( mxGetField(mem, 0, "num_steps") );
    mwSize nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    double *Sx = mxGetPr( mxGetField(mem, 0, "Sx"));
    double *Su = mxGetPr( mxGetField(mem, 0, "Su"));
    
    double *Jac_x = Jac[0];
    double *Jac_u = Jac[1];   
    
    sim_erk_workspace *workspace = (sim_erk_workspace *) sim_erk_cast_workspace(mem, forw_sens, work_);
    
    double *xt = workspace->xt;
    double *K = workspace->K;
    double *dKx = workspace->dKx;
    double *dKu = workspace->dKu;
    double *jacX_t = workspace->jacX_t;
    double *jacU_t = workspace->jacU_t;

    mwSize jx = nx*nx;
    mwSize ju = nx*nu;
    
    double *vde_in[5];
    double *vde_out[2];
    double *ode_in[3];
    double *ode_out[1];
     
    memcpy(&xn[0], &x[0], nx*sizeof(double));
    
    set_zeros(num_stages*nx, K);
    
    if (forw_sens){
        memcpy(&Jac_x[0], &Sx[0], nx*nx*sizeof(double));
        memcpy(&Jac_u[0], &Su[0], nx*nu*sizeof(double));
    }
      
    ode_in[1] = in[1];
    ode_in[2] = in[2];
    
    if (forw_sens){
        vde_in[1] = in[1];
        vde_in[2] = in[2];
    }
    
    for (istep = 0; istep < num_steps; istep++) {
        
        for (s = 0; s < num_stages; s++) {
            memcpy(&xt[0], &xn[0], nx*sizeof(double));
            
            if (forw_sens){
                memcpy(&jacX_t[0], &Jac_x[0], nx*nx*sizeof(double));
                memcpy(&jacU_t[0], &Jac_u[0], nx*nu*sizeof(double));
            }
            for (j = 0; j < s; j++){
                a = A[j * num_stages + s];
                if (a!=0){
                    a *= h;                   
                    daxpy(&nx, &a, K+j*nx, &one_i, xt, &one_i);
                    
                    if (forw_sens){
                        daxpy(&jx, &a, dKx+j*jx, &one_i, jacX_t, &one_i);
                        daxpy(&ju, &a, dKu+j*ju, &one_i, jacU_t, &one_i);
                    }
                                 
                }
            }
            
            ode_in[0] = xt;
            ode_out[0] = K+s*nx;
            f_Fun(ode_in,ode_out);         
                 
            if(forw_sens){
                vde_in[0] = xt;
                vde_in[3] = jacX_t;
                vde_in[4] = jacU_t;
                vde_out[0] = dKx+s*jx;
                vde_out[1] = dKu+s*ju;
                vde_Fun(vde_in, vde_out);
            }
                      
        }
        
        for (s = 0; s < num_stages; s++){
            b = h * B[s];
            daxpy(&nx, &b, K+s*nx, &one_i, xn, &one_i);
            
            if(forw_sens){
                daxpy(&jx, &b, dKx+s*jx, &one_i, Jac_x, &one_i);
                daxpy(&ju, &b, dKu+s*ju, &one_i, Jac_u, &one_i);
            }
        }
        
    }
    
    return 0;      
}

int sim_irk_calculate_workspace_size(const mxArray *mem, bool forw_sens)
{
    mwSize nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    mwSize num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    mwSize ldG = num_stages*nx;

    int size = sizeof(sim_irk_workspace);
    
    size += 4 * sizeof(double *); // impl_ode_in
    size += 4 * sizeof(double *); // impl_ode_out
     
    size += nx * sizeof(double); // xt
    size += 2 *ldG * sizeof(double); // K, rG
    size += ldG*ldG * sizeof(double); // JGK
    
    if (forw_sens){        
        size += ldG*nx * sizeof(double); // JGx
        size += ldG*nu * sizeof(double); // JGu
        size += ldG*nx * sizeof(double); // JKx
        size += ldG*nu * sizeof(double); // JKu
    }
      
    size += nx*nx*sizeof(double); //impl_ode_out[1]
    size += nx*nu*sizeof(double); //impl_ode_out[2]
    size += nx*nx*sizeof(double); //impl_ode_out[3]
    
    size += ldG * sizeof(mwIndex); // IPIV
 
    return size;
}

void *sim_irk_cast_workspace(const mxArray *mem, bool forw_sens, void *raw_memory){
    mwSize nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    mwSize num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    mwSize ldG = num_stages*nx;
    
    char *c_ptr = (char *)raw_memory;
    
    sim_irk_workspace *workspace = (sim_irk_workspace *) c_ptr;
    c_ptr += sizeof(sim_irk_workspace);
    
    workspace->impl_ode_in = (double **)c_ptr;
    c_ptr += 4 * sizeof(double *);
    
    workspace->impl_ode_out = (double **)c_ptr;
    c_ptr += 4 * sizeof(double *);
    
    workspace->xt = (double *)c_ptr;
    c_ptr += nx * sizeof(double);
    
    workspace->K = (double *)c_ptr;
    c_ptr += ldG * sizeof(double);
    
    workspace->rG = (double *)c_ptr;
    c_ptr += ldG * sizeof(double);
    
    workspace->JGK = (double *)c_ptr;
    c_ptr += ldG * ldG * sizeof(double);
    
    if (forw_sens){
        workspace->JGx = (double *)c_ptr;
        c_ptr += ldG*nx * sizeof(double);

        workspace->JGu = (double *)c_ptr;
        c_ptr += ldG*nu * sizeof(double);

        workspace->JKx = (double *)c_ptr;
        c_ptr += ldG*nx * sizeof(double);

        workspace->JKu = (double *)c_ptr;
        c_ptr += ldG*nu * sizeof(double);
    }
    
    workspace->impl_ode_out[1] = (double *)c_ptr;
    c_ptr += nx*nx * sizeof(double);
    
    workspace->impl_ode_out[2] = (double *)c_ptr;
    c_ptr += nx*nu * sizeof(double);
    
    workspace->impl_ode_out[3] = (double *)c_ptr;
    c_ptr += nx*nx * sizeof(double);
    
    workspace->IPIV = (mwIndex *)c_ptr;
    c_ptr += ldG * sizeof(mwIndex);
    
    assert((char*)raw_memory + sim_irk_calculate_workspace_size(mem, forw_sens) >= c_ptr);

    return (void *)workspace;

}

int sim_irk(double **in, double **out, double **Jac, const mxArray *mem, bool forw_sens, void *work_){    
    int istep, i, j, k, iter;
    double a,b;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
      
    double *x = in[0];
    double *xn = out[0];
   
    double *Jac_x = Jac[0];
    double *Jac_u = Jac[1];
    
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
    
    sim_irk_workspace *workspace = (sim_irk_workspace *) sim_irk_cast_workspace(mem, forw_sens, work_);
    
    double *xt = workspace->xt;
    double *K = workspace->K;
    double *rG = workspace->rG;
    double *JGK = workspace->JGK;
    double *JGx = workspace->JGx;
    double *JGu = workspace->JGu;
    double *JKx = workspace->JKx;
    double *JKu = workspace->JKu;
    double **impl_ode_in = workspace->impl_ode_in;
    double **impl_ode_out = workspace->impl_ode_out;
    mwIndex *IPIV = workspace->IPIV;
       
    mwSize jx = nx*nx;
    mwSize ju = nx*nu;
    mwSize ldG = num_stages*nx;
    
    mwIndex INFO = 0;
             
    memcpy(&xn[0], &x[0], nx*sizeof(double));
    
    if (forw_sens){
        memcpy(&Jac_x[0], &Sx[0], nx*nx*sizeof(double));
        memcpy(&Jac_u[0], &Su[0], nx*nu*sizeof(double));
    }
    
    impl_ode_in[1] = in[1]; // u
    impl_ode_in[2] = in[2]; // p
    
    set_zeros(ldG, K);
        
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
        if (forw_sens){
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
        }
       
        // update xn
        for (i = 0; i < num_stages; i++){
            b = h * B[i];
            daxpy(&nx, &b, K+i*nx, &one_i, xn, &one_i);
        }
        
    }
    
    return 0;
    
}


