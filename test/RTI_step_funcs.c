#include "mex.h"
#include "stdlib.h"
#include "string.h"

#include "RTI_step_common.h"
#include "RTI_step_funcs.h"
#include "casadi_wrapper.h"
#include "mpc_common.h"
#include "blas.h"

void qp_generation(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *Cu, double *CxN, 
        double *gx, double *gu, double *a, double *ds0, double *lc, double *uc, double *lb_du, double *ub_du, 
        double *x0, double *z, double *xN, double *y, double *yN, double *od, double *W, double *WN, double *lb, double *ub,
        double *lbN, double *ubN, double *lbu, double *ubu, int lin_obj,
        rti_step_dims *dim, rti_step_workspace *work)
{
    size_t nx = dim->nx;
    size_t nu = dim->nu;
    size_t np = dim->np;
    size_t ny = dim->ny;
    size_t nyN = dim->nyN;
    size_t nc = dim->nc;
    size_t ncN = dim->ncN;   
    size_t N = dim->N;
    
    double **Jac = work->Jac;
    double *Jac_N = work->Jac_N;
    
    size_t i=0,j=0;
    size_t nz = nx+nu;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - z[i];
    
    double *Sens[2];    
    double *Cons[2];
      
    double *vec_in[4];
    double *vec_out[2];    
    vec_in[3] = W;
    
    double *ode_in[3];
        
    for(i=0;i<N;i++){
        vec_in[0] = z+i*nz;
        vec_in[1] = od+i*np;
        vec_in[2] = y+i*ny;
        
        // control bounds
        for (j=0;j<nu;j++){
            lb_du[i*nu+j] = lbu[i*nu+j]-z[i*nz+nx+j];
            ub_du[i*nu+j] = ubu[i*nu+j]-z[i*nz+nx+j];
        }
        
        // integration      
        vec_out[0] = a+i*nx;
        Sens[0] = A + i*nx*nx;
        Sens[1] = B + i*nx*nu;
        
        F_Fun(vec_in, vec_out);
        D_Fun(vec_in, Sens);
        
        if (i < N-1){
            for (j=0;j<nx;j++)
                vec_out[0][j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                vec_out[0][j] -= xN[j];
        }

        
        // Hessian
        if (!lin_obj){
            Ji_Fun(vec_in, Jac);
            dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, Jac[0], &ny, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
            dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, Jac[0], &ny, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
            dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, Jac[1], &ny, Jac[1], &ny, &zero, R+i*nu*nu, &nu);          
        }
        
        // gradient
        vec_out[0] = gx+i*nx;
        vec_out[1] = gu+i*nu;
        gi_Fun(vec_in, vec_out);
        
        // constraint residual
        if (nc>0){  
            vec_in[0]=z+i*nz;
            vec_in[1]=z+i*nz+nx;
            vec_in[2]=od+i*np; 
            vec_out[0] = lc + i*nc;
            path_con_Fun(vec_in, vec_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - vec_out[0][j];
                vec_out[0][j] = lb[i*nc+j] - vec_out[0][j];            
            }
        
            // constraint Jacobian
            Cons[0] = Cx+i*nc*nx;
            Cons[1] = Cu+i*nc*nu;
            Ci_Fun(vec_in, Cons);
        }
    }
    
    // terminal data
    vec_in[0] = xN;
    vec_in[1] = od+N*np;
    vec_in[2] = yN;
    vec_in[3] = WN;
    
    if (!lin_obj){
        JN_Fun(vec_in, Jac_N);
        dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, Jac_N, &nyN, Jac_N, &nyN, &zero, Q+N*nx*nx, &nx);
    }
    
    vec_out[0] = gx+N*nx;
    gN_Fun(vec_in, vec_out);

    if (ncN>0){
        vec_out[0] = lc + N*nc;
        path_con_N_Fun(vec_in, vec_out);
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ubN[j] - vec_out[0][j];
            vec_out[0][j] = lbN[j] - vec_out[0][j];            
        }

        CN_Fun(vec_in, CxN);
    }
}

void condensing(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *Cu, double *CxN, 
        double *gx, double *gu, double *a, double *ds0, double *lc, double *uc,
        double *G, double *Hc, double *gc, double *Cc, double *lcc, double *ucc,
        int iter, bool cond_save,
        rti_step_dims *dim, rti_step_workspace *work)
{
    size_t nx = dim->nx;
    size_t nu = dim->nu;
    size_t np = dim->np;
    size_t ny = dim->ny;
    size_t nyN = dim->nyN;
    size_t nc = dim->nc;
    size_t ncN = dim->ncN;   
    size_t N = dim->N;
    
    double *L = work->L;
    double *w_vec = work->w_vec;
    double *w_vec_dup = work->w_vec_dup;
    double *W_mat = work->W_mat;
    double *W_mat_dup = work->W_mat_dup;
    double *Hi = work->Hi;
    double *Cci = work->Cci;
    double *CcN = work->CcN;
    
    size_t i=0,j=0;
    size_t nz = nx+nu;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
    
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
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, Q+N*nx*nx, &nx, G+(i*N+N-1)*nx*nu, &nx, &zero, W_mat, &nx);
            for(j=N-1;j>i;j--){        
                dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, S+j*nx*nu, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, Hi, &nu);                     
                dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+j*nx*nu, &nx, W_mat, &nx, &one_d, Hi, &nu);
                Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);

                Block_Fill_Trans(nu, nu, Hi, Hc, i*nu, j*nu, N*nu);

                dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, W_mat, &nx, &zero, W_mat_dup, &nx); 
                dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, Q+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &one_d, W_mat_dup, &nx); 
                memcpy(W_mat,W_mat_dup,nx*nu*sizeof(double));
            }
            memcpy(Hi,R+i*nu*nu,nu*nu*sizeof(double));
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+i*nx*nu, &nx, W_mat, &nx, &one_d, Hi, &nu);
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
                dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, CxN, &ncN, G+(i*N+N-1)*nx*nu, &nx, &zero, CcN, &ncN);
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
    memcpy(w_vec,gx+N*nx,nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,Q+N*nx*nx,&nx,L+N*nx,&one_i,&one_d,w_vec,&one_i);
    for(i=N-1;i>0;i--){
        memcpy(gc+i*nu,gu+i*nu,nu*sizeof(double));
        dgemv(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,L+i*nx,&one_i,&one_d,gc+i*nu,&one_i);
        dgemv(Trans,&nx,&nu,&one_d,B+i*nx*nu,&nx,w_vec,&one_i,&one_d,gc+i*nu,&one_i);
         
        memcpy(w_vec_dup,gx+i*nx,nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,L+i*nx,&one_i,&one_d,w_vec_dup,&one_i);
        dgemv(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,w_vec,&one_i,&one_d,w_vec_dup,&one_i);
        memcpy(w_vec,w_vec_dup,nx*sizeof(double));
    }   
    memcpy(gc,gu,nu*sizeof(double));
    dgemv(Trans,&nx,&nu,&one_d,S,&nx,L,&one_i,&one_d,gc,&one_i);
    dgemv(Trans,&nx,&nu,&one_d,B,&nx,w_vec,&one_i,&one_d,gc,&one_i);
   
    /* Compute cc */
    if (nc>0){                    
        for(i=0;i<N;i++){
            memcpy(lcc+i*nc,lc+i*nc,nc*sizeof(double));
            dgemv(nTrans,&nc,&nx,&minus_one_d,Cx+i*nc*nx,&nc,L+i*nx,&one_i,&one_d,lcc+i*nc,&one_i); 

            memcpy(ucc+i*nc,uc+i*nc,nc*sizeof(double));
            dgemv(nTrans,&nc,&nx,&minus_one_d,Cx+i*nc*nx,&nc,L+i*nx,&one_i,&one_d,ucc+i*nc,&one_i);
        }        
    }   
    
    /* Compute ccN */
    if (ncN>0){          
        memcpy(lcc+N*nc,lc+N*nc,ncN*sizeof(double));
        dgemv(nTrans,&ncN,&nx,&minus_one_d,CxN,&ncN,L+N*nx,&one_i,&one_d,lcc+N*nc,&one_i);
        memcpy(ucc+N*nc,uc+N*nc,ncN*sizeof(double));
        dgemv(nTrans,&ncN,&nx,&minus_one_d,CxN,&ncN,L+N*nx,&one_i,&one_d,ucc+N*nc,&one_i);
    }
}

void recover(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *CxN,
        double *gx, double *a, double *ds0, double *du, double *multipliers,
        double *dz, double *dxN, double *lambda_new, double *mu_new, double *muN_new, double *mu_u_new,
        rti_step_dims *dim)
{
    size_t nx = dim->nx;
    size_t nu = dim->nu;
    size_t nc = dim->nc;
    size_t ncN = dim->ncN;   
    size_t N = dim->N;
    
    size_t nz = nx+nu;
    
    memcpy(mu_u_new, multipliers, N*nu*sizeof(double));
    memcpy(mu_new, multipliers+N*nu, N*nc*sizeof(double));
    memcpy(muN_new, multipliers+N*nu+N*nc, ncN*sizeof(double));
    
    int i;        
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
    
    memcpy(dz, ds0, nx*sizeof(double)); 
    
    for (i=0;i<N;i++){        
        if (i<N-1){
            memcpy(dz+(i+1)*nz, a+i*nx, nx*sizeof(double));          
            dgemv(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,dz+i*nz,&one_i,&one_d,dz+(i+1)*nz,&one_i);
            dgemv(nTrans,&nx,&nu,&one_d,B+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,dz+(i+1)*nz,&one_i);
        }else{
            memcpy(dxN, a+i*nx, nx*sizeof(double));
            dgemv(nTrans,&nx,&nx,&one_d,A+i*nx*nx,&nx,dz+i*nz,&one_i,&one_d,dxN,&one_i);
            dgemv(nTrans,&nx,&nu,&one_d,B+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,dxN,&one_i);
        }
                
        memcpy(dz+i*nz+nx, du+i*nu, nu*sizeof(double));
    }
    
    memcpy(lambda_new+N*nx, gx+N*nx, nx*sizeof(double));
    dgemv(nTrans,&nx,&nx,&one_d,Q+N*nx*nx,&nx,dxN,&one_i,&one_d,lambda_new+N*nx,&one_i);
    
    if (ncN>0){
        dgemv(Trans,&ncN,&nx,&one_d,CxN,&ncN,muN_new,&one_i,&one_d,lambda_new+N*nx,&one_i);
    }
    for (i=N-1;i>-1;i--){
        memcpy(lambda_new+i*nx,gx+i*nx, nx*sizeof(double));
        dgemv(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,dz+i*nz,&one_i,&one_d,lambda_new+i*nx,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,S+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,lambda_new+i*nx,&one_i);
        dgemv(Trans,&nx,&nx,&one_d,A+i*nx*nx,&nx,lambda_new+(i+1)*nx,&one_i,&one_d,lambda_new+i*nx,&one_i);
        
        if (nc>0)
            dgemv(Trans,&nc,&nx,&one_d,Cx+i*nc*nx,&nc,mu_new+i*nc,&one_i,&one_d,lambda_new+i*nx,&one_i);       
    }
}

void line_search(double *dz, double *dxN, double *lambda_new, double *mu_new, double *muN_new, 
        double *z, double *xN, double *lambda, double *mu, double *muN,
        rti_step_dims *dim)
{
    size_t nx = dim->nx;
    size_t nu = dim->nu;    
    size_t nc = dim->nc;
    size_t ncN = dim->ncN;   
    size_t N = dim->N;
       
    int i;
    for (i=0;i<N*(nx+nu);i++)
        z[i] += dz[i];
    for (i=0;i<nx;i++)
        xN[i] += dxN[i];
    
    memcpy(lambda, lambda_new, (N+1)*nx*sizeof(double));        
    if (nc>0)        
        memcpy(mu, mu_new, N*nc*sizeof(double));        
    if (ncN>0)
        memcpy(muN, muN_new, ncN*sizeof(double));
    
}
