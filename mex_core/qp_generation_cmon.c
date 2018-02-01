

#include "mex.h"
#include "string.h"

#include "casadi_wrapper.h"
#include "mpc_common.h"

#include "blas.h"

static double *Jac[2];
static double *Jac_N;
static double *num_pri;
static double *den_pri;
static double *num_dual;
static double *den_dual;

static bool mem_alloc = false;

void exitFcn_sim(){
    if (mem_alloc){
        mxFree(Jac[0]);
        mxFree(Jac[1]);
        mxFree(Jac_N);
        mxFree(num_pri);
        mxFree(den_pri);
        mxFree(num_dual);
        mxFree(den_dual);
    }
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *z = mxGetPr( mxGetField(prhs[0], 0, "z") );
    double *xN = mxGetPr( mxGetField(prhs[0], 0, "xN") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *QN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub") );
    double *lbN = mxGetPr( mxGetField(prhs[0], 0, "lbN") );
    double *ubN = mxGetPr( mxGetField(prhs[0], 0, "ubN") );
    double *x0 = mxGetPr( mxGetField(prhs[0], 0, "x0") );
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ubu") );
    double *lambda = mxGetPr( mxGetField(prhs[0], 0, "lambda") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu") );
    double *muN = mxGetPr( mxGetField(prhs[0], 0, "muN") );
    double *mu_u = mxGetPr( mxGetField(prhs[0], 0, "mu_u") );
    
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    mwSize nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
    
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    mwSize nw = N*nz+nx;
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    mwSignedIndex one_i = 1;
      
    double *Qh = mxGetPr( mxGetField(prhs[2], 0, "Q_h") );
    double *S = mxGetPr( mxGetField(prhs[2], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[2], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[2], 0, "A_sens") );
    double *B = mxGetPr( mxGetField(prhs[2], 0, "B_sens") );
    double *Cx = mxGetPr( mxGetField(prhs[2], 0, "Cx") );
    double *Cu = mxGetPr( mxGetField(prhs[2], 0, "Cu") );
    double *gx = mxGetPr( mxGetField(prhs[2], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[2], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[2], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[2], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[2], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[2], 0, "uc") );
    double *lb_du = mxGetPr( mxGetField(prhs[2], 0, "lb_du") );
    double *ub_du = mxGetPr( mxGetField(prhs[2], 0, "ub_du") );
    double *CxN = mxGetPr( mxGetField(prhs[2], 0, "CxN") );
    
    double *F_old = mxGetPr( mxGetField(prhs[2], 0, "F_old") );
    double *CMON_pri = mxGetPr( mxGetField(prhs[2], 0, "CMON_pri") );
    double *CMON_dual = mxGetPr( mxGetField(prhs[2], 0, "CMON_dual") );
    double *q_pri = mxGetPr( mxGetField(prhs[2], 0, "dz") );
    double *q_dual = mxGetPr( mxGetField(prhs[2], 0, "q_dual") );
    double threshold_pri = mxGetScalar( mxGetField(prhs[2], 0, "threshold_pri") );
    double threshold_dual = mxGetScalar( mxGetField(prhs[2], 0, "threshold_dual") );
    double *perc = mxGetPr( mxGetField(prhs[2], 0, "perc") );
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - z[i];
    
    // allocate memory
    double *Sens[2];    
    double *Cons[2];
      
//     double *vec_in[4];
    double *vec_in[6];    
//     double *vec_out[2];
    double *vec_out[3];
    vec_in[3] = Q;
    double *tmp = (double *)mxMalloc(nz * sizeof(double));
    vec_out[2] = (double *)mxMalloc(nz * sizeof(double));   
    double *ode_in[3];

    if (!mem_alloc){
        
        Jac[0] = (double *) mxCalloc(ny*nx, sizeof(double));
        mexMakeMemoryPersistent(Jac[0]); 
        Jac[1] = (double *) mxCalloc(ny*nu, sizeof(double));
        mexMakeMemoryPersistent(Jac[1]); 
        Jac_N = (double *) mxCalloc(nyN*nx, sizeof(double));
        mexMakeMemoryPersistent(Jac_N);
        
        num_pri = (double *) mxCalloc(nx, sizeof(double));
        mexMakeMemoryPersistent(num_pri);
        den_pri = (double *) mxCalloc(nx, sizeof(double));
        mexMakeMemoryPersistent(den_pri);
        
        num_dual = (double *) mxCalloc(nz, sizeof(double));
        mexMakeMemoryPersistent(num_dual);
        den_dual = (double *) mxCalloc(nz, sizeof(double));
        mexMakeMemoryPersistent(den_dual);
        
        mem_alloc=true;
        mexAtExit(exitFcn_sim);
    }
    
    double num_norm_pri, den_norm_pri, num_norm_dual, den_norm_dual;
    int num_updated = 0;
    
    // start loop
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
        F_Fun(vec_in, vec_out);
        
        // sensitivity
        Sens[0] = A + i*nx*nx;
        Sens[1] = B + i*nx*nu;        
        // primal CMON              
        dgemv(nTrans, &nx, &nx, &one_d, Sens[0], &nx, q_pri+i*nz, &one_i, &zero, den_pri, &one_i);            
        dgemv(nTrans, &nx, &nu, &one_d, Sens[1], &nx, q_pri+i*nz+nx, &one_i, &one_d, den_pri, &one_i);            
        den_norm_pri = dnrm2(&nx, den_pri, &one_i);                                  
        if (den_norm_pri>0){           
            memcpy(&num_pri[0], &a[i*nx], nx*sizeof(double));
            daxpy(&nx, &minus_one_d, F_old+i*nx, &one_i, num_pri, &one_i);
//             dgemv(nTrans, &nx, &nx, &minus_one_d, A+i*nx*nx, &nx, q+i*nz, &one_i, &one_d, num_pri, &one_i);
//             dgemv(nTrans, &nx, &nu, &minus_one_d, B+i*nx*nu, &nx, q+i*nz+nx, &one_i, &one_d, num_pri, &one_i);
            daxpy(&nx, &minus_one_d, den_pri, &one_i, num_pri, &one_i);
            num_norm_pri=dnrm2(&nx, num_pri, &one_i);                
            CMON_pri[i] = num_norm_pri/den_norm_pri;
        }else {
            CMON_pri[i] = 10000000;
        }
          
        // dual CMON
        dgemv(Trans, &nx, &nx, &one_d, Sens[0], &nx, q_dual+(i+1)*nx, &one_i, &zero, den_dual, &one_i);
        dgemv(Trans, &nx, &nu, &one_d, Sens[1], &nx, q_dual+(i+1)*nx, &one_i, &zero, den_dual+nx, &one_i);
        den_norm_dual = dnrm2(&nz, den_dual, &one_i);
        if (den_norm_dual>0){   
            vec_in[4] = q_dual+(i+1)*nx;
            vec_in[5] = mu+i*nc;
            vec_out[0] = tmp;
            vec_out[1] = num_dual;
            adj_Fun(vec_in, vec_out);
            daxpy(&nz, &minus_one_d, den_dual, &one_i, num_dual, &one_i);
            num_norm_dual=dnrm2(&nz, num_dual, &one_i); 
            CMON_dual[i] = num_norm_dual/den_norm_dual;
        }else {
            CMON_dual[i] = 10000000;
        }
        
        if (CMON_pri[i]>threshold_pri || CMON_dual[i]>threshold_dual){
            num_updated++;
            D_Fun(vec_in, Sens);
//             memcpy(&F_old[i*nx], &a[i*nx], nx*sizeof(double));
//             set_zeros(nz,q+i*nz);
        }
        memcpy(&F_old[i*nx], &a[i*nx], nx*sizeof(double));
        
        // integration residual
        if (i < N-1){
            for (j=0;j<nx;j++)
                a[i*nx+j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                a[i*nx+j] -= xN[j];
        }
               
        // Hessian
        Ji_Fun(vec_in, Jac);
        dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, Jac[0], &ny, Jac[0], &ny, &zero, Qh+i*nx*nx, &nx);
        dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, Jac[0], &ny, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
        dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, Jac[1], &ny, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
        
        // gradient
        vec_out[0] = gx+i*nx;
        vec_out[1] = gu+i*nu;
        gi_Fun(vec_in, vec_out);
        
//         vec_in[4] = lambda+(i+1)*nx;
//         vec_in[5] = mu+i*nc;         
//         vec_out[0] = dL+i*nz;
//         adj_Fun(vec_in, vec_out);      
//         for (j=0;j<nx;j++)
//             gx[i*nx+j] =vec_out[0][j] + vec_out[1][j] - lambda[i*nx+j];        
//         for (j=0;j<nu;j++)
//             gu[i*nu+j] =vec_out[0][nx+j] + vec_out[1][nx+j] + mu_u[i*nu+j];       
//         if (nc>0){
//             for (j=0;j<nx;j++)
//                 gx[i*nx+j] += vec_out[2][j];
//             for (j=0;j<nu;j++)
//                 gu[i*nu+j] += vec_out[2][nx+j];
//         }
               
        // constraint residual
        if (nc>0){        
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
    vec_in[3] = QN;
    
    JN_Fun(vec_in, Jac_N);
    dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, Jac_N, &nyN, Jac_N, &nyN, &zero, Qh+N*nx*nx, &nx);
    
    vec_out[0] = gx+N*nx;
    gN_Fun(vec_in, vec_out);
    
//     vec_in[4] = muN;    
//     vec_out[0] = gx+N*nx;
//     adjN_Fun(vec_in, vec_out);
//     if (ncN>0){
//         for (j=0;j<nx;j++)
//             gx[N*nx+j] += vec_out[1][j];
//     }

    if (ncN>0){
        vec_out[0] = lc + N*nc;
        path_con_N_Fun(vec_in, vec_out);
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ubN[j] - vec_out[0][j];
            vec_out[0][j] = lbN[j] - vec_out[0][j];            
        }

        CN_Fun(vec_in, CxN);
    }
    
    perc[0] = 100*num_updated/N;
    
    mxFree(tmp);
    mxFree(vec_out[2]);
    
}