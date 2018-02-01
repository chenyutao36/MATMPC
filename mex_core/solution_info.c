#include "mex.h"
#include "string.h"
#include <stdbool.h>

#include "casadi_wrapper.h"
#include "sim.h"

#include "lapack.h"
#include "blas.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

static double *vec_out[3];
static double *L = NULL;
static double *eq_res_vec = NULL;
static double *lu, *uu;
static void *workspace = NULL;
static bool mem_alloc_info = false;

void exitFcn_info(){
    if (mem_alloc_info){
        mxFree(vec_out[1]);
        mxFree(vec_out[2]);
        mxFree(L);
        mxFree(eq_res_vec);
        mxFree(lu);
        mxFree(uu);
        if (workspace!=NULL)
            mxFree(workspace);
    }
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{    
    double *z = mxGetPr( mxGetField(prhs[0], 0, "z") );
    double *xN = mxGetPr( mxGetField(prhs[0], 0, "xN") );
    double *lambda = mxGetPr( mxGetField(prhs[0], 0, "lambda") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu") );
    double *muN = mxGetPr( mxGetField(prhs[0], 0, "muN") );
    double *mu_u = mxGetPr( mxGetField(prhs[0], 0, "mu_u") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *QN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub") );
    double *lbN = mxGetPr( mxGetField(prhs[0], 0, "lbN") );
    double *ubN = mxGetPr( mxGetField(prhs[0], 0, "ubN") );
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ubu") );
    
    double *ds0 = mxGetPr( mxGetField(prhs[2], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[2], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[2], 0, "uc") );
     
    mwSize nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    mwSize np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    mwSize nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    mwSize nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
    
    sim_opts opts;
    opts.forw_sens = false;
      
    mwIndex i=0,j=0;
    mwSize nz = nx+nu;
    mwSize nw = N*nz+nx;
    mwSize neq = (N+1)*nx;
    mwSize nineq = N*nc+ncN;
    char *nTrans = "N", *Trans="T", *Norm="O";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
    
    double *vec_in[6];
    double *ode_in[3];
    double *Sens[2];
    
    if (!mem_alloc_info){
        vec_out[1] = (double *)mxMalloc(nz * sizeof(double));
        mexMakeMemoryPersistent(vec_out[1]);
        vec_out[2] = (double *)mxMalloc(nz * sizeof(double));
        mexMakeMemoryPersistent(vec_out[2]);     
        L = (double *)mxMalloc( nw * sizeof(double));
        mexMakeMemoryPersistent(L);
        
        eq_res_vec = (double *)mxMalloc( neq * sizeof(double));
        mexMakeMemoryPersistent(eq_res_vec);
        
        lu = (double *)mxMalloc( N*nu * sizeof(double));
        mexMakeMemoryPersistent(lu);
        uu = (double *)mxMalloc( N*nu * sizeof(double));
        mexMakeMemoryPersistent(uu);
        
        int size;
        switch(sim_method){
            case 0:
                size = 0;
                break;
            case 1:
                size = sim_erk_calculate_workspace_size(prhs[2],&opts);
                break;
            case 2:
                size = sim_irk_calculate_workspace_size(prhs[2],&opts);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
         
        }   
        
        if (size > 0){
            workspace = mxMalloc(size);
            mexMakeMemoryPersistent(workspace);  
        }
        
        mem_alloc_info = true;     
        mexAtExit(exitFcn_info);
    }
    
    vec_in[3] = Q;
    memcpy(eq_res_vec, ds0, nx*sizeof(double));
    
    double *work;
    double KKT=0, eq_res=0, ineq_res=0;
    
    for (i=0;i<N;i++){
        vec_in[0] = z+i*nz;
        vec_in[1] = od+i*np;
        vec_in[2] = y+i*ny;
        vec_in[4] = lambda+(i+1)*nx;
        vec_in[5] = mu+i*nc;
           
        vec_out[0] = L+i*nz;
        adj_Fun(vec_in, vec_out);
        
        for (j=0;j<nx;j++){
            vec_out[1][j] -= lambda[i*nx+j];
        }
        
        daxpy(&nz, &one_d, vec_out[1], &one_i, vec_out[0], &one_i);
        daxpy(&nu, &one_d, mu_u+i*nu, &one_i, vec_out[0]+nx, &one_i);
        if (nc>0)
            daxpy(&nz, &one_d, vec_out[2], &one_i, vec_out[0], &one_i);
        
        vec_out[0] = eq_res_vec+(i+1)*nx;        
        switch(sim_method){
            case 0:
                F_Fun(vec_in, vec_out);
                break;
            case 1:               
                ode_in[0]=z+i*nz;
                ode_in[1]=z+i*nz+nx;
                ode_in[2]=od+i*np;
                sim_erk(ode_in, vec_out, Sens, prhs[2], &opts, workspace);
                break;
            case 2:                         
                ode_in[0]=z+i*nz;
                ode_in[1]=z+i*nz+nx;
                ode_in[2]=od+i*np;
                sim_irk(ode_in, vec_out, Sens, prhs[2], &opts, workspace);
                break;
            default :
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
        
        if (i < N-1){
            for (j=0;j<nx;j++)
                vec_out[0][j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                vec_out[0][j] -= xN[j];
        }
        
        for (j=0;j<nu;j++){
            lu[i*nu+j] = lbu[i*nu+j] - z[i*nz+nx+j];
            uu[i*nu+j] = ubu[i*nu+j] - z[i*nz+nx+j];
        }
        
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
        }
    }
    vec_in[0] = xN;
    vec_in[1] = od+N*np;
    vec_in[2] = yN;
    vec_in[3] = QN;
    vec_in[4] = muN;
    
    vec_out[0] = L+N*nz;
    adjN_Fun(vec_in, vec_out);
    
    if (ncN>0)
        daxpy(&nx, &one_d, vec_out[1], &one_i, vec_out[0], &one_i);
    
    if (ncN>0){
        vec_out[0] = lc + N*nc;
        path_con_N_Fun(vec_in, vec_out); 
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ubN[j] - vec_out[0][j];
            vec_out[0][j] = lbN[j] - vec_out[0][j];            
        }
    }
         
    eq_res = dlange(Norm, &neq, &one_i, eq_res_vec, &one_i, work);
    KKT = dlange(Norm, &nw, &one_i, L, &one_i, work);
    
    for (i=0;i<N*nu;i++)
        ineq_res += MIN(uu[i],0) + MAX(lu[i],0);  
    for (i=0;i<nineq;i++)
        ineq_res += MIN(uc[i],0) + MAX(lc[i],0);
    
    plhs[0] = mxCreateDoubleScalar(eq_res); // eq_res
    plhs[1] = mxCreateDoubleScalar(ineq_res); // ineq_res
    plhs[2] = mxCreateDoubleScalar(KKT); // KKT
    
}