#include "mex.h"
#include "string.h"

#include "sim.h"
#include "erk.h"
#include "irk.h"
#include "casadi_wrapper.h"

#include "blas.h"
#include "lapack.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

static double *eq_res_vec = NULL;
static double *z_new = NULL, *xN_new = NULL;

static sim_opts *opts = NULL;
static sim_in *in = NULL;
static sim_out *out = NULL;
static sim_erk_workspace *erk_workspace = NULL;
static sim_irk_workspace *irk_workspace = NULL;
static bool mem_alloc = false;

void exitFcn(){
    if (mem_alloc){
        mxFree(eq_res_vec);
        mxFree(z_new);
        mxFree(xN_new);
    }
    if (erk_workspace!=NULL)
        sim_erk_workspace_free(opts, erk_workspace);
    if (irk_workspace!=NULL)
        sim_irk_workspace_free(opts, irk_workspace);
    if (opts!=NULL)
        sim_opts_free(opts);
    if (in!=NULL)
        sim_in_free(in);
    if (out!=NULL)
        sim_out_free(out);    
}

double eval_cons_res(double *z, double *xN, double *od, double *ds0, double *lb, double *ub, double *lc, double *uc,
                   double *lbN, double *ubN, double *lbu, double *ubu, size_t nx, size_t nu, size_t nc, size_t ncN,
                   size_t N, size_t np, double *eq_res_vec, int sim_method, sim_opts *opts, sim_in *in, sim_out *out,
                   sim_erk_workspace *erk_workspace, sim_irk_workspace *irk_workspace)
{
    mwIndex i=0,j=0;
    
    mwSize nz = nx+nu;
    mwSize nw = N*nz;
    mwSize neq = (N+1)*nx;
    mwSize nineq = N*nc+ncN;
    
    double *casadi_in[3];
    double *casadi_out[1];
    double *work;
    double eq_res=0, ineq_res=0, cons_res;
    
    char *nTrans = "N", *Trans="T", *Norm="O";
    double one_d = 1.0, zero = 0.0;
    mwSignedIndex one_i = 1;
    
    double *lu = (double *)mxMalloc( N*nu * sizeof(double));        
    double *uu = (double *)mxMalloc( N*nu * sizeof(double));
    
    memcpy(eq_res_vec, ds0, nx*sizeof(double));
           
    for (i=0;i<N;i++){      
        switch(sim_method){
            case 0:
                casadi_in[0]=z+i*nz;
                casadi_out[0] = eq_res_vec+(i+1)*nx;
                F_Fun(casadi_in, casadi_out);
                break;
            case 1:               
                in->x = z+i*nz;
                in->u = z+i*nz+nx;
                in->p = od+i*np;
                out->xn = eq_res_vec+(i+1)*nx;
                sim_erk(in, out, opts, erk_workspace);
                break;
            case 2:  
                in->x = z+i*nz;
                in->u = z+i*nz+nx;
                in->p = od+i*np;
                out->xn = eq_res_vec+(i+1)*nx;
                sim_irk(in, out, opts, irk_workspace);
                break;
            default :
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
        
        if (i < N-1){
            for (j=0;j<nx;j++)
                eq_res_vec[(i+1)*nx+j] -= z[(i+1)*nz+j];
        }else{
            for (j=0;j<nx;j++)
                eq_res_vec[(i+1)*nx+j] -= xN[j];
        }
        
        for (j=0;j<nu;j++){
            lu[i*nu+j] = lbu[i*nu+j] - z[i*nz+nx+j];
            uu[i*nu+j] = ubu[i*nu+j] - z[i*nz+nx+j];
        }
            
        if (nc>0){
            casadi_in[0]=z+i*nz;
            casadi_in[1]=z+i*nz+nx;
            casadi_in[2]=od+i*np;
            casadi_out[0] = lc + i*nc;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - casadi_out[0][j];
                casadi_out[0][j] = lb[i*nc+j] - casadi_out[0][j];            
            }
        }          
    }
        
    if (ncN>0){
        casadi_in[0] = xN;
        casadi_in[1] = od+N*np;
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out);
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ubN[j] - casadi_out[0][j];
            casadi_out[0][j] = lbN[j] - casadi_out[0][j];            
        }
    }
        
    eq_res = dlange(Norm, &neq, &one_i, eq_res_vec, &one_i, work);
            
    for (i=0;i<N*nu;i++)
        ineq_res += MIN(uu[i],0) + MAX(lu[i],0);  
    for (i=0;i<nineq;i++)
        ineq_res += MIN(uc[i],0) + MAX(lc[i],0);
        
    cons_res = eq_res + ineq_res;
    
    mxFree(lu);
    mxFree(uu);

    return cons_res;
}

double eval_curv(double *Q, double *S, double *R, double *dz, double *dxN,
                 size_t nx, size_t nu, size_t N)
{
    mwIndex i;
    mwSize nz = nx+nu;
    double curv = 0.0;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    mwSignedIndex one_i = 1; 
    
    double *tmp = (double *)mxMalloc(nz*sizeof(double));
    
    for(i=0; i<N; i++)
    {
        dgemv(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,dz+i*nz,&one_i,&zero,tmp,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,S+i*nx*nu,&nx,dz+i*nz+nx,&one_i,&one_d,tmp,&one_i);
        curv += ddot(&nx, dz+i*nz, &one_i, tmp, &one_i);
        
        dgemv(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,dz+i*nz,&one_i,&zero,tmp+nx,&one_i);
        dgemv(nTrans,&nu,&nu,&one_d,R+i*nu*nu,&nu,dz+i*nz+nx,&one_i,&one_d,tmp+nx,&one_i);
        curv += ddot(&nu, dz+i*nz+nx, &one_i, tmp+nx, &one_i);
    }
    
    dgemv(nTrans,&nx,&nx,&one_d,Q+N*nx*nx,&nx,dxN,&one_i,&zero,tmp,&one_i);
    curv += ddot(&nx, dxN, &one_i, tmp, &one_i);
    
    mxFree(tmp);
    
    return curv;
}

double eval_grad(double *gx, double *gu, double *dz, double *dxN,
                 size_t nx, size_t nu, size_t N)
{
    mwIndex i;
    mwSize nz = nx+nu;
    double grad = 0.0;   
    
    mwSignedIndex one_i = 1; 
    
    for(i=0; i<N; i++)
    {
        grad += ddot(&nx, dz+i*nz, &one_i, gx+i*nx, &one_i);
        grad += ddot(&nu, dz+i*nz+nx, &one_i, gu+i*nu, &one_i);
    }
    
    grad += ddot(&nx, dxN, &one_i, gx+N*nx, &one_i);
 
    return grad;
}

double eval_obj(double *z, double *xN, double *od, double *y, double *yN, double *W, double *WN,
                size_t nx, size_t nu, size_t np, size_t ny, size_t N)
{
    mwIndex i;
    mwSize nz = nx+nu;
    double *in[4];
    double *out[1];
    double obj=0.0;
    
    in[3] = W;
    out[0] = (double *) mxMalloc(sizeof(double));
    for (i=0; i<N; i++){
        in[0] = z+i*nz;
        in[1] = od+i*np;
        in[2] = y+i*ny;
        
        obji_Fun(in, out);
        obj += *out[0];
    }
    in[0] = xN;
    in[1] = od+N*np;
    in[2] = yN;
    in[3] = WN;
    
    objN_Fun(in, out);
    obj += *out[0];
    
    mxFree(out[0]);
    
    return obj;
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *z = mxGetPr( mxGetField(prhs[1], 0, "z") );
    double *xN = mxGetPr( mxGetField(prhs[1], 0, "xN") );
    double *lambda = mxGetPr( mxGetField(prhs[1], 0, "lambda") );
    double *mu = mxGetPr( mxGetField(prhs[1], 0, "mu") );
    double *muN = mxGetPr( mxGetField(prhs[1], 0, "muN") );
    double *mu_u = mxGetPr( mxGetField(prhs[1], 0, "mu_u") );
    double *od = mxGetPr( mxGetField(prhs[1], 0, "od") );
    double *lb = mxGetPr( mxGetField(prhs[1], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[1], 0, "ub") );
    double *lbN = mxGetPr( mxGetField(prhs[1], 0, "lbN") );
    double *ubN = mxGetPr( mxGetField(prhs[1], 0, "ubN") );
    double *y = mxGetPr( mxGetField(prhs[1], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[1], 0, "yN") );
    double *W = mxGetPr( mxGetField(prhs[1], 0, "W") );
    double *WN = mxGetPr( mxGetField(prhs[1], 0, "WN") );
    double *lbu = mxGetPr( mxGetField(prhs[1], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[1], 0, "ubu") );   
    
    mwSize nx = mxGetScalar( mxGetField(prhs[2], 0, "nx") );
    mwSize nu = mxGetScalar( mxGetField(prhs[2], 0, "nu") );
    mwSize nc = mxGetScalar( mxGetField(prhs[2], 0, "nc") );
    mwSize ncN = mxGetScalar( mxGetField(prhs[2], 0, "ncN") );
    mwSize N = mxGetScalar( mxGetField(prhs[2], 0, "N") );    
    mwSize np = mxGetScalar( mxGetField(prhs[2], 0, "np") ); if(np==0) np++;
    mwSize ny = mxGetScalar( mxGetField(prhs[2], 0, "ny") );
            
    mwSize nz = nx+nu;
    mwSize nw = N*nz;
    mwSize neq = (N+1)*nx;
    mwSize nineq = N*nc;
    mwSize nbu = N*nu;
    
    double one_d = 1.0;
    mwSignedIndex one_i = 1;
    
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[0], 0, "gu") );    
    double *dz = mxGetPr( mxGetField(prhs[0], 0, "dz") );
    double *dxN = mxGetPr( mxGetField(prhs[0], 0, "dxN") );
    double *lambda_new = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu_new = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *muN_new = mxGetPr( mxGetField(prhs[0], 0, "muN_new") );
    double *mu_u_new = mxGetPr( mxGetField(prhs[0], 0, "mu_u_new") );
    double *lc = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[0], 0, "uc") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
          
    double rho = mxGetScalar( mxGetField(prhs[0], 0, "rho") );
    double eta = mxGetScalar( mxGetField(prhs[0], 0, "eta") );
    double tau = mxGetScalar( mxGetField(prhs[0], 0, "tau") );
    double mu_safty = mxGetScalar( mxGetField(prhs[0], 0, "mu_safty") );
    double *mu_merit = mxGetPr( mxGetField(prhs[0], 0, "mu_merit") );
    
    mwIndex i=0,j=0;
    
    int sim_method = mxGetScalar( mxGetField(prhs[0], 0, "sim_method") );
    mwSize sqp_maxit = mxGetScalar( mxGetField(prhs[0], 0, "sqp_maxit") );   
       
    if (!mem_alloc){       
        eq_res_vec = (double *)mxMalloc( neq * sizeof(double));
        mexMakeMemoryPersistent(eq_res_vec);
        
        z_new = (double *)mxMalloc( nw * sizeof(double));
        mexMakeMemoryPersistent(z_new);
        xN_new = (double *)mxMalloc( nx * sizeof(double));
        mexMakeMemoryPersistent(xN_new);
              
        switch(sim_method){
            case 0:
                break;
            case 1:
                opts = sim_opts_create(prhs[0]);
                opts->forw_sens = false;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                erk_workspace = sim_erk_workspace_create(opts);               
                sim_erk_workspace_init(opts, prhs[0], erk_workspace);  
                break;
            case 2:
                opts = sim_opts_create(prhs[0]);
                opts->forw_sens = false;
                in = sim_in_create(opts);              
                out = sim_out_create(opts);                
                irk_workspace = sim_irk_workspace_create(opts);               
                sim_irk_workspace_init(opts, prhs[0], irk_workspace);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
         
        } 
                
        mem_alloc = true;     
        mexAtExit(exitFcn);
    }
    
    
    // backtracking
    double cons_res;
    double sigma = 0.0, pd, grad, mu_lb=0, obj, obj_new, dir_grad;
    double alpha = 1.0;
    bool newpoint = false;
    if (sqp_maxit > 1 ){
        cons_res = eval_cons_res(z, xN, od, ds0, lb, ub, lc, uc,
                                 lbN, ubN, lbu, ubu, nx, nu, nc, ncN,
                                 N, np, eq_res_vec, sim_method, opts, in, out,
                                 erk_workspace, irk_workspace);
              
        pd = eval_curv(Q, S, R, dz, dxN, nx, nu, N);
        grad = eval_grad(gx, gu, dz, dxN, nx, nu, N);
                
        if (pd>0)
            sigma = 0.5;
        
        if (cons_res>0)
            mu_lb=(grad+sigma*pd)/(1-rho)/cons_res;
                
        if (mu_merit[0]<mu_lb)
            mu_merit[0]=mu_lb*mu_safty;
        
        obj = eval_obj(z, xN, od, y, yN, W, WN,
                       nx, nu, np, ny, N) + mu_merit[0]*cons_res;
                
        dir_grad = grad - mu_merit[0] * cons_res;
                        
        while (!newpoint && alpha > 1E-4){
            memcpy(z_new, z, nw*sizeof(double));
            memcpy(xN_new, xN, nx*sizeof(double));
            
            daxpy(&nw, &alpha, dz, &one_i, z_new, &one_i); 
            daxpy(&nx, &alpha, dxN, &one_i, xN_new, &one_i);
            cons_res = eval_cons_res(z_new, xN_new, od, ds0, lb, ub, lc, uc,
                                     lbN, ubN, lbu, ubu, nx, nu, nc, ncN,
                                     N, np, eq_res_vec, sim_method, opts, in, out,
                                     erk_workspace, irk_workspace);
                        
            obj_new = eval_obj(z_new, xN_new, od, y, yN, W, WN,
                               nx, nu, np, ny, N) + mu_merit[0]*cons_res;
                        
            if (obj_new <= obj + eta*alpha*dir_grad)
                newpoint = true;
            else
                alpha *= tau;
        }
        
    }
               
    // update
    double inc = 1.0 - alpha;
     
    daxpy(&nw, &alpha, dz, &one_i, z, &one_i); 
    daxpy(&nx, &alpha, dxN, &one_i, xN, &one_i);
    
    for (i=0;i<neq;i++)
        lambda[i] *= inc;
    daxpy(&neq, &alpha, lambda_new, &one_i, lambda, &one_i);
    
    for (i=0;i<nbu;i++)
        mu_u[i] *= inc;
    daxpy(&nbu, &alpha, mu_u_new, &one_i, mu_u, &one_i);
    
    if (nc>0){        
        for (i=0;i<nineq;i++)
            mu[i] *= inc;
        daxpy(&nineq, &alpha, mu_new, &one_i, mu, &one_i);                    
    }
    
    if (ncN>0){                   
        for (i=0;i<ncN;i++)
            muN[i] *= inc;
        daxpy(&ncN, &alpha, muN_new, &one_i, muN, &one_i);
    }
}
