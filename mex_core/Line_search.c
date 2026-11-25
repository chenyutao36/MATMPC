#include "mex.h"
#include "string.h"

#include "sim_workspace.h"
#include "casadi_wrapper.h"

#include "blas.h"
#include "lapack.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

static double *eq_res_vec = NULL;
static double *x_new = NULL, *u_new = NULL;

static sim_workspace integrator;
static bool mem_alloc = false;
static size_t workspace_shape[2];

static void free_scratch(void){
    if (eq_res_vec != NULL){
        mxFree(eq_res_vec);
        eq_res_vec = NULL;
    }
    if (x_new != NULL){
        mxFree(x_new);
        x_new = NULL;
    }
    if (u_new != NULL){
        mxFree(u_new);
        u_new = NULL;
    }
    mem_alloc = false;
}

static void exitFcn(void){
    free_scratch();
    sim_workspace_free(&integrator);
}

double eval_cons_res(double *x, double *u, double *od, double *z, double *x0, double *lb, double *ub, double *lc, double *uc,
                   double *lbx, double *ubx, double *lbu, double *ubu, size_t nx, size_t nu, size_t nz, size_t nc, size_t ncN,
                   size_t N, size_t np, size_t nbx, double *nbx_idx, double *eq_res_vec, const sim_workspace *sim)
{
    int i=0,j=0;
    sim_in *in = sim->in;
    sim_out *out = sim->out;
    
    size_t neq = (N+1)*nx;
    size_t nineq = N*nc+ncN;
    
    double *casadi_in[4];
    double *casadi_out[1];
    double *work;
    double eq_res=0, ineq_res=0, cons_res;
    
    char *nTrans = "N", *Trans="T", *Norm="O";
    double one_d = 1.0, zero = 0.0;
    size_t one_i = 1;
    int idx;
    
    double *lu = (double *)mxMalloc( N*nu * sizeof(double));        
    double *uu = (double *)mxMalloc( N*nu * sizeof(double));
    double *lx = (double *)mxMalloc( N*nbx * sizeof(double));        
    double *ux = (double *)mxMalloc( N*nbx * sizeof(double));
    double *zn = (double *)mxMalloc( nz * sizeof(double));
    
    for (j=0;j<nx;j++)
        eq_res_vec[j] = x0[j] - x[j];
           
    for (i=0;i<N;i++){      
        switch(sim->method){
            case 1:               
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                out->xn = eq_res_vec+(i+1)*nx;
                sim_erk(in, out, sim->opts, sim->erk);
                break;
            case 2:  
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                out->xn = eq_res_vec+(i+1)*nx;
                sim_irk_ode(in, out, sim->opts, sim->irk_ode);
                break;
            case 3:  
                in->x = x+i*nx;
                in->u = u+i*nu;
                in->p = od+i*np;
                in->z = z+i*nz;
                out->xn = eq_res_vec+(i+1)*nx;
                out->zn = zn;
                sim_irk_dae(in, out, sim->opts, sim->irk_dae);
                break;
            default :
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
                
        for (j=0;j<nx;j++)
            eq_res_vec[(i+1)*nx+j] -= x[(i+1)*nx+j];
        
        
        for (j=0;j<nbx;j++){
            idx = (int)nbx_idx[j]-1;
            lx[i*nbx+j] = lbx[i*nbx+j] - x[(i+1)*nx+idx];
            ux[i*nbx+j] = ubx[i*nbx+j] - x[(i+1)*nx+idx];
        }
        
        for (j=0;j<nu;j++){
            lu[i*nu+j] = lbu[i*nu+j] - u[i*nu+j];
            uu[i*nu+j] = ubu[i*nu+j] - u[i*nu+j];
        }
            
        if (nc>0){
            casadi_in[0]=x+i*nx;
            casadi_in[1]=u+i*nu;
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
        casadi_in[0] = x+N*nx;
        casadi_in[1] = od+N*np;
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out);
        for (j=0;j<ncN;j++){
            uc[i*nc+j] = ub[N*nc+j] - casadi_out[0][j];
            casadi_out[0][j] = lb[N*nc+j] - casadi_out[0][j];            
        }
    }
        
    eq_res = dlange(Norm, &neq, &one_i, eq_res_vec, &neq, work);
            
    for (i=0;i<N*nu;i++)
        ineq_res += MAX(-1*uu[i],0) + MAX(lu[i],0);
    for (i=0;i<N*nbx;i++)
        ineq_res += MAX(-1*ux[i],0) + MAX(lx[i],0);
    for (i=0;i<nineq;i++)
        ineq_res += MAX(-1*uc[i],0) + MAX(lc[i],0);
        
    cons_res = eq_res + ineq_res;
    
    mxFree(lu);
    mxFree(uu);
    mxFree(lx);
    mxFree(ux);
    mxFree(zn);

    return cons_res;
}

double eval_curv(double *Q, double *S, double *R, double *dx, double *du,
                 size_t nx, size_t nu, size_t N)
{
    int i;
    double curv = 0.0;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1; 
    
    double *tmp = (double *)mxMalloc((nx+nu)*sizeof(double));
    
    for(i=0; i<N; i++)
    {
        dgemv(nTrans,&nx,&nx,&one_d,Q+i*nx*nx,&nx,dx+i*nx,&one_i,&zero,tmp,&one_i);
        dgemv(nTrans,&nx,&nu,&one_d,S+i*nx*nu,&nx,du+i*nu,&one_i,&one_d,tmp,&one_i);
        curv += ddot(&nx, dx+i*nx, &one_i, tmp, &one_i);
        
        dgemv(Trans,&nx,&nu,&one_d,S+i*nx*nu,&nx,dx+i*nx,&one_i,&zero,tmp+nx,&one_i);
        dgemv(nTrans,&nu,&nu,&one_d,R+i*nu*nu,&nu,du+i*nu,&one_i,&one_d,tmp+nx,&one_i);
        curv += ddot(&nu, du+i*nu, &one_i, tmp+nx, &one_i);
    }
    
    dgemv(nTrans,&nx,&nx,&one_d,Q+N*nx*nx,&nx,dx+N*nx,&one_i,&zero,tmp,&one_i);
    curv += ddot(&nx, dx+N*nx, &one_i, tmp, &one_i);
    
    mxFree(tmp);
    
    return curv;
}

double eval_grad(double *gx, double *gu, double *dx, double *du,
                 size_t nx, size_t nu, size_t N)
{
    int i;
    
    double grad = 0.0;   
    
    size_t one_i = 1; 
    
    for(i=0; i<N; i++)
    {
        grad += ddot(&nx, dx+i*nx, &one_i, gx+i*nx, &one_i);
        grad += ddot(&nu, du+i*nu, &one_i, gu+i*nu, &one_i);
    }
    
    grad += ddot(&nx, dx+N*nx, &one_i, gx+N*nx, &one_i);
 
    return grad;
}

double eval_obj(double *x, double *u, double *od, double *y, double *yN, double *W, double *WN,
                size_t nx, size_t nu, size_t np, size_t ny, size_t N)
{
    int i;
    
    double *casadi_in[5];
    double *casadi_out[1];
    double obj=0.0;
    
    
    casadi_out[0] = (double *) mxMalloc(sizeof(double));
    for (i=0; i<N; i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = od+i*np;
        casadi_in[3] = y+i*ny;
        casadi_in[4] = W+i*ny;
        
        obji_Fun(casadi_in, casadi_out);
        obj += *casadi_out[0];
    }
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    
    objN_Fun(casadi_in, casadi_out);
    obj += *casadi_out[0];
    
    mxFree(casadi_out[0]);
    
    return obj;
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *x = mxGetPr( mxGetField(prhs[1], 0, "x") );
    double *u = mxGetPr( mxGetField(prhs[1], 0, "u") );
    // double *z = mxGetPr( mxGetField(prhs[1], 0, "z") );
    double *lambda = mxGetPr( mxGetField(prhs[1], 0, "lambda") );
    double *mu = mxGetPr( mxGetField(prhs[1], 0, "mu") );
    double *mu_u = mxGetPr( mxGetField(prhs[1], 0, "mu_u") );
    double *mu_x = mxGetPr( mxGetField(prhs[1], 0, "mu_x") );
    double *od = mxGetPr( mxGetField(prhs[1], 0, "od") );
    double *lb = mxGetPr( mxGetField(prhs[1], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[1], 0, "ub") );
    double *y = mxGetPr( mxGetField(prhs[1], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[1], 0, "yN") );
    double *W = mxGetPr( mxGetField(prhs[1], 0, "W") );
    double *WN = mxGetPr( mxGetField(prhs[1], 0, "WN") );
    double *lbu = mxGetPr( mxGetField(prhs[1], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[1], 0, "ubu") );  
    double *lbx = mxGetPr( mxGetField(prhs[1], 0, "lbx") );
    double *ubx = mxGetPr( mxGetField(prhs[1], 0, "ubx") ); 
    double *x0 = mxGetPr( mxGetField(prhs[1], 0, "x0") );
    
    size_t nx = mxGetScalar( mxGetField(prhs[2], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[2], 0, "nu") );
    size_t nz = mxGetScalar( mxGetField(prhs[2], 0, "nz") );
    size_t nc = mxGetScalar( mxGetField(prhs[2], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[2], 0, "ncN") );
    size_t N = mxGetScalar( mxGetField(prhs[2], 0, "N") );    
    size_t np = mxGetScalar( mxGetField(prhs[2], 0, "np") ); if(np==0) np++;
    size_t ny = mxGetScalar( mxGetField(prhs[2], 0, "ny") );
    size_t nbx = mxGetScalar( mxGetField(prhs[2], 0, "nbx") );
    double *nbx_idx = mxGetPr( mxGetField(prhs[2], 0, "nbx_idx") );
            
    size_t neq = (N+1)*nx;
    size_t nineq = N*nc+ncN;
    size_t nx_t = nx*(N+1);
    size_t nu_t = nu*N;
    size_t nbu_t = N*nu;
    size_t nbx_t = N*nbx;
    
    double one_d = 1.0;
    size_t one_i = 1;
    
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *gx = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[0], 0, "gu") );    
    double *dx = mxGetPr( mxGetField(prhs[0], 0, "dx") );
    double *du = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *lambda_new = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu_new = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *mu_u_new = mxGetPr( mxGetField(prhs[0], 0, "mu_u_new") );
    double *mu_x_new = mxGetPr( mxGetField(prhs[0], 0, "mu_x_new") );
    double *lc = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[0], 0, "uc") );
    double *a = mxGetPr( mxGetField(prhs[0], 0, "a") );
          
    double rho = mxGetScalar( mxGetField(prhs[0], 0, "rho") );
    double eta = mxGetScalar( mxGetField(prhs[0], 0, "eta") );
    double tau = mxGetScalar( mxGetField(prhs[0], 0, "tau") );
    double mu_safty = mxGetScalar( mxGetField(prhs[0], 0, "mu_safty") );
    double *mu_merit = mxGetPr( mxGetField(prhs[0], 0, "mu_merit") );
    double *alpha = mxGetPr( mxGetField(prhs[0], 0, "alpha") );
    
    double *q_dual = mxGetPr( mxGetField(prhs[0], 0, "q_dual") );          
    double *dmu = mxGetPr( mxGetField(prhs[0], 0, "dmu") );

    double *z = mxGetPr( mxGetField(prhs[0], 0, "z_out") );
    
    
    int i=0,j=0;
    
    int sim_method = mxGetScalar( mxGetField(prhs[0], 0, "sim_method") );
    int sqp_maxit = mxGetScalar( mxGetField(prhs[0], 0, "sqp_maxit") ); 
    int sqp_it = mxGetScalar( mxGetField(prhs[0], 0, "sqp_it") ); 
       
    if (integrator.opts == NULL)
        mexAtExit(exitFcn);
    sim_workspace_prepare(&integrator, prhs[0], sim_method, nx, nu, nz, false, false);

    const size_t shape[] = { nx_t, nu_t };
    if (!mem_alloc || memcmp(shape, workspace_shape, sizeof(shape)) != 0){
        free_scratch();
        eq_res_vec = (double *)mxMalloc( neq * sizeof(double));
        mexMakeMemoryPersistent(eq_res_vec);
        
        x_new = (double *)mxMalloc( nx_t * sizeof(double));
        mexMakeMemoryPersistent(x_new);
        u_new = (double *)mxMalloc( nu_t * sizeof(double));
        mexMakeMemoryPersistent(u_new);
              
                
        memcpy(workspace_shape, shape, sizeof(shape));
        mem_alloc = true;
    }
    
    
    // backtracking
    double cons_res;
    double sigma, pd, grad, mu_lb=0, obj, obj_new, dir_grad, obj_tmp;
    int newpoint = 0;
    alpha[0] = 1;
    if (sqp_maxit > 1 && sqp_it > 0){
        cons_res = eval_cons_res(x, u, od, z, x0, lb, ub, lc, uc,
                                 lbx, ubx, lbu, ubu, nx, nu, nz, nc, ncN,
                                 N, np, nbx, nbx_idx, eq_res_vec, &integrator);
              
        pd = eval_curv(Q, S, R, dx, du, nx, nu, N);
        grad = eval_grad(gx, gu, dx, du, nx, nu, N);
                
        if (pd>0)
            sigma = 0.5;
        else
            sigma = 0.0;
        
//         if (cons_res>0)
        mu_lb=(grad+sigma*pd)/(1-rho)/cons_res;
                
        if (mu_merit[0]<mu_lb)
            mu_merit[0]=mu_lb*mu_safty;
        
        obj = eval_obj(x, u, od, y, yN, W, WN, nx, nu, np, ny, N);
        obj += mu_merit[0]*cons_res;
                
        dir_grad = grad - mu_merit[0] * cons_res;
        
//         mexPrintf("\nl:%f  pd:%f  mu_lb:%f  phi:%f  D:%f\n",cons_res,pd,mu_lb,obj,dir_grad);
                                
        while (newpoint!=1 && alpha[0] > 1E-8){
            memcpy(x_new, x, nx_t*sizeof(double));
            memcpy(u_new, u, nu_t*sizeof(double));
            
            daxpy(&nx_t, alpha, dx, &one_i, x_new, &one_i); 
            daxpy(&nu_t, alpha, du, &one_i, u_new, &one_i);
            cons_res = eval_cons_res(x_new, u_new, od, z, x0, lb, ub, lc, uc,
                                     lbx, ubx, lbu, ubu, nx, nu, nz, nc, ncN,
                                     N, np, nbx, nbx_idx, eq_res_vec, &integrator);
                        
            obj_new = eval_obj(x_new, u_new, od, y, yN, W, WN, nx, nu, np, ny, N);
            obj_new += mu_merit[0]*cons_res;
            
            // dir_grad = grad - mu_merit[0] * cons_res;
            
            obj_tmp = obj + eta*alpha[0]*dir_grad;
            
//             mexPrintf("l:%f  phi:%f  phi_rhs:%f\n",cons_res,obj_new,obj_tmp);
                        
            if (obj_new <= obj_tmp)
                newpoint = 1;
            else
                alpha[0] *= tau;
        }
                
    }
               
    // update     
    daxpy(&nx_t, alpha, dx, &one_i, x, &one_i);
    
    daxpy(&nu_t, alpha, du, &one_i, u, &one_i);
    
    for (i=0;i<neq;i++)
        q_dual[i] = alpha[0]*(lambda_new[i] - lambda[i]);
    daxpy(&neq, &one_d, q_dual, &one_i, lambda, &one_i);
    
    for (i=0;i<nbu_t;i++)
        dmu[i] = alpha[0]*(mu_u_new[i] - mu_u[i]);
    daxpy(&nbu_t, &one_d, dmu, &one_i, mu_u, &one_i);
    
    if (nbx_t>0){
        for (i=0;i<nbx_t;i++)
            dmu[nbu_t+i] = alpha[0]*(mu_x_new[i] - mu_x[i]);
        daxpy(&nbx_t, &one_d, dmu+nbu_t, &one_i, mu_x, &one_i);

    }
    
    if (nineq>0){        
        for (i=0;i<nineq;i++)
            dmu[nbu_t+nbx_t+i] = alpha[0]*(mu_new[i] - mu[i]);
        daxpy(&nineq, &one_d, dmu+nbu_t+nbx_t, &one_i, mu, &one_i);                
    }
    
}
