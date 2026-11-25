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

static double eval_cons_res(double *x, double *u, double *od, double *z, double *x0, double *lb, double *ub, double *lc, double *uc,
                   double *lbx, double *ubx, double *lbu, double *ubu, size_t nx, size_t nu, size_t nz, size_t nc, size_t ncN,
                   size_t N, size_t np, size_t nbx, double *nbx_idx, double *eq_res_vec, const sim_workspace *sim, const double *grid)
{
    int i=0,j=0;
    sim_in *in = sim->in;
    sim_out *out = sim->out;
    
    size_t neq = (N+1)*nx;
    size_t nineq = N*nc+ncN;
    
    double *casadi_in[4];
    double *casadi_out[1];
    double *work = NULL;
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
        sim_set_stage_step(sim->opts, grid, i);
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

static double eval_curv(double *Q, double *S, double *R, double *dx, double *du,
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

static double eval_grad(double *gx, double *gu, double *dx, double *du,
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

static double eval_obj(double *x, double *u, double *od, double *y, double *yN, double *W, double *WN,
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3 || nlhs < 1 || nlhs > 2) {
        mexErrMsgIdAndTxt("MATMPC:LineSearch:Arguments",
                         "Call [alpha, meritWeight] = Line_search(mem, input, settings).");
    }

    const mxArray *mem = prhs[0];
    const mxArray *input = prhs[1];
    const mxArray *settings = prhs[2];
    const double current_merit = mxGetScalar(mxGetField(mem, 0, "mu_merit"));
    double accepted_merit = current_merit;
    double alpha = 1.0;
    const int sqp_maxit = mxGetScalar(mxGetField(mem, 0, "sqp_maxit"));
    const int sqp_it = mxGetScalar(mxGetField(mem, 0, "sqp_it"));

    /* RTI and the first SQP iteration retain their full step. */
    if (sqp_maxit > 1 && sqp_it > 0) {
        double *x = mxGetPr(mxGetField(input, 0, "x"));
        double *u = mxGetPr(mxGetField(input, 0, "u"));
        double *od = mxGetPr(mxGetField(input, 0, "od"));
        double *lb = mxGetPr(mxGetField(input, 0, "lb"));
        double *ub = mxGetPr(mxGetField(input, 0, "ub"));
        double *y = mxGetPr(mxGetField(input, 0, "y"));
        double *yN = mxGetPr(mxGetField(input, 0, "yN"));
        double *W = mxGetPr(mxGetField(input, 0, "W"));
        double *WN = mxGetPr(mxGetField(input, 0, "WN"));
        double *lbu = mxGetPr(mxGetField(input, 0, "lbu"));
        double *ubu = mxGetPr(mxGetField(input, 0, "ubu"));
        double *lbx = mxGetPr(mxGetField(input, 0, "lbx"));
        double *ubx = mxGetPr(mxGetField(input, 0, "ubx"));
        double *x0 = mxGetPr(mxGetField(input, 0, "x0"));

        size_t nx = mxGetScalar(mxGetField(settings, 0, "nx"));
        size_t nu = mxGetScalar(mxGetField(settings, 0, "nu"));
        size_t nz = mxGetScalar(mxGetField(settings, 0, "nz"));
        size_t nc = mxGetScalar(mxGetField(settings, 0, "nc"));
        size_t ncN = mxGetScalar(mxGetField(settings, 0, "ncN"));
        size_t N = mxGetScalar(mxGetField(settings, 0, "N"));
        size_t np = mxGetScalar(mxGetField(settings, 0, "np"));
        if (np == 0)
            np = 1;
        size_t ny = mxGetScalar(mxGetField(settings, 0, "ny"));
        size_t nbx = mxGetScalar(mxGetField(settings, 0, "nbx"));
        double *nbx_idx = mxGetPr(mxGetField(settings, 0, "nbx_idx"));
        size_t nx_t = nx * (N + 1);
        size_t nu_t = nu * N;
        size_t nineq = N * nc + ncN;
        size_t one = 1;

        double *Q = mxGetPr(mxGetField(mem, 0, "Q"));
        double *S = mxGetPr(mxGetField(mem, 0, "S"));
        double *R = mxGetPr(mxGetField(mem, 0, "R"));
        double *gx = mxGetPr(mxGetField(mem, 0, "gx"));
        double *gu = mxGetPr(mxGetField(mem, 0, "gu"));
        double *dx = mxGetPr(mxGetField(mem, 0, "dx"));
        double *du = mxGetPr(mxGetField(mem, 0, "du"));
        double *z = mxGetPr(mxGetField(mem, 0, "z_out"));
        const double rho = mxGetScalar(mxGetField(mem, 0, "rho"));
        const double eta = mxGetScalar(mxGetField(mem, 0, "eta"));
        const double tau = mxGetScalar(mxGetField(mem, 0, "tau"));
        const double mu_safety = mxGetScalar(mxGetField(mem, 0, "mu_safty"));
        const int sim_method = mxGetScalar(mxGetField(mem, 0, "sim_method"));

        if (!(tau > 0.0 && tau < 1.0)) {
            mexErrMsgIdAndTxt("MATMPC:LineSearch:Damping", "Line search damping tau must be between zero and one.");
        }
        if (integrator.opts == NULL)
            mexAtExit(exitFcn);
        sim_workspace_prepare(&integrator, mem, sim_method, nx, nu, nz, false, false);
        const double *grid = sim_grid_indices(mem, N, false);

        const size_t shape[] = {nx_t, nu_t};
        if (!mem_alloc || memcmp(shape, workspace_shape, sizeof(shape)) != 0) {
            free_scratch();
            eq_res_vec = mxMalloc(nx_t * sizeof(double));
            mexMakeMemoryPersistent(eq_res_vec);
            x_new = mxMalloc(nx_t * sizeof(double));
            mexMakeMemoryPersistent(x_new);
            u_new = mxMalloc(nu_t * sizeof(double));
            mexMakeMemoryPersistent(u_new);
            memcpy(workspace_shape, shape, sizeof(shape));
            mem_alloc = true;
        }

        double *lc = mxMalloc(nineq * sizeof(double));
        double *uc = mxMalloc(nineq * sizeof(double));
        const double residual = eval_cons_res(x, u, od, z, x0, lb, ub, lc, uc, lbx, ubx, lbu, ubu,
                                              nx, nu, nz, nc, ncN, N, np, nbx, nbx_idx, eq_res_vec, &integrator, grid);
        const double curvature = eval_curv(Q, S, R, dx, du, nx, nu, N);
        const double gradient = eval_grad(gx, gu, dx, du, nx, nu, N);
        const double penalty_denominator = (1.0 - rho) * residual;
        double merit = current_merit;
        if (penalty_denominator > 0.0) {
            const double lower_bound = (gradient + 0.5 * MAX(curvature, 0.0)) / penalty_denominator;
            if (lower_bound > merit)
                merit = lower_bound * mu_safety;
        }

        const double objective = eval_obj(x, u, od, y, yN, W, WN, nx, nu, np, ny, N) + merit * residual;
        const double direction = gradient - merit * residual;
        alpha = 0.0;
        if (mxIsFinite(objective) && mxIsFinite(direction)) {
            for (double trial = 1.0; trial > 1e-8; trial *= tau) {
                memcpy(x_new, x, nx_t * sizeof(double));
                memcpy(u_new, u, nu_t * sizeof(double));
                daxpy(&nx_t, &trial, dx, &one, x_new, &one);
                daxpy(&nu_t, &trial, du, &one, u_new, &one);

                const double trial_residual = eval_cons_res(x_new, u_new, od, z, x0, lb, ub, lc, uc, lbx, ubx, lbu, ubu,
                                                            nx, nu, nz, nc, ncN, N, np, nbx, nbx_idx, eq_res_vec, &integrator, grid);
                const double trial_objective = eval_obj(x_new, u_new, od, y, yN, W, WN, nx, nu, np, ny, N) +
                                               merit * trial_residual;
                if (mxIsFinite(trial_objective) && trial_objective <= objective + eta * trial * direction) {
                    alpha = trial;
                    accepted_merit = merit;
                    break;
                }
            }
        }

        mxFree(lc);
        mxFree(uc);
    }

    plhs[0] = mxCreateDoubleScalar(alpha);
    if (nlhs > 1)
        plhs[1] = mxCreateDoubleScalar(accepted_merit);
}
