
#include "mex.h"
#include "string.h"
#include "stdbool.h"

#include "RTI_step_common.h"
#include "RTI_step_funcs.h"
#include "mpc_common.h"

#include <math.h>               /* isinf */
#include "qpsolver_dense.h"

static void *work = NULL;
static double *lb_qore = NULL;
static double *ub_qore = NULL;
static double *x_qore = NULL;
static QoreProblemDense *problem[1] = {0}; 

void exitFcn(){
    
    if (work!=NULL)
        mxFree(work);
    
    mxFree(lb_qore);
    mxFree(ub_qore);
    mxFree(x_qore);
    
    if (problem!=NULL)
        QPDenseFree(&problem[0]);
    
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    /* from input */
    double *x = mxGetPr( mxGetField(prhs[0], 0, "x") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "u") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *W = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *WN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
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
    
    /* from settings */
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    int ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    int nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    int nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    
    /* from memory */ 
    double *Q = mxGetPr( mxGetField(prhs[2], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[2], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[2], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[2], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[2], 0, "B") );
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
    double *CN = mxGetPr( mxGetField(prhs[2], 0, "CN") );
    int lin_obj = mxGetScalar( mxGetField(prhs[2], 0, "lin_obj") );    
    
    double *G = mxGetPr( mxGetField(prhs[2], 0, "G")  );
    double *Hc = mxGetPr( mxGetField(prhs[2], 0, "Hc")  );
    double *gc = mxGetPr( mxGetField(prhs[2], 0, "gc")  );
//     double *Cc = mxGetPr( mxGetField(prhs[2], 0, "Cc")  );
    double *Cc = mxGetPr( mxGetField(prhs[2], 0, "Cc_qore")  );
    double *lcc = mxGetPr( mxGetField(prhs[2], 0, "lcc")  );
    double *ucc = mxGetPr( mxGetField(prhs[2], 0, "ucc")  );
    int iter = mxGetScalar( mxGetField(prhs[2], 0, "iter") );
    int hot_start = mxGetScalar( mxGetField(prhs[2], 0, "hot_start") );
    
    double *dx = mxGetPr( mxGetField(prhs[2], 0, "dx") );
    double *du = mxGetPr( mxGetField(prhs[2], 0, "du") );

    bool cond_save = (hot_start==1) && (lin_obj==1) ;
    
    rti_step_dims dim;
    dim.nx = nx;
    dim.nu = nu;
    dim.np = np;
    dim.ny = ny;
    dim.nyN = nyN;
    dim.nc = nc;
    dim.ncN = ncN;   
    dim.N = N;
    
    int size = rti_step_calculate_workspace_size(&dim);
    int err;
    if (work==NULL){
        work = mxMalloc(size);
        mexMakeMemoryPersistent(work); 
        
        lb_qore = mxMalloc((N*nu+N*nc+ncN)*sizeof(double));
        ub_qore = mxMalloc((N*nu+N*nc+ncN)*sizeof(double));
        mexMakeMemoryPersistent(lb_qore); 
        mexMakeMemoryPersistent(ub_qore); 
        
        x_qore = mxMalloc((N*nu+N*nc+ncN)*sizeof(double));
        mexMakeMemoryPersistent(x_qore); 
        
        err = QPDenseNew(&problem[0], N*nu, N*nc+ncN);
        err = QPDenseSetInt(problem[0], "prtfreq", -1);
        
        mexAtExit(exitFcn);
    }
    
    rti_step_workspace *workspace = (rti_step_workspace *) rti_step_cast_workspace(&dim, work);
    
    qp_generation(Q, S, R, A, B, Cx, Cu, CN, 
        gx, gu, a, ds0, lc, uc, lb_du, ub_du, 
        x0, x, u, y, yN, od, W, WN, lb, ub,
        lbN, ubN, lbu, ubu, lin_obj,
        &dim, workspace);
    
    condensing(Q, S, R, A, B, Cx, Cu, CN, 
        gx, gu, a, ds0, lc, uc, 
        G, Hc, gc, Cc, lcc, ucc,
        iter, cond_save,
        &dim, workspace);
    
//     mxArray *qpoases_in[10];
//     mxArray *qpoases_out[7];
//     double *du;
//     if (iter == 1){
//         qpoases_in[0] = mxCreateString("i");        
//         qpoases_in[1] = mxGetField(prhs[2], 0, "Hc");
//         qpoases_in[2] = mxGetField(prhs[2], 0, "gc");
//         qpoases_in[3] = mxGetField(prhs[2], 0, "Cc");
//         qpoases_in[4] = mxGetField(prhs[2], 0, "lb_du");
//         qpoases_in[5] = mxGetField(prhs[2], 0, "ub_du");
//         qpoases_in[6] = mxGetField(prhs[2], 0, "lcc");
//         qpoases_in[7] = mxGetField(prhs[2], 0, "ucc");
//         qpoases_in[8] = mxGetField(prhs[2], 0, "qpoases_opt");
//         
//         mexCallMATLAB(7, qpoases_out, 9, qpoases_in, "qpOASES_sequence");
//         
//         mxSetField(prhs[2], 0, "warm_start", qpoases_out[0]);       
//         du = mxGetPr(qpoases_out[1]);        
//     }else{
//         if (!hot_start){
//             qpoases_in[0] = mxCreateString("m");
//             qpoases_in[1] = mxGetField(prhs[2], 0, "warm_start");
//             qpoases_in[2] = mxGetField(prhs[2], 0, "Hc");
//             qpoases_in[3] = mxGetField(prhs[2], 0, "gc");
//             qpoases_in[4] = mxGetField(prhs[2], 0, "Cc");
//             qpoases_in[5] = mxGetField(prhs[2], 0, "lb_du");
//             qpoases_in[6] = mxGetField(prhs[2], 0, "ub_du");
//             qpoases_in[7] = mxGetField(prhs[2], 0, "lcc");
//             qpoases_in[8] = mxGetField(prhs[2], 0, "ucc");       
//             qpoases_in[9] = mxGetField(prhs[2], 0, "qpoases_opt");
//             mexCallMATLAB(6, qpoases_out, 10, qpoases_in, "qpOASES_sequence");
//         }else{
//             qpoases_in[0] = mxCreateString("h");
//             qpoases_in[1] = mxGetField(prhs[2], 0, "warm_start");
//             qpoases_in[2] = mxGetField(prhs[2], 0, "gc");
//             qpoases_in[3] = mxGetField(prhs[2], 0, "lb_du");
//             qpoases_in[4] = mxGetField(prhs[2], 0, "ub_du");
//             qpoases_in[5] = mxGetField(prhs[2], 0, "lcc");
//             qpoases_in[6] = mxGetField(prhs[2], 0, "ucc");       
//             qpoases_in[7] = mxGetField(prhs[2], 0, "qpoases_opt");
//             mexCallMATLAB(6, qpoases_out, 8, qpoases_in, "qpOASES_sequence");
//         }                        
//         du = mxGetPr(qpoases_out[0]);
//     }
    
    memcpy(lb_qore, lb_du, N*nu*sizeof(double));
    memcpy(lb_qore+N*nu, lcc, (N*nc+ncN)*sizeof(double));
    memcpy(ub_qore, ub_du, N*nu*sizeof(double));
    memcpy(ub_qore+N*nu, ucc, (N*nc+ncN)*sizeof(double));
    if (iter==1){
        err = QPDenseSetData(problem[0], N*nu, N*nc+ncN, Cc, Hc);
        err = QPDenseOptimize(problem[0], lb_qore, ub_qore, gc, NULL, NULL);
    }else{
        if (!hot_start){
            //         err = QPDenseSetInt(problem[0], "warmstrategy", 1);
            err = QPDenseUpdateMatrices(problem[0], N*nu, N*nc+ncN, Cc, Hc);
        }
        err = QPDenseOptimize(problem[0], lb_qore, ub_qore, gc, NULL, NULL);
    }   
    err = QPDenseGetDblVector(problem[0], "primalsol", x_qore);
    memcpy(du, x_qore, N*nu*sizeof(double));
    
    recover(A, B, a, ds0, dx, du, &dim);
    
    line_search(dx, du, x, u, &dim);
}