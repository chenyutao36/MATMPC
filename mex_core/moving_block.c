/*                
 *       min     ( 1/2* U_new'* Hc_r* U_new  + gc_r'* U_new ) 
 *       U_new
 *       s.t.   lcx <= Ccx_r * U_new <= ucx
 *              lcc <= Ccg_r * U_new <= ucg
 *              lb_du <= T * U_new <= ub_du
 *
 *
 *      where U= T* U_new,          T is a (N x r) matrix, with rank(T)=r, 
                                    r<N, composed by N one elements.   
 *      Hc_r = T'*Hc*T 
 *      gc_r = T'*gc
 *      Ccx_r = Ccx*T
 *      Ccg_r = Ccg*T  
 */


#include "mex.h"
#include "string.h"

#include "mpc_common.h"

#include "blas.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{  

    double  *Hc, *gc, *Ccg, *Ccx; 
       
    Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc")  );
    gc = mxGetPr( mxGetField(prhs[0], 0, "gc")  );
    Ccg = mxGetPr( mxGetField(prhs[0], 0, "Ccg")  );
    Ccx = mxGetPr( mxGetField(prhs[0], 0, "Ccx")  );
    
    double  *Hc_r, *gc_r, *Ccg_r, *Ccx_r;
    Hc_r = mxGetPr( mxGetField(prhs[0], 0, "Hc_r")  );
    gc_r = mxGetPr( mxGetField(prhs[0], 0, "gc_r")  );
    Ccg_r = mxGetPr( mxGetField(prhs[0], 0, "Ccg_r")  );
    Ccx_r = mxGetPr( mxGetField(prhs[0], 0, "Ccx_r")  );
    
    double *T = mxGetPr( mxGetField(prhs[0], 0, "T")  );
    
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    size_t r = mxGetScalar( mxGetField(prhs[0], 0, "r") );
    
    double *temp = mxCalloc(r*nu*N*nu,sizeof(double));
    
    size_t newU = r*nu;
    size_t oldU = N*nu;
    size_t ncx = N*nbx;
    size_t ncg = N*nc+ncN;
    
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, zero = 0.0, minus_one = -1.0;
    size_t one_i = 1;
    
    // Hc_r
    dgemm(Trans, nTrans, &newU, &oldU, &oldU, &one_d, T, &oldU, Hc, &oldU, &zero, temp, &newU); 
    dgemm(nTrans, nTrans, &newU, &newU, &oldU, &one_d, temp, &newU, T, &oldU, &zero, Hc_r, &newU);
    
    //gc_r
    dgemv(Trans,&oldU,&newU,&one_d,T,&oldU,gc,&one_i,&zero,gc_r,&one_i);
    
    //Ccx_r
    if (ncx>0)
        dgemm(nTrans, nTrans, &ncx, &newU, &oldU, &one_d, Ccx, &ncx, T, &oldU, &zero, Ccx_r, &ncx);
    
    //Ccg_r
    if (ncg>0)
        dgemm(nTrans, nTrans, &ncg, &newU, &oldU, &one_d, Ccg, &ncg, T, &oldU, &zero, Ccg_r, &ncg);
    
    mxFree(temp);
}