
#include "mex.h"
#include <stdbool.h>
#include "string.h"

#include "mpc_common.h"
#include "condense_affine.h"

#include "blas.h"

static double *G = NULL, *W_mat = NULL, *affine_workspace = NULL;
static double *Hi = NULL, *Cci = NULL, *Ccxi = NULL, *CcN = NULL;
static bool mem_alloc = false;
static size_t workspace_shape[6];

static void exitFcn(void){
    if (mem_alloc){
        mxFree(G);
        mxFree(W_mat);
        mxFree(affine_workspace);
        mxFree(Hi);
        mxFree(Cci);   
        mxFree(CcN);
        mxFree(Ccxi);   
    }
    mem_alloc = false;
}

static const mxArray *affine_array(const mxArray *mem, const char *name, size_t count)
{
    const mxArray *array = mxGetField(mem, 0, name);

    if (array == NULL || !mxIsDouble(array) || mxIsComplex(array) || mxIsSparse(array) ||
        mxGetNumberOfElements(array) != count) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:AffineBuffer",
                         "mem.%s must contain the expected number of full real double elements.", name);
    }
    return array;
}

static matmpc_condense_buffer affine_input(const mxArray *mem, const char *name, size_t count)
{
    matmpc_condense_buffer result = {mxGetPr(affine_array(mem, name, count)), count};
    return result;
}

static matmpc_condense_output affine_output(const mxArray *mem, const char *name, size_t count)
{
    matmpc_condense_output result = {mxGetPr(affine_array(mem, name, count)), count};
    return result;
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{    
    if (nlhs != 0 || nrhs != 2 || !mxIsStruct(prhs[0]) || !mxIsStruct(prhs[1]) ||
        mxGetNumberOfElements(prhs[0]) != 1 || mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:Arguments", "Use Condensing(mem, settings) with scalar structures.");
    }

    /*Inputs*/
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *Cgx = mxGetPr( mxGetField(prhs[0], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[0], 0, "Cgu") );
    double *CgN = mxGetPr( mxGetField(prhs[0], 0, "CgN") );
        
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") );
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
            
    /*Outputs*/
    double *Hc, *Ccg, *Ccx;
       
    Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc")  );
    Ccg = mxGetPr( mxGetField(prhs[0], 0, "Ccg")  );
    Ccx = mxGetPr( mxGetField(prhs[0], 0, "Ccx")  );
    const matmpc_condense_dimensions dimensions = {nx, nu, N, nbx, nc, ncN};
    matmpc_condense_affine_sizes sizes;
    if (matmpc_condense_affine_get_sizes(&dimensions, &sizes) != MATMPC_CONDENSE_OK) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:AffineDimensions", "The affine workspace dimensions are invalid.");
    }
    const mxArray *mem = prhs[0];
    const matmpc_condense_affine_inputs affine = {
        affine_input(mem, "Q", sizes.q), affine_input(mem, "S", sizes.s),
        affine_input(mem, "A", sizes.a), affine_input(mem, "B", sizes.s),
        affine_input(mem, "Cx", sizes.cx), affine_input(mem, "Cgx", sizes.cgx),
        affine_input(mem, "CgN", sizes.cgn), affine_input(mem, "gx", sizes.states),
        affine_input(mem, "gu", sizes.controls), affine_input(mem, "a", sizes.defect),
        affine_input(mem, "ds0", nx), affine_input(mem, "lc", sizes.general),
        affine_input(mem, "uc", sizes.general), affine_input(mem, "lb_dx", sizes.state),
        affine_input(mem, "ub_dx", sizes.state)
    };
    const matmpc_condense_affine_outputs shifted = {
        affine_output(mem, "gc", sizes.controls),
        affine_output(mem, "lcc", sizes.general), affine_output(mem, "ucc", sizes.general),
        affine_output(mem, "lxc", sizes.state), affine_output(mem, "uxc", sizes.state)
    };
       
    /*Allocate memory*/
    int i=0,j=0;
     
    const size_t shape[] = { nx, nu, nc, ncN, nbx, N };
    if (!mem_alloc || memcmp(shape, workspace_shape, sizeof(shape)) != 0){
        exitFcn();
        mexAtExit(exitFcn);
        G = (double *)mxMalloc(nx*N*N*nu * sizeof(double));
        W_mat = (double *)mxMalloc(nx*N*N*nu * sizeof(double));
        affine_workspace = (double *)mxMalloc(sizes.workspace * sizeof(double));
        Hi = (double *)mxMalloc(nu*nu * sizeof(double));
        Cci = (double *)mxMalloc(nc*nu * sizeof(double));
        CcN = (double *)mxMalloc(ncN*nu * sizeof(double));
        Ccxi = (double *)mxMalloc(nbx*nu * sizeof(double));

        mexMakeMemoryPersistent(G);
        mexMakeMemoryPersistent(W_mat);
        mexMakeMemoryPersistent(affine_workspace);
        mexMakeMemoryPersistent(Hi);
        mexMakeMemoryPersistent(Cci);
        mexMakeMemoryPersistent(CcN);
        mexMakeMemoryPersistent(Ccxi);
        memcpy(workspace_shape, shape, sizeof(shape));
        mem_alloc = true;
    }
   
    char *nTrans = "N", *Trans="T", *SIDE = "L", *UPLO = "L";
    double one_d = 1.0, zero = 0.0;
              
    /*Start the loop*/
        
    /* compute G */
    for(i=0;i<N;i++){
        memcpy(G+(i*N+i)*nx*nu, B+i*nx*nu, nx*nu*sizeof(double));
        for (j=i+1;j<N;j++){
            dgemm(nTrans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, G+(i*N+j)*nx*nu, &nx);
        }
    }
    
    /* Compute Hc */
    for(i=0;i<N;i++){
        dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+N*nx*nx, &nx, G+(i*N+N-1)*nx*nu, &nx, &zero, W_mat+(i*N+N-1)*nx*nu, &nx);
        for(j=N-1;j>i;j--){        
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, S+j*nx*nu, &nx, G+(i*N+j-1)*nx*nu, &nx, &zero, Hi, &nu);                     
            dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+j*nx*nu, &nx, W_mat+(i*N+j)*nx*nu, &nx, &one_d, Hi, &nu);
            Block_Fill(nu, nu, Hi, Hc, j*nu, i*nu, N*nu);
            Block_Fill_Trans(nu, nu, Hi, Hc, i*nu, j*nu, N*nu);
            dgemm(Trans, nTrans, &nx, &nu, &nx, &one_d, A+j*nx*nx, &nx, W_mat+(i*N+j)*nx*nu, &nx, &zero, W_mat+(i*N+j-1)*nx*nu, &nx); 
            dsymm(SIDE, UPLO, &nx, &nu, &one_d, Q+j*nx*nx, &nx, G+(i*N+j-1)*nx*nu, &nx, &one_d, W_mat+(i*N+j-1)*nx*nu, &nx);
        }
        memcpy(Hi,R+i*nu*nu,nu*nu*sizeof(double));
        dgemm(Trans, nTrans, &nu, &nu, &nx, &one_d, B+i*nx*nu, &nx, W_mat+(i*N+i)*nx*nu, &nx, &one_d, Hi, &nu);
        Block_Fill(nu, nu, Hi, Hc, i*nu, i*nu, N*nu);      
    }
    
    /* Compute Cc */
    if (nc>0){         
        for(i=0;i<N;i++){
            Block_Fill(nc, nu, Cgu+i*nc*nu, Ccg, i*nc, i*nu, N*nc+ncN);
            for(j=i+1;j<N;j++){   
                dgemm(nTrans, nTrans, &nc, &nu, &nx, &one_d, Cgx+j*nc*nx, &nc, G+(i*N+j-1)*nx*nu, &nx, &zero, Cci, &nc);
                Block_Fill(nc, nu, Cci, Ccg, j*nc, i*nu, N*nc+ncN);
            }    
        }  
    }
    
    /* Compute CcN */
    if (ncN>0){          
        for(i=0;i<N;i++){                 
            dgemm(nTrans, nTrans, &ncN, &nu, &nx, &one_d, CgN, &ncN, G+(i*N+N-1)*nx*nu, &nx, &zero, CcN, &ncN);
            Block_Fill(ncN, nu, CcN, Ccg, N*nc, i*nu, N*nc+ncN);
        }
    }
    
    /* Compute Ccx */
    if (nbx>0){         
        for(i=0;i<N;i++){
            for(j=i+1;j<=N;j++){   
                dgemm(nTrans, nTrans, &nbx, &nu, &nx, &one_d, Cx, &nbx, G+(i*N+j-1)*nx*nu, &nx, &zero, Ccxi, &nbx);
//                    Block_Fill(nbx, nu, Ccxi, Ccx, j*nbx, i*nu, (N+1)*nbx);
                Block_Fill(nbx, nu, Ccxi, Ccx, (j-1)*nbx, i*nu, N*nbx);
            }    
        }  
    }
         
    const matmpc_condense_output scratch = {affine_workspace, sizes.workspace};
    if (matmpc_condense_affine(&dimensions, &affine, &shifted, scratch) != MATMPC_CONDENSE_OK) {
        mexErrMsgIdAndTxt("MATMPC:Condensing:AffineBuffer", "Unable to condense the affine terms.");
    }
}
