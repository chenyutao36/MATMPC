
#include "mex.h"
#include "string.h"

#include "sim_workspace.h"
#include "casadi_wrapper.h"
#include "mpc_common.h"
#include "qp_gauss_newton.h"

#include <math.h>
#include <stdint.h>

#include "blas.h"

static sim_workspace integrator;
static bool mem_alloc = false;
static size_t workspace_shape[4];
static double *native_qp_workspace;
static size_t native_qp_workspace_count;
static double *Hes[1];
static double *HesN[1];
static double *Jac[2];
static double *JacN[1];
static double *temp[3];

static void free_scratch(void){
    if (Hes[0] != NULL){
        mxFree(Hes[0]);
        Hes[0] = NULL;
    }
    if (HesN[0] != NULL){
        mxFree(HesN[0]);
        HesN[0] = NULL;
    }
    if (Jac[0] != NULL){
        mxFree(Jac[0]);
        Jac[0] = NULL;
    }
    if (Jac[1] != NULL){
        mxFree(Jac[1]);
        Jac[1] = NULL;
    }
    if (JacN[0] != NULL){
        mxFree(JacN[0]);
        JacN[0] = NULL;
    }
    if (temp[0] != NULL){
        mxFree(temp[0]);
        temp[0] = NULL;
    }
    if (temp[1] != NULL){
        mxFree(temp[1]);
        temp[1] = NULL;
    }
    if (temp[2] != NULL){
        mxFree(temp[2]);
        temp[2] = NULL;
    }
    mem_alloc = false;
}

static void exitFcn(void){
    if (native_qp_workspace != NULL)
        mxFree(native_qp_workspace);
    free_scratch();
    sim_workspace_free(&integrator);
}

static double scalar_field(const mxArray *value, const char *name)
{
    const mxArray *field = mxGetField(value, 0, name);
    if (field == NULL || !mxIsNumeric(field) || mxIsComplex(field) || mxIsSparse(field) ||
        mxGetNumberOfElements(field) != 1 || !mxIsFinite(mxGetScalar(field))) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Configuration", "Invalid QP configuration field '%s'.", name);
    }
    return mxGetScalar(field);
}

static size_t dimension(const mxArray *value, const char *name)
{
    const double scalar = scalar_field(value, name);
    if (scalar < 0.0 || scalar > (double)(PTRDIFF_MAX / sizeof(double)) || floor(scalar) != scalar)
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Configuration", "Invalid QP dimension '%s'.", name);
    return (size_t)scalar;
}

static matmpc_erk_view input_field(const mxArray *value, const char *name, size_t rows, size_t columns, bool allow_vector)
{
    const size_t count = rows * columns;
    if (count == 0) {
        const matmpc_erk_view empty = {NULL, 0};
        return empty;
    }
    const mxArray *field = mxGetField(value, 0, name);
    if (field == NULL || !mxIsDouble(field) || mxIsComplex(field) || mxIsSparse(field) ||
        mxGetNumberOfElements(field) != count || mxGetNumberOfDimensions(field) != 2 ||
        !((mxGetM(field) == rows && mxGetN(field) == columns) ||
          (allow_vector && (mxGetM(field) == 1 || mxGetN(field) == 1)))) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Input", "Invalid real double array for '%s'.", name);
    }
    const matmpc_erk_view result = {mxGetPr(field), count};
    return result;
}

static matmpc_erk_buffer output_field(mxArray *value, const char *name, size_t rows, size_t columns)
{
    mxArray *matrix = mxCreateDoubleMatrix(rows, columns, mxREAL);
    mxSetField(value, 0, name, matrix);
    const matmpc_erk_buffer result = {mxGetPr(matrix), rows * columns};
    return result;
}

/* CasADi's legacy wrapper argument arrays are mutable C pointers; the model only reads them. */
static int qp_rhs(void *context, const double *x, const double *u, const double *p, double *rate)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p, NULL};
    double *out[] = {rate};
    f_Fun(in, out);
    return 0;
}

static int qp_forward(void *context, const double *x, const double *u, const double *p,
                      const double *Sx, const double *Su, double *dSx, double *dSu)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p, (double *)Sx, (double *)Su, NULL};
    double *out[] = {dSx, dSu};
    vde_Fun(in, out);
    return 0;
}

static int qp_stage_cost(void *context, const double *x, const double *u, const double *p,
                         const double *reference, const double *weight,
                         double *Jx, double *Ju, double *gx, double *gu)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p, (double *)reference, (double *)weight};
    double *jacobian[] = {Jx, Ju};
    double *gradient[] = {gx, gu};
    Ji_Fun(in, jacobian);
    gi_Fun(in, gradient);
    return 0;
}

static int qp_terminal_cost(void *context, const double *x, const double *p,
                            const double *reference, const double *weight, double *J, double *gx)
{
    (void)context;
    double *in[] = {(double *)x, (double *)p, (double *)reference, (double *)weight};
    double *jacobian[] = {J};
    double *gradient[] = {gx};
    JN_Fun(in, jacobian);
    gN_Fun(in, gradient);
    return 0;
}

static int qp_path_constraint(void *context, const double *x, const double *u, const double *p,
                              double *value, double *Cx, double *Cu)
{
    (void)context;
    double *in[] = {(double *)x, (double *)u, (double *)p};
    double *constraint[] = {value};
    double *jacobian[] = {Cx, Cu};
    path_con_Fun(in, constraint);
    Ci_Fun(in, jacobian);
    return 0;
}

static int qp_terminal_constraint(void *context, const double *x, const double *p, double *value, double *C)
{
    (void)context;
    double *in[] = {(double *)x, (double *)p};
    double *constraint[] = {value};
    double *jacobian[] = {C};
    path_con_N_Fun(in, constraint);
    CN_Fun(in, jacobian);
    return 0;
}

static int qp_irk_stage(void *context, double step, const double *x, const double *u, const double *p,
                        const double *z, double *xn, double *zn, double *A, double *B)
{
    sim_workspace *work = context;
    sim_in in = {(double *)x, (double *)u, (double *)p, (double *)z, NULL};
    sim_out out = {xn, zn, A, B, NULL};
    work->opts->h = step;
    if (work->method == 3)
        return sim_irk_dae(&in, &out, work->opts, work->irk_dae);
    return sim_irk_ode(&in, &out, work->opts, work->irk_ode);
}

static void check_model_dimensions(const matmpc_qp_config *config)
{
    size_t dimensions[4];
    const int dynamics = matmpc_model_dimensions(dimensions);
    if (dynamics < 0 || (dynamics > 0 && (dimensions[0] != config->integration.shape.nx ||
            dimensions[1] != config->integration.shape.nu || dimensions[2] != config->np || dimensions[3] != config->nz))) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:ModelDimensions", "QP dimensions do not match the compiled dynamics.");
    }
    const int costs = matmpc_qp_model_dimensions(dimensions);
    if (costs < 0 || (costs > 0 && (dimensions[0] != config->ny || dimensions[1] != config->ny_terminal ||
            dimensions[2] != config->nc || dimensions[3] != config->nc_terminal))) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:ModelDimensions", "QP dimensions do not match the compiled objective and constraints.");
    }
}

static mxArray *generate_explicit_qp(const mxArray *input, const mxArray *settings, const mxArray *mem,
                                     const mxArray *block_boundaries)
{
    if (!mxIsStruct(input) || !mxIsStruct(settings) || !mxIsStruct(mem) ||
        mxGetNumberOfElements(input) != 1 || mxGetNumberOfElements(settings) != 1 || mxGetNumberOfElements(mem) != 1) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Configuration", "input, settings and mem must be scalar structures.");
    }

    const size_t method = dimension(mem, "sim_method");
    const size_t stages = dimension(mem, "num_stages");
    const size_t nz = dimension(settings, "nz");
    if (dimension(mem, "hessian") != 0 || (nz != 0) != (method == 3) ||
        method < 1 || method > 3 || (method == 1 && stages != 4)) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Mode", "Explicit QP output requires Gauss-Newton and a compatible ERK or IRK integrator.");
    }

    matmpc_qp_config config = {0};
    const size_t nx = dimension(settings, "nx");
    const size_t nu = dimension(settings, "nu");
    const size_t N = dimension(settings, "N");
    config.stage_grid.data = sim_grid_indices(mem, N, false);
    if (config.stage_grid.data != NULL)
        config.stage_grid.count = N + 1;
    config.integration.shape.nx = nx;
    config.nz = nz;
    config.integration.shape.nu = nu;
    config.integration.shape.stages = stages;
    config.integration.shape.steps = dimension(mem, "num_steps");
    config.integration.shape.forward = true;
    config.integration.step = scalar_field(mem, "h");
    if (method != 1) {
        config.integrate_stage = qp_irk_stage;
        config.integration_context = &integrator;
    }
    config.horizon = N;
    config.np = dimension(settings, "np");
    config.ny = dimension(settings, "ny");
    config.ny_terminal = dimension(settings, "nyN");
    config.nc = dimension(settings, "nc");
    config.nc_terminal = dimension(settings, "ncN");
    config.nbx = dimension(settings, "nbx");
    config.diagonal_floor = scalar_field(mem, "reg");
    config.levenberg_marquardt = mpc_levenberg_marquardt(mem);
    size_t blocks = 0;
    if (block_boundaries != NULL) {
        if (config.stage_grid.count != 0) {
            mexErrMsgIdAndTxt("MATMPC:IncompatibleGrid", "Move blocking requires a uniform shooting grid.");
        }
        blocks = dimension(mem, "r");
        if (blocks == 0 || blocks > N || !mxIsDouble(block_boundaries) || mxIsComplex(block_boundaries) ||
            mxIsSparse(block_boundaries) || mxGetNumberOfDimensions(block_boundaries) != 2 ||
            (mxGetM(block_boundaries) != 1 && mxGetN(block_boundaries) != 1) ||
            mxGetNumberOfElements(block_boundaries) != blocks + 1) {
            mexErrMsgIdAndTxt("MATMPC:QPGeneration:ControlBlocks", "Expected r + 1 real double control-block boundaries.");
        }
        config.control_blocks.data = mxGetPr(block_boundaries);
        config.control_blocks.count = blocks + 1;
    }

    check_model_dimensions(&config);
    size_t count;
    if (dimension(mem, "nx") != nx || dimension(mem, "nu") != nu || dimension(mem, "nz") != nz ||
        matmpc_qp_workspace_size(&config, &count) != MATMPC_QP_OK) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Configuration", "Invalid QP workspace dimensions.");
    }
    if (method == 1) {
        config.integration.a = input_field(mem, "A_tab", 4, 4, false);
        config.integration.b = input_field(mem, "B_tab", 1, 4, true);
        config.state_seed = input_field(mem, "Sx", nx, nx, false);
        config.control_seed = input_field(mem, "Su", nx, nu, false);
    } else {
        if (integrator.opts == NULL)
            mexAtExit(exitFcn);
        sim_workspace_prepare(&integrator, mem, (int)method, nx, nu, nz, true, false);
    }

    const size_t ny = config.ny, nyN = config.ny_terminal;
    const size_t nc = config.nc, ncN = config.nc_terminal;
    const size_t nbx = config.nbx;
    const size_t path_rows = N * nc + ncN;
    const matmpc_qp_input problem = {
        input_field(input, "x0", nx, 1, true), input_field(input, "x", nx, N + 1, false),
        input_field(input, "u", nu, N, false), input_field(input, "od", config.np, N + 1, false),
        input_field(input, "y", ny, N, false), input_field(input, "yN", nyN, 1, true),
        input_field(input, "W", ny, N, false), input_field(input, "WN", nyN, 1, true),
        input_field(input, "lbu", nu, N, true), input_field(input, "ubu", nu, N, true),
        input_field(input, "lbx", nbx, N, true), input_field(input, "ubx", nbx, N, true),
        input_field(input, "lb", path_rows, 1, true), input_field(input, "ub", path_rows, 1, true),
        input_field(input, "z", nz, N, false)
    };
    if (nbx != 0) {
        const matmpc_erk_view indices = input_field(settings, "nbx_idx", nbx, 1, true);
        size_t *native_indices = mxMalloc(nbx * sizeof(size_t));
        for (size_t j = 0; j < nbx; ++j) {
            const double index = indices.data[j];
            if (!(index >= 1.0 && index <= (double)nx) || floor(index) != index)
                mexErrMsgIdAndTxt("MATMPC:QPGeneration:StateBounds", "State-bound indices must lie within the model's state vector.");
            native_indices[j] = (size_t)index - 1;
        }
        config.bounded_states = native_indices;
    }

    const char *fields[] = {"A", "B", "a", "ds0", "Q", "S", "R", "gx", "gu", "Cgx", "Cgu", "CgN", "lc", "uc",
                            "lb_du", "ub_du", "lb_dx", "ub_dx", "R_block", "gu_block", "z_out"};
    int field_count = (int)(sizeof(fields) / sizeof(fields[0])) - 1 - (blocks == 0 ? 2 : 0);
    if (nz != 0)
        fields[field_count++] = "z_out";
    mxArray *arrays = mxCreateStructMatrix(1, 1, field_count, fields);
    const matmpc_erk_buffer empty = {NULL, 0};
    const matmpc_qp_output output = {
        output_field(arrays, "A", nx, N * nx), output_field(arrays, "B", nx, N * nu),
        output_field(arrays, "a", nx, N), output_field(arrays, "ds0", nx, 1),
        output_field(arrays, "Q", nx, (N + 1) * nx), output_field(arrays, "S", nx, N * nu),
        output_field(arrays, "R", nu, N * nu), output_field(arrays, "gx", nx, N + 1), output_field(arrays, "gu", nu, N),
        output_field(arrays, "Cgx", nc, N * nx), output_field(arrays, "Cgu", nc, N * nu),
        output_field(arrays, "CgN", ncN, nx), output_field(arrays, "lc", path_rows, 1), output_field(arrays, "uc", path_rows, 1),
        output_field(arrays, "lb_du", N * nu, 1), output_field(arrays, "ub_du", N * nu, 1),
        output_field(arrays, "lb_dx", N * nbx, 1), output_field(arrays, "ub_dx", N * nbx, 1),
        blocks != 0 ? output_field(arrays, "R_block", nu, blocks * nu) : empty,
        blocks != 0 ? output_field(arrays, "gu_block", nu, blocks) : empty,
        nz != 0 ? output_field(arrays, "z_out", nz, N) : empty
    };

    if (native_qp_workspace_count != count) {
        double *next = mxMalloc(count * sizeof(double));
        mexAtExit(exitFcn);
        mexMakeMemoryPersistent(next);
        if (native_qp_workspace != NULL)
            mxFree(native_qp_workspace);
        native_qp_workspace = next;
        native_qp_workspace_count = count;
    }
    const matmpc_qp_model model = {
        {NULL, qp_rhs, qp_forward, NULL}, NULL, qp_stage_cost, qp_terminal_cost, qp_path_constraint, qp_terminal_constraint
    };
    const matmpc_erk_buffer workspace = {native_qp_workspace, native_qp_workspace_count};
    const matmpc_qp_status status = matmpc_qp_gauss_newton(&config, &model, &problem, &output, workspace);
    if (config.bounded_states != NULL)
        mxFree((void *)config.bounded_states);
    if (status != MATMPC_QP_OK)
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Assembly", "Native QP assembly failed (status %d).", (int)status);
    return arrays;
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    if ((nrhs != 3 && (nrhs != 4 || nlhs != 1)) || nlhs > 1) {
        mexErrMsgIdAndTxt("MATMPC:QPGeneration:Arguments", "Call qp = qp_generation(input, settings, mem) with optional control-block boundaries.");
    }
    if (nlhs == 1) {
        plhs[0] = generate_explicit_qp(prhs[0], prhs[1], prhs[2], nrhs == 4 ? prhs[3] : NULL);
        return;
    }
    /* Retained for legacy no-output calls and the other integration/Hessian modes. */

    double *x = mxGetPr( mxGetField(prhs[0], 0, "x") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "u") );
    double *z = mxGetPr( mxGetField(prhs[0], 0, "z") );
    double *y = mxGetPr( mxGetField(prhs[0], 0, "y") );
    double *yN = mxGetPr( mxGetField(prhs[0], 0, "yN") );
    double *od = mxGetPr( mxGetField(prhs[0], 0, "od") );
    double *W = mxGetPr( mxGetField(prhs[0], 0, "W") );
    double *WN = mxGetPr( mxGetField(prhs[0], 0, "WN") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub") );
    double *x0 = mxGetPr( mxGetField(prhs[0], 0, "x0") );
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lbu") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ubu") );
    double *lbx = mxGetPr( mxGetField(prhs[0], 0, "lbx") );
    double *ubx = mxGetPr( mxGetField(prhs[0], 0, "ubx") );
    
    size_t nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    size_t nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    size_t nz = mxGetScalar( mxGetField(prhs[1], 0, "nz") );
    size_t np = mxGetScalar( mxGetField(prhs[1], 0, "np") ); if(np==0) np++;
    size_t ny = mxGetScalar( mxGetField(prhs[1], 0, "ny") );
    size_t nyN = mxGetScalar( mxGetField(prhs[1], 0, "nyN") );
    size_t nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    size_t ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    size_t nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    double *nbx_idx = mxGetPr( mxGetField(prhs[1], 0, "nbx_idx") );
    size_t N = mxGetScalar( mxGetField(prhs[1], 0, "N") );    
    int sim_method = mxGetScalar( mxGetField(prhs[2], 0, "sim_method") );
    
    int i=0,j=0;
    char *nTrans = "N", *Trans="T", *UPLO="L";
    double one_d = 1.0, zero = 0.0, minus_one_d = -1.0;
    size_t one_i = 1;
    int idx;
      
    double *Q = mxGetPr( mxGetField(prhs[2], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[2], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[2], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[2], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[2], 0, "B") );
    double *Cx = mxGetPr( mxGetField(prhs[2], 0, "Cx") );
    double *Cgx = mxGetPr( mxGetField(prhs[2], 0, "Cgx") );
    double *Cgu = mxGetPr( mxGetField(prhs[2], 0, "Cgu") );
    double *CgN = mxGetPr( mxGetField(prhs[2], 0, "CgN") );
    double *gx = mxGetPr( mxGetField(prhs[2], 0, "gx") );
    double *gu = mxGetPr( mxGetField(prhs[2], 0, "gu") );   
    double *a = mxGetPr( mxGetField(prhs[2], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[2], 0, "ds0") );
    double *lc = mxGetPr( mxGetField(prhs[2], 0, "lc") );
    double *uc = mxGetPr( mxGetField(prhs[2], 0, "uc") );
    double *lb_du = mxGetPr( mxGetField(prhs[2], 0, "lb_du") );
    double *ub_du = mxGetPr( mxGetField(prhs[2], 0, "ub_du") );
    double *lb_dx = mxGetPr( mxGetField(prhs[2], 0, "lb_dx") );
    double *ub_dx = mxGetPr( mxGetField(prhs[2], 0, "ub_dx") );
    double *z_out = mxGetPr( mxGetField(prhs[2], 0, "z_out") );
    
    double reg = mxGetScalar( mxGetField(prhs[2], 0, "reg") );
    const double lm = mpc_levenberg_marquardt(prhs[2]);
    int hessian_type = mxGetScalar( mxGetField(prhs[2], 0, "hessian") );
    
    for (i=0;i<nx;i++)
        ds0[i] = x0[i] - x[i];
    
    // allocate memory    
    double *Cons[2];
      
    double *casadi_in[5];
    double *casadi_out[2];    
          
    if (integrator.opts == NULL)
        mexAtExit(exitFcn);
    sim_workspace_prepare(&integrator, prhs[2], sim_method, nx, nu, nz, true, false);

    const size_t shape[] = { nx, nu, ny, nyN };
    if (!mem_alloc || memcmp(shape, workspace_shape, sizeof(shape)) != 0){
        free_scratch();
        
        Hes[0] = (double *) mxMalloc(ny*ny * sizeof(double));	
        mexMakeMemoryPersistent(Hes[0]); 	
        HesN[0] = (double *) mxMalloc(nyN*nyN * sizeof(double));	
        mexMakeMemoryPersistent(HesN[0]); 	
                       
        Jac[0] = (double *) mxMalloc(ny*nx * sizeof(double));	
        mexMakeMemoryPersistent(Jac[0]); 	
        Jac[1] = (double *) mxMalloc(ny*nu * sizeof(double));	
        mexMakeMemoryPersistent(Jac[1]); 	
        JacN[0] = (double *) mxMalloc(nyN*nx * sizeof(double));	
        mexMakeMemoryPersistent(JacN[0]);       
        
        temp[0] = (double *) mxMalloc(ny*nx * sizeof(double));	
        mexMakeMemoryPersistent(temp[0]); 	
        temp[1] = (double *) mxMalloc(ny*nu * sizeof(double));	
        mexMakeMemoryPersistent(temp[1]);        
        temp[2] = (double *) mxMalloc(nyN*nx * sizeof(double));	
        mexMakeMemoryPersistent(temp[2]); 
        
        memcpy(workspace_shape, shape, sizeof(shape));
        mem_alloc = true;
    }
    
    // start loop
    for(i=0;i<N;i++){
        casadi_in[0] = x+i*nx;
        casadi_in[1] = u+i*nu;
        casadi_in[2] = od+i*np;
        casadi_in[3] = y+i*ny;
        casadi_in[4] = W+i*ny;
        
        // control bounds
        for (j=0;j<nu;j++){
            lb_du[i*nu+j] = lbu[i*nu+j]-u[i*nu+j];
            ub_du[i*nu+j] = ubu[i*nu+j]-u[i*nu+j];
        }
        
        // state bounds
        for (j=0;j<nbx;j++){
            idx = (int)nbx_idx[j]-1;
            lb_dx[i*nbx+j] = lbx[i*nbx+j]-x[(i+1)*nx+idx];
            ub_dx[i*nbx+j] = ubx[i*nbx+j]-x[(i+1)*nx+idx];
        }
        
        // integration                      
        switch(sim_method){
            case 1:
                integrator.in->x = x+i*nx;
                integrator.in->u = u+i*nu;
                integrator.in->p = od+i*np;
                integrator.out->xn = a+i*nx;
                integrator.out->Sx = A + i*nx*nx;
                integrator.out->Su = B + i*nx*nu;
                sim_erk(integrator.in, integrator.out, integrator.opts, integrator.erk);
                break;
            case 2:
                integrator.in->x = x+i*nx;
                integrator.in->u = u+i*nu;
                integrator.in->p = od+i*np;
                integrator.in->z = z+i*nz;
                integrator.out->xn = a+i*nx;
                integrator.out->Sx = A + i*nx*nx;
                integrator.out->Su = B + i*nx*nu;
                sim_irk_ode(integrator.in, integrator.out, integrator.opts, integrator.irk_ode);
                break;
            case 3:
                integrator.in->x = x+i*nx;
                integrator.in->u = u+i*nu;
                integrator.in->p = od+i*np;
                integrator.in->z = z+i*nz;
                integrator.out->xn = a+i*nx;
                integrator.out->Sx = A + i*nx*nx;
                integrator.out->Su = B + i*nx*nu;
                integrator.out->zn = z_out + i*nz;
                sim_irk_dae(integrator.in, integrator.out, integrator.opts, integrator.irk_dae);
                break;
            default:
                mexErrMsgTxt("Please choose a supported integrator");
                break;
        }
       
        // equality residual        
        for (j=0;j<nx;j++)
            a[i*nx+j] -= x[(i+1)*nx+j];
       
        // Hessian
        Ji_Fun(casadi_in, Jac);
        switch(hessian_type){                
            case 0:                   
                dgemm(Trans, nTrans, &nx, &nx, &ny, &one_d, Jac[0], &ny, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
                dgemm(Trans, nTrans, &nx, &nu, &ny, &one_d, Jac[0], &ny, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
                dgemm(Trans, nTrans, &nu, &nu, &ny, &one_d, Jac[1], &ny, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                break;
            
            case 1:
                Hi_Fun(casadi_in, Hes);
                dgemm(Trans, nTrans, &nx, &ny, &ny, &one_d, Jac[0], &ny, Hes[0], &ny, &zero, temp[0], &nx);
                dgemm(nTrans, nTrans, &nx, &nx, &ny, &one_d, temp[0], &nx, Jac[0], &ny, &zero, Q+i*nx*nx, &nx);
                
                dgemm(nTrans, nTrans, &nx, &nu, &ny, &one_d, temp[0], &nx, Jac[1], &ny, &zero, S+i*nx*nu, &nx);
                
                dgemm(Trans, nTrans, &nu, &ny, &ny, &one_d, Jac[1], &ny, Hes[0], &ny, &zero, temp[1], &nu);
                dgemm(nTrans, nTrans, &nu, &nu, &ny, &one_d, temp[1], &nu, Jac[1], &ny, &zero, R+i*nu*nu, &nu);
                
                break;
            default:
                mexErrMsgTxt("Please choose a supported Hessian type");
                break;
                            
        }
        mpc_regularize_hessian_diagonal(nx, Q+i*nx*nx, reg, lm);
        mpc_regularize_hessian_diagonal(nu, R+i*nu*nu, reg, lm);
        
        // gradient
        casadi_out[0] = gx+i*nx;
        casadi_out[1] = gu+i*nu;
        gi_Fun(casadi_in, casadi_out);
                
        // constraint residual
        if (nc>0){  
            casadi_out[0] = lc + i*nc;
            path_con_Fun(casadi_in, casadi_out);
            for (j=0;j<nc;j++){
                uc[i*nc+j] = ub[i*nc+j] - casadi_out[0][j];
                casadi_out[0][j] = lb[i*nc+j] - casadi_out[0][j];            
            }
        
            // constraint Jacobian
            Cons[0] = Cgx+i*nc*nx;
            Cons[1] = Cgu+i*nc*nu;
            Ci_Fun(casadi_in, Cons);
        }
    }
    
    // terminal data
    casadi_in[0] = x+N*nx;
    casadi_in[1] = od+N*np;
    casadi_in[2] = yN;
    casadi_in[3] = WN;
    
    JN_Fun(casadi_in, JacN);
    switch(hessian_type){            
        case 0:
            dgemm(Trans, nTrans, &nx, &nx, &nyN, &one_d, JacN[0], &nyN, JacN[0], &nyN, &zero, Q+N*nx*nx, &nx);
            break;
        case 1:
            HN_Fun(casadi_in, HesN);
            dgemm(Trans, nTrans, &nx, &nyN, &nyN, &one_d, JacN[0], &nyN, HesN[0], &nyN, &zero, temp[2], &nx);
            dgemm(nTrans, nTrans, &nx, &nx, &nyN, &one_d, temp[2], &nx, JacN[0], &nyN, &zero, Q+N*nx*nx, &nx);
            break;
        default:
            mexErrMsgTxt("Please choose a supported Hessian type");
            break;
            
    }
    mpc_regularize_hessian_diagonal(nx, Q+N*nx*nx, reg, lm);
        
    casadi_out[0] = gx+N*nx;
    gN_Fun(casadi_in, casadi_out);

    if (ncN>0){
        casadi_out[0] = lc + N*nc;
        path_con_N_Fun(casadi_in, casadi_out);
        for (j=0;j<ncN;j++){
            uc[N*nc+j] = ub[N*nc+j] - casadi_out[0][j];
            casadi_out[0][j] = lb[N*nc+j] - casadi_out[0][j];            
        }

        CN_Fun(casadi_in, &CgN);
    }
    
}
