#include "irk_ode.h"

#include "casadi_wrapper.h"

sim_irk_ode_workspace *sim_irk_ode_workspace_create(sim_opts *opts)
{
    const matmpc_erk_shape shape = {
        opts->nx, opts->nu, opts->num_stages, opts->num_steps, opts->forw_sens_flag, opts->adj_sens_flag
    };
    size_t doubles, pivots;
    if (matmpc_irk_ode_workspace_size(&shape, &doubles, &pivots) != MATMPC_IRK_OK)
        mexErrMsgIdAndTxt("MATMPC:IRK:Workspace", "Invalid IRK ODE workspace dimensions.");
    return sim_irk_mex_create(opts, doubles, pivots);
}

void sim_irk_ode_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_ode_workspace *work)
{
    sim_irk_mex_init(opts, mem, work);
}

static int residual(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, double *r)
{
    (void)context;
    double *args[] = {(double *)x, (double *)u, (double *)p, (double *)xdot, NULL};
    impl_f_Fun(args, &r);
    return 0;
}

static int jacobians(void *context, const double *x, const double *u, const double *p,
                     const double *xdot, double *Fx, double *Fu, double *Fxdot)
{
    (void)context;
    double *args[] = {(double *)x, (double *)u, (double *)p, (double *)xdot, NULL};
    impl_jac_f_x_Fun(args, &Fx);
    impl_jac_f_xdot_Fun(args, &Fxdot);
    if (Fu != NULL)
        impl_jac_f_u_Fun(args, &Fu);
    return 0;
}

int sim_irk_ode(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_ode_workspace *work)
{
    const size_t nx = opts->nx, nu = opts->nu;
    const matmpc_irk_config config = sim_irk_mex_config(opts, work);
    const matmpc_irk_ode_model model = {NULL, residual, jacobians};
    const matmpc_erk_input input = {
        {in->x, nx}, {in->u, nu}, in->p, {work->Sx, nx * nx}, {work->Su, nx * nu},
        {in->lambda, opts->adj_sens_flag ? nx : 0}
    };
    const matmpc_erk_output output = {
        {out->xn, nx}, {out->Sx, opts->forw_sens_flag ? nx * nx : 0},
        {out->Su, opts->forw_sens_flag ? nx * nu : 0}, {out->adj_sens, opts->adj_sens_flag ? nx + nu : 0}
    };
    sim_irk_mex_check_status(matmpc_irk_ode_integrate(&config, &model, &input, &output, work->scratch), "IRK ODE");
    return 0;
}

void sim_irk_ode_workspace_free(sim_opts *opts, sim_irk_ode_workspace *work)
{
    (void)opts;
    sim_irk_mex_free(work);
}
