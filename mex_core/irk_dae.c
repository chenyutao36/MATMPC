#include "irk_dae.h"

#include "casadi_wrapper.h"
#include "irk_dae_kernel.h"

sim_irk_dae_workspace *sim_irk_dae_workspace_create(sim_opts *opts)
{
    const matmpc_erk_shape shape = {
        opts->nx, opts->nu, opts->num_stages, opts->num_steps, opts->forw_sens_flag, opts->adj_sens_flag
    };
    size_t doubles, pivots;
    if (matmpc_irk_dae_workspace_size(&shape, opts->nz, &doubles, &pivots) != MATMPC_IRK_OK)
        mexErrMsgIdAndTxt("MATMPC:IRK:Workspace", "Invalid IRK DAE workspace dimensions.");
    return sim_irk_mex_create(opts, doubles, pivots);
}

void sim_irk_dae_workspace_init(sim_opts *opts, const mxArray *mem, sim_irk_dae_workspace *work)
{
    sim_irk_mex_init(opts, mem, work);
}

static int residual(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, const double *z, double *f, double *g)
{
    (void)context;
    double *args[] = {(double *)x, (double *)u, (double *)p, (double *)xdot, (double *)z};
    impl_f_Fun(args, &f);
    g_Fun(args, &g);
    return 0;
}

static int jacobian(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, const double *z, double *Fx, double *Fu, double *Fxdot,
                    double *Fz, double *Gx, double *Gu, double *Gz)
{
    (void)context;
    double *args[] = {(double *)x, (double *)u, (double *)p, (double *)xdot, (double *)z};
    impl_jac_f_x_Fun(args, &Fx);
    impl_jac_f_xdot_Fun(args, &Fxdot);
    impl_jac_f_z_Fun(args, &Fz);
    impl_jac_g_x_Fun(args, &Gx);
    impl_jac_g_z_Fun(args, &Gz);
    if (Fu != NULL)
        impl_jac_f_u_Fun(args, &Fu);
    if (Gu != NULL)
        impl_jac_g_u_Fun(args, &Gu);
    return 0;
}

int sim_irk_dae(sim_in *in, sim_out *out, sim_opts *opts, sim_irk_dae_workspace *work)
{
    const size_t nx = opts->nx, nu = opts->nu, nz = opts->nz;
    const matmpc_irk_config config = sim_irk_mex_config(opts, work);
    const matmpc_irk_dae_model model = {NULL, residual, jacobian};
    const matmpc_irk_dae_input input = {
        {{in->x, nx}, {in->u, nu}, in->p, {work->Sx, nx * nx}, {work->Su, nx * nu},
         {in->lambda, opts->adj_sens_flag ? nx : 0}},
        {in->z, nz}
    };
    const matmpc_irk_dae_output output = {
        {{out->xn, nx}, {out->Sx, opts->forw_sens_flag ? nx * nx : 0},
         {out->Su, opts->forw_sens_flag ? nx * nu : 0}, {out->adj_sens, opts->adj_sens_flag ? nx + nu : 0}},
        {out->zn, nz}
    };
    sim_irk_mex_check_status(matmpc_irk_dae_integrate(&config, nz, &model, &input, &output, work->scratch), "IRK DAE");
    return 0;
}

void sim_irk_dae_workspace_free(sim_opts *opts, sim_irk_dae_workspace *work)
{
    (void)opts;
    sim_irk_mex_free(work);
}
