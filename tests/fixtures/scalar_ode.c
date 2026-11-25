#include "casadi_wrapper.h"
#include <math.h>

/* Analytic model x' = u. For DAE linkage, the algebraic equation is z = 0. */
void f_Fun(double **in, double **out) { out[0][0] = in[1][0]; }
void vde_Fun(double **in, double **out) { out[0][0] = 0.0; out[1][0] = 1.0; }
void adj_ERK_Fun(double **in, double **out) { out[0][0] = 0.0; out[0][1] = in[3][0]; }
void impl_f_Fun(double **in, double **out) { out[0][0] = in[3][0] - in[1][0]; }
void g_Fun(double **in, double **out) { out[0][0] = in[4][0]; }
void impl_jac_f_x_Fun(double **in, double **out) { out[0][0] = 0.0; }
void impl_jac_f_u_Fun(double **in, double **out) { out[0][0] = -1.0; }
void impl_jac_f_xdot_Fun(double **in, double **out) { out[0][0] = 1.0; }
void impl_jac_f_z_Fun(double **in, double **out) { out[0][0] = 0.0; }
void impl_jac_g_x_Fun(double **in, double **out) { out[0][0] = 0.0; }
void impl_jac_g_u_Fun(double **in, double **out) { out[0][0] = 0.0; }
void impl_jac_g_z_Fun(double **in, double **out) { out[0][0] = 1.0; }

/* Half-weighted squared control error, with a terminal state error. */
void Ji_Fun(double **in, double **out) { out[0][0] = 0.0; out[1][0] = sqrt(in[4][0]); }
void Hi_Fun(double **in, double **out) { out[0][0] = 1.0; }
void gi_Fun(double **in, double **out) { out[0][0] = 0.0; out[1][0] = in[4][0] * (in[1][0] - in[3][0]); }
void JN_Fun(double **in, double **out) { out[0][0] = sqrt(in[3][0]); }
void HN_Fun(double **in, double **out) { out[0][0] = 1.0; }
void gN_Fun(double **in, double **out) { out[0][0] = in[3][0] * (in[0][0] - in[2][0]); }
void adj_Fun(double **in, double **out)
{
    out[0][0] = 0.0;
    out[0][1] = in[4][0] * (in[1][0] - in[3][0]);
    out[1][0] = 0.0;
    out[1][1] = 0.0;
}
void adjN_Fun(double **in, double **out) { gN_Fun(in, out); out[1][0] = 0.0; }
void obji_Fun(double **in, double **out)
{
    const double error = in[1][0] - in[3][0];
    out[0][0] = 0.5 * in[4][0] * error * error;
}
void objN_Fun(double **in, double **out)
{
    const double error = in[0][0] - in[2][0];
    out[0][0] = 0.5 * in[3][0] * error * error;
}

/* No path constraints are enabled in the scalar test model. */
void path_con_Fun(double **in, double **out) { }
void Ci_Fun(double **in, double **out) { }
void path_con_N_Fun(double **in, double **out) { }
void CN_Fun(double **in, double **out) { }
