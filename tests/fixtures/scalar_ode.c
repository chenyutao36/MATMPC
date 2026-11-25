#include "casadi_wrapper.h"

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
