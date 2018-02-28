#ifndef CASADI_WRAPPER_H_
#define CASADI_WRAPPER_H_

void f_Fun(double **in, double **out);
void vde_Fun(double **in, double **out);
void impl_f_Fun(double **in, double **out);
void impl_jac_x_Fun(double **in, double **out);
void impl_jac_u_Fun(double **in, double **out);
void impl_jac_xdot_Fun(double **in, double **out);
void F_Fun(double **in, double **out);
void D_Fun(double **in, double **out);
void Ji_Fun(double **in, double **out);
void gi_Fun(double **in, double **out);
void path_con_Fun(double **in, double **out);
void Ci_Fun(double **in, double **out);
void JN_Fun(double **in, double *out);
void gN_Fun(double **in, double **out);
void path_con_N_Fun(double **in, double **out);
void CN_Fun(double **in, double **out);
void adj_Fun(double **in, double **out);
void adjN_Fun(double **in, double **out);
void obji_Fun(double **in, double **out);
void objN_Fun(double **in, double **out);

#endif
