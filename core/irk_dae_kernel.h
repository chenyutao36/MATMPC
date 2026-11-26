#ifndef MATMPC_IRK_DAE_KERNEL_H
#define MATMPC_IRK_DAE_KERNEL_H

#include "irk_kernel.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    void *context;
    /* F(x,u,p,xdot,z) = 0, G(x,u,p,z) = 0. Callbacks return zero on success. */
    int (*residual)(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, const double *z, double *f, double *g);
    /* Dense column-major Jacobians. Fu and Gu are NULL when unnecessary. */
    int (*jacobian)(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, const double *z, double *Fx, double *Fu,
                    double *Fxdot, double *Fz, double *Gx, double *Gu, double *Gz);
} matmpc_irk_dae_model;

typedef struct {
    matmpc_erk_input differential;
    matmpc_erk_view z; /* Initial algebraic guess, repeated at each stage. */
} matmpc_irk_dae_input;

typedef struct {
    matmpc_erk_output differential;
    matmpc_erk_buffer z; /* Algebraic state of the last collocation stage. */
} matmpc_irk_dae_output;

/* Returns capacities in double and size_t elements. nx, nz, stages, and steps are positive. */
matmpc_irk_status matmpc_irk_dae_workspace_size(const matmpc_erk_shape *shape, size_t nz,
                                               size_t *double_count, size_t *pivot_count);

/*
 * No allocations or retained state. Buffers must not overlap and may be reused
 * after any result. Outputs are committed after a successful interval.
 * Optional sensitivity_weights is nx-by-((nx+nz)*stages), scaled by step/base_step.
 * An empty view derives step*b(i) for the differential components of stage i.
 * Forward uses the supplied seeds. Adjoint differentiates with respect to the
 * interval input x and u with the supplied terminal differential multiplier.
 */
matmpc_irk_status matmpc_irk_dae_integrate(const matmpc_irk_config *config, size_t nz,
                                          const matmpc_irk_dae_model *model,
                                          const matmpc_irk_dae_input *input,
                                          const matmpc_irk_dae_output *output,
                                          matmpc_irk_workspace workspace);

#ifdef __cplusplus
}
#endif

#endif
