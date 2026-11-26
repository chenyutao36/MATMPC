#ifndef MATMPC_IRK_KERNEL_H
#define MATMPC_IRK_KERNEL_H

#include "erk_kernel.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MATMPC_IRK_OK,
    MATMPC_IRK_INVALID_ARGUMENT,
    MATMPC_IRK_SIZE_OVERFLOW,
    MATMPC_IRK_WORKSPACE_TOO_SMALL,
    MATMPC_IRK_MODEL_ERROR,
    MATMPC_IRK_SINGULAR,
    MATMPC_IRK_NONFINITE,
    MATMPC_IRK_NONCONVERGENCE
} matmpc_irk_status;

typedef struct {
    matmpc_erk_buffer doubles;
    size_t *pivots;
    size_t pivot_count;
} matmpc_irk_workspace;

typedef struct {
    matmpc_erk_shape shape;
    double step;         /* Duration of each of shape.steps substeps. */
    matmpc_erk_view a;   /* Square, column-major Butcher tableau. */
    matmpc_erk_view b;
    size_t newton_iterations;
    double newton_tolerance; /* Applied to the norm of the final Newton correction. */
    /* Optional column-major final-state sensitivity weights at base_step.
       An empty view derives the diagonal stage weights from step * b. */
    matmpc_erk_view sensitivity_weights;
    double base_step;
} matmpc_irk_config;

typedef struct {
    void *context;
    /* F(x,u,p,xdot) = 0. Callbacks return zero on success. */
    int (*residual)(void *context, const double *x, const double *u, const double *p,
                    const double *xdot, double *r);
    /* Column-major Jacobians. Fu is NULL when the control Jacobian is unnecessary. */
    int (*jacobians)(void *context, const double *x, const double *u, const double *p,
                     const double *xdot, double *Fx, double *Fu, double *Fxdot);
} matmpc_irk_ode_model;

/* Returns the capacity of each caller-owned workspace partition in elements. */
matmpc_irk_status matmpc_irk_ode_workspace_size(const matmpc_erk_shape *shape,
                                               size_t *double_count, size_t *pivot_count);

/* No allocations are performed. Buffers must not overlap. Discard results after an error.
   Optional ODE sensitivity_weights is nx-by-(nx * stages). It affects derivatives only.
   The adjoint differentiates the interval at its input x and u, independently of the seeds. */
matmpc_irk_status matmpc_irk_ode_integrate(const matmpc_irk_config *config,
                                         const matmpc_irk_ode_model *model,
                                         const matmpc_erk_input *input,
                                         const matmpc_erk_output *output,
                                         matmpc_irk_workspace workspace);

#ifdef __cplusplus
}
#endif

#endif
