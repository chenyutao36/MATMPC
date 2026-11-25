#ifndef MATMPC_QP_GAUSS_NEWTON_H
#define MATMPC_QP_GAUSS_NEWTON_H

#include "erk_kernel.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MATMPC_QP_OK,
    MATMPC_QP_INVALID_ARGUMENT,
    MATMPC_QP_SIZE_OVERFLOW,
    MATMPC_QP_WORKSPACE_TOO_SMALL,
    MATMPC_QP_MODEL_ERROR
} matmpc_qp_status;

typedef struct {
    matmpc_erk_config integration; /* Four-stage explicit integration with forward sensitivities. */
    size_t horizon, np, ny, ny_terminal, nc, nc_terminal;
    const size_t *bounded_states; /* Zero-based indices of length nbx. */
    size_t nbx;
    matmpc_erk_view state_seed, control_seed;
    double diagonal_floor;
    /* N+1 boundary indices or an empty view for a uniform grid. Each interval scales
     * the base microstep without changing the integration step count. */
    matmpc_erk_view stage_grid;
    matmpc_erk_view control_blocks; /* Uniform grid only: integer stage boundaries from 0 to horizon. */
    double levenberg_marquardt; /* Unscaled diagonal shift after the legacy floor. Zero preserves the legacy QP. */
} matmpc_qp_config;

/* Callbacks return zero on success. Jacobians are column-major; objective Jacobians
 * already include the weights used by the Gauss-Newton approximation. */
typedef struct {
    matmpc_erk_model dynamics;
    void *context;
    int (*stage_cost)(void *context, const double *x, const double *u, const double *p,
                      const double *reference, const double *weight,
                      double *Jx, double *Ju, double *gx, double *gu);
    int (*terminal_cost)(void *context, const double *x, const double *p,
                         const double *reference, const double *weight, double *J, double *gx);
    int (*path_constraint)(void *context, const double *x, const double *u, const double *p,
                           double *value, double *Cx, double *Cu);
    int (*terminal_constraint)(void *context, const double *x, const double *p,
                               double *value, double *C);
} matmpc_qp_model;

/* Stage matrices are concatenated in time. All views carry exact element counts. */
typedef struct {
    matmpc_erk_view x0, x, u, parameters; /* x and parameters include the terminal stage. */
    matmpc_erk_view reference, reference_terminal, weight, weight_terminal;
    matmpc_erk_view lower_u, upper_u, lower_x, upper_x;
    matmpc_erk_view lower_path, upper_path; /* Path bounds followed by terminal bounds. */
} matmpc_qp_input;

typedef struct {
    matmpc_erk_buffer A, B, a, ds0;
    matmpc_erk_buffer Q, S, R, gx, gu;
    matmpc_erk_buffer Cgx, Cgu, CgN, lc, uc;
    matmpc_erk_buffer lb_du, ub_du, lb_dx, ub_dx;
    matmpc_erk_buffer R_block, gu_block; /* With blocks, stage R stays raw and the diagonal floor applies once per block. */
} matmpc_qp_output;

/* Workspace size in doubles. No allocations or backend initialization occur in this kernel. */
matmpc_qp_status matmpc_qp_workspace_size(const matmpc_qp_config *config, size_t *count);

/* Writable buffers must be separate from inputs and one another. Discard outputs on error. */
matmpc_qp_status matmpc_qp_gauss_newton(const matmpc_qp_config *config, const matmpc_qp_model *model,
                                        const matmpc_qp_input *input, const matmpc_qp_output *output,
                                        matmpc_erk_buffer workspace);

#ifdef __cplusplus
}
#endif

#endif
