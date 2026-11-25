#ifndef MATMPC_ERK_KERNEL_H
#define MATMPC_ERK_KERNEL_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MATMPC_ERK_OK,
    MATMPC_ERK_INVALID_ARGUMENT,
    MATMPC_ERK_SIZE_OVERFLOW,
    MATMPC_ERK_WORKSPACE_TOO_SMALL,
    MATMPC_ERK_MODEL_ERROR
} matmpc_erk_status;

typedef struct {
    const double *data;
    size_t count;
} matmpc_erk_view;

typedef struct {
    double *data;
    size_t count;
} matmpc_erk_buffer;

typedef struct {
    size_t nx;
    size_t nu;
    size_t stages;
    size_t steps;
    bool forward;
    bool adjoint;
} matmpc_erk_shape;

typedef struct {
    matmpc_erk_shape shape;
    double step;
    matmpc_erk_view a; /* Column-major, strictly lower triangular. */
    matmpc_erk_view b;
} matmpc_erk_config;

/* Callbacks return zero on success. Parameters are passed through unchanged. */
typedef struct {
    void *context;
    int (*rhs)(void *context, const double *x, const double *u, const double *p, double *rate);
    /* Column-major tangent rates: dSx = fx*Sx, dSu = fx*Su + fu. */
    int (*forward)(void *context, const double *x, const double *u, const double *p,
                   const double *Sx, const double *Su, double *dSx, double *dSu);
    /* Product with the transposed RHS Jacobian, ordered [fx' * lambda; fu' * lambda]. */
    int (*adjoint)(void *context, const double *x, const double *u, const double *p,
                   const double *lambda, double *product);
} matmpc_erk_model;

typedef struct {
    matmpc_erk_view x;
    matmpc_erk_view u;
    const double *p;
    matmpc_erk_view Sx;     /* nx-by-nx seed, used only in forward mode. */
    matmpc_erk_view Su;     /* nx-by-nu seed, used only in forward mode. */
    matmpc_erk_view lambda; /* nx terminal adjoint, used only in adjoint mode. */
} matmpc_erk_input;

typedef struct {
    matmpc_erk_buffer x;
    matmpc_erk_buffer Sx;
    matmpc_erk_buffer Su;
    matmpc_erk_buffer adjoint; /* nx state derivatives followed by nu control derivatives. */
} matmpc_erk_output;

/* Size is returned in double elements. The caller owns the reusable workspace and all results. */
matmpc_erk_status matmpc_erk_workspace_size(const matmpc_erk_shape *shape, size_t *count);

/* Buffers must not overlap. Discard output after an error. No allocations are performed. */
matmpc_erk_status matmpc_erk_integrate(const matmpc_erk_config *config, const matmpc_erk_model *model,
                                     const matmpc_erk_input *input, const matmpc_erk_output *output,
                                     matmpc_erk_buffer workspace);

#ifdef __cplusplus
}
#endif

#endif
