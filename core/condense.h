#ifndef MATMPC_CONDENSE_H
#define MATMPC_CONDENSE_H

#include <stddef.h>

typedef enum {
    MATMPC_CONDENSE_OK = 0,
    MATMPC_CONDENSE_INVALID_DIMENSIONS,
    MATMPC_CONDENSE_INVALID_BUFFER
} matmpc_condense_status;

typedef struct {
    size_t nx;
    size_t nu;
    size_t horizon;
    size_t state_bounds;
    size_t path_constraints;
    size_t terminal_constraints;
} matmpc_condense_dimensions;

typedef struct {
    const double *data;
    size_t count;
} matmpc_condense_buffer;

typedef struct {
    double *data;
    size_t count;
} matmpc_condense_output;

#endif
