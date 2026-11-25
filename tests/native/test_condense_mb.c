#include "condense_mb.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void compare(const double *actual, const double *expected, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (!(fabs(actual[i] - expected[i]) < 1e-10)) {
            fprintf(stderr, "At %lu expected %.17g, got %.17g\n", (unsigned long)i, expected[i], actual[i]);
            exit(EXIT_FAILURE);
        }
    }
}

static void test_shared_input_cross_term(void)
{
    const matmpc_condense_dimensions dims = {1, 1, 2, 0, 0, 0};
    const size_t partition[] = {0, 2};
    const size_t invalid[] = {1, 2};
    const double q[] = {0, 1, 0}, s[] = {0, 1}, a[] = {0, 0}, b[] = {1, 0};
    const double gx[] = {0, 0, 0}, defect[] = {0, 0}, zero = 0, r = 1;
    const matmpc_condense_mb_inputs input = {
        {{q, 3}, {s, 2}, {a, 2}, {b, 2}, {NULL, 0}, {NULL, 0}, {NULL, 0},
         {gx, 3}, {&zero, 1}, {defect, 2}, {&zero, 1}, {NULL, 0}, {NULL, 0}, {NULL, 0}, {NULL, 0}},
        {&r, 1}, {NULL, 0}
    };
    double hc, gc, workspace[5];
    const matmpc_condense_mb_outputs output = {
        {{&hc, 1}, {NULL, 0}, {NULL, 0}}, {{&gc, 1}, {NULL, 0}, {NULL, 0}, {NULL, 0}, {NULL, 0}}
    };
    const matmpc_condense_output scratch = {workspace, 5};

    /* x1 = u0 = u1 = v, so (x1 + u1)^2/2 = 2*v^2. */
    require(matmpc_condense_mb(&dims, partition, 2, &input, &output, scratch) == MATMPC_CONDENSE_OK,
            "Single-block condensing failed");
    require(hc == 4 && gc == 0, "The shared-input Hessian must include both mixed terms");
    require(matmpc_condense_mb(&dims, invalid, 2, &input, &output, scratch) == MATMPC_CONDENSE_INVALID_DIMENSIONS &&
            hc == 4 && gc == 0, "Reject invalid partitions before writing outputs");
    const matmpc_condense_output short_scratch = {workspace, 4};
    require(matmpc_condense_mb(&dims, partition, 2, &input, &output, short_scratch) == MATMPC_CONDENSE_INVALID_BUFFER,
            "Reject undersized workspace");
}

static void test_matches_full_projection(void)
{
    const matmpc_condense_dimensions dims = {2, 2, 4, 2, 1, 2};
    const size_t partitions[][3] = {{0, 1, 4}, {0, 3, 4}, {0, 2, 4}, {0, 1, 4}};
    const double q[] = {2, 1, 99, 3, 4, -1, 101, 2, 6, 2, -91, 5, 1, 0, -92, 2, 2, 1, -93, 3};
    const double s[] = {1, -1, 2, 0, 0, 3, 0, 2, 1, -1, 2, 0, 2, -2, 0, 3};
    const double r[] = {1, 0.5, 0.5, 1, 2, -1, -1, 3, 2, 0, 0, 1, 2, 1, 1, 4};
    const double a[] = {1, 0, 2, -1, 0.5, 1, 0, 2, 1, 0.5, 1, -1, 0.5, -2, 1.5, 1};
    const double b[] = {1, -1, 0, 3, 2, 0, 0, 1, 2, 0, -2, 1, 1.5, 0.5, -1, 2};
    const double cx[] = {2, 0.5, -1, 3}, cgx[] = {1, 2, -1, 1.5, 0.5, 2, 1, -3};
    const double cgu[] = {1, 2, 0, -2, 0.5, 1, -1, 5}, cgn[] = {1, 3, -2, 0.5};
    const double gx[] = {1, -1, -2, 3, 4, 1, 2, -1, 1, -3}, gu[] = {0.25, -0.5, 2, 1, -2, 3, 1, 2};
    const double defect[] = {0.5, -0.5, 1, -2, 3, 1, -1, 2}, ds0[] = {1, -2};
    const double lc[] = {-10, -5, -8, -6, -2, -9}, uc[] = {10, 8, 11, 12, 6, 13};
    const double lbx[] = {-4, -7, -3, -11, -3, -5, -6, -7}, ubx[] = {5, 10, 6, 7, 3, 8, 5, 7};
    const matmpc_condense_affine_inputs affine = {
        {q, 20}, {s, 16}, {a, 16}, {b, 16}, {cx, 4}, {cgx, 8}, {cgn, 4},
        {gx, 10}, {gu, 8}, {defect, 8}, {ds0, 2}, {lc, 6}, {uc, 6}, {lbx, 8}, {ubx, 8}
    };
    const matmpc_condense_matrix_inputs matrices = {
        affine.q, affine.s, {r, 16}, affine.a, affine.b, affine.cx, affine.cgx, {cgu, 8}, affine.cgn
    };
    double hc[64], ccx[64], ccg[48], gc[8], lo[6], hi[6], xlo[8], xhi[8], scratch[132];
    const matmpc_condense_matrix_outputs full_matrix = {{hc, 64}, {ccx, 64}, {ccg, 48}};
    const matmpc_condense_affine_outputs full_affine = {{gc, 8}, {lo, 6}, {hi, 6}, {xlo, 8}, {xhi, 8}};
    const matmpc_condense_output full_scratch = {scratch, 132};
    double hmb[16], xcmb[32], gcmb[24], gmb[4], lmb[6], umb[6], xlmb[8], xumb[8], buffer[24];
    const matmpc_condense_mb_outputs blocked = {
        {{hmb, 16}, {xcmb, 32}, {gcmb, 24}}, {{gmb, 4}, {lmb, 6}, {umb, 6}, {xlmb, 8}, {xumb, 8}}
    };
    const matmpc_condense_output blocked_scratch = {buffer, 24};
    matmpc_condense_mb_sizes sizes;

    require(matmpc_condense_matrix(&dims, &matrices, &full_matrix, full_scratch) == MATMPC_CONDENSE_OK &&
            matmpc_condense_affine(&dims, &affine, &full_affine, full_scratch) == MATMPC_CONDENSE_OK,
            "Reference full condensing failed");
    require(matmpc_condense_mb_get_sizes(&dims, 2, &sizes) == MATMPC_CONDENSE_OK && sizes.workspace == 24,
            "Scratch must fit one response column and two costate blocks");

    for (size_t partition = 0; partition < sizeof(partitions) / sizeof(partitions[0]); ++partition) {
        const size_t split = partitions[partition][1];
        double rb[8] = {0}, gb[4] = {0};
        double expected_h[16] = {0}, expected_x[32] = {0}, expected_g[24] = {0}, expected_c[4] = {0};
        matmpc_condense_mb_inputs input = {affine, {rb, 8}, {cgu, 8}};
        input.affine.gu.data = gb;
        input.affine.gu.count = 4;
        for (size_t stage = 0; stage < 4; ++stage) {
            const size_t block = stage < split ? 0 : 1;
            for (size_t i = 0; i < 4; ++i) {
                rb[4 * block + i] += r[4 * stage + i];
            }
            for (size_t i = 0; i < 2; ++i) {
                gb[2 * block + i] += gu[2 * stage + i];
            }
        }
        /* Independently sum stage rows and columns according to T. */
        for (size_t col = 0; col < 8; ++col) {
            const size_t reduced_col = (col / 2 < split ? 0 : 2) + col % 2;
            expected_c[reduced_col] += gc[col];
            for (size_t row = 0; row < 8; ++row) {
                const size_t reduced_row = (row / 2 < split ? 0 : 2) + row % 2;
                expected_h[4 * reduced_col + reduced_row] += hc[8 * col + row];
                expected_x[8 * reduced_col + row] += ccx[8 * col + row];
            }
            for (size_t row = 0; row < 6; ++row) {
                expected_g[6 * reduced_col + row] += ccg[6 * col + row];
            }
        }

        require(matmpc_condense_mb(&dims, partitions[partition], 3, &input, &blocked, blocked_scratch) == MATMPC_CONDENSE_OK,
                "Move-blocking condensing failed");
        compare(hmb, expected_h, 16);
        compare(gmb, expected_c, 4);
        compare(xcmb, expected_x, 32);
        compare(gcmb, expected_g, 24);
        compare(lmb, lo, 6);
        compare(umb, hi, 6);
        compare(xlmb, xlo, 8);
        compare(xumb, xhi, 8);
    }
}

int main(void)
{
    test_shared_input_cross_term();
    test_matches_full_projection();
    return EXIT_SUCCESS;
}
