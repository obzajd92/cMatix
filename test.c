#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Ritz.c"

#define WATERMARK_SIZE 8192
#define WATERMARK_BYTE 0xAA

void run_profiled_algorithm(int benchmark_mode, int *err_code) {
    volatile unsigned char stack_buffer[WATERMARK_SIZE];
    for (int i = 0; i < WATERMARK_SIZE; i++) stack_buffer[i] = WATERMARK_BYTE;

    reset_heap_counters();

    int n = 150;
    int nnz = n + 2 * (n - 1);
    CSRMatrix *A = create_csr_matrix(n, nnz);
    if (!A) { *err_code = -1; return; }

    int idx = 0;
    for (int i = 0; i < n; i++) {
        A->row_ptr[i] = idx;
        if (i > 0) { A->col_ind[idx] = i - 1; A->values[idx] = -1.0; idx++; }
        A->col_ind[idx] = i; A->values[idx] = 2.0; idx++;
        if (i < n - 1) { A->col_ind[idx] = i + 1; A->values[idx] = -1.0; idx++; }
    }
    A->row_ptr[n] = idx;

    int max_basis = 25;
    double **Q = (double**)malloc((max_basis + 1) * sizeof(double*));
    double **H = (double**)malloc(max_basis * sizeof(double*));
    double *w = (double*)calloc(n, sizeof(double));

    if (!Q || !H || !w) {
        free(Q); free(H); free(w); free_csr_matrix(A);
        *err_code = -2; return;
    }
    for (int i = 0; i <= max_basis; i++) Q[i] = (double*)calloc(n, sizeof(double));
    for (int i = 0; i < max_basis; i++) H[i] = (double*)calloc(max_basis, sizeof(double));

    for (int i = 0; i < n; i++) Q[i] = 1.0 / sqrt(n);

    int iterations = benchmark_mode ? 500 : 1;
    for (int run = 0; run < iterations; run++) {
        for (int j = 0; j < max_basis; j++) {
            if (spmv(A, Q[j], w) != 0) { *err_code = -3; break; }
            for (int i = 0; i <= j; i++) {
                H[i][j] = dot_product(Q[i], w, n);
                for (int k = 0; k < n; k++) w[k] -= H[i][j] * Q[i][k];
            }
            double h_next = sqrt(dot_product(w, w, n));
            if (j + 1 < max_basis && h_next > 1e-12) {
                for (int k = 0; k < n; k++) Q[j + 1][k] = w[k] / h_next;
            }
        }
    }

    int unused_bytes = 0;
    while (unused_bytes < WATERMARK_SIZE && stack_buffer[unused_bytes] == WATERMARK_BYTE) {
        unused_bytes++;
    }
    int peak_stack_used = WATERMARK_SIZE - unused_bytes;
    size_t peak_heap_used = get_peak_heap();

    if (benchmark_mode) {
        printf("PEAK_STACK_BYTES|%d\n", peak_stack_used);
        printf("PEAK_HEAP_BYTES|%zu\n", peak_heap_used);
    }

    for (int i = 0; i <= max_basis; i++) free(Q[i]); free(Q);
    for (int i = 0; i < max_basis; i++) free(H[i]); free(H);
    free(w); free_csr_matrix(A);
}

int main(int argc, char **argv) {
    int benchmark_mode = 0;
    char *flag_label = "unknown";
    
    if (argc > 1 && strcmp(argv, "--benchmark") == 0) {
        benchmark_mode = 1;
        if (argc > 2) flag_label = argv[2];
    }

    int err_code = 0;
    if (benchmark_mode) {
        clock_t start = clock();
        run_profiled_algorithm(1, &err_code);
        clock_t end = clock();
        
        if (err_code != 0) {
            printf("BENCHMARK_ERROR|%s|%d\n", flag_label, err_code);
            return 1;
        }
        double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("BENCHMARK_TIME|%s|%7.4f\n", flag_label, time_spent);
    } else {
        printf("\n=== RUNNING SECURITY TESTING BARRIER MATRIX ===\n");
        run_profiled_algorithm(0, &err_code);
        if (err_code == 0) {
            printf("[STATUS] Execution tracking completed cleanly. [PASSED]\n");
        } else {
            printf("[STATUS] Internal tracking framework crashed! Code: %d [FAILED]\n", err_code);
        }
        printf("====================================================\n\n");
    }
    return 0;
}
