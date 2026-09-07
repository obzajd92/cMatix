#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPSILON 1e-12
#define DEFLATION_TOL 1e-10

// --- HEAP TRACKING FRAMEWORK INTERNALS ---
static size_t current_heap_allocated = 0;
static size_t peak_heap_allocated = 0;

void* tracked_malloc(size_t size) {
    if (size == 0) return NULL;
    // Prefix each block with its size metadata
    size_t total_size = size + sizeof(size_t);
    void *ptr = malloc(total_size);
    if (!ptr) return NULL;

    *(size_t*)ptr = size;
    current_heap_allocated += size;
    if (current_heap_allocated > peak_heap_allocated) {
        peak_heap_allocated = current_heap_allocated;
    }
    return (void*)((char*)ptr + sizeof(size_t));
}

void* tracked_calloc(size_t num, size_t size) {
    size_t total_bytes = num * size;
    void *ptr = tracked_malloc(total_bytes);
    if (ptr) {
        memset(ptr, 0, total_bytes);
    }
    return ptr;
}

void tracked_free(void *ptr) {
    if (!ptr) return;
    void *base_ptr = (void*)((char*)ptr - sizeof(size_t));
    size_t size = *(size_t*)base_ptr;
    current_heap_allocated -= size;
    free(base_ptr);
}

size_t get_peak_heap() {
    return peak_heap_allocated;
}

void reset_heap_counters() {
    current_heap_allocated = 0;
    peak_heap_allocated = 0;
}

// Override standard memory allocation routines for all downstream code paths
#define malloc(x) tracked_malloc(x)
#define calloc(x,y) tracked_calloc(x,y)
#define free(x) tracked_free(x)

// --- Data Structures for Sparse Systems (CSR) ---
typedef struct {
    int n;          
    int nnz;        
    double *values; 
    int *col_ind;   
    int *row_ptr;   
} CSRMatrix;

CSRMatrix* create_csr_matrix(int n, int nnz) {
    if (n <= 0 || nnz <= 0) return NULL;
    CSRMatrix *M = (CSRMatrix*)malloc(sizeof(CSRMatrix));
    if (!M) return NULL;
    
    M->n = n; 
    M->nnz = nnz;
    M->values = (double*)calloc(nnz, sizeof(double));
    M->col_ind = (int*)calloc(nnz, sizeof(int));
    M->row_ptr = (int*)calloc(n + 1, sizeof(int));
    
    if (!M->values || !M->col_ind || !M->row_ptr) {
        if (M->values) free(M->values);
        if (M->col_ind) free(M->col_ind);
        if (M->row_ptr) free(M->row_ptr);
        free(M);
        return NULL;
    }
    return M;
}

void free_csr_matrix(CSRMatrix *M) {
    if (M) {
        free(M->values); 
        free(M->col_ind); 
        free(M->row_ptr); 
        free(M);
    }
}

int spmv(CSRMatrix *A, const double *x, double *y) {
    if (!A || !x || !y) return -1;
    for (int i = 0; i < A->n; i++) {
        y[i] = 0.0;
        int start = A->row_ptr[i];
        int end = A->row_ptr[i + 1];
        if (start < 0 || end > A->nnz || start > end) return -2;
        for (int k = start; k < end; k++) {
            y[i] += A->values[k] * x[A->col_ind[k]];
        }
    }
    return 0;
}

double dot_product(const double *u, const double *v, int n) {
    if (!u || !v || n <= 0) return 0.0;
    double dot = 0.0;
    for (int i = 0; i < n; i++) dot += u[i] * v[i];
    return dot;
}

int solve_hessenberg_with_deflation(double **H, int m, double *real, double *imag) {
    if (!H || m <= 0 || !real || !imag) return -1;

    double **T = (double**)malloc(m * sizeof(double*));
    if (!T) return -2;
    for (int i = 0; i < m; i++) {
        T[i] = (double*)calloc(m, sizeof(double));
        if (!T[i]) {
            for (int k = 0; k < i; k++) free(T[k]);
            free(T); return -2;
        }
        for (int j = 0; j < m; j++) T[i][j] = H[i][j];
    }

    int max_iters = 300, n = m;
    while (n > 0) {
        int iter = 0;
        while (iter < max_iters) {
            int split_idx = n - 1;
            while (split_idx > 0) {
                double thresh = DEFLATION_TOL * (fabs(T[split_idx - 1][split_idx - 1]) + fabs(T[split_idx][split_idx]));
                if (fabs(T[split_idx][split_idx - 1]) <= thresh) break;
                split_idx--;
            }

            if (split_idx == n - 1) {
                real[n - 1] = T[n - 1][n - 1]; imag[n - 1] = 0.0;
                n--; break;
            }

            double shift = T[n - 1][n - 1];
            for (int i = split_idx; i < n; i++) T[i][i] -= shift;

            double *c = (double*)calloc(m, sizeof(double));
            double *s = (double*)calloc(m, sizeof(double));
            if (!c || !s) {
                free(c); free(s);
                for (int k = 0; k < m; k++) free(T[k]); free(T); return -3;
            }

            for (int i = split_idx; i < n - 1; i++) {
                double r = hypot(T[i][i], T[i + 1][i]);
                c[i] = (r > EPSILON) ? T[i][i] / r : 1.0;
                s[i] = (r > EPSILON) ? T[i + 1][i] / r : 0.0;
                for (int j = i; j < n; j++) {
                    double t1 = T[i][j], t2 = T[i + 1][j];
                    T[i][j] = c[i] * t1 + s[i] * t2;
                    T[i + 1][j] = -s[i] * t1 + c[i] * t2;
                }
            }
            for (int i = split_idx; i < n - 1; i++) {
                for (int j = 0; j <= i + 1; j++) {
                    double t1 = T[j][i], t2 = T[j][i + 1];
                    T[j][i] = c[i] * t1 + s[i] * t2;
                    T[j][i + 1] = -s[i] * t1 + c[i] * t2;
                }
            }
            for (int i = split_idx; i < n; i++) T[i][i] += shift;
            free(c); free(s); iter++;
        }
        if (iter == max_iters) {
            real[n - 1] = T[n - 1][n - 1]; imag[n - 1] = 0.0; n--;
        }
    }
    for (int i = 0; i < m; i++) free(T[i]); free(T);
    return 0;
}
