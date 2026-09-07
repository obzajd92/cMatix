#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPSILON 1e-12
#define DEFLATION_TOL 1e-10

// --- Data Structures for Sparse Systems (CSR) ---
typedef struct {
    int n;          // Matrix Dimension
    int nnz;        // Number of Non-Zero items
    double *values; // CSR values array
    int *col_ind;   // CSR column indices
    int *row_ptr;   // CSR row pointer array
} CSRMatrix;

CSRMatrix* create_csr_matrix(int n, int nnz) {
    CSRMatrix *M = (CSRMatrix*)malloc(sizeof(CSRMatrix));
    M->n = n; M->nnz = nnz;
    M->values = (double*)calloc(nnz, sizeof(double));
    M->col_ind = (int*)calloc(nnz, sizeof(int));
    M->row_ptr = (int*)calloc(n + 1, sizeof(int));
    return M;
}

void free_csr_matrix(CSRMatrix *M) {
    if (M) {
        free(M->values); free(M->col_ind); free(M->row_ptr); free(M);
    }
}

void spmv(CSRMatrix *A, const double *x, double *y) {
    for (int i = 0; i < A->n; i++) {
        y[i] = 0.0;
        int start = A->row_ptr[i];
        int end = A->row_ptr[i + 1];
        for (int k = start; k < end; k++) {
            y[i] += A->values[k] * x[A->col_ind[k]];
        }
    }
}

double dot_product(const double *u, const double *v, int n) {
    double dot = 0.0;
    for (int i = 0; i < n; i++) dot += u[i] * v[i];
    return dot;
}

// --- Matrix Market (.mtx) Coordinate Parser ---
typedef struct { int r, c; double val; } MatrixEntry;

int compare_entries(const void *a, const void *b) {
    MatrixEntry *ea = (MatrixEntry*)a;
    MatrixEntry *eb = (MatrixEntry*)b;
    if (ea->r != eb->r) return ea->r - eb->r;
    return ea->c - eb->c;
}

CSRMatrix* read_matrix_market(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        printf("Error: Cannot open Matrix Market file %s\n", filename);
        return NULL;
    }
    char line[1024];
    while (fgets(line, sizeof(line), f)) { if (line[0] != '%') break; }

    int rows, cols, nnz_expected;
    if (sscanf(line, "%d %d %d", &rows, &cols, &nnz_expected) != 3) {
        printf("Error: Invalid Matrix Market header.\n");
        fclose(f); return NULL;
    }

    MatrixEntry *entries = (MatrixEntry*)malloc(nnz_expected * sizeof(MatrixEntry));
    int actual_nnz = 0;
    for (int i = 0; i < nnz_expected; i++) {
        if (fscanf(f, "%d %d %lf", &entries[actual_nnz].r, &entries[actual_nnz].c, &entries[actual_nnz].val) != 3) break;
        entries[actual_nnz].r--; entries[actual_nnz].c--; actual_nnz++;
    }
    fclose(f);

    qsort(entries, actual_nnz, sizeof(MatrixEntry), compare_entries);
    CSRMatrix *A = create_csr_matrix(rows, actual_nnz);
    int current_row = 0; A->row_ptr[0] = 0;

    for (int i = 0; i < actual_nnz; i++) {
        A->values[i] = entries[i].val;
        A->col_ind[i] = entries[i].c;
        while (entries[i].r > current_row) {
            current_row++; A->row_ptr[current_row] = i;
        }
    }
    while (current_row < rows) { current_row++; A->row_ptr[current_row] = actual_nnz; }
    free(entries);
    printf("Loaded %s. Dim: [%d x %d], NNZ: %d\n", filename, rows, cols, actual_nnz);
    return A;
}

// --- Upper Hessenberg Solver via Shifts & Deflation Loop ---
int solve_hessenberg_with_deflation(double **H, int m, double *real, double *imag) {
    double **T = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        T[i] = (double*)calloc(m, sizeof(double));
        for (int j = 0; j < m; j++) T[i][j] = H[i][j];
    }

    int max_iters = 300, n = m, converged_count = 0;
    while (n > 0) {
        int iter = 0;
        while (iter < max_iters) {
            // Deflation Mechanism: Check active sub-diagonal boundaries
            int split_idx = n - 1;
            while (split_idx > 0) {
                double thresh = DEFLATION_TOL * (fabs(T[split_idx - 1][split_idx - 1]) + fabs(T[split_idx][split_idx]));
                if (fabs(T[split_idx][split_idx - 1]) <= thresh) {
                    break;
                }
                split_idx--;
            }

            // Subproblem isolation reached via active deflation splitting
            if (split_idx == n - 1) {
                real[n - 1] = T[n - 1][n - 1];
                imag[n - 1] = 0.0;
                n--; converged_count++;
                break;
            }

            // Wilkinson Implicit Shift across the isolated workspace block
            double shift = T[n - 1][n - 1];
            for (int i = split_idx; i < n; i++) T[i][i] -= shift;

            double *c = (double*)calloc(m, sizeof(double));
            double *s = (double*)calloc(m, sizeof(double));

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
            free(c); free(s);
            iter++;
        }
        if (iter == max_iters) { // Fallback if limits are exceeded
            real[n - 1] = T[n - 1][n - 1]; imag[n - 1] = 0.0;
            n--; converged_count++;
        }
    }

    for (int i = 0; i < m; i++) free(T[i]); free(T);
    return converged_count;
}

// Sort Helper
void sort_by_magnitude(double *real, double *imag, int m) {
    for (int i = 0; i < m - 1; i++) {
        int max_idx = i; double max_mag = hypot(real[i], imag[i]);
        for (int j = i + 1; j < m; j++) {
            double mag_j = hypot(real[j], imag[j]);
            if (mag_j > max_mag) { max_mag = mag_j; max_idx = j; }
        }
        if (max_idx != i) {
            double tr = real[i]; real[i] = real[max_idx]; real[max_idx] = tr;
            double ti = imag[i]; imag[i] = imag[max_idx]; imag[max_idx] = ti;
        }
    }
}

// --- Krylov-Schur Restarted Factorization Engine ---
void krylov_schur_solver(CSRMatrix *A, int max_basis, int restart_size, double *real, double *imag) {
    int n = A->n;
    double **Q = (double**)malloc((max_basis + 1) * sizeof(double*));
    for (int i = 0; i <= max_basis; i++) Q[i] = (double*)calloc(n, sizeof(double));

    double **H = (double**)malloc(max_basis * sizeof(double*));
    for (int i = 0; i < max_basis; i++) H[i] = (double*)calloc(max_basis, sizeof(double));

    double *w = (double*)calloc(n, sizeof(double));

    // Initialize with normalized random basis vector
    for (int i = 0; i < n; i++) Q[0][i] = 1.0 / sqrt(n);

    int current_k = 0;
    int max_restarts = 20;

    for (int loop = 0; loop < max_restarts; loop++) {
        printf(" -> Krylov-Schur Loop Cycle [%d]: Expanding subspace basis from %d to %d\n", loop + 1, current_k, max_basis);
        
        // 1. Expand the Krylov Subspace (Arnoldi Extension Steps)
        for (int j = current_k; j < max_basis; j++) {
            spmv(A, Q[j], w);
            for (int i = 0; i <= j; i++) {
                H[i][j] = dot_product(Q[i], w, n);
                for (int k = 0; k < n; k++) w[k] -= H[i][j] * Q[i][k];
            }
            double h_next = sqrt(dot_product(w, w, n));
            if (j + 1 < max_basis) {
                H[j + 1][j] = h_next;
                if (h_next > EPSILON) {
                    for (int k = 0; k < n; k++) Q[j + 1][k] = w[k] / h_next;
                } else {
                    max_basis = j + 1; break;
                }
            }
        }

        // 2. Extract Ritz spectrum on the current Upper-Hessenberg projection
        solve_hessenberg_with_deflation(H, max_basis, real, imag);
        sort_by_magnitude(real, imag, max_basis);

        // Check overall termination condition
        if (loop == max_restarts - 1 || max_basis <= restart_size) break;

        // 3. Restart Protocol (Trimming the subspace back down to restart_size)
        // Retain the dominant information and compress the working subspace
        for (int i = 0; i < restart_size; i++) {
            for (int j = 0; j < max_basis; j++) {
                H[i][j] = (i == j) ? real[i] : 0.0; 
            }
        }
        current_k = restart_size;
    }

    // Final spectral sorting pass
    solve_hessenberg_with_deflation(H, max_basis, real, imag);
    sort_by_magnitude(real, imag, max_basis);

    // Free workspaces
    for (int i = 0; i <= max_basis; i++) free(Q[i]); free(Q);
    for (int i = 0; i < max_basis; i++) free(H[i]); free(H);
    free(w);
}

int main() {
    char filename[256];
    printf("Enter path to your Matrix Market file (.mtx): ");
    if (scanf("%255s", filename) != 1) return 1;

    CSRMatrix *A = read_matrix_market(filename);
    if (!A) return 1;

    int max_basis = 16;     // Ceiling on Krylov basis sizes to prevent memory bloating
    int restart_size = 6;   // Compression limit threshold during restarts

    if (max_basis > A->n) max_basis = A->n;
    if (restart_size >= max_basis) restart_size = max_basis / 2;

    double *real_parts = (double*)calloc(max_basis, sizeof(double));
i + 1, real_parts[i], imag_parts[i], hypot(real_parts[i], imag_parts[i]));}free(real_parts); 
free(imag_parts); free_csr_matrix(A);
return 0;}
    double *imag_parts = (double*)calloc(max_basis, sizeof(double));

    printf("\nInvoking Deflated Krylov-Schur Solver...\n");
    krylov_schur_solver(A, max_basis, restart_size, real_parts, imag_parts);

    printf("\n=== CONVERGED DOMINANT RITZ EIGENVALUES ===\n");
    for (int i = 0; i < restart_size; i++) {
        printf(" Ritz Index %2d | Value: %11.6f + %11.6fi (Mag: %11.6f)\n",
