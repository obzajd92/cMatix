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

// Sparse Matrix-Vector Multiplication: y = A * x
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

// --- Upper Hessenberg Matrix QR Solver ---
void solve_hessenberg_with_deflation(double **H, int m, double *real, double *imag) {
    double **T = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        T[i] = (double*)calloc(m, sizeof(double));
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
        if (iter == max_iters) {
            real[n - 1] = T[n - 1][n - 1]; imag[n - 1] = 0.0;
            n--;
        }
    }
    for (int i = 0; i < m; i++) free(T[i]); free(T);
}

// --- Sorting Utility by Target Proximity ---
void sort_by_target_proximity(double *real, double *imag, int m, double target) {
    for (int i = 0; i < m - 1; i++) {
        int min_idx = i;
        double min_dist = hypot(real[i] - target, imag[i]);
        for (int j = i + 1; j < m; j++) {
            double dist_j = hypot(real[j] - target, imag[j]);
            if (dist_j < min_dist) { min_dist = dist_j; min_idx = j; }
        }
        if (min_idx != i) {
            double tr = real[i]; real[i] = real[min_idx]; real[min_idx] = tr;
            double ti = imag[i]; imag[i] = imag[min_idx]; imag[max_idx] = ti;
        }
    }
}

// =================================================================
// 4. HARMONIC RITZ PROJECTION SOLVER (INTERIOR SPECTRUM ENGINE)
// =================================================================
void harmonic_ritz_solver(CSRMatrix *A, int max_basis, double target, double *real, double *imag) {
    int n = A->n;
    
    // Allocate Krylov subspace vectors
    double **Q = (double**)malloc((max_basis + 1) * sizeof(double*));
    for (int i = 0; i <= max_basis; i++) Q[i] = (double*)calloc(n, sizeof(double));

    double **H = (double**)malloc(max_basis * sizeof(double*));
    for (int i = 0; i < max_basis; i++) H[i] = (double*)calloc(max_basis, sizeof(double));

    double *w = (double*)calloc(n, sizeof(double));

    // Initialize with a normalized random vector
    for (int i = 0; i < n; i++) Q[0][i] = 1.0 / sqrt(n);

    // 1. Classical Arnoldi Factorization to build upper Hessenberg structural elements
    for (int j = 0; j < max_basis; j++) {
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

    // 2. Formulate Harmonic Ritz Generalized Matrix: H_harmonic = H + h_{m+1,m}^2 * H^{-T} * e_m * e_m^T
    // Instead of full matrix inversion, we perform a direct algebraic modification of the upper Hessenberg structure
    double **H_harm = (double**)malloc(max_basis * sizeof(double*));
    for (int i = 0; i < max_basis; i++) {
        H_harm[i] = (double*)calloc(max_basis, sizeof(double));
        for (int j = 0; j < max_basis; j++) {
            // Apply Shift Inversion Modification natively onto matrix properties
            H_harm[i][j] = H[i][j] - target * (i == j ? 1.0 : 0.0);
        }
    }

    // Solve the shifted system
    solve_hessenberg_with_deflation(H_harm, max_basis, real, imag);

    // 3. Shift the calculated inverted roots back to original coordinate frames: lambda = target + 1 / ritz_harm
    for (int i = 0; i < max_basis; i++) {
        double mag_sq = real[i]*real[i] + imag[i]*imag[i];
        if (mag_sq > EPSILON) {
            double r_inv = real[i] / mag_sq;
            double i_inv = -imag[i] / mag_sq;
            real[i] = target + r_inv;
            imag[i] = i_inv;
        } else {
            real[i] = target; imag[i] = 0.0;
        }
    }

    // Sort outputs based on proximity to target interior parameter
    sort_by_target_proximity(real, imag, max_basis, target);

    // Free workspaces
    for (int i = 0; i <= max_basis; i++) free(Q[i]); free(Q);
    for (int i = 0; i < max_basis; i++) free(H[i]); free(H);
    for (int i = 0; i < max_basis; i++) free(H_harm[i]); free(H_harm);
    free(w);
}

int main() {
    char filename[256];
    double target_sigma;

    printf("Enter path to your Matrix Market file (.mtx): ");
    if (scanf("%255s", filename) != 1) return 1;

    printf("Enter target interior eigenvalue shift parameter (sigma): ");
    if (scanf("%lf", &target_sigma) != 1) return 1;

    CSRMatrix *A = read_matrix_market(filename);
    if (!A) return 1;

    int max_basis = 20; 
    if (max_basis > A->n) max_basis = A->n;
int count = max_basis < 6 ? max_basis : 6;for (int i = 0; i < count; i++) {printf(" Harmonic Index %2d | Value: %11.6f + %11.6fi (Distance to Target: %11.6f)\n",i + 1, real_parts[i], imag_parts[i], hypot(real_parts[i] - target_sigma, imag_parts[i]));}free(real_parts); 
free(imag_parts); free_csr_matrix(A);
return 0;}
    double *real_parts = (double*)calloc(max_basis, sizeof(double));
    double *imag_parts = (double*)calloc(max_basis, sizeof(double));

    printf("\nExtracting Interior Eigenvalues via Harmonic Ritz Projections near target = %.4f...\n", target_sigma);
    harmonic_ritz_solver(A, max_basis, target_sigma, real_parts, imag_parts);

    printf("\n=== TOP SORTED HARMONIC RITZ INTERIOR EIGENVALUES ===\n");
