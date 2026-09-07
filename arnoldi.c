#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPSILON 1e-12

// ==========================================
// 1. DATA STRUCTURES FOR SPARSE SYSTEMS (CSR)
// ==========================================
typedef struct {
    int n;          // Matrix Dimension
    int nnz;        // Number of Non-Zero items
    double *values; // CSR values array
    int *col_ind;   // CSR column indices
    int *row_ptr;   // CSR row pointer array
} CSRMatrix;

CSRMatrix* create_csr_matrix(int n, int nnz) {
    CSRMatrix *M = (CSRMatrix*)malloc(sizeof(CSRMatrix));
    M->n = n;
    M->nnz = nnz;
    M->values = (double*)calloc(nnz, sizeof(double));
    M->col_ind = (int*)calloc(nnz, sizeof(int));
    M->row_ptr = (int*)calloc(n + 1, sizeof(int));
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

// ==========================================
// 2. READ MATRIX MARKET (.MTX) COORDINATE FORMAT
// ==========================================
typedef struct {
    int r, c;
    double val;
} MatrixEntry;

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
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
    }

    int rows, cols, nnz_expected;
    if (sscanf(line, "%d %d %d", &rows, &cols, &nnz_expected) != 3) {
        printf("Error: Invalid Matrix Market size header declaration.\n");
        fclose(f);
        return NULL;
    }

    MatrixEntry *entries = (MatrixEntry*)malloc(nnz_expected * sizeof(MatrixEntry));
    int actual_nnz = 0;

    for (int i = 0; i < nnz_expected; i++) {
        if (fscanf(f, "%d %d %lf", &entries[actual_nnz].r, &entries[actual_nnz].c, &entries[actual_nnz].val) != 3) {
            break;
        }
        entries[actual_nnz].r--; // Convert 1-based indexing to 0-based
        entries[actual_nnz].c--;
        actual_nnz++;
    }
    fclose(f);

    qsort(entries, actual_nnz, sizeof(MatrixEntry), compare_entries);

    CSRMatrix *A = create_csr_matrix(rows, actual_nnz);
    int current_row = 0;
    A->row_ptr[0] = 0;

    for (int i = 0; i < actual_nnz; i++) {
        A->values[i] = entries[i].val;
        A->col_ind[i] = entries[i].c;
        while (entries[i].r > current_row) {
            current_row++;
            A->row_ptr[current_row] = i;
        }
    }
    while (current_row < rows) {
        current_row++;
        A->row_ptr[current_row] = actual_nnz;
    }

    free(entries);
    printf("Loaded %s. Size: [%d x %d], Non-Zeros: %d\n", filename, rows, cols, actual_nnz);
    return A;
}

// ==========================================
// 3. ARNOLDI FACTORIZATION METHOD
// ==========================================
void arnoldi_factorization(CSRMatrix *A, int m_steps, double **Q, double **H) {
    int n = A->n;
    double *w = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) Q[0][i] = 1.0 / sqrt(n);

    for (int j = 0; j < m_steps; j++) {
        spmv(A, Q[j], w);

        for (int i = 0; i <= j; i++) {
            H[i][j] = dot_product(Q[i], w, n);
            for (int k = 0; k < n; k++) w[k] -= H[i][j] * Q[i][k];
        }

        double h_next = sqrt(dot_product(w, w, n));
        if (j + 1 < m_steps) {
            H[j + 1][j] = h_next;
            if (h_next > EPSILON) {
                for (int k = 0; k < n; k++) Q[j + 1][k] = w[k] / h_next;
            } else {
                printf("Happy breakdown reached at Arnoldi step %d\n", j);
                free(w);
                return;
            }
        }
    }
    free(w);
}

// QR algorithm for Upper Hessenberg Matrices
void solve_hessenberg_eigenvalues(double **H, int m, double *real_parts, double *imag_parts) {
    double **T = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        T[i] = (double*)calloc(m, sizeof(double));
        for (int j = 0; j < m; j++) T[i][j] = H[i][j];
    }

    int max_iters = 500, n = m;
    while (n > 0) {
        int iter = 0;
        while (iter < max_iters) {
            int converged = 1;
            if (n > 1 && fabs(T[n - 1][n - 2]) > EPSILON) converged = 0;

            if (converged) {
                real_parts[n - 1] = T[n - 1][n - 1];
                imag_parts[n - 1] = 0.0;
                n--;
                break;
            }

            double shift = T[n - 1][n - 1];
            for (int i = 0; i < n; i++) T[i][i] -= shift;

            double *c = (double*)calloc(n, sizeof(double));
            double *s = (double*)calloc(n, sizeof(double));

            for (int i = 0; i < n - 1; i++) {
                double r = hypot(T[i][i], T[i + 1][i]);
                c[i] = (r > EPSILON) ? T[i][i] / r : 1.0;
                s[i] = (r > EPSILON) ? T[i + 1][i] / r : 0.0;

                for (int j = i; j < n; j++) {
                    double t1 = T[i][j], t2 = T[i + 1][j];
                    T[i][j] = c[i] * t1 + s[i] * t2;
                    T[i + 1][j] = -s[i] * t1 + c[i] * t2;
                }
            }

            for (int i = 0; i < n - 1; i++) {
                for (int j = 0; j <= i + 1; j++) {
                    double t1 = T[j][i], t2 = T[j][i + 1];
                    T[j][i] = c[i] * t1 + s[i] * t2;
                    T[j][i + 1] = -s[i] * t1 + c[i] * t2;
                }
            }

            for (int i = 0; i < n; i++) T[i][i] += shift;
            free(c); free(s);
            iter++;
        }
        if (iter == max_iters) {
            real_parts[n - 1] = T[n - 1][n - 1];
            imag_parts[n - 1] = 0.0;
            n--;
        }
    }
    for (int i = 0; i < m; i++) free(T[i]); free(T);
}

void sort_complex_by_magnitude(double *real, double *imag, int m) {
    for (int i = 0; i < m - 1; i++) {
        int max_idx = i;
        double mag_max = hypot(real[i], imag[i]);
        for (int j = i + 1; j < m; j++) {
            double mag_j = hypot(real[j], imag[j]);
            if (mag_j > mag_max) {
                mag_max = mag_j;
                max_idx = j;
            }
        }
        if (max_idx != i) {
            double tr = real[i]; real[i] = real[max_idx]; real[max_idx] = tr;
            double ti = imag[i]; imag[i] = imag[max_idx]; imag[max_idx] = ti;
        }
    }
}

// ==========================================
// 4. EXPORT ENGINE (CSV & JSON)
// ==========================================
void export_arnoldi_results(double *real, double *imag, int m) {
    int choice;
    printf("\nExport System? (1=CSV, 2=JSON, 3=Both, 0=No): ");
    if (scanf("%d", &choice) != 1 || choice <= 0) return;

    if (choice == 1 || choice == 3) {
        FILE *f = fopen("arnoldi_output.csv", "w");
        if (f) {
            fprintf(f, "Index,Real,Imaginary,Magnitude\n");
            for (int i = 0; i < m; i++) {
                fprintf(f, "%d,%.6f,%.6f,%.6f\n", i + 1, real[i], imag[i], hypot(real[i], imag[i]));
            }
            fclose(f);
            printf("Saved 'arnoldi_output.csv'\n");
        }
    }
    if (choice == 2 || choice == 3) {
        FILE *f = fopen("arnoldi_output.json", "w");
        if (f) {
            fprintf(f, "{\n  \"arnoldi_eigenvalues\": [\n");
            for (int i = 0; i < m; i++) {
                fprintf(f, "    {\n      \"index\": %d,\n      \"real\": %.6f,\n      \"imaginary\": %.6f,\n      \"magnitude\": %.6f\n    }%s\n",
                        i + 1, real[i], imag[i], hypot(real[i], imag[i]), (i == m - 1) ? "" : ",");
            }
            fprintf(f, "  ]\n}\n");
            fclose(f);
            printf("Saved 'arnoldi_output.json'\n");
        }
    }
}

int main() {
    char filename[256];
    printf("Enter path to your Matrix Market file (.mtx): ");
    if (scanf("%255s", filename) != 1) return 1;

    CSRMatrix *A = read_matrix_market(filename);
    if (!A) return 1;

    int m_steps = 15; 
    if (m_steps > A->n) m_steps = A->n;

    double **Q = (double**)malloc((m_steps + 1) * sizeof(double*));
    for (int i = 0; i <= m_steps; i++) Q[i] = (double*)calloc(A->n, sizeof(double));

    double **H = (double**)malloc(m_steps * sizeof(double*));
    for (int i = 0; i < m_steps; i++) H[i] = (double*)calloc(m_steps, sizeof(double));

    printf("Running Arnoldi Iteration...\n");
    arnoldi_factorization(A, m_steps, Q, H);

    double *real_parts = (double*)calloc(m_steps, sizeof(double));
    double *imag_parts = (double*)calloc(m_steps, sizeof(double));

    printf("Solving Hessenberg system via Givens QR...\n");
    solve_hessenberg_eigenvalues(H, m_steps, real_parts, imag_parts);

    sort_complex_by_magnitude(real_parts, imag_parts, m_steps);

    printf("\n=== TOP SORTED ARNOLDI RITZ EIGENVALUES (BY MAGNITUDE) ===\n");
    for (int i = 0; i < 5 && i < m_steps; i++) {
        printf(" Ritz Index %2d | Value: %9.4f + %9.4fi (Mag: %9.4f)\n",
               i + 1, real_parts[i], imag_parts[i], hypot(real_parts[i], imag_parts[i]));
    }

    export_arnoldi_results(real_parts, imag_parts, m_steps);

    for (int i = 0; i <= m_steps; i++) free(Q[i]); free(Q);
for (int i = 0; i < m_steps; i++) free(H[i]); 
free(H);
free(real_parts); 
free(imag_parts);
free_csr_matrix(A);
return 0;
}