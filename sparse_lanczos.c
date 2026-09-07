#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1e-9

// ==========================================
// 1. DATA STRUCTURES FOR SPARSE SYSTEMS
// ==========================================
typedef struct {
    int n;          // Matrix Dimension
    int nnz;        // Number of Non-Zero items
    double *values; // CSR values array
    int *col_ind;   // CSR column indices
    int *row_ptr;   // CSR row pointer array
} CSRMatrix;

// Allocate memory for the CSR structure
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

// Helper vector methods
double dot_product(const double *u, const double *v, int n) {
    double dot = 0.0;
    for (int i = 0; i = 0.0 ? r : -r));

            double s = 1.0, c = 1.0, p = 0.0;
            for (int i = m_idx - 1; i >= l; i--) {
                double f = s * e[i];
                double b = c * e[i];
                r = hypot(f, g);
                e[i + 1] = r;
                if (r == 0.0) {
                    d[i + 1] -= p;
                    e[m_idx] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;

                // Accumulate transformation to update eigenvectors matrix
                for (int k = 0; k < m; k++) {
                    double tmp = evecs_tri[k][i + 1];
                    evecs_tri[k][i + 1] = s * evecs_tri[k][i] + c * tmp;
                    evecs_tri[k][i] = c * evecs_tri[k][i] - s * tmp;
                }
            }
            if (r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l] = g;
            e[m_idx] = 0.0;
        }
    }

    for (int i = 0; i < m; i++) evals[i] = d[i];
    free(d); free(e);
}

// ==========================================
// 3. SORTING ALGORITHM (BY MAGNITUDE)
// ==========================================
void sort_by_magnitude(double *evals, double **evecs, int n, int k) {
    for (int i = 0; i < k - 1; i++) {
        int max_idx = i;
        for (int j = i + 1; j < k; j++) {
            if (fabs(evals[j]) > fabs(evals[max_idx])) {
                max_idx = j;
            }
        }
        if (max_idx != i) {
            // Swap Eigenvalues
            double temp_val = evals[i];
            evals[i] = evals[max_idx];
            evals[max_idx] = temp_val;
            
            // Swap Column Eigenvectors
            for (int r = 0; r < n; r++) {
                double temp_vec = evecs[r][i];
                evecs[r][i] = evecs[r][max_idx];
                evecs[r][max_idx] = temp_vec;
            }
        }
    }
}

// ==========================================
// 4. LARGE SPARSE LANCZOS EIGEN ENGINE
// ==========================================
void sparse_lanczos_solver(CSRMatrix *A, int m_steps, double *final_evals, double **final_evecs) {
    int n = A->n;
    if (m_steps > n) m_steps = n;

    // Allocate Krylov subspace basis array storage blocks
    double **Q = (double**)malloc((m_steps + 1) * sizeof(double*));
    for (int i = 0; i <= m_steps; i++) Q[i] = (double*)calloc(n, sizeof(double));

    double *alpha = (double*)calloc(m_steps, sizeof(double));
    double *beta = (double*)calloc(m_steps, sizeof(double));
    double *w = (double*)calloc(n, sizeof(double));

    // 1. Initial random vector normalized to unit length
    for (int i = 0; i < n; i++) Q[0][i] = 1.0 / sqrt(n);

    // 2. Lanczos Iteration Loop
    int converged_steps = m_steps;
    for (int j = 0; j < m_steps; j++) {
        spmv(A, Q[j], w); // w = A * q_j

        if (j > 0) {
            for (int i = 0; i < n; i++) w[i] -= beta[j - 1] * Q[j - 1][i];
        }

        alpha[j] = dot_product(w, Q[j], n);
        for (int i = 0; i < n; i++) w[i] -= alpha[j] * Q[j][i];

        // Gram-Schmidt Reorthogonalization to maintain numerical precision
        for (int k = 0; k <= j; k++) {
            double dot = dot_product(w, Q[k], n);
            for (int i = 0; i < n; i++) w[i] -= dot * Q[k][i];
        }

        double next_beta = sqrt(dot_product(w, w, n));
        if (next_beta < EPSILON) {
            converged_steps = j + 1;
            break;
        }

        if (j < m_steps - 1) {
            beta[j] = next_beta;
            for (int i = 0; i < n; i++) Q[j + 1][i] = w[i] / beta[j];
        }
    }

    // 3. Solve smaller Tridiagonal Subsystem
    double **evecs_tri = (double**)malloc(converged_steps * sizeof(double*));
    for (int i = 0; i < converged_steps; i++) evecs_tri[i] = (double*)calloc(converged_steps, sizeof(double));
    
    double *tri_evals = (double*)calloc(converged_steps, sizeof(double));
    solve_tridiagonal_eigen(alpha, beta, converged_steps, tri_evals, evecs_tri);

    // 4. Compute original space Eigenvectors: V = Q * V_tri
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < converged_steps; j++) {
            final_evecs[i][j] = 0.0;
            for (int k = 0; k < converged_steps; k++) {
                final_evecs[i][j] += Q[k][i] * evecs_tri[k][j];
            }
        }
    }
    for (int j = 0; j < converged_steps; j++) final_evals[j] = tri_evals[j];

    // 5. Apply Magnitude Sorting Options
    sort_by_magnitude(final_evals, final_evecs, n, converged_steps);

    // Free workspaces
    for (int i = 0; i <= m_steps; i++) free(Q[i]); free(Q);
    for (int i = 0; i < converged_steps; i++) free(evecs_tri[i]); free(evecs_tri);
    free(alpha); free(beta); free(w); free(tri_evals);
}

// ==========================================
// 5. EXPORT UTILITIES (CSV & JSON)
// ==========================================
void export_results(double *evals, double **evecs, int n, int k) {
    int choice;
    printf("Export Output Configuration? (1=CSV, 2=JSON, 3=Both, 0=No): ");
    if (scanf("%d", &choice) != 1 || choice <= 0) return;

    if (choice == 1 || choice == 3) {
        FILE *f = fopen("sparse_output.csv", "w");
        if (f) {
            fprintf(f, "Eigen-Index,Eigenvalue");
            for (int i = 0; i < n; i++) fprintf(f, ",Component_%d", i + 1);
            fprintf(f, "\n");
            for (int j = 0; j < k; j++) {
                fprintf(f, "%d,%.6f", j + 1, evals[j]);
                for (int i = 0; i < n; i++) fprintf(f, ",%.6f", evecs[i][j]);
                fprintf(f, "\n");
            }
            fclose(f);
            printf("Saved sparse_output.csv\n");
        }
    }
    if (choice == 2 || choice == 3) {
        FILE *f = fopen("sparse_output.json", "w");
        if (f) {
            fprintf(f, "{\n  \"sorted_eigen_system\": [\n");
            for (int j = 0; j < k; j++) {
                fprintf(f, "    {\n      \"index\": %d,\n      \"eigenvalue\": %.6f,\n      \"eigenvector\": [", j + 1, evals[j]);
                for (int i = 0; i < n; i++) fprintf(f, "%.6f%s", evecs[i][j], (i == n - 1) ? "" : ", ");
                fprintf(f, "]\n    }%s\n", (j == k - 1) ? "" : ",");
            }
            fprintf(f, "  ]\n}\n");
            fclose(f);
            printf("Saved sparse_output.json\n");
        }
    }
}

// ==========================================
// MAIN TESTING FRAMEWORK
// ==========================================
int main() {
    // Generate a sparse tri-diagonal matrix size N=1000 dynamically to avoid dense allocation
    int n = 1000; 
    int nnz = n + 2 * (n - 1); // 3 non-zeros per inner row layout
    CSRMatrix *A = create_csr_matrix(n, nnz);

    int idx = 0;
    for (int i = 0; i < n; i++) {
        A->row_ptr[i] = idx;
        if (i > 0) {
            A->col_ind[idx] = i - 1;
            A->values[idx] = -1.0;
            idx++;
        }
        A->col_ind[idx] = i;
A->values[idx] = 2.0; // Dominant diagonal elementsidx++;if (i < n - 1) {A->col_ind[idx] = i + 1;A->values[idx] = -1.0;idx++;}}A->row_ptr[n] = idx;int steps = 15; // Track top 15 dominant eigenvaluesdouble evals = (double)calloc(steps, sizeof(double));double evecs = (double)malloc(n * sizeof(double*));for (int i = 0; i < n; i++) evecs[i] = (double*)calloc(steps, sizeof(double));printf("Executing Lanczos Sparse Eigendecomposition Solver on N = %d Matrix...\n", n);sparse_lanczos_solver(A, steps, evals, evecs);printf("\n=== TOP SORTED DOMINANT EIGENVALUES (BY MAGNITUDE) ===\n");for (int i = 0; i < 5; i++) {printf(" Index %2d | Eigenvalue: %10.6f\n", i + 1, evals[i]);}export_results(evals, evecs, n, steps);// Cleanupfree(evals);for (int i = 0; i < n; i++) free(evecs[i]); free(evecs);free_csr_matrix(A);return 0;}