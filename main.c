#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1e-9

// --- Dynamic Memory Helpers ---
double** allocate_matrix(int rows, int cols) {
    double **M = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        M[i] = (double*)calloc(cols, sizeof(double));
    }
    return M;
}

void free_matrix(double **M, int rows) {
    for (int i = 0; i < rows; i++) free(M[i]);
    free(M);
}

void print_matrix(const char *name, double **M, int rows, int cols) {
    printf("--- %s ---\n", name);
    for (int i = 0; i < rows; i++) {
        printf("  [");
        for (int j = 0; j < cols; j++) {
            printf(" %10.4f ", M[i][j]);
        }
        printf("]\n");
    }
    printf("\n");
}

int is_symmetric(double **M, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (fabs(M[i][j] - M[j][i]) > 1e-6) return 0; 
        }
    }
    return 1; 
}

void matrix_multiply(double **A, double **B, double **C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < n; k++) C[i][j] += A[i][k] * B[k][j];
        }
    }
}

void matrix_transpose(double **A, double **B, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) B[j][i] = A[i][j];
    }
}

void draw_decay_chart(double *errors, int total_iters) {
    printf("\n      === LOGARITHMIC OFF-DIAGONAL DECAY VISUALIZATION ===\n");
    printf(" Iter | Error Level (Logarithmic Scale Lower to the Left)\n");
    printf("------|---------------------------------------------------\n");
    for (int i = 0; i < total_iters; i++) {
        double err = errors[i];
        if (err < 1e-15) err = 1e-15; 
        double log_val = log10(err);
        int bar_length = (int)((log_val + 10.0) * 3.5); 
        if (bar_length < 1) bar_length = 1;
        if (bar_length > 45) bar_length = 45;
        printf(" %3d  | ", i + 1);
        for (int b = 0; b  max_val) {
                    max_val = fabs(D[i][j]);
                    p = i; q = j;
                }
            }
        }

        error_history[iter - 1] = max_val;
        for (int i = 0; i < n; i++) {
            eigenvalue_history[iter - 1][i] = non_symmetric_flag ? sqrt(fmax(0.0, D[i][i])) : D[i][i];
            for(int j = 0; j < n; j++) eigenvector_history[iter - 1][i][j] = V[i][j];
        }
        actual_iters = iter;

        if (max_val < EPSILON) {
            printf("Convergence reached at iteration %d!\n", iter - 1);
            break;
        }

        double theta = (D[q][q] - D[p][p]) / (2.0 * D[p][q]);
        double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(1.0 + theta * theta));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = c * t;

        double d_pp = D[p][p];
        double d_qq = D[q][q];
        D[p][p] = c * c * d_pp - 2.0 * s * c * D[p][q] + s * s * d_qq;
        D[q][q] = s * s * d_pp + 2.0 * s * c * D[p][q] + c * c * d_qq;
        D[p][q] = 0.0;
        D[q][p] = 0.0;

        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                double d_ip = D[i][p];
                double d_iq = D[i][q];
                D[i][p] = c * d_ip - s * d_iq;
                D[p][i] = D[i][p];
                D[i][q] = s * d_ip + c * d_iq;
                D[q][i] = D[i][q];
            }
        }

        // Simultaneously rotate and update complete coordinate spaces inside V
        for (int i = 0; i < n; i++) {
            double v_ip = V[i][p];
            double v_iq = V[i][q];
            V[i][p] = c * v_ip - s * v_iq;
            V[i][q] = s * v_ip + c * v_iq;
        }
    }
    
    draw_decay_chart(error_history, actual_iters);
    print_matrix("Final Computed Eigenvectors (As Columns)", V, n, n);

    // Prompt user to export logs to standard data files
    int choice;
    printf("Export logs? (1=CSV, 2=JSON, 3=Both, 0=No): ");
    if (scanf("%d", &choice) == 1 && choice > 0) {
        if (choice == 1 || choice == 3) {
            FILE *f = fopen("tracking_history.csv", "w");
            if (f) {
                fprintf(f, "Iteration,Max_OffDiag_Error");
                for (int i = 0; i < n; i++) fprintf(f, ",Eigenvalue_%d", i + 1);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) fprintf(f, ",Eigenvector_V[%d][%d]", i, j);
                fprintf(f, "\n");

                for (int it = 0; it < actual_iters; it++) {
                    fprintf(f, "%d,%.6e", it + 1, error_history[it]);
                    for (int i = 0; i < n; i++) fprintf(f, ",%.6f", eigenvalue_history[it][i]);
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++) fprintf(f, ",%.6f", eigenvector_history[it][i][j]);
                    fprintf(f, "\n");
                }
                fclose(f);
                printf("Written to 'tracking_history.csv'\n");
            }
        }
        if (choice == 2 || choice == 3) {
            FILE *f = fopen("tracking_history.json", "w");
            if (f) {
                fprintf(f, "{\n  \"iterations\": [\n");
                for (int it = 0; it < actual_iters; it++) {
                    fprintf(f, "    {\n      \"step\": %d,\n", it + 1);
                    fprintf(f, "      \"max_off_diagonal_error\": %.6e,\n", error_history[it]);
                    fprintf(f, "      \"eigenvalues\": [");
                    for (int i = 0; i < n; i++) fprintf(f, "%.6f%s", eigenvalue_history[it][i], (i == n - 1) ? "" : ", ");
                    fprintf(f, "],\n      \"eigenvectors\": [\n");
                    for (int i = 0; i < n; i++) {
                        fprintf(f, "        [");
                        for (int j = 0; j < n; j++) fprintf(f, "%.6f%s", eigenvector_history[it][i][j], (j == n - 1) ? "" : ", ");
                        fprintf(f, "]%s\n", (i == n - 1) ? "" : ",");
                    }
                    fprintf(f, "      ]\n    }%s\n", (it == actual_iters - 1) ? "" : ",");
                }
                fprintf(f, "  ]\n}\n");
                fclose(f);
                printf("Written to 'tracking_history.json'\n");
            }
        }
    }

    free(error_history);
    free_matrix(eigenvalue_history, max_iter);
    for(int i=0; i<max_iter; i++) free_matrix(eigenvector_history[i], n);
    free(eigenvector_history);
    free_matrix(D, n);
    free_matrix(V, n);
}

int main() {
    int n = 3;
    double **A = allocate_matrix(n, n);
    // Hardcoded test configuration for a 3x3 scenario matrix A
    A[0][0] = 2.0; A[0][1] = 1.0; A[0][2] = 0.0;
    A[1][0] = 1.0; A[1][1] = 3.0; A[1][2] = 1.0;
    A[2][0] = 0.0; A[2][1] = 1.0; A[2][2] = 4.0;

    print_matrix("Input Scenario Matrix A", A, n, n);
    track_eigenvalues_and_vectors(A, n);

    free_matrix(A, n);
    return 0;
}
