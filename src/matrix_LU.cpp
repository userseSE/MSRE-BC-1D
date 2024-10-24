#include "matrix_LU.hpp"

// Function to perform forward substitution for Ly = b
void forward_substitution(double L[N][N], double y[N], double b[N]) {
    for (int i = 0; i < N; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }
}

// Function to perform back substitution for Ux = y
void back_substitution(double U[N][N], double x[N], double y[N]) {
    for (int i = N - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < N; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

void lu_decompose(double A[N][N], double L[N][N], double U[N][N]) {
    // Initialize L to be an identity matrix and U to be a zero matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                L[i][j] = 1.0;  // Diagonal of L is 1
            } else {
                L[i][j] = 0.0;
            }
            U[i][j] = 0.0;
        }
    }

    // Perform the LU decomposition
    for (int k = 0; k < N; ++k) {
        // Upper Triangular U
        for (int j = k; j < N; ++j) {
            double sum = 0.0;
            for (int p = 0; p < k; ++p) {
                sum += L[k][p] * U[p][j];
            }
            U[k][j] = A[k][j] - sum;
        }

        // Lower Triangular L
        for (int i = k + 1; i < N; ++i) {
            double sum = 0.0;
            for (int p = 0; p < k; ++p) {
                sum += L[i][p] * U[p][k];
            }
            L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }
}

// Function to invert matrix A using LU decomposition
void invert_LU(double A[N][N]){
    double L[N][N] = {0};  // Lower triangular matrix
    double U[N][N] = {0};  // Upper triangular matrix
    double inv_A[N][N] = {0};  // Inverse of A

    // Perform LU decomposition of A (you should have LU decomposition in matrix_LU.hpp)
    lu_decompose(A, L, U);  // Assumes lu_decompose(A, L, U) decomposes A into L and U

    // Temporary vectors for forward and back substitution
    double b[N];  
    double y[N];
    double x[N];

    // For each column of the identity matrix
    for (int col = 0; col < N; ++col) {
        // Set b to be the current column of the identity matrix
        for (int i = 0; i < N; ++i) {
            b[i] = (i == col) ? 1.0 : 0.0;
        }

        // Solve Ly = b using forward substitution
        forward_substitution(L, y, b);

        // Solve Ux = y using back substitution
        back_substitution(U, x, y);

        // Copy result to the corresponding column of inv_A
        for (int i = 0; i < N; ++i) {
            inv_A[i][col] = x[i];
        }
    }

    // Copy inv_A back to A
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = inv_A[i][j];
        }
    }
}
