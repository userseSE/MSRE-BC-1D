#include "matrix_LU.hpp"
#include "parameters.hpp"
// #include <cstring>
// #include <stdexcept>
// #include <cmath>
#include <Eigen/Dense>

void invert_LU(float A[N][N]) {
    // Step 1: Convert the input array to an Eigen matrix (try using dynamic-size MatrixXd)
    Eigen::MatrixXd eigen_A(N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            eigen_A(i, j) = A[i][j];
        }
    }

    // Step 2: Perform LU decomposition and calculate the inverse
    Eigen::MatrixXd inv_A = eigen_A.inverse();

    // Step 3: Copy the result back to the original array
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = inv_A(i, j);
        }
    }
}
// // Function to convert CSR matrix to Eigen sparse matrix
// void csr_multiply(const CSRMatrix &A, const float *x, float *result) {
//     memset(result, 0, N * sizeof(float)); // Initialize result vector to zero
//     for (int i = 0; i < N; ++i) {
//         for (int j = A.row_pointers[i]; j < A.row_pointers[i + 1]; ++j) {
//             result[i] += A.values[j] * x[A.col_indices[j]];
//         }
//     }
// }

// // Jacobi method for inverting sparse matrices in CSR format
// void jacobi_inverse_sparse(const CSRMatrix &A, CSRMatrix &A_inv) {
//     float x[N];           // Solution vector for each column of the inverse
//     float b[N];           // Identity column vector
//     float r[N];           // Residual
//     const int max_iter = 1000;   // Maximum iterations
//     const float tolerance = 1e-6; // Tolerance for convergence
    
//     // Iterate over each column to compute each column of the inverse
//     for (int k = 0; k < N; ++k) {
//         // Initialize b as the k-th column of the identity matrix
//         memset(b, 0, N * sizeof(float));
//         b[k] = 1.0;

//         // Initial guess x = 0
//         memset(x, 0, N * sizeof(float));

//         // Jacobi iteration
//         for (int iter = 0; iter < max_iter; ++iter) {
//             bool converged = true;

//             for (int i = 0; i < N; ++i) {
//                 float sum = 0.0;
//                 for (int j = A.row_pointers[i]; j < A.row_pointers[i + 1]; ++j) {
//                     if (A.col_indices[j] != i) {
//                         sum += A.values[j] * x[A.col_indices[j]];
//                     }
//                 }
//                 float new_xi = (b[i] - sum) / A.values[A.row_pointers[i]]; // A(i, i) is the diagonal
//                 if (std::fabs(new_xi - x[i]) > tolerance) {
//                     converged = false;
//                 }
//                 x[i] = new_xi;
//             }

//             if (converged) break;
//         }

//         // Copy the result x into the corresponding column of A_inv
//         for (int i = 0; i < N; ++i) {
//             A_inv.values[A_inv.row_pointers[i] + k] = x[i];
//         }
//     }
// }

// void inverse_sparse(const CSRMatrix &A_csr, CSRMatrix &A_inv) {
//     // Initialize A_inv with the structure of A_csr
//     std::memcpy(A_inv.row_pointers, A_csr.row_pointers, (N + 1) * sizeof(int));
//     std::memcpy(A_inv.col_indices, A_csr.col_indices, MAX_NON_ZERO * sizeof(int));

//     // Compute inverse column by column
//     jacobi_inverse_sparse(A_csr, A_inv);
// }