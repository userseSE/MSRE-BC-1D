#include <Eigen/Dense>
#include "matrix_LU.hpp"

void invert_LU(double A[N][N]) {
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
