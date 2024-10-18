#include "PartialPivLU.hpp"
#include <stdexcept>
#include <iostream>

// Constructor
PartialPivLU::PartialPivLU() : is_initialized(false) {}

// Perform LU decomposition with partial pivoting
bool PartialPivLU::compute(const Eigen::MatrixXd& A) {
    int N = A.rows();
    lu = A;  // Copy A to lu for in-place modification
    perm.resize(N);  // Permutation vector to track row swaps

    for (int i = 0; i < N; ++i) {
        perm[i] = i;  // Initialize permutation vector
    }

    for (int k = 0; k < N; ++k) {
        // Find the pivot: largest absolute value in column k
        int max_row = k;
        for (int i = k + 1; i < N; ++i) {
            if (std::abs(lu(i, k)) > std::abs(lu(max_row, k))) {
                max_row = i;
            }
        }

        // Check for singularity
        if (lu(max_row, k) == 0) {
            std::cerr << "Matrix is singular!" << std::endl;
            return false;
        }

        // Swap rows in LU matrix and track the permutation
        if (max_row != k) {
            lu.row(k).swap(lu.row(max_row));
            std::swap(perm[k], perm[max_row]);
        }

        // Perform Gaussian elimination
        for (int i = k + 1; i < N; ++i) {
            lu(i, k) /= lu(k, k);  // Compute the multiplier
            for (int j = k + 1; j < N; ++j) {
                lu(i, j) -= lu(i, k) * lu(k, j);
            }
        }
    }

    is_initialized = true;
    return true;
}

// Solving the system Ax = b using the LU decomposition
Eigen::VectorXd PartialPivLU::solve(const Eigen::VectorXd& b) {
    if (!is_initialized) {
        throw std::runtime_error("LU decomposition has not been initialized.");
    }

    int N = b.size();
    Eigen::VectorXd x = b;

    // Apply the permutation to the right-hand side
    for (int i = 0; i < N; ++i) {
        x[i] = b[perm[i]];
    }

    // Forward substitution to solve Ly = Pb
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < i; ++j) {
            x[i] -= lu(i, j) * x[j];
        }
    }

    // Backward substitution to solve Ux = y
    for (int i = N - 1; i >= 0; --i) {
        for (int j = i + 1; j < N; ++j) {
            x[i] -= lu(i, j) * x[j];
        }
        x[i] /= lu(i, i);
    }

    return x;
}

// Get the LU matrix
Eigen::MatrixXd PartialPivLU::matrixLU() const {
    if (!is_initialized) {
        throw std::runtime_error("LU decomposition has not been initialized.");
    }
    return lu;
}

// Get the permutation vector
std::vector<int> PartialPivLU::permutationP() const {
    if (!is_initialized) {
        throw std::runtime_error("LU decomposition has not been initialized.");
    }
    return perm;
}
