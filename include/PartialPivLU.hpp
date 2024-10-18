#ifndef PARTIAL_PIV_LU_HPP
#define PARTIAL_PIV_LU_HPP

#include <Eigen/Dense>
#include <vector>

class PartialPivLU {
public:
    // Constructor
    PartialPivLU();

    // Perform LU decomposition with partial pivoting
    bool compute(const Eigen::MatrixXd& A);

    // Solving the system Ax = b using the LU decomposition
    Eigen::VectorXd solve(const Eigen::VectorXd& b);

    // Get the LU matrix
    Eigen::MatrixXd matrixLU() const;

    // Get the permutation vector
    std::vector<int> permutationP() const;

private:
    Eigen::MatrixXd lu;  // LU decomposition matrix
    std::vector<int> perm;  // Permutation vector for partial pivoting
    bool is_initialized;
};

#endif // PARTIAL_PIV_LU_HPP
