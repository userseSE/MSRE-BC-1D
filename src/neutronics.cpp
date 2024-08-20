#include "neutronics.hpp"
#include "parameters.hpp"
#include "ode_solver.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

// Include necessary linear algebra headers
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::SparseLU;

std::pair<std::vector<double>, std::vector<double>> neutronics(const std::vector<double>& y_n, 
                                                               double rho, int step) {
    const double Keff = 1.0 / (1.0 - rho);

    // Finite difference matrix for the second derivative using Crank-Nicolson method
    VectorXd main_diag = -2 * VectorXd::Ones(N);
    VectorXd off_diag = VectorXd::Ones(N - 1);
    MatrixXd D2 = MatrixXd::Zero(N, N);
    // Set the main diagonal
    D2.diagonal() = main_diag;
    // Set the superdiagonal (above the main diagonal)
    D2.diagonal(1) = off_diag;
    // Set the subdiagonal (below the main diagonal)
    D2.diagonal(-1) = off_diag;
    // Divide by dz^2
    D2 /= (dz * dz);

    // Apply Dirichlet boundary conditions for zero flux at boundaries
    D2.row(0).setZero();
    D2.row(N - 1).setZero();
    D2(0, 0) = 1.0;
    D2(N - 1, N - 1) = 1.0;

    // Convert to sparse matrix format
    SparseMatrix<double> D2_sparse = D2.sparseView();
    SparseMatrix<double> I(N, N);
    I.setIdentity();
    
    SparseMatrix<double> A = I - 0.5 * dt * V * D * D2_sparse;
    SparseMatrix<double> B = I + 0.5 * dt * V * D * D2_sparse;
    
    // Function to solve the ODE system
    auto pde_to_ode_neutronics = [&](double t, const std::vector<double>& y) {
        VectorXd phi = Eigen::Map<const VectorXd>(y.data(), N);
        VectorXd lambda_ci = VectorXd::Zero(N);
        
        for (int i = 0; i < 6; ++i) {
            lambda_ci += lambda_i[i] * Eigen::Map<const VectorXd>(&y[(i+1)*N], N);
        }

        VectorXd rhs_phi = B * phi + dt * V * ((-sigma_a + (1.0 - Beta) * nu_sigma_f / Keff) * phi + lambda_ci);
        // std::cout << "Matrix B size: " << B.rows() << "x" << B.cols() << std::endl;
        // std::cout << "Vector phi size: " << phi.size() << std::endl;
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(A);
        VectorXd phi_new = solver.solve(rhs_phi);
        
        VectorXd dphi_dt = (phi_new - phi) / dt;
        
        std::vector<double> dci_dt(6 * N);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < N; ++j) {
                dci_dt[i * N + j] = beta[i] * (nu_sigma_f / Keff) * phi[j] - lambda_i[i] * y[(i+1)*N + j];
            }
        }
        
        std::vector<double> result(dphi_dt.data(), dphi_dt.data() + dphi_dt.size());
        result.insert(result.end(), dci_dt.begin(), dci_dt.end());

        return result;
    };
    
    // Initial condition vector
    std::vector<double> y0;
    if (step == 0) {
        y0.reserve(7 * N);
        y0.insert(y0.end(), phi_0, phi_0 + N);
        for (int i = 0; i < 6; ++i) {
            y0.insert(y0.end(), c0, c0 + N);
        }
    } else {
        y0 = y_n;
    }
    // std::cout << y0.back() << std::endl;
    // Solve the system of ODEs
    std::vector<double> solution_y_n = ode_solver(y0, {}, pde_to_ode_neutronics);
    // std::cout << solution_y_n.size() << std::endl;
    // Extract the solution at the last time step
    std::vector<double> phi(solution_y_n.begin(), solution_y_n.begin() + N);
    std::vector<double> q_prime = phi; // Assuming q_prime is phi for now

    return {solution_y_n, q_prime};
}
