#include "neutronics.hpp"
#include "ode_solver.hpp"
#include "parameters.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

// Include necessary linear algebra headers
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/numeric/odeint.hpp>

// using Eigen::BiCGSTAB; // Include the BiCGSTAB solver
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using namespace boost::numeric::odeint;

typedef Eigen::VectorXd state_type;

std::pair<Eigen::VectorXd, Eigen::VectorXd>
neutronics(const state_type &y_n, const std::vector<double> &rho, int step, double nu_sigma_f) {
    VectorXd Keff = VectorXd::Zero(N);
    std::transform(rho.data(), rho.data() + N, Keff.data(),
                 [](double r) { return (1.0 / (1.0 - r)); });
    // std::cout << "Keff[N/2]: " << Keff[N / 2] << std::endl;
    // std::cout << "Keff_avg: " << Keff.mean() << std::endl;

    // Finite difference matrix for the second derivative using Crank-Nicolson method
    VectorXd main_diag = (-2 / (dz * dz)) * VectorXd::Ones(N);
    VectorXd off_diag = (1 / (dz * dz)) * VectorXd::Ones(N - 1);
    MatrixXd D2 = MatrixXd::Zero(N, N);

    // Set the main diagonal
    D2.diagonal() = main_diag;

    // Set the superdiagonal (above the main diagonal)
    D2.diagonal(1) = off_diag;

    // Set the subdiagonal (below the main diagonal)
    D2.diagonal(-1) = off_diag;

    // Apply Dirichlet boundary conditions for zero flux at boundaries
    D2.row(0).setZero();
    D2.row(N - 1).setZero();
    D2(0, 0) = 1.0 / (dz*dz);
    D2(N - 1, N - 1) = 1.0 / (dz*dz);
    // // Reflective boundary at z = 0 (forward difference)
    // D2(0, 0) = -1 / dz;
    // D2(0, 1) = 1 / dz;

    // // Reflective boundary at z = L (backward difference)
    // D2(N - 1, N - 2) = -1 / dz;
    // D2(N - 1, N - 1) = 1 / dz;

    // Convert to sparse matrix format (this is effectively a CSC matrix)
    SparseMatrix<double> D2_sparse = D2.sparseView();

    SparseMatrix<double> I(N, N);
    I.setIdentity();

    // Define matrices A and B for Crank-Nicolson method using transposed D2_csc
    SparseMatrix<double> A = I - 0.5 * dt * V * D * D2_sparse;
    // Solve the system
    Eigen::SparseLU<SparseMatrix<double>> solver;
    solver.compute(A);
    SparseMatrix<double> B = I + 0.5 * dt * V * D * D2_sparse;

    // ODE system function compatible with odeint
    std::function<void(double, const VectorXd &, VectorXd &)> pde_to_ode_neutronics = 
    [&](double t, const VectorXd &y, VectorXd &dydt) {
        VectorXd phi = y.head(N);
        VectorXd lambda_ci = VectorXd::Zero(N);

        // Compute delayed neutron precursor contributions
        for (int i = 0; i < 6; ++i) {
            lambda_ci += lambda_i[i] * y.segment((i + 1) * N, N);
        }


        // Right-hand side for the Crank-Nicolson method
        Eigen::VectorXd rho_eigen = Eigen::VectorXd::Map(rho.data(), rho.size());

        VectorXd rhs_phi = B * phi + dt * V * ((-sigma_a + (1.0 - Beta + rho_eigen.array()) * nu_sigma_f)
                .matrix()
                .cwiseProduct(phi) + lambda_ci);

        
        VectorXd phi_new = solver.solve(rhs_phi);

        // Calculate time derivative of phi
        VectorXd dphi_dt = (phi_new - phi) / dt;

        // Compute dci/dt
        VectorXd dci_dt(6 * N);
        for (int i = 0; i < 6; ++i) {
            dci_dt.segment(i * N, N) = beta[i] * (nu_sigma_f) * phi - lambda_i[i] * y.segment((i + 1) * N, N);
        }

        // Populate dydt with the derivatives
        dydt.head(N) = dphi_dt;
        dydt.tail(6 * N) = dci_dt;
    };

    // Initial condition vector
    state_type y0(7 * N);
    if (step == 0) {
        y0 << phi_0, c1, c2, c3, c4, c5, c6;
    } else {
        y0 = y_n;
    }

    // Solve the system of ODEs
    state_type solution_y_n = ode_solver(y0, pde_to_ode_neutronics);  // Pass t0, t_end, dt

    // Extract the solution at the last time step
    VectorXd phi = solution_y_n.head(N);
    VectorXd q_prime = phi; // Assuming q_prime is phi for now

    return {solution_y_n, q_prime};
}
