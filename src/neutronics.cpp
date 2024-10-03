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

using Eigen::BiCGSTAB; // Include the BiCGSTAB solver
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using namespace boost::numeric::odeint;

std::pair<std::vector<std::vector<double>>, std::vector<double>> neutronics(const std::vector<std::vector<double>> &y_n, const std::vector<double> &rho,
           int step) {

  VectorXd Keff = VectorXd::Zero(N);
  std::transform(rho.begin(), rho.end(), Keff.data(),
                 [](double r) { return (1.0 / (1.0 - r)); });
  std::cout << "Keff[N/2]: " << Keff[N / 2] << std::endl;
  std::cout << "Keff_avg: " << Keff.mean() << std::endl;

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
  D2(0, 0) = 1.0;
  D2(N - 1, N - 1) = 1.0;

  // Convert to sparse matrix format (effectively a CSC matrix)
  SparseMatrix<double> D2_sparse = D2.sparseView();
  SparseMatrix<double> I(N, N);
  I.setIdentity();

  // Define matrices A and B for Crank-Nicolson method
  SparseMatrix<double> A = I - 0.5 * dt * V * D * D2_sparse;
  SparseMatrix<double> B = I + 0.5 * dt * V * D * D2_sparse;

  // ODE system function compatible with odeint
  auto pde_to_ode_neutronics = [&](const std::vector<std::vector<double>> &y,
                                   std::vector<std::vector<double>> &dydt, double t) {
    Eigen::Map<const VectorXd> phi(y[0].data(), N); // Map the first row (neutron flux)
    VectorXd lambda_ci = VectorXd::Zero(N);

    // Compute delayed neutron precursor contributions
    for (int i = 0; i < 6; ++i) {
      lambda_ci += lambda_i[i] * Eigen::Map<const VectorXd>(y[i + 1].data(), N);
    }

    // Right-hand side for Crank-Nicolson method
    VectorXd rhs_phi =
        B * phi + dt * V * ((-sigma_a + (1.0 - Beta) * nu_sigma_f * Keff.array()).matrix().cwiseProduct(phi) +
                            lambda_ci);

    BiCGSTAB<SparseMatrix<double>> solver;
    solver.compute(A);
    VectorXd phi_new = solver.solve(rhs_phi);

    // Calculate time derivative of phi
    VectorXd dphi_dt = (phi_new - phi) / dt;

    // Populate dydt[0] with dphi_dt (neutron flux)
    std::copy(dphi_dt.data(), dphi_dt.data() + N, dydt[0].begin());

    // Populate dydt for delayed neutron precursors
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < N; ++j) {
        dydt[i + 1][j] = beta[i] * (nu_sigma_f * Keff[j]) * phi[j] - (lambda_i[i] * y[i + 1][j]);
      }
    }
  };

  // Initial condition vector (2D)
  std::vector<std::vector<double>> y0(7, std::vector<double>(N));
  if (step == 0) {
    // Initialize neutron flux and delayed neutron precursors in separate rows
    std::copy(neutronics_initial_conditions[0].begin(), neutronics_initial_conditions[0].end(), y0[0].begin());
    std::copy(neutronics_initial_conditions[1].begin(), neutronics_initial_conditions[1].end(), y0[1].begin());
    std::copy(neutronics_initial_conditions[2].begin(), neutronics_initial_conditions[2].end(), y0[2].begin());
    std::copy(neutronics_initial_conditions[3].begin(), neutronics_initial_conditions[3].end(), y0[3].begin());
    std::copy(neutronics_initial_conditions[4].begin(), neutronics_initial_conditions[4].end(), y0[4].begin());
    std::copy(neutronics_initial_conditions[5].begin(), neutronics_initial_conditions[5].end(), y0[5].begin());
    std::copy(neutronics_initial_conditions[6].begin(), neutronics_initial_conditions[6].end(), y0[6].begin());
  } else {
    y0 = y_n;
  }

  // Solve the system of ODEs using 2D vectors
  std::vector<std::vector<double>> solution_y_n = ode_solver(y0, pde_to_ode_neutronics);

  // Extract the solution at the last time step (neutron flux in the first row)
  std::vector<double> phi = solution_y_n[0];
  std::vector<double> q_prime = phi; // Assuming q_prime is phi for now

  return {solution_y_n, q_prime};
}

