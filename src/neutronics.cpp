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

std::pair<std::vector<double>, std::vector<double>>
neutronics(const std::vector<double> &y_n, const std::vector<double> &rho,
           int step) {
  VectorXd Keff = VectorXd::Zero(N);
  std::transform(rho.begin(), rho.end(), Keff.data(),
                 [](double r) { return (1.0 /(1.0 - r)); });
  std::cout << "Keff[N/2]: " << Keff[N / 2] << std::endl;
  std::cout<<"Keff_avg: "<<Keff.mean()<<std::endl;

  // Finite difference matrix for the second derivative using Crank-Nicolson
  // method
  VectorXd main_diag = (-2 / (dz * dz)) * VectorXd::Ones(N);
  VectorXd off_diag = (1 / (dz * dz)) * VectorXd::Ones(N - 1);
  MatrixXd D2 = MatrixXd::Zero(N, N);

  // Set the main diagonal
  D2.diagonal() = main_diag;

  // Set the superdiagonal (above the main diagonal)
  D2.diagonal(1) = off_diag;

  // Set the subdiagonal (below the main diagonal)
  D2.diagonal(-1) = off_diag;

  // Divide by dz^2
  //   D2 /= (dz * dz);

  // Apply Dirichlet boundary conditions for zero flux at boundaries
  D2.row(0).setZero();
  D2.row(N - 1).setZero();
  D2(0, 0) = 1.0;
  D2(N - 1, N - 1) = 1.0;

  // Convert to sparse matrix format (this is effectively a CSC matrix)
  SparseMatrix<double> D2_sparse = D2.sparseView();

  // If you need the matrix in CSC format explicitly, transpose it
  SparseMatrix<double> D2_csc = D2_sparse.transpose();

  SparseMatrix<double> I(N, N);
  I.setIdentity();

  // Define matrices A and B for Crank-Nicolson method using transposed D2_csc
  SparseMatrix<double> A = I - 0.5 * dt * V * D * D2_sparse;
  SparseMatrix<double> B = I + 0.5 * dt * V * D * D2_sparse;

  // ODE system function compatible with odeint
  auto pde_to_ode_neutronics = [&](const std::vector<double> &y,
                                   std::vector<double> &dydt, double t) {
    Eigen::Map<const VectorXd> phi(y.data(), N);
    VectorXd lambda_ci = VectorXd::Zero(N);

    // Compute delayed neutron precursor contributions
    for (int i = 0; i < 6; ++i) {
      lambda_ci += lambda_i[i] * Eigen::Map<const VectorXd>(&y[(i + 1) * N], N);
    }

    // Right-hand side for the Crank-Nicolson method
    VectorXd rhs_phi =
        B * phi + dt * V *
                      ((-sigma_a + (1.0 - Beta) * nu_sigma_f * Keff.array())
                           .matrix()
                           .cwiseProduct(phi) +
                       lambda_ci);

    BiCGSTAB<SparseMatrix<double>> solver;
    solver.compute(A);
    // solver.setTolerance(1e-15); // Stricter tolerance
    VectorXd phi_new = solver.solve(rhs_phi);

    // Calculate time derivative of phi
    VectorXd dphi_dt = (phi_new - phi) / dt;

    std::vector<double> dci_dt(6 * N);
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < N; ++j) {
        dci_dt[i * N + j] = beta[i] * (nu_sigma_f * Keff[j]) * phi[j] -
                            (lambda_i[i] * y[(i + 1) * N + j]);
      }
    }

    // std::cout << "dci_dt: " << dci_dt[0] << std::endl;

    // Populate dydt with the derivatives
    std::copy(dphi_dt.data(), dphi_dt.data() + dphi_dt.size(), dydt.begin());
    std::copy(dci_dt.begin(), dci_dt.end(), dydt.begin() + N);
  };

  // Initial condition vector
  std::vector<double> y0;
  if (step == 0) {
    y0.reserve(7 * N);
    y0.insert(y0.end(), phi_0, phi_0 + N);
    y0.insert(y0.end(), c1, c1 + N);
    y0.insert(y0.end(), c2, c2 + N);
    y0.insert(y0.end(), c3, c3 + N);
    y0.insert(y0.end(), c4, c4 + N);
    y0.insert(y0.end(), c5, c5 + N);
    y0.insert(y0.end(), c6, c6 + N);
  } else {
    y0 = y_n;
  }

  // Solve the system of ODEs
  std::vector<double> solution_y_n = ode_solver(y0, pde_to_ode_neutronics);

  // Extract the solution at the last time step
  std::vector<double> phi(solution_y_n.begin(), solution_y_n.begin() + N);
  std::vector<double> q_prime = phi; // Assuming q_prime is phi for now

  return {solution_y_n, q_prime};
}

