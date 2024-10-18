#include "neutronics.hpp"
#include "ode_solver.hpp"
#include "parameters.hpp"
#include "PartialPivLU.hpp"
#include <iostream>

// Include necessary linear algebra headers
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// typedef Eigen::VectorXd state_type;

Eigen::VectorXd neutronics(const state_type &y_n,
                           const std::vector<double> &rho, int step,
                           Parameters &params) {
  VectorXd Keff = VectorXd::Zero(N);
  std::transform(rho.data(), rho.data() + N, Keff.data(),
                 [](double r) { return (1.0 / (1.0 + r)); });

  double dz = params.dz;
  double dt = params.dt;
  double V1 = params.V1;
  double V2 = params.V2;
  double D1 = params.D1;
  double D2 = params.D2;
  double sigma_a1 = params.sigma_a1;
  double sigma_a2 = params.sigma_a2;
  double nu_sigma_f1 = params.nu_sigma_f1;
  double nu_sigma_f2 = params.nu_sigma_f2;
  double sigma_s12 = params.sigma_s12;
  double Beta = params.Beta;
  VectorXd lambda_i = VectorXd::Map(params.lambda_i.data(), params.lambda_i.size());
  VectorXd beta = VectorXd::Map(params.beta.data(), params.beta.size());

  // Finite difference matrix for the second derivative using Crank-Nicolson
  // method
  VectorXd main_diag = (-2 / (dz * dz)) * VectorXd::Ones(N);
  VectorXd off_diag = (1 / (dz * dz)) * VectorXd::Ones(N - 1);
  MatrixXd D3 = MatrixXd::Zero(N, N);

  // Set the main diagonal
  D3.diagonal() = main_diag;

  // Set the superdiagonal (above the main diagonal)
  D3.diagonal(1) = off_diag;

  // Set the subdiagonal (below the main diagonal)
  D3.diagonal(-1) = off_diag;

  // Apply Dirichlet boundary conditions for zero flux at boundaries
  D3.row(0).setZero();
  D3.row(N - 1).setZero();
  D3(0, 0) = 1.0 / (dz * dz);
  D3(N - 1, N - 1) = 1.0 / (dz * dz);

  MatrixXd I = MatrixXd::Identity(N, N);     // Identity matrix
  MatrixXd A1 = I - 0.5 * dt * V1 * D1 * D3; // No need for sparseView() here
  MatrixXd A2 = I - 0.5 * dt * V2 * D2 * D3;
  MatrixXd B1 = I + 0.5 * dt * V1 * D1 * D3;
  MatrixXd B2 = I + 0.5 * dt * V2 * D2 * D3;

  PartialPivLU solver1;
  PartialPivLU solver2;
  solver1.compute(A1);
  solver2.compute(A2);

  // ODE system function compatible with odeint
  std::function<void(double, const VectorXd &, VectorXd &)>
      pde_to_ode_neutronics = [&](double t, const VectorXd &y, VectorXd &dydt) {
        // Split y into fast flux (phi1), thermal flux (phi2), and delayed
        // neutron precursors
        VectorXd phi1 = y.head(N);       // Fast group
        VectorXd phi2 = y.segment(N, N); // Thermal group
        VectorXd lambda_ci = VectorXd::Zero(N);

        // Compute delayed neutron precursor contributions
        for (int i = 0; i < 6; ++i) {
          lambda_ci += lambda_i[i] * y.segment((i + 2) * N, N);
        }

        // Right-hand side for fast group
        VectorXd rhs_phi1 = B1 * phi1 + dt * V1 * ((-sigma_a1 + (1.0 - Beta) * ((nu_sigma_f1 + nu_sigma_f2) / Keff.array())).matrix()
                     .cwiseProduct(phi1) + lambda_ci);
        // VectorXd rhs_phi1 = B1 * phi1 + dt * V1 * ((-sigma_a1 + (1.0 - Beta)
        // * (nu_sigma_f1 / Keff.array())).matrix().cwiseProduct(phi1) +
        // lambda_ci); Right-hand side for thermal group
        VectorXd rhs_phi2 =
            B2 * phi2 + dt * V2 * ((-sigma_a2 * phi2) + (sigma_s12 * phi1));

        // Solve for the new fluxes
        VectorXd phi1_new = solver1.solve(rhs_phi1);
        VectorXd phi2_new = solver2.solve(rhs_phi2);

        // Calculate time derivative of phi1 and phi2
        VectorXd dphi1_dt = (phi1_new - phi1) / dt;
        VectorXd dphi2_dt = (phi2_new - phi2) / dt;

        // phi2 = VectorXd::Zero(N);

        // Compute dci/dt for delayed neutron precursors
        VectorXd dci_dt(6 * N);
        for (int i = 0; i < 6; ++i) {
          dci_dt.segment(i * N, N) =
              beta[i] * (nu_sigma_f1 + nu_sigma_f2) * (phi1 + phi2) -
              lambda_i[i] * y.segment((i + 2) * N, N);
        }

        // Populate dydt with the derivatives
        dydt.head(N) = dphi1_dt;
        dydt.segment(N, N) = dphi2_dt;
        dydt.tail(6 * N) = dci_dt;
      };

  // Initial condition vector
  state_type y0(8 *
                N); // Updated to hold both phi1, phi2, and 6 precursor groups
  if (step == 0) {
    y0 << params.phi1_0, params.phi2_0, params.c1, params.c2, params.c3,
        params.c4, params.c5, params.c6;
  } else {
    y0 = y_n;
  }

  // Solve the system of ODEs
  state_type solution_y_n =
      ode_solver(y0, pde_to_ode_neutronics, step); // Pass t0, t_end, dt

  // Extract the solution at the last time step
  VectorXd phi1 = solution_y_n.head(N);
  VectorXd phi2 = solution_y_n.segment(N, N);
  VectorXd q_prime = phi1 + phi2; // Combine fluxes for q_prime

  return {solution_y_n};
}
