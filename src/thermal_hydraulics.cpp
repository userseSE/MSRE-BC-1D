#include "thermal_hydraulics.hpp"
#include "ode_solver.hpp"
#include "parameters.hpp"
// #include <algorithm>
// #include <cmath>
#include <iostream>
// #include <vector>

// Include necessary linear algebra headers
#include <Eigen/Dense>
// #include <Eigen/Sparse>
// #include <Eigen/SparseLU>
// #include <boost/numeric/odeint.hpp>

using Eigen::MatrixXd;
// using Eigen::SparseMatrix;
using Eigen::VectorXd;
// using namespace boost::numeric::odeint;

// Discretize the spatial domain
// SparseMatrix<double> AT(N, N);
MatrixXd AT = MatrixXd::Zero(N, N); // Initialize dense matrix with zeros

void initialize_AT(Parameters &params) {
  double dz = params.dz;
  AT = MatrixXd::Zero(N, N);

  for (int i = 0; i < N; ++i) {
    if (i == 0) {
      // Set boundary conditions at the first row
      AT(0, 0) = 1.0 / (dz * dz);
    } else if (i == N - 1) {
      // Set boundary conditions at the last row
      AT(N - 1, N - 1) = 1.0 / (dz * dz);
    } else {
      // Inner part of the matrix
      AT(i, i) = -2.0 / (dz * dz);
      if (i < N - 1) {
        AT(i, i + 1) = 1.0 / (dz * dz);
      }
      if (i > 0) {
        AT(i, i - 1) = 1.0 / (dz * dz);
      }
    }
  }
}

VectorXd thermal_hydraulics(VectorXd &y_th, const VectorXd &q_prime,
                            double Ts_core_0, int step, Parameters &params) {

  // Ensure AT is initialized
  static bool initialized = false;
  if (!initialized) {
    initialize_AT(params);
    initialized = true;
  }

  // Set boundary conditions
  y_th[0] = Ts_core_0;
  y_th[N - 1] = params.fixed_boundary_sL;
  y_th[2 * N - 1] = params.fixed_boundary_gL;

  // Define the ODE system for the thermal-hydraulics model compatible with
  // odeint
  std::function<void(double, const VectorXd &, VectorXd &)> pde_to_ode_th =
      [&](double t, const VectorXd &y, VectorXd &dydt) {
        // Split the input state vector y into temperature_fuel and
        // temperature_graphite
        VectorXd temperature_fuel = y.head(N);
        VectorXd temperature_graphite = y.tail(N);

        // Compute the time derivatives for temperature_fuel and
        // temperature_graphite
        VectorXd temperature_fuel_dt =
            params.a_th * (AT * temperature_fuel) +
            params.b_th * (temperature_graphite - temperature_fuel) +
            params.d_th * q_prime;

        VectorXd temperature_graphite_dt =
            params.c_th * (temperature_fuel - temperature_graphite) +
            params.e_th * q_prime;

        // // Apply time-varying boundary conditions
        temperature_fuel_dt[0] = params.bc_s0 - temperature_fuel[0];
        // temperature_fuel_dt[N - 1] = params.bc_sL - temperature_fuel[N - 1];
        temperature_graphite_dt[0] = params.bc_g0 - temperature_graphite[0];
        // temperature_graphite_dt[N - 1] = params.bc_gL -
        // temperature_graphite[N - 1];
        temperature_fuel_dt[N - 1] = 0;
        temperature_graphite_dt[N - 1] = 0;

        // Copy the derivatives to the dydt vector
        dydt.head(N) = temperature_fuel_dt;
        dydt.tail(N) = temperature_graphite_dt;
      };

  // Initial condition vector
  VectorXd y0(2 * N);
  if (step == 0) {
    y0.head(N) = params.initialS;
    y0.tail(N) = params.initialG;
  } else {
    y0 = y_th;
  }

  // Solve the ODE system
  VectorXd solution_y_th = ode_solver(y0, pde_to_ode_th, step);

  return solution_y_th;
}
