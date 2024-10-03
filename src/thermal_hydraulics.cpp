#include "thermal_hydraulics.hpp"
#include "ode_solver.hpp"
#include "parameters.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

// Include necessary linear algebra headers
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <boost/numeric/odeint.hpp>

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using namespace boost::numeric::odeint;

const double a_th = Vc;
const double b_th = U_gs / (Ms * c_p_s);
const double c_th = U_sg / (Mg * c_p_g);
const double d_th = L * gamma / (Ms * c_p_s);
const double e_th = L * (1 - gamma) / (Mg * c_p_g);

// Discretize the spatial domain
SparseMatrix<double> AT(N, N);

void initialize_AT() {

AT = SparseMatrix<double>(N, N);
std::vector<Eigen::Triplet<double>> tripletList;
tripletList.reserve(3 * N - 2);

for (int i = 0; i < N; ++i) {
  if (i == 0) {
    // Set boundary conditions at the first row
    tripletList.emplace_back(0, 0, 1.0 / (dz * dz));
  } else if (i == N - 1) {
    // Set boundary conditions at the last row
    tripletList.emplace_back(N - 1, N - 1, 1.0 / (dz * dz));
  } else {
    // Inner part of the matrix
    tripletList.emplace_back(i, i, -2.0 / (dz * dz));
    if (i < N - 1) {
      tripletList.emplace_back(i, i + 1, 1.0 / (dz * dz));
    }
    if (i > 0) {
      tripletList.emplace_back(i, i - 1, 1.0 / (dz * dz));
    }
  }
}

// Now set the sparse matrix with the triplet list
AT.setFromTriplets(tripletList.begin(), tripletList.end());
}

std::vector<double> thermal_hydraulics(std::vector<double> &y_th,
                                       const std::vector<double> &q_prime,
                                       double Ts_core_0, int step) {
  // Ensure AT is initialized
  static bool initialized = false;
  if (!initialized) {
    initialize_AT();
    initialized = true;
  }

  // Set boundary conditions
  if (step == 0) {
    y_th[0] = Ts_core_0;
  } else {
    y_th[0] = Ts_core_0;
  }

  // Define the ODE system for the thermal-hydraulics model compatible with
  // odeint
  auto pde_to_ode_th = [&](const std::vector<double> &y,
                           std::vector<double> &dydt, double t) {
    // Split the input state vector y into temperature_fuel and
    // temperature_graphite
    VectorXd temperature_fuel = Eigen::Map<const VectorXd>(y.data(), N);
    VectorXd temperature_graphite = Eigen::Map<const VectorXd>(y.data() + N, N);

    // Compute the time derivatives for temperature_fuel and
    // temperature_graphite
    VectorXd temperature_fuel_dt =
        a_th * (AT * temperature_fuel) +
        b_th * (temperature_graphite - temperature_fuel) +
        d_th * VectorXd::Map(q_prime.data(), N);

    VectorXd temperature_graphite_dt =
        c_th * (temperature_fuel - temperature_graphite) +
        e_th * VectorXd::Map(q_prime.data(), N);

    // Apply time-varying boundary conditions
    temperature_fuel_dt[0] = bc_s0 - temperature_fuel[0];
    temperature_fuel_dt[N - 1] = bc_sL - temperature_fuel[N - 1];
    temperature_graphite_dt[0] = bc_g0 - temperature_graphite[0];
    temperature_graphite_dt[N - 1] = bc_gL - temperature_graphite[N - 1];

    // Copy the derivatives to the dydt vector
    std::copy(temperature_fuel_dt.data(), temperature_fuel_dt.data() + N,
              dydt.begin());
    std::copy(temperature_graphite_dt.data(),
              temperature_graphite_dt.data() + N, dydt.begin() + N);
  };

  // Initial condition vector
  std::vector<double> y0(2 * N);
  if (step == 0) {
    std::copy(initialS.data(), initialS.data() + N, y0.begin());
    std::copy(initialG.data(), initialG.data() + N, y0.begin() + N);
  } else {
    y0 = y_th;
  }

  // Solve the ODE system
  // std::vector<double> solution_y_th = ode_solver(y0, pde_to_ode_th);
  std::vector<double> solution_y_th = y0;
  solution_y_th = ode_solver(y0, pde_to_ode_th);

  return solution_y_th;
}