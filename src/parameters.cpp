#include <Eigen/Dense>
#include <array>
// #include <numeric>
// #include <vector>
// #include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "parameters.hpp"

void initialize_neutronics(Parameters& params) {
    double lambda_sum = std::accumulate(params.lambda_i.begin(), params.lambda_i.end(), 0.0);

    for (int i = 0; i < N; ++i) {
        params.phi1_0[i] = 1.0e13;  // Initial neutron flux
        params.phi2_0[i] = 0.5e13;  // Initial neutron flux
        double x = i * params.dz;  // Compute the spatial position
        params.c1[i] = ((params.beta[0] * (params.nu_sigma_f1)) / (params.lambda_i[0]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c2[i] = ((params.beta[1] * (params.nu_sigma_f1)) / (params.lambda_i[1]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c3[i] = ((params.beta[2] * (params.nu_sigma_f1)) / (params.lambda_i[2]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c4[i] = ((params.beta[3] * (params.nu_sigma_f1)) / (params.lambda_i[3]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c5[i] = ((params.beta[4] * (params.nu_sigma_f1)) / (params.lambda_i[4]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c6[i] = ((params.beta[5] * (params.nu_sigma_f1)) / (params.lambda_i[5]/6)) * (params.phi1_0[i]+params.phi2_0[i]);  
    }
    // Finite difference matrix for the second derivative using Crank-Nicolson
  // method
  VectorXd main_diag = (-2 / (params.dz * params.dz)) * VectorXd::Ones(N);
  VectorXd off_diag = (1 / (params.dz * params.dz)) * VectorXd::Ones(N - 1);
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
  D3(0, 0) = 1.0 / (params.dz * params.dz);
  D3(N - 1, N - 1) = 1.0 / (params.dz * params.dz);

  MatrixXd I = MatrixXd::Identity(N, N);     // Identity matrix
  params.A1 = I - 0.5 * params.dt * params.V1 * params.D1 * D3; 
  params.A2 = I - 0.5 * params.dt * params.V2 * params.D2 * D3;
  params.B1 = I + 0.5 * params.dt * params.V1 * params.D1 * D3;
  params.B2 = I + 0.5 * params.dt * params.V2 * params.D2 * D3;
  
  params.A1 = params.A1.inverse();
  params.A2 = params.A2.inverse();
}

void initialize_thermal_hydraulics(Parameters& params) {
    for (int i = 0; i < N; ++i) {
        double position = static_cast<double>(i) * params.L / (N - 1);
        params.initialS[i] = params.bc_s0 + (params.bc_sL - params.bc_s0) * ((0.5 + 0.5 * std::sin(M_PI * (position) / (params.L * 2))) * 0.8);
        params.initialG[i] = params.bc_g0 + (params.bc_gL - params.bc_g0) * ((0.5 + 0.5 * std::sin(M_PI * (position) / (params.L * 2))) * 1.05);
    }
}

void initialize_heat_exchanger_1(Parameters& params) {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * params.L_HX / (Nx - 1);
        params.u_init[i] = params.u_L + (params.u_L - params.u_H) * (0.5 + 0.5 * std::sin(M_PI * (position / params.L_HX)));
        params.v_init[i] = params.v_L + (params.v_L - params.v_H) * (0.5 + 0.5 * std::sin(M_PI * (position / params.L_HX)) * 1.05);
    }
}

void initialize_heat_exchanger_2(Parameters& params) {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * params.L_HX2 / (Nx - 1);
        params.u2_init[i] = params.u2_L + (params.u2_H - params.u2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / params.L_HX2)));
        params.v2_init[i] = params.v2_L + (params.v2_H - params.v2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / params.L_HX2)) * 1.05);
    }
}

void initialize_parameters(Parameters& params) {
    initialize_neutronics(params);
    initialize_thermal_hydraulics(params);
    initialize_heat_exchanger_1(params);
    initialize_heat_exchanger_2(params);
}
