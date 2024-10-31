#include <cmath>
// #include <iostream>

#include "matrix_LU.hpp"
#include "parameters.hpp"
#include "csr_matrix.hpp"

void initialize_neutronics(Parameters &params) {

  for (int i = 0; i < N; ++i) {
    params.phi1_0[i] = 1.0e13f; // Initial neutron flux
    params.phi2_0[i] = 0.5e13f; // Initial neutron flux
    params.c1[i] = ((params.beta[0] * (params.nu_sigma_f1)) / (params.lambda_i[0] / 6)) * (params.phi1_0[i] + params.phi2_0[i]);
    params.c2[i] = ((params.beta[1] * (params.nu_sigma_f1)) / (params.lambda_i[1] / 6)) * (params.phi1_0[i] + params.phi2_0[i]);
    params.c3[i] = ((params.beta[2] * (params.nu_sigma_f1)) / (params.lambda_i[2] / 6)) * (params.phi1_0[i] + params.phi2_0[i]);
    params.c4[i] = ((params.beta[3] * (params.nu_sigma_f1)) / (params.lambda_i[3] / 6)) * (params.phi1_0[i] + params.phi2_0[i]);
    params.c5[i] = ((params.beta[4] * (params.nu_sigma_f1)) / (params.lambda_i[4] / 6)) * (params.phi1_0[i] + params.phi2_0[i]);
    params.c6[i] = ((params.beta[5] * (params.nu_sigma_f1)) / (params.lambda_i[5] / 6)) * (params.phi1_0[i] + params.phi2_0[i]);
  }

  // Finite difference matrix for the second derivative using Crank-Nicolson
  float D3[N][N];
  // Initialize D3 to zero and set the main diagonal and off-diagonals
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      D3[i][j] = 0.0; // Initialize all elements to zero
    }
    D3[i][i] = -2.0 / (params.dz * params.dz); // Set main diagonal
  }
  // Now set the superdiagonal and subdiagonal separately
  for (int i = 0; i < N - 1; ++i) {
    D3[i][i + 1] = 1.0 / (params.dz * params.dz); // Set superdiagonal
    D3[i + 1][i] = 1.0 / (params.dz * params.dz); // Set subdiagonal
  }
  // Apply Dirichlet boundary conditions for zero flux at boundaries
  for (int j = 0; j < N; ++j) {
    D3[0][j] = 0.0;
    D3[N - 1][j] = 0.0;
  }
  D3[0][0] = 1.0 / (params.dz * params.dz);
  D3[N - 1][N - 1] = 1.0 / (params.dz * params.dz);
  float I[N][N];
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      I[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
  // Compute params.A1, params.A2, params.B1, and params.B2
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      params.A1[i][j] = I[i][j] - 0.5 * params.dt * params.V1 * params.D1 * D3[i][j];
      params.A2[i][j] = I[i][j] - 0.5 * params.dt * params.V2 * params.D2 * D3[i][j];
      params.B1[i][j] = I[i][j] + 0.5 * params.dt * params.V1 * params.D1 * D3[i][j];
      params.B2[i][j] = I[i][j] + 0.5 * params.dt * params.V2 * params.D2 * D3[i][j];
    }
  }
  invert_LU(params.A1);
  invert_LU(params.A2);
  sparse_csr(params.B1, params.B1_csr);
  sparse_csr(params.B2, params.B2_csr);
}

// Simple sine approximation function (using a few terms from the Taylor series)
float approx_sin(float x) {
  // Sin(x) ≈ x - x^3 / 6 + x^5 / 120 (valid for small x)
  float x2 = x * x;
  return x - (x2 * x) / 6.0 + (x2 * x2 * x) / 120.0;
}

void initialize_thermal_hydraulics(Parameters &params) {
  for (int i = 0; i < N; ++i) {
    float position = static_cast<float>(i) * params.L / (N - 1);
    params.initialS[i] = params.bc_s0 + (params.bc_sL - params.bc_s0) * ((0.5 + 0.5 * approx_sin(M_PI * position / (params.L * 2))) * 0.8);
    params.initialG[i] = params.bc_g0 + (params.bc_gL - params.bc_g0) * ((0.5 + 0.5 * approx_sin(M_PI * position / (params.L * 2))) * 1.05);
  }
  for (int i = 0; i < N; ++i) {
    if (i == 0) {
      // Set boundary conditions at the first row
      params.AT[0][0] = 1.0 / (params.dz * params.dz);
    } else if (i == N - 1) {
      // Set boundary conditions at the last row
      params.AT[N - 1][N - 1] = 1.0 / (params.dz * params.dz);
    } else {
      // Inner part of the matrix
      params.AT[i][i] = -2.0 / (params.dz * params.dz);
      if (i < N - 1) {
        params.AT[i][i + 1] = 1.0 / (params.dz * params.dz);
      }
      if (i > 0) {
        params.AT[i][i - 1] = 1.0 / (params.dz * params.dz);
      }
    }
  }
  sparse_csr(params.AT, params.AT_csr);
}

void initialize_heat_exchanger_1(Parameters &params) {
  for (int i = 0; i < Nx; ++i) {
    float position = static_cast<float>(i) * params.L_HX / (Nx - 1);
    params.u_init[i] = params.u_L + (params.u_L - params.u_H) * (0.5 + 0.5 * approx_sin(M_PI * (position / params.L_HX)));
    params.v_init[i] = params.v_L + (params.v_L - params.v_H) * (0.5 + 0.5 * approx_sin(M_PI * (position / params.L_HX)) * 1.05);
  }
  for (int i = 0; i < Nx; ++i) {
    params.A_HX[i][i] = -2.0 / (params.dx * params.dx);
    if (i < Nx - 1) {
      params.A_HX[i][i + 1] = 1.0 / (params.dx * params.dx);
    }
    if (i > 0) {
      params.A_HX[i][i - 1] = 1.0 / (params.dx * params.dx);
    }
  }
  params.A_HX[0][0] = -1.0 / (params.dx * params.dx);
  params.A_HX[Nx - 1][Nx - 1] = -1.0 / (params.dx * params.dx);
  sparse_csr(params.A_HX, params.A_HX_csr);
}

void initialize_heat_exchanger_2(Parameters &params) {
  for (int i = 0; i < Nx; ++i) {
    float position = static_cast<float>(i) * params.L_HX2 / (Nx - 1);
    params.u2_init[i] = params.u2_L + (params.u2_H - params.u2_L) * (0.5 + 0.7 * approx_sin(M_PI * (position / params.L_HX2)));
    params.v2_init[i] = params.v2_L + (params.v2_H - params.v2_L) * (0.5 + 0.7 * approx_sin(M_PI * (position / params.L_HX2)) * 1.05);
  }
  for (int i = 0; i < Nx; ++i) {
    params.A_HX2[i][i] = -2.0 / (params.dx * params.dx);
    if (i < Nx - 1) {
      params.A_HX2[i][i + 1] = 1.0 / (params.dx * params.dx);
    }
    if (i > 0) {
      params.A_HX2[i][i - 1] = 1.0 / (params.dx * params.dx);
    }
  }
  params.A_HX2[0][0] = -1.0 / (params.dx * params.dx);
  params.A_HX2[Nx - 1][Nx - 1] = -1.0 / (params.dx * params.dx);
  sparse_csr(params.A_HX2, params.A_HX2_csr);
}

void initialize_parameters(Parameters &params) {
  initialize_neutronics(params);
  initialize_thermal_hydraulics(params);
  initialize_heat_exchanger_1(params);
  initialize_heat_exchanger_2(params);
}
