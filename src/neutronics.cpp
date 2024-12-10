#include <cmath>
#include <iostream>

#include "neutronics.hpp"
#include "ode_solver_neutronics.hpp"
#include "parameters.hpp"

// Standalone function definition
void pde_to_ode_neutronics(const float y[length_neutr], float dydt[length_neutr], Param_Neutronics &params_neutr, const float Keff[N]) {
  // std::cout << "Neutronics PDE to ODE solver called" << std::endl;
  float lambda_ci[N]; // For summing contributions
  // Initialize lambda_ci to zero
  for (int i = 0; i < N; ++i) {
    lambda_ci[i] = 0.0;
  }
  // Delayed neutron precursor contributions
  for (int i = 0; i < 6; ++i) {
    const float *ci = &y[(i + 2) * N]; // Access ci for each group
    for (int j = 0; j < N; ++j) {
      lambda_ci[j] += params_neutr.lambda_i[i] * ci[j]; // Sum contributions to lambda_ci
    }
  }
  float rhs_phi1[N];
  float rhs_phi2[N];
  for (int i = 0; i < N; ++i) {
    rhs_phi1[i] = 0.0;
    rhs_phi2[i] = 0.0;
    // Perform sparse matrix-vector multiplication for B1_csr * y
    for (int idx = params_neutr.B1_csr.row_pointers[i]; idx < params_neutr.B1_csr.row_pointers[i + 1]; ++idx) {
      int j = params_neutr.B1_csr.col_indices[idx];
      rhs_phi1[i] += params_neutr.B1_csr.values[idx] * y[j];
      rhs_phi2[i] += params_neutr.B2_csr.values[idx] * y[N + j];
    }
    // Add the remaining terms as before
    rhs_phi1[i] += params_neutr.dt * params_neutr.V1 * (((-params_neutr.sigma_a1 + (1.0 - params_neutr.Beta) * ((params_neutr.nu_sigma_f1 + params_neutr.nu_sigma_f2) / Keff[i])) * y[i]) + lambda_ci[i]);
    rhs_phi2[i] += params_neutr.dt * params_neutr.V2 * (-params_neutr.sigma_a2 * y[N + i] + params_neutr.sigma_s12 * y[i]);
  }
  // float rhs_phi2[N];
  // for (int i = 0; i < N; ++i) {
  //   rhs_phi2[i] = 0.0;
  //   // Perform sparse matrix-vector multiplication for B2_csr * y[N:]
  //   for (int idx = params_neutr.B2_csr.row_pointers[i]; idx < params_neutr.B2_csr.row_pointers[i + 1]; ++idx) {
  //     int j = params_neutr.B2_csr.col_indices[idx];
  //     rhs_phi2[i] += params_neutr.B2_csr.values[idx] * y[N + j];
  //   }
  //   // Add the remaining terms as before
  //   rhs_phi2[i] += params_neutr.dt * params_neutr.V2 * (-params_neutr.sigma_a2 * y[N + i] + params_neutr.sigma_s12 * y[i]);
  // }
  float phi1_new[N];
  float phi2_new[N];
  for (int i = 0; i < N; ++i) {
    phi1_new[i] = 0.0;
    phi2_new[i] = 0.0;
    for (int j = 0; j < N; ++j) {
      phi1_new[i] += params_neutr.A1[i][j] * rhs_phi1[j]; 
      phi2_new[i] += params_neutr.A2[i][j] * rhs_phi2[j]; 
    }
  }
  // Calculate time derivative of phi1 and phi2
  for (int i = 0; i < N; ++i) {
    dydt[i] = (phi1_new[i] - y[i]) / params_neutr.dt;         // dphi1/dt
    dydt[N + i] = (phi2_new[i] - y[N + i]) / params_neutr.dt; // dphi2/dt
  }

  // Compute dci/dt for delayed neutron precursors
  for (int i = 0; i < 6; ++i) {
    // float *dci_dt = &dydt[(i + 2) * N];
    for (int j = 0; j < N; ++j) {
      dydt[(i + 2) * N + j] = params_neutr.beta[i] * (params_neutr.nu_sigma_f1 + params_neutr.nu_sigma_f2) * (y[j] + y[N + j]) - params_neutr.lambda_i[i] * y[(i + 2) * N + j];
    }
  }
}

// Modified neutronics function using the function pointer
void neutronics(float y_n[length_neutr], const float rho[N], int step, Param_Neutronics& params, float min_step_neutr) {
  // Param_Neutronics params;
  float Keff[N] = {0};
  for (int i = 0; i < N; ++i) {
    Keff[i] = (1.0 / (1.0 + rho[i]));
  }

  // Initialize the initial condition vector
  int index = 0;
  if (step == 0) {
    // Concatenate phi1_0, phi2_0, and c1 through c6
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.phi1_0[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.phi2_0[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.c1[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.c2[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.c3[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.c4[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.c5[i];
    }
    for (int i = 0; i < N; ++i) {
      y_n[index++] = params.c6[i];
    }
  }

  // Function pointer to the ODE system
  OdeFuncPointer ode_func = pde_to_ode_neutronics;

  // Solve the system of ODEs using ode_solver
  ode_solver_neutr(y_n, ode_func, step, params, Keff, min_step_neutr); // Pass t0, t_end, dt and necessary parameters
}
