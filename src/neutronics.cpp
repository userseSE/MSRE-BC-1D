#include <cmath>

#include "neutronics.hpp"
#include "ode_solver_neutronics.hpp"
#include "parameters.hpp"

// Define the function pointer type
// typedef void (*OdeFuncPointer)(float, const float[length_neutr],
// float[length_neutr], Parameters &params, const float Keff[N]);

// Standalone function definition
void pde_to_ode_neutronics(float t, const float y[length_neutr],
                           float dydt[length_neutr], Param_Neutronics &params,
                           const float Keff[N]) {
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
      lambda_ci[j] +=
          params.lambda_i[i] * ci[j]; // Sum contributions to lambda_ci
    }
  }

  // Compute the right-hand side for phi1 (fast group)
  // float rhs_phi1[N];
  // for (int i = 0; i < N; ++i) {
  //     rhs_phi1[i] = 0.0;
  //     for (int j = 0; j < N; ++j) {
  //         rhs_phi1[i] += params.B1[i][j] * y[j]; // Matrix multiplication for
  //         B1 * phi1
  //     }
  //     rhs_phi1[i] += params.dt * params.V1 * (((-params.sigma_a1 + (1.0 -
  //     params.Beta) * ((params.nu_sigma_f1 + params.nu_sigma_f2) / Keff[i])) *
  //     y[i]) + lambda_ci[i]);
  //     // std::cout<<params.dt * params.V1 * (-params.sigma_a1 + (1.0 -
  //     params.Beta) * ((params.nu_sigma_f1 + params.nu_sigma_f2) /
  //     Keff[i]))<<std::endl;
  // }
  float rhs_phi1[N];
  for (int i = 0; i < N; ++i) {
    rhs_phi1[i] = 0.0;
    // Perform sparse matrix-vector multiplication for B1_csr * y
    for (int idx = params.B1_csr.row_pointers[i]; idx < params.B1_csr.row_pointers[i + 1]; ++idx) {
      int j = params.B1_csr.col_indices[idx];
      rhs_phi1[i] += params.B1_csr.values[idx] * y[j];
    }
    // Add the remaining terms as before
    rhs_phi1[i] += params.dt * params.V1 * (((-params.sigma_a1 + (1.0 - params.Beta) *
               ((params.nu_sigma_f1 + params.nu_sigma_f2) / Keff[i])) * y[i]) + lambda_ci[i]);
  }

  // for (int i = 0; i < N; ++i) {
  //     for (int j = 0; j < N; ++j) {
  //         std::cout  << params.B1[i][j]<< " ";
  //     }
  //     std::cout << std::endl;
  // }
  // for (int i = 0; i < N; ++i) {
  //     std::cout << "rhs_phi1[" << i << "]: " << rhs_phi1[i] << std::endl;
  // }

  // Compute the right-hand side for phi2 (thermal group)
    // float rhs_phi2[N];
    // for (int i = 0; i < N; ++i) {
    //   rhs_phi2[i] = 0.0;
    //   for (int j = 0; j < N; ++j) {
    //     rhs_phi2[i] +=
    //         params.B2[i][j] * y[N + j]; // Matrix multiplication for B2 * phi2
    //   }
    //   rhs_phi2[i] += params.dt * params.V2 *
    //                  (-params.sigma_a2 * y[N + i] + params.sigma_s12 * y[i]);
    // }
  float rhs_phi2[N];
  for (int i = 0; i < N; ++i) {
    rhs_phi2[i] = 0.0;

    // Perform sparse matrix-vector multiplication for B2_csr * y[N:]
    for (int idx = params.B2_csr.row_pointers[i]; idx < params.B2_csr.row_pointers[i + 1]; ++idx) {
      int j = params.B2_csr.col_indices[idx];
      rhs_phi2[i] += params.B2_csr.values[idx] * y[N + j];
    }
    // Add the remaining terms as before
    rhs_phi2[i] += params.dt * params.V2 *
                   (-params.sigma_a2 * y[N + i] + params.sigma_s12 * y[i]);
  }

  // for (int i = 0; i < N; ++i) {
  //     std::cout << "rhs_phi2[" << i << "]: " << rhs_phi2[i] << std::endl;
  // }

  // Solve for the new fluxes (phi1_new and phi2_new)
    float phi1_new[N];
    float phi2_new[N];
    for (int i = 0; i < N; ++i) {
      phi1_new[i] = 0.0;
      phi2_new[i] = 0.0;
      for (int j = 0; j < N; ++j) {
        phi1_new[i] += params.A1[i][j] * rhs_phi1[j]; // Solve A1 * rhs_phi1
        phi2_new[i] += params.A2[i][j] * rhs_phi2[j]; // Solve A2 * rhs_phi2
      }
    }
  // float phi1_new[N];
  // float phi2_new[N];
  // for (int i = 0; i < N; ++i) {
  //   phi1_new[i] = 0.0;
  //   phi2_new[i] = 0.0;
  //   // Sparse matrix-vector multiplication for A1_csr * rhs_phi1
  //   for (int idx = params.A1_csr_inv.row_pointers[i]; idx < params.A1_csr_inv.row_pointers[i + 1]; ++idx) {
  //     int j = params.A1_csr_inv.col_indices[idx];
  //     phi1_new[i] += params.A1_csr_inv.values[idx] * rhs_phi1[j];
  //   }
  //   // Sparse matrix-vector multiplication for A2_csr * rhs_phi2
  //   for (int idx = params.A2_csr_inv.row_pointers[i]; idx < params.A2_csr_inv.row_pointers[i + 1]; ++idx) {
  //     int j = params.A2_csr_inv.col_indices[idx];
  //     phi2_new[i] += params.A2_csr_inv.values[idx] * rhs_phi2[j];
  //   }
  // }

  // for (int i = 0; i < N; ++i) {
  //     for (int j = 0; j < N; ++j) {
  //         std::cout  << params.A1[i][j]<< " ";
  //     }
  //     std::cout << std::endl;
  // }

  // Calculate time derivative of phi1 and phi2
  for (int i = 0; i < N; ++i) {
    dydt[i] = (phi1_new[i] - y[i]) / params.dt;         // dphi1/dt
    dydt[N + i] = (phi2_new[i] - y[N + i]) / params.dt; // dphi2/dt
  }

  // Compute dci/dt for delayed neutron precursors
  for (int i = 0; i < 6; ++i) {
    // float *dci_dt = &dydt[(i + 2) * N];
    for (int j = 0; j < N; ++j) {
      dydt[(i + 2) * N + j] = params.beta[i] *
                                  (params.nu_sigma_f1 + params.nu_sigma_f2) *
                                  (y[j] + y[N + j]) -
                              params.lambda_i[i] * y[(i + 2) * N + j];
    }
  }
  // for (int i = 0 ; i< N; ++i){
  //     std::cout << "dydt[" << i << "]: " << dydt[i] << std::endl;
  // }
  // std::cout << "Neutronics PDE to ODE solver finished" << std::endl;
}

// Modified neutronics function using the function pointer
void neutronics(float y_n[length_neutr], const float rho[N], int step,
                Param_Neutronics &params) {
  // std::cout << "Neutronics is called for step " << step << std::endl;
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
  ode_solver_neutr(y_n, ode_func, step, params,
                   Keff); // Pass t0, t_end, dt and necessary parameters
}
