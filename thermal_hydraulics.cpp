#include "thermal_hydraulics.hpp"
#include "ode_solver_th.hpp"
#include "parameters.hpp"

void pde_to_ode_th(float t, const float y[length_th], float dydt[length_th], Param_Thermal &params, const float q_prime[N]) {
  for (int i = 0; i < N; ++i) {
#pragma HLS UNROLL
    // Start with the matrix-vector multiplication part
    dydt[i] = 0.0;
    for (int idx = params.AT_csr.row_pointers[i]; idx < params.AT_csr.row_pointers[i + 1]; ++idx) {
#pragma HLS UNROLL
    	int j = params.AT_csr.col_indices[idx];
        dydt[i] += params.a_th * params.AT_csr.values[idx] * y[j];
    }
    // Add the temperature difference term
    dydt[i] += params.b_th * (y[N+i] - y[i]);
    // Add the heat source term
    dydt[i] += params.d_th * q_prime[i];
    dydt[N + i] = params.c_th * (y[i] - y[N + i]) + params.e_th * q_prime[i];
}

  // Apply time-varying boundary conditions
  dydt[0] = params.bc_s0 - y[0];
  dydt[N] = params.bc_g0 - y[N];
  dydt[N - 1] = 0;
  dydt[N * 2 - 1] = 0;
};

void thermal_hydraulics(float y_th[length_th], const float y_n[length_neutr],
                        float Ts_core_0, int step) {
  Param_Thermal params;

  float q_prime[N];
  for (int i = 0; i < N; ++i) {
#pragma HLS UNROLL
      q_prime[i] = (y_n[i] + y_n[N + i]) * params.sigma_f * params.A / params.flux_to_power;
  }
  // Set boundary conditions
  y_th[0] = Ts_core_0;
  y_th[N - 1] = params.fixed_boundary_sL;
  y_th[2 * N - 1] = params.fixed_boundary_gL;

  // Initial condition vector
  float y0[2 * N];
  if (step == 0) {
    for (int i = 0; i < N; ++i) {
#pragma HLS UNROLL
      y_th[i] = params.initialS[i];
    }
    for (int i = 0; i < N; ++i) {
#pragma HLS UNROLL
      y_th[N + i] = params.initialG[i];
    }
  }

  // Solve the ODE system
  ode_solver_th(y_th, pde_to_ode_th, step, params, q_prime);
}
