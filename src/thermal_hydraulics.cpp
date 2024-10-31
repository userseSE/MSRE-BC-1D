#include "thermal_hydraulics.hpp"
#include "ode_solver_th.hpp"
#include "parameters.hpp"
#include <iostream>

void pde_to_ode_th(double t, const double y[length_th], double dydt[length_th],
                   Parameters &params, const double q_prime[N]) {
  for (int i = 0; i < N; ++i) {
    // Start with the matrix-vector multiplication part
    dydt[i] = 0;
    for (int j = 0; j < N; ++j) {
        dydt[i] += params.a_th * params.AT[i][j] * y[j];
    }
    // Add the temperature difference term
    dydt[i] += params.b_th * (y[N+i] - y[i]);
    // Add the heat source term
    dydt[i] += params.d_th * q_prime[i];
    // std::cout << "dydt[" << i << "]: " << dydt[i] << std::endl;
}
  // temperature_graphite_dt
  for (int i = 0; i < N; ++i) {
    dydt[N + i] = params.c_th * (y[i] - y[N + i]) + params.e_th * q_prime[i];
  }

  // Apply time-varying boundary conditions
  dydt[0] = params.bc_s0 - y[0];
  dydt[N] = params.bc_g0 - y[N];
  dydt[N - 1] = 0;
  dydt[N * 2 - 1] = 0;
  // for(int i=0;i<N*2;i++){
  //   std::cout<<"dydt["<<i<<"]"<<dydt[i]<<std::endl;
  // }
};

void thermal_hydraulics(double y_th[length_th], const double q_prime[N],
                        double Ts_core_0, int step, Parameters &params) {
  std::cout << "thermal_hydraulics called" << std::endl;
  // Set boundary conditions
  y_th[0] = Ts_core_0;
  y_th[N - 1] = params.fixed_boundary_sL;
  y_th[2 * N - 1] = params.fixed_boundary_gL;

  // Initial condition vector
  double y0[2 * N];
  if (step == 0) {
    for (int i = 0; i < N; ++i) {
      y_th[i] = params.initialS[i];
    }
    for (int i = 0; i < N; ++i) {
      y_th[N + i] = params.initialG[i];
    }
  }

  // Solve the ODE system
  ode_solver_th(y_th, pde_to_ode_th, step, params, q_prime);
}
