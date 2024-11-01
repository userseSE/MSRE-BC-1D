#include "HX1.hpp"
#include "HX2.hpp"
// #include "data_saving.hpp"
#include "neutronics.hpp"
#include "parameters.hpp"
#include "reactivity.hpp"
#include "thermal_hydraulics.hpp"

#include <iostream>

void run_simulation(int simulation_id) {
  Param_Neutronics params_neutr;
  Param_Thermal params_th;
  Param_HX1 params_hx1;
  Param_HX2 params_hx2;
  Param_React params_react;
  // initialize_neutronics(params_neutr);
  // initialize_thermal_hydraulics(params_th);
  // initialize_heat_exchanger_1(params_hx1);
  // initialize_heat_exchanger_2(params_hx2);
  // initialize_reactivity(params_react);
  // Initialization
  float rho[N];
  float y_n[length_neutr];
  float q_prime[N];
  float y_th[length_th];
  float y_hx1[length_hx];
  float y_hx2[length_hx];

  float Tss_HX2_0 = 0.0;
  float Ts_HX1_0 = 0.0;
  float Tss_HX1_0 = 0.0;
  float Tsss_pp_0 = 0.0;
  float Tsss_HX2_0 = 0.0;

  // VectorXd buffer_hx_c, buffer_c_hx, buffer_r_hx, buffer_hx_r,
  //     buffer_r_pp, buffer_pp_r;
  // Arrays for time-dependent data (replacing VectorXd)
  float rho_matrix[time_span] = {0};  // Initialize to 0
  float phi_middle_matrix[time_span * 2] = {0};
  float ci_middle_matrix[time_span * 6] = {0};  // 2D array for ci
  float temperature_fuel_middle_matrix[time_span] = {0};
  
  // Arrays for state variables (replacing VectorXd)
  float temperature_fuel[N] = {0};  // Fuel temperature
  float temperature_graphite[N] = {0};  // Graphite temperature
  
  // Arrays for heat exchanger temperatures (replacing VectorXd)
  float Ts_HX1[Nx] = {0};  // Heat exchanger 1 temperatures
  float Tss_HX1[Nx] = {0};  // Secondary side of heat exchanger 1
  float Tss_HX2[Nx] = {0};  // Heat exchanger 2 temperatures
  float Tsss_HX2[Nx] = {0};  // Secondary side of heat exchanger 2

  for (int step = 0; step < time_span; ++step) {
    std::cout << "Simulation " << simulation_id << " - Time step: " << step << std::endl;
    neutronics(y_n, rho, step, params_neutr);
    std::cout << "y_n[100]: " << y_n[100] << std::endl;
    // for (int i = 0; i < N; ++i) {
    //   q_prime[i] = (y_n[i] + y_n[N + i]) * params_neutr.sigma_f * params_neutr.A / params_neutr.flux_to_power;
    // }
    // std::cout << "q_prime[N/2]: " << q_prime[N/2] << std::endl;
    phi_middle_matrix[step] = y_n[N / 2];
    phi_middle_matrix[step + time_span] = y_n[N / 2 + N];
    for (int i = 0; i < 6; ++i) {
      ci_middle_matrix[step + i*time_span] = y_n[(i * N + (i + 1) * N) / 2 + 2 * N];
    }

    thermal_hydraulics(y_th, y_n,Ts_HX1_0, step, params_th);
    std::cout << "fuel_th[N/2]: " << y_th[N/4] << std::endl;
    float Ts_core_L = y_th[N - 1];
    temperature_fuel_middle_matrix[step] = y_th[N / 4];

    HX1(y_hx1, Ts_core_L, Tss_HX1_0, step, params_hx1);
    Ts_HX1_0 = y_hx1[0];
    float Tss_HX1_L = y_hx1[Nx - 1];


    HX2(y_hx2, Tss_HX1_L, step, params_hx2);
    Tss_HX2_0 = y_hx2[0];
    float Tsss_HX2_L = y_hx2[Nx - 1];

    reactivity(y_th, y_th + N, step, params_react, rho);
    std::cout << "rho[N/2]: " << rho[N/2] << std::endl;
    for (int i = 0; i < N; i++){
      rho_matrix[step] += rho[i];
    }
    rho_matrix[step] = (rho_matrix[step]/N) * 1e5;
  }

  // Save results for this simulation in a specific folder
  std::string folder = "results/simulation_" + std::to_string(simulation_id);
  // save_results(rho_matrix, phi_middle_matrix, ci_middle_matrix, temperature_fuel_middle_matrix, folder);
  for (int i = 0; i < N; ++i) {
    rho[i] = rho[i] * 1e5;
  }
  // save_spacial_results(rho, y_n, y_th, y_hx1, y_hx2, folder);
}

int main() {
  int simulation_id = 1;
  run_simulation(simulation_id);
  return 0;
}