#include "HX1.hpp"
#include "HX2.hpp"
#include "data_saving.hpp"
#include "neutronics.hpp"
#include "parameters.hpp"
#include "reactivity.hpp"
#include "thermal_hydraulics.hpp"

#include <ctime>
#include <iostream>
#include <numeric>

void run_simulation(int simulation_id) {
  Parameters params;
  initialize_parameters(params);
  // Initialization
  double rho[N];
  double y_n[length_neutr];

  double q_prime[N];
  double y_th[length_th];
  double y_hx1[length_hx];
  double y_hx2[length_hx];

  double Tss_HX2_0 = 0.0;
  double Ts_HX1_0 = 0.0;
  double Tss_HX1_0 = 0.0;
  double Tsss_pp_0 = 0.0;
  double Tsss_HX2_0 = 0.0;

  // VectorXd buffer_hx_c, buffer_c_hx, buffer_r_hx, buffer_hx_r,
  //     buffer_r_pp, buffer_pp_r;
  // Arrays for time-dependent data (replacing VectorXd)
  double rho_matrix[time_span] = {0};  // Initialize to 0
  double phi_middle_matrix[time_span * 2] = {0};
  double ci_middle_matrix[time_span * 6] = {0};  // 2D array for ci
  double temperature_fuel_middle_matrix[time_span] = {0};
  
  // Arrays for state variables (replacing VectorXd)
  double temperature_fuel[N] = {0};  // Fuel temperature
  double temperature_graphite[N] = {0};  // Graphite temperature
  
  // Arrays for heat exchanger temperatures (replacing VectorXd)
  double Ts_HX1[Nx] = {0};  // Heat exchanger 1 temperatures
  double Tss_HX1[Nx] = {0};  // Secondary side of heat exchanger 1
  double Tss_HX2[Nx] = {0};  // Heat exchanger 2 temperatures
  double Tsss_HX2[Nx] = {0};  // Secondary side of heat exchanger 2

  for (int step = 0; step < time_span; ++step) {
    std::cout << "Simulation " << simulation_id << " - Time step: " << step << std::endl;
    neutronics(y_n, rho, step, params);
    std::cout << "y_n[100]: " << y_n[100] << std::endl;
    for (int i = 0; i < N; ++i) {
      q_prime[i] = (y_n[i] + y_n[N + i]) * params.sigma_f * params.A / params.flux_to_power;
    }
    std::cout << "q_prime[N/2]: " << q_prime[N/2] << std::endl;
    phi_middle_matrix[step] = y_n[N / 2];
    phi_middle_matrix[step + time_span] = y_n[N / 2 + N];
    for (int i = 0; i < 6; ++i) {
      ci_middle_matrix[step + i*time_span] = y_n[(i * N + (i + 1) * N) / 2 + 2 * N];
    }

    // double Ts_core_0 = transport_delay(Ts_HX1_0, params.tau_hx_c, params.Ts_in,
    //                                    buffer_hx_c, step);
    thermal_hydraulics(y_th, q_prime,Ts_HX1_0, step, params);
    std::cout << "fuel_th[N/2]: " << y_th[N/4] << std::endl;
    double Ts_core_L = y_th[N - 1];
    temperature_fuel_middle_matrix[step] = y_th[N / 4];

    // double Ts_HX1_L = transport_delay(Ts_core_L, params.tau_c_hx, params.Ts_out,
    //                                   buffer_c_hx, step);
    // Tss_HX1_0 = transport_delay(Tss_HX2_0, params.tau_r_hx, params.Tss_in,
    //                             buffer_r_hx, step);

    HX1(y_hx1, Ts_core_L, Tss_HX1_0, step, params);
    // Ts_HX1 = y_hx1.head(Nx);
    // Tss_HX1 = y_hx1.tail(Nx);
    Ts_HX1_0 = y_hx1[0];
    double Tss_HX1_L = y_hx1[Nx - 1];

    // double Tss_HX2_L = transport_delay(Tss_HX1_L, params.tau_hx_r,
    //                                    params.Ts_out, buffer_hx_r, step);
    // Tsss_HX2_0 = transport_delay(Tsss_pp_0, tau_pp_r, Tsss_in, buffer_pp_r,
    // step);

    HX2(y_hx2, Tss_HX1_L, step, params);
    // Tss_HX2 = y_hx2.head(Nx);
    // Tsss_HX2 = y_hx2.tail(Nx);
    Tss_HX2_0 = y_hx2[0];
    double Tsss_HX2_L = y_hx2[Nx - 1];

    // double Tsss_pp_L = transport_delay(Tsss_HX2_L, params.tau_r_pp,
    //                                    params.Tsss_out, buffer_r_pp, step);
    // double y_pp = power_plant_temp(Tsss_pp_L, step);
    // Tsss_pp_0 = y_pp;
    // for (int i = 0; i < N; ++i) {
    //   std::cout <<y_th[i] << " ";
    // }
    // std::cout << std::endl;
    reactivity(y_th, y_th + N, step, time_span, params, rho);
    std::cout << "rho[N/2]: " << rho[N/2] << std::endl;
    rho_matrix[step] = (std::accumulate(rho, rho + N, 0.0) / N) * 1e5;
    // rho_matrix[step] = rho[N / 2] * 1e5;
  }

  // Save results for this simulation in a specific folder
  std::string folder = "results/simulation_" + std::to_string(simulation_id);
  save_results(rho_matrix, phi_middle_matrix,
               ci_middle_matrix, temperature_fuel_middle_matrix, folder);
  for (int i = 0; i < N; ++i) {
    rho[i] = rho[i] * 1e5;
  }
  save_spacial_results(rho, y_n, y_th, y_hx1, y_hx2, folder);
}

int main() {
  int simulation_id = 1;
  run_simulation(simulation_id);
  return 0;
}