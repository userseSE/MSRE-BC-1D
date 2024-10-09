#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include "HX1.hpp"
#include "HX2.hpp"
#include "data_saving.hpp"
#include "neutronics.hpp"
#include "parameters.hpp"
#include "power_plant.hpp"
#include "reactivity.hpp"
#include "thermal_hydraulics.hpp"
#include "transport_delay.hpp"

using Eigen::VectorXd;

std::mutex parameter_mutex;

void run_simulation(double nu_sigma_f_value, double rho_0_value_input, int simulation_id) {
  {
    // Lock the mutex to safely modify global parameters
    std::lock_guard<std::mutex> guard(parameter_mutex);
    nu_sigma_f = nu_sigma_f_value;
    std::cout << "Simulation ID: " << simulation_id << " - Assigned nu_sigma_f: " << nu_sigma_f << std::endl;
    rho_0_value = rho_0_value_input;
    std::cout << "Simulation ID: " << simulation_id << " - Assigned rho_0_value: " << rho_0_value << std::endl;
  }

  int time_span = 20000;
  double rho_insertion = 0.0; // pcm

  // Initialization
  VectorXd rho = VectorXd::Zero(N);
  initialize_parameters();
  VectorXd y_n = VectorXd::Zero(7 * N);
  VectorXd q_prime = VectorXd::Zero(N);
  VectorXd y_th = VectorXd::Zero(2 * N);
  VectorXd y_hx1 = VectorXd::Zero(2 * Nx);
  VectorXd y_hx2 = VectorXd::Zero(2 * Nx);

  double Tss_HX2_0 = 0.0;
  double Ts_HX1_0 = 0.0;
  double Tss_HX1_0 = 0.0;
  double Tsss_pp_0 = 0.0;
  double Tsss_HX2_0 = 0.0;

  Eigen::VectorXd buffer_hx_c, buffer_c_hx, buffer_r_hx, buffer_hx_r,
      buffer_r_pp, buffer_pp_r;

  VectorXd rho_matrix = VectorXd::Zero(time_span);
  VectorXd phi_middle_matrix = VectorXd::Zero(time_span);
  Eigen::MatrixXd ci_middle_matrix = Eigen::MatrixXd::Zero(time_span, 6);
  VectorXd temperature_fuel_middle_matrix = VectorXd::Zero(time_span);

  VectorXd phi = VectorXd::Zero(N);
  VectorXd ci = VectorXd::Zero(6 * N);
  VectorXd temperature_fuel = VectorXd::Zero(N);
  VectorXd temperature_graphite = VectorXd::Zero(N);
  VectorXd Ts_HX1 = VectorXd::Zero(Nx);
  VectorXd Tss_HX1 = VectorXd::Zero(Nx);
  VectorXd Tss_HX2 = VectorXd::Zero(Nx);
  VectorXd Tsss_HX2 = VectorXd::Zero(Nx);

  for (int step = 0; step < time_span; ++step) {
    std::cout << "Simulation " << simulation_id << " - Time step: " << step << std::endl;
    std::vector<double> rho_vec(rho.data(), rho.data() + rho.size());
    std::tie(y_n, phi) = neutronics(y_n, rho_vec, step, nu_sigma_f);    // phi: n*cm−2*s−1
    std::cout << "phi[0]: " << phi[0] << std::endl;
    std::cout << "phi[100]: " << phi[100] << std::endl;
    for (int i = 0; i < N; ++i) {
      q_prime[i] = phi[i] * sigma_f * A / flux_to_power;
    }
    std::cout << "q_prime[0]: " << q_prime[0] << std::endl;
    std::cout << "q_prime[N-1]: " << q_prime[N - 1] << std::endl;
    // phi = y_n.head(N);
    ci = y_n.tail(6 * N);

    phi_middle_matrix[step] = phi[N / 2];
    for (int i = 0; i < 6; ++i) {
      ci_middle_matrix(step, i) = ci[(i * N + (i + 1) * N) / 3];
    }

    double Ts_core_0 =
        transport_delay(Ts_HX1_0, tau_hx_c, Ts_in, buffer_hx_c, step);
    y_th = thermal_hydraulics(y_th, q_prime, Ts_core_0, step);
    std::cout << "fuel_th[N/2]: " << y_th[N/4] << std::endl;
    temperature_fuel = y_th.head(N);
    temperature_graphite = y_th.tail(N);

    double Ts_core_L = y_th[N - 1];
    temperature_fuel_middle_matrix[step] = temperature_fuel[N / 3];

    double Ts_HX1_L =
        transport_delay(Ts_core_L, tau_c_hx, Ts_out, buffer_c_hx, step);
    Tss_HX1_0 = transport_delay(Tss_HX2_0, tau_r_hx, Tss_in, buffer_r_hx, step);

    y_hx1 = HX1(y_hx1, Ts_HX1_L, Tss_HX1_0, step);
    Ts_HX1 = y_hx1.head(Nx);
    Tss_HX1 = y_hx1.tail(Nx);

    Ts_HX1_0 = Ts_HX1[0];
    double Tss_HX1_L = Tss_HX1[Nx - 1];

    double Tss_HX2_L =
        transport_delay(Tss_HX1_L, tau_hx_r, Ts_out, buffer_hx_r, step);
    Tsss_HX2_0 =
        transport_delay(Tsss_pp_0, tau_pp_r, Tsss_in, buffer_pp_r, step);

    y_hx2 = HX2(y_hx2, Tss_HX2_L, Tsss_HX2_0, step);
    Tss_HX2 = y_hx2.head(Nx);
    Tsss_HX2 = y_hx2.tail(Nx);

    Tss_HX2_0 = Tss_HX2[0];
    double Tsss_HX2_L = Tsss_HX2[Nx - 1];

    double Tsss_pp_L =
        transport_delay(Tsss_HX2_L, tau_r_pp, Tsss_out, buffer_r_pp, step);
    double y_pp = power_plant_temp(Tsss_pp_L, step);
    Tsss_pp_0 = y_pp;

    rho = reactivity(initialS, initialG, temperature_fuel, temperature_graphite,
                     step, time_span, rho_insertion, rho_0_value);
    rho_matrix[step] = rho.sum();
  }

  // Save results for this simulation in a specific folder
  std::string folder = "../results/simulation_" + std::to_string(simulation_id);
  save_results(rho_matrix, phi_middle_matrix, ci_middle_matrix, temperature_fuel_middle_matrix, folder);
  save_spacial_results(phi, ci, temperature_fuel, temperature_graphite, Ts_HX1, Tss_HX1, Tss_HX2, Tsss_HX2, folder);
}

// int main() {
//   // Define parameter ranges
//   std::vector<double> nu_sigma_f_values;
//   int num_slices_nu_sigma_f = 10; // Define the number of slices
//   double nu_sigma_f_start = 0.0044120;
//   double nu_sigma_f_end = 0.0044130;
//   for (int i = 0; i < num_slices_nu_sigma_f; ++i) {
//     nu_sigma_f_values.push_back(nu_sigma_f_start + i * (nu_sigma_f_end - nu_sigma_f_start) / (num_slices_nu_sigma_f - 1));
//   }

//   double rho_0_value = 0.0;

//   std::vector<std::thread> threads;
//   int simulation_id = 0;

//   // Create threads for each parameter combination
//   for (double nu_sigma_f : nu_sigma_f_values) {
//     int current_simulation_id = simulation_id++;
//     threads.emplace_back([nu_sigma_f, rho_0_value, current_simulation_id]() {
//       run_simulation(nu_sigma_f, rho_0_value, current_simulation_id);
//     });
//   }

//   // Join all threads
//   for (auto &thread : threads) {
//     if (thread.joinable()) {
//       thread.join();
//     }
//   }

//   return 0;
// }

int main(){
    double nu_sigma_f = 0.0044121797;    //0.004412224, 0.004411798 for -9.6e-5, 0.004412258 for N=500, 0.0044121796-0
    double rho_0_value = 0.096e-5;       //-9.6e-5 to -10e-5
    int simulation_id = 0;
    run_simulation(nu_sigma_f, rho_0_value, simulation_id);
    return 0;
}