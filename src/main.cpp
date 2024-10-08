#include <iostream>
#include <numeric>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

#include "parameters.hpp"
#include "reactivity.hpp"
#include "neutronics.hpp"
#include "thermal_hydraulics.hpp"
#include "HX1.hpp"
#include "HX2.hpp"
#include "transport_delay.hpp"
#include "power_plant.hpp"
#include "data_saving.hpp"

using Eigen::VectorXd;

int main() {
    int time_span = 3000;

    double rho_insertion = 0.0;  // pcm

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

    Eigen::VectorXd buffer_hx_c, buffer_c_hx, buffer_r_hx, buffer_hx_r, buffer_r_pp, buffer_pp_r;
    
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
        std::cout << "Time step: " << step << std::endl;
        std::vector<double> rho_vec(rho.data(), rho.data() + rho.size());
        std::tie(y_n, q_prime) = neutronics(y_n, rho_vec, step);
        std::cout << "phi[0]: " << q_prime[0] << std::endl;
        std::cout << "phi[N-1]: " << q_prime[N-1] << std::endl;
        for (int i = 0; i < N; ++i) {
            q_prime[i] = q_prime[i] * sigma_f * A / flux_to_power;
        }
        std::cout << "q_prime[0]: " << q_prime[0] << std::endl;
        std::cout << "q_prime[N-1]: " << q_prime[N-1] << std::endl;
        phi = y_n.head(N);
        ci = y_n.tail(6 * N);

        phi_middle_matrix[step] = phi[N / 2];
        for (int i = 0; i < 6; ++i) {
            ci_middle_matrix(step, i) = ci[(i * N + (i + 1) * N) / 3];
        }

        double Ts_core_0 = transport_delay(Ts_HX1_0, tau_hx_c, Ts_in, buffer_hx_c, step);
        std::cout << "y_th: " << y_th[0] << std::endl;
        std::cout << "y_th_N: " << y_th[2 * N - 1] << std::endl;
        y_th = thermal_hydraulics(y_th, q_prime, Ts_core_0, step);

        temperature_fuel = y_th.head(N);
        temperature_graphite = y_th.tail(N);

        double Ts_core_L = y_th[N - 1];
        temperature_fuel_middle_matrix[step] = temperature_fuel[N / 3];

        double Ts_HX1_L = transport_delay(Ts_core_L, tau_c_hx, Ts_out, buffer_c_hx, step);
        Tss_HX1_0 = transport_delay(Tss_HX2_0, tau_r_hx, Tss_in, buffer_r_hx, step);

        y_hx1 = HX1(y_hx1, Ts_HX1_L, Tss_HX1_0, step);
        Ts_HX1 = y_hx1.head(Nx);
        Tss_HX1 = y_hx1.tail(Nx);

        Ts_HX1_0 = Ts_HX1[0];
        double Tss_HX1_L = Tss_HX1[Nx - 1];

        double Tss_HX2_L = transport_delay(Tss_HX1_L, tau_hx_r, Tss_out, buffer_hx_r, step);
        Tsss_HX2_0 = transport_delay(Tsss_pp_0, tau_pp_r, Tsss_in, buffer_pp_r, step);

        y_hx2 = HX2(y_hx2, Tss_HX2_L, Tsss_HX2_0, step);
        Tss_HX2 = y_hx2.head(Nx);
        Tsss_HX2 = y_hx2.tail(Nx);

        Tss_HX2_0 = Tss_HX2[0];
        double Tsss_HX2_L = Tsss_HX2[Nx - 1];

        double Tsss_pp_L = transport_delay(Tsss_HX2_L, tau_r_pp, Tsss_out, buffer_r_pp, step);
        double y_pp = power_plant_temp(Tsss_pp_L, step);
        Tsss_pp_0 = y_pp;

        rho = reactivity(initialS, initialG, temperature_fuel, temperature_graphite, step, time_span, rho_insertion);
        rho_matrix[step] = rho.sum();
    }

    // Plotting
    save_results(rho_matrix, phi_middle_matrix, ci_middle_matrix, temperature_fuel_middle_matrix);
    save_spacial_results(phi, ci, temperature_fuel, temperature_graphite, Ts_HX1, Tss_HX1, Tss_HX2, Tsss_HX2);

    return 0;
}
