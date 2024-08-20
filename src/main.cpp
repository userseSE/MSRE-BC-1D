// File: src/main.cpp
#include <iostream>
#include <vector>
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

int main() {
    int time_span = 100;

    double rho_insertion = 0.0;  // pcm

    // Initialization
    double rho = 0.0;
    
    std::vector<double> y_n(7 * N, 0.0);
    std::vector<double> q_prime(N, 0.0);
    std::vector<double> y_th(2 * N, 0.0);
    std::vector<double> y_hx1(2 * Nx, 0.0);
    std::vector<double> y_hx2(2 * Nx, 0.0);
    
    double Tss_HX2_0 = 0.0;
    double Ts_HX1_0 = 0.0;
    double Tss_HX1_0 = 0.0;
    double Tsss_pp_0 = 0.0;
    double Tsss_HX2_0=0.0;
    
    std::vector<double> buffer_hx_c, buffer_c_hx, buffer_r_hx, buffer_hx_r, buffer_r_pp, buffer_pp_r;
    
    std::vector<double> rho_matrix(time_span, 0.0);
    std::vector<double> phi_middle_matrix(time_span, 0.0);
    std::vector<std::vector<double>> ci_middle_matrix(time_span, std::vector<double>(6, 0.0));
    std::vector<double> temperature_fuel_middle_matrix(time_span, 0.0);

    std::vector<double> phi(N, 0.0);
    std::vector<double> ci(6 * N, 0.0);
    std::vector<double> temperature_fuel(N, 0.0);
    std::vector<double> temperature_graphite(N, 0.0);
    std::vector<double> Ts_HX1(Nx, 0.0);
    std::vector<double> Tss_HX1(Nx, 0.0);
    std::vector<double> Tss_HX2(Nx, 0.0);
    std::vector<double> Tsss_HX2(Nx, 0.0);

    for (int step = 0; step < time_span; ++step) {
        std::cout << "Time step: " << step << std::endl;
        
        std::tie(y_n, q_prime) = neutronics(y_n, rho, step);
        
        phi=std::vector<double>(y_n.begin(), y_n.begin() + N);
        ci=std::vector<double>(y_n.begin() + N, y_n.end());

        phi_middle_matrix[step] = phi[N / 2];
        for (int i = 0; i < 6; ++i) {
            ci_middle_matrix[step][i] = ci[(i * N + (i + 1) * N) / 2];
        }
        // std::cout<<"test n"<<std::endl;
        double Ts_core_0 = transport_delay(Ts_HX1_0, tau_hx_c, Ts_in, buffer_hx_c, step);
        // std::cout<<"test delay"<<std::endl;
        y_th = thermal_hydraulics(y_th, q_prime, Ts_core_0, step);
        // std::cout<<"test th"<<std::endl;
        temperature_fuel=std::vector<double>(y_th.begin(), y_th.begin() + N);
        temperature_graphite=std::vector<double>(y_th.begin() + N, y_th.end());

        double Ts_core_L = y_th.back();
        temperature_fuel_middle_matrix[step] = temperature_fuel[N / 2];

        double Ts_HX1_L = transport_delay(Ts_core_L, tau_c_hx, Ts_out, buffer_c_hx, step);
        Tss_HX1_0 = transport_delay(Tss_HX2_0, tau_r_hx, Tss_in, buffer_r_hx, step);

        y_hx1 = HX1(y_hx1, Ts_HX1_L, Tss_HX1_0, step);
        Ts_HX1=std::vector<double> (y_hx1.begin(), y_hx1.begin() + Nx);
        Tss_HX1=std::vector<double> (y_hx1.begin() + Nx, y_hx1.end());

        Ts_HX1_0 = Ts_HX1.front();
        double Tss_HX1_L = Tss_HX1.back();

        double Tss_HX2_L = transport_delay(Tss_HX1_L, tau_hx_r, Tss_out, buffer_hx_r, step);
        Tsss_HX2_0 = transport_delay(Tsss_pp_0, tau_pp_r, Tsss_in, buffer_pp_r, step);

        y_hx2 = HX2(y_hx2, Tss_HX2_L, Tsss_HX2_0, step);
        Tss_HX2 = std::vector<double> (y_hx2.begin(), y_hx2.begin() + Nx);
        Tsss_HX2 = std::vector<double> (y_hx2.begin() + Nx, y_hx2.end());

        Tss_HX2_0 = Tss_HX2.front();
        double Tsss_HX2_L = Tsss_HX2.back();

        double Tsss_pp_L = transport_delay(Tsss_HX2_L, tau_r_pp, Tsss_out, buffer_r_pp, step);
        double y_pp = power_plant_temp(Tsss_pp_L, step);
        Tsss_pp_0 = y_pp;

        rho = reactivity(temperature_fuel, temperature_graphite, step, time_span, rho_insertion);
        rho_matrix[step] = rho;
    }

    // Plotting
    // std::cout << "Plotting results..." << std::endl;
    save_results(rho_matrix, phi_middle_matrix, ci_middle_matrix, temperature_fuel_middle_matrix);
    std::cout << "Results have been plotted." << std::endl;
    save_spacial_results(phi, ci, temperature_fuel, temperature_graphite, Ts_HX1, Tss_HX1, Tss_HX2, Tsss_HX2);
    return 0;
}

