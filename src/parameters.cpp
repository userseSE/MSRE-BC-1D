#include <Eigen/Dense>
#include <array>
#include <numeric>
#include <vector>
#include <cmath>

#include "parameters.hpp"

void initialize_neutronics(Parameters& params) {
    double lambda_sum = std::accumulate(params.lambda_i.begin(), params.lambda_i.end(), 0.0);

    for (int i = 0; i < N; ++i) {
        params.phi1_0[i] = 1.0e13;  // Initial neutron flux
        params.phi2_0[i] = 0.4e13;  // Initial neutron flux
        double x = i * params.dz;  // Compute the spatial position
        params.c1[i] = ((params.beta[0] * (params.nu_sigma_f1)) / (params.lambda_i[0]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c2[i] = ((params.beta[1] * (params.nu_sigma_f1)) / (params.lambda_i[1]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c3[i] = ((params.beta[2] * (params.nu_sigma_f1)) / (params.lambda_i[2]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c4[i] = ((params.beta[3] * (params.nu_sigma_f1)) / (params.lambda_i[3]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c5[i] = ((params.beta[4] * (params.nu_sigma_f1)) / (params.lambda_i[4]/6)) * (params.phi1_0[i]+params.phi2_0[i]);
        params.c6[i] = ((params.beta[5] * (params.nu_sigma_f1)) / (params.lambda_i[5]/6)) * (params.phi1_0[i]+params.phi2_0[i]);  
    }
}

void initialize_thermal_hydraulics(Parameters& params) {
    for (int i = 0; i < N; ++i) {
        double position = static_cast<double>(i) * params.L / (N - 1);
        params.initialS[i] = params.bc_s0 + (params.bc_sL - params.bc_s0) * ((0.5 + 0.5 * std::sin(M_PI * (position) / (params.L * 2))) * 0.8);
        params.initialG[i] = params.bc_g0 + (params.bc_gL - params.bc_g0) * ((0.5 + 0.5 * std::sin(M_PI * (position) / (params.L * 2))) * 1.05);
    }
}

void initialize_heat_exchanger_1(Parameters& params) {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * params.L_HX / (Nx - 1);
        params.u_init[i] = params.u_L + (params.u_L - params.u_H) * (0.5 + 0.5 * std::sin(M_PI * (position / params.L_HX)));
        params.v_init[i] = params.v_L + (params.v_L - params.v_H) * (0.5 + 0.5 * std::sin(M_PI * (position / params.L_HX)) * 1.05);
    }
}

void initialize_heat_exchanger_2(Parameters& params) {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * params.L_HX2 / (Nx - 1);
        params.u2_init[i] = params.u2_L + (params.u2_H - params.u2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / params.L_HX2)));
        params.v2_init[i] = params.v2_L + (params.v2_H - params.v2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / params.L_HX2)) * 1.05);
    }
}

void initialize_parameters(Parameters& params) {
    initialize_neutronics(params);
    initialize_thermal_hydraulics(params);
    initialize_heat_exchanger_1(params);
    initialize_heat_exchanger_2(params);
}
