#include "parameters.hpp"
#include <cmath>    // for sin function
#include <numeric>  // for accumulate function
#include <algorithm> // for fill_n
#include <vector>
#include <iostream>

// Neutronics
const double dt = 0.3;
const double L = 172;
const double dz = L / (N - 1);
const double V = 4e7;
const double D = 0.390016;
const double sigma_a = 0.0835;
const double nu_sigma_f = 3.33029e-2;
const std::array<double, 6> beta = {0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088};
const double Beta = std::accumulate(beta.begin(), beta.end(), 0.0);
const double delta = Beta * nu_sigma_f;
const std::array<double, 6> lambda_i = {0.0126, 0.0337, 0.139, 0.325, 1.13, 2.5};
double phi_0[N];
double c1[N];
double c2[N];
double c3[N];
double c4[N];
double c5[N];
double c6[N];

void initialize_neutronics() {
    // std::fill_n(phi_0, N, 522654);  // Initialize phi_0 with 522654
    double lambda_sum = std::accumulate(lambda_i.begin(), lambda_i.end(), 0.0);

    for (int i = 0; i < N; ++i) {
        phi_0[i] = 522654 * sin(M_PI * i / (N - 1));
        c1[i] = (beta[0] * nu_sigma_f / lambda_i[0]) * phi_0[i];
        c2[i] = (beta[1] * nu_sigma_f / lambda_i[1]) * phi_0[i];
        c3[i] = (beta[2] * nu_sigma_f / lambda_i[2]) * phi_0[i];
        c4[i] = (beta[3] * nu_sigma_f / lambda_i[3]) * phi_0[i];
        c5[i] = (beta[4] * nu_sigma_f / lambda_i[4]) * phi_0[i];
        c6[i] = (beta[5] * nu_sigma_f / lambda_i[5]) * phi_0[i];
    }
}

// Thermal-Hydraulics
const double c_p_s = 1983;
const double Vc = 0.2;
const double Ms = 1448;
const double Mg = 3687;
const double gamma = 0.93;
const double U = 36000;
const double c_p_g = 1757;
const double bc_s0 = 910;
const double bc_sL = 958.15;
const double bc_g0 = 920;
const double bc_gL = 968.71;
double initialS[N];
double initialG[N];

void initialize_thermal_hydraulics() {
    for (int i = 0; i < N; ++i) {
        double position = static_cast<double>(i) * L / (N - 1);
        initialS[i] = bc_s0 + (bc_sL - bc_s0) * (0.5 + 0.5 * std::sin(M_PI * (position - 0.3 * L) / L));
        initialG[i] = bc_g0 + (bc_gL - bc_g0) * (0.5 + 0.5 * std::sin(M_PI * (position - 0.3 * L) / L));
    }
}

// Heat Exchanger 1
const double L_HX = 2;
const double dx = L / (N - 1);
const double V_he_s = 75.7 * 1e-3 / 23.6;
const double V_he_ss = 53.6 * 1e-3 / 23.6;
const double U_hx = 82800;
const double M_he_s = 342;
const double M_he_ss = 117;
const double c_p_ss = 2416;
// Initial conditions
const double u_L = 908.15;
const double u_H = 958;
const double v_L = 824.85;
const double v_H = 866.45;
double u_init[Nx];
double v_init[Nx];

void initialize_heat_exchanger_1() {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * L_HX / (Nx - 1);
        u_init[i] = u_L + (u_L - u_H) * (0.5 + 0.5 * std::sin(M_PI * (position / L_HX - 0.5)));
        v_init[i] = v_L + (v_L - v_H) * (0.5 + 0.5 * std::sin(M_PI * (position / L_HX - 0.5)));
    }
}

// Heat Exchanger 2
const double L_HX2 = 2;
const double V_he2_s = 53.6 * 1e-3 / 23.6;
const double V_he2_ss = 33.6 * 1e-3 / 23.6;
const double U2_hx = 82800;
const double M_he2_s = 117;
const double M_he2_ss = 100;
const double c_p_sss = 2416;
// Initial conditions
const double u2_L = 824;
const double u2_H = 866;
const double v2_L = 744;
const double v2_H = 786;
double u2_init[Nx];
double v2_init[Nx];

void initialize_heat_exchanger_2() {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * L_HX2 / (Nx - 1);
        u2_init[i] = u2_L + (u2_H - u2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / L_HX2 - 0.5)));
        v2_init[i] = v2_L + (v2_H - v2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / L_HX2 - 0.5)));
    }
}

// Reactivity
const double alpha_f = -5.904E-5;
const double alpha_g = -6.624E-5;
const double tau_l = 16.73;
const double tau_c = 8.46;

// Transport Delays
// const double tau_hx_c = 9;
// const double tau_c_hx = 4;
// const double tau_hx_r = 5;
// const double tau_r_hx = 8;
// const double tau_r_pp = 10;
// const double tau_pp_r = 10;

const double tau_hx_c = 5;
const double tau_c_hx = 2;
const double tau_hx_r = 3;
const double tau_r_hx = 4;
const double tau_r_pp = 5;
const double tau_pp_r = 5;

// Initial Conditions
const double Ts_in = bc_s0;
const double Ts_out = bc_sL;
const double Tss_in = v_init[0];
const double Tss_out = v_init[Nx - 1];
const double Tsss_in = v2_init[0];
const double Tsss_out = v2_init[Nx - 1];

// Initialization function to call all individual initializations
void initialize_parameters() {
    initialize_neutronics();
    initialize_thermal_hydraulics();
    initialize_heat_exchanger_1();
    initialize_heat_exchanger_2();
}
