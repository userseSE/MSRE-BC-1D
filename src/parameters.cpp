#include <Eigen/Dense>
#include <array>
#include <numeric>
#include <vector>
#include <cmath>

#include "parameters.hpp"

// TODO: try modify sigma_a and nu_sigma_f, and modify initial conditions figures
// Neutronics
double dt = 0.001;
const double L = 172;
const double dz = L / (N - 1);
const double A = 4094;
const double flux_to_power = 3.12e10;
const double V = 1.103497 * 1e7;
const double D = 0.96343 * 7;
const double sigma_a = 0.002161939172413793;    // 0.002161939172413793, smaller-explode, 0.002340 diverge, 0.002341 converge to 0
double nu_sigma_f = 0.0044122185;     // 0.004417172-0.00441718, 0.004411764705882353
const double sigma_f = 0.004411764705882353 / 2.41;
const std::array<double, 6> beta = {0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088};
const double Beta = std::accumulate(beta.begin(), beta.end(), 0.0);
const std::array<double, 6> lambda_i = {0.0126, 0.0337, 0.139, 0.325, 1.13, 2.5};

double t0 = 0;
double t1 = 1;

Eigen::VectorXd phi_0 = Eigen::VectorXd::Zero(N);  // Initial neutron flux
Eigen::VectorXd c1 = Eigen::VectorXd::Zero(N);
Eigen::VectorXd c2 = Eigen::VectorXd::Zero(N);
Eigen::VectorXd c3 = Eigen::VectorXd::Zero(N);
Eigen::VectorXd c4 = Eigen::VectorXd::Zero(N);
Eigen::VectorXd c5 = Eigen::VectorXd::Zero(N);
Eigen::VectorXd c6 = Eigen::VectorXd::Zero(N);

void initialize_neutronics() {
    double lambda_sum = std::accumulate(lambda_i.begin(), lambda_i.end(), 0.0);

    for (int i = 0; i < N; ++i) {
        phi_0[i] = 1.0e13;  // Initial neutron flux
        double x = i * dz;  // Compute the spatial position
        // phi_0[i] = 1e13 * std::sin(M_PI * x / L);  // Sine wave initialization, scaled by 1e13
        c1[i] = ((beta[0] * nu_sigma_f) / (lambda_i[0]/6)) * phi_0[i];
        c2[i] = ((beta[1] * nu_sigma_f) / (lambda_i[1]/6)) * phi_0[i];
        c3[i] = ((beta[2] * nu_sigma_f) / (lambda_i[2]/6)) * phi_0[i];
        c4[i] = ((beta[3] * nu_sigma_f) / (lambda_i[3]/6)) * phi_0[i];
        c5[i] = ((beta[4] * nu_sigma_f) / (lambda_i[4]/6)) * phi_0[i];
        c6[i] = ((beta[5] * nu_sigma_f) / (lambda_i[5]/6)) * phi_0[i];  
    }
}

// Thermal-Hydraulics
const double c_p_s = 2090;
const double Vc = 20;   
const double Ms = 1448;
const double Mg = 3687;
const double gamma = 0.93;
const double U_sg = 36000;
const double U_gs = 36000;
const double c_p_g = 1757;

// const double bc_s0 = 700;
// const double bc_sL = 800;
// const double bc_g0 = 700;
// const double bc_gL = 850;
const double bc_s0 = 500;
const double bc_sL = 600;
const double bc_g0 = 500;
const double bc_gL = 650;

Eigen::VectorXd initialS = Eigen::VectorXd::Zero(N);  // Initial temperature in salt
Eigen::VectorXd initialG = Eigen::VectorXd::Zero(N);  // Initial temperature in graphite

void initialize_thermal_hydraulics() {
    for (int i = 0; i < N; ++i) {
        double position = static_cast<double>(i) * L / (N - 1);
        initialS[i] = bc_s0 + (bc_sL - bc_s0) * ((0.5 + 0.5 * std::sin(M_PI * (position) / (L * 2))) * 0.8);
        initialG[i] = bc_g0 + (bc_gL - bc_g0) * ((0.5 + 0.5 * std::sin(M_PI * (position) / (L * 2))) * 1.05);
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

const double u_L = 700;
const double u_H = 850;
const double v_L = 700;
const double v_H = 800;

// const double u_L = 40;
// const double u_H = 40;
// const double v_L = 42;
// const double v_H = 46;

Eigen::VectorXd u_init = Eigen::VectorXd::Zero(Nx);  // Initial temperature profile for fluid 1
Eigen::VectorXd v_init = Eigen::VectorXd::Zero(Nx);  // Initial temperature profile for fluid 2

void initialize_heat_exchanger_1() {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * L_HX / (Nx - 1);
        u_init[i] = u_L + (u_L - u_H) * (0.5 + 0.5 * std::sin(M_PI * (position / L_HX)));
        v_init[i] = v_L + (v_L - v_H) * (0.5 + 0.5 * std::sin(M_PI * (position / L_HX)) * 1.05);
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

const double u2_L = 300;
const double u2_H = 400;
const double v2_L = 200;
const double v2_H = 250;

// const double u2_L = 1;
// const double u2_H = 1;
// const double v2_L = 1;
// const double v2_H = 1;

Eigen::VectorXd u2_init = Eigen::VectorXd::Zero(Nx);  // Initial temperature profile for fluid 1 in exchanger 2
Eigen::VectorXd v2_init = Eigen::VectorXd::Zero(Nx);  // Initial temperature profile for fluid 2 in exchanger 2

void initialize_heat_exchanger_2() {
    for (int i = 0; i < Nx; ++i) {
        double position = static_cast<double>(i) * L_HX2 / (Nx - 1);
        u2_init[i] = u2_L + (u2_H - u2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / L_HX2)));
        v2_init[i] = v2_L + (v2_H - v2_L) * (0.5 + 0.7 * std::sin(M_PI * (position / L_HX2)) * 1.05);
    }
}

// Reactivity
const double alpha_f =  - 5.904e-5;
const double alpha_g =  - 6.624e-5;
const double tau_l = 16.73;
const double tau_c = 8.46;
const double max_rho_change = 1e-4;
double rho_0_value = 0.0;

// Transport Delays
const double tau_hx_c = 9;
const double tau_c_hx = 4;
const double tau_hx_r = 5;
const double tau_r_hx = 8;
const double tau_r_pp = 10;
const double tau_pp_r = 10;

// const double tau_hx_c = 4;
// const double tau_c_hx = 4;
// const double tau_hx_r = 5;
// const double tau_r_hx = 8;
// const double tau_r_pp = 6;
// const double tau_pp_r = 6;

// Initial Conditions
const double Ts_in = bc_s0;
const double Ts_out = bc_sL;
const double Tss_in = v_L;
const double Tss_out = v_H;
const double Tsss_in = v2_L;
const double Tsss_out = v2_H;

// Initialization function to call all individual initializations
void initialize_parameters() {
    initialize_neutronics();
    initialize_thermal_hydraulics();
    initialize_heat_exchanger_1();
    initialize_heat_exchanger_2();
}
