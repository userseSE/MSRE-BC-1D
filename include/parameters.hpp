#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <Eigen/Dense>
#include <array>
#include <numeric>
#include <vector>

// Spatial discretization constants
constexpr int N = 200;  // spatial discretization
constexpr int Nx = N;   // spatial discretization

struct Parameters {
    // seperate submodules into different structs
    // Neutronics
    static constexpr double dt = 0.01;
    static constexpr double L = 172;
    static constexpr double A = 4094;
    static constexpr double flux_to_power = 3.12e10;
    static constexpr double dz = L / (N - 1);
    static constexpr double V = 1.103497 * 1e7;
    static constexpr double D = 0.96343 * 7;
    static constexpr double sigma_a = 0.002161939172413793;
    static constexpr double nu_sigma_f = 0.0044120763;
    static constexpr double sigma_f = 0.004411764705882353 / 2.41;
    static constexpr std::array<double, 6> beta = {0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088};
    static constexpr double Beta = std::accumulate(beta.begin(), beta.end(), 0.0);
    // static constexpr double Beta = 0.000228 + 0.000788 + 0.000664 + 0.000736 + 0.000136 + 0.000088;
    static constexpr std::array<double, 6> lambda_i = {0.0126, 0.0337, 0.139, 0.325, 1.13, 2.5};

    static constexpr double t0 = 0;
    static constexpr double t1 = 1;

    // TODO: init put into another struct, initial conditions set as const
    Eigen::VectorXd phi_0 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c1 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c2 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c3 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c4 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c5 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c6 = Eigen::VectorXd::Zero(N);

    // Thermal-Hydraulics
    static constexpr double c_p_s = 2090;
    static constexpr double Vc = 20;
    static constexpr double Ms = 1448;
    static constexpr double Mg = 3687;
    static constexpr double gamma = 0.93;
    static constexpr double U_sg = 36000;
    static constexpr double U_gs = 36000;
    static constexpr double c_p_g = 1757;
    static constexpr double bc_s0 = 800;
    static constexpr double bc_sL = 900;
    static constexpr double bc_g0 = 850;
    static constexpr double bc_gL = 950;
    static constexpr double a_th = Vc;
    static constexpr double b_th = U_gs / (Ms * c_p_s);
    static constexpr double c_th = U_sg / (Mg * c_p_g);
    static constexpr double d_th = L * gamma / (Ms * c_p_s);
    static constexpr double e_th = L * (1 - gamma) / (Mg * c_p_g);

    Eigen::VectorXd initialS = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd initialG = Eigen::VectorXd::Zero(N);

    // Heat Exchanger 1
    static constexpr double L_HX = 2;
    static constexpr double dx = L / (N - 1);
    static constexpr double V_he_s = 171.2;
    static constexpr double V_he_ss = 105.7;
    static constexpr double U_hx = 82800;
    static constexpr double M_he_s = 342;
    static constexpr double M_he_ss = 117;
    static constexpr double c_p_ss = 2416;
    static constexpr double u_L = 900;
    static constexpr double u_H = 950;
    static constexpr double v_L = 850;
    static constexpr double v_H = 900;
    static constexpr double C1_1 = V_he_s;
    static constexpr double C2_1 = -U_hx / (M_he_s * c_p_s);
    static constexpr double C3_1 = V_he_ss;
    static constexpr double C4_1 = U_hx / (M_he_ss * c_p_ss);

    Eigen::VectorXd u_init = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd v_init = Eigen::VectorXd::Zero(N);

    // TODO: cosimulation - generate and use 3 sets of temperatures from brayton cycle, min, max, nominal in simulink
    // Heat Exchanger 2
    static constexpr double L_HX2 = 2;
    static constexpr double V_he2_s = 105.7;
    static constexpr double V_he2_ss = 45.3;
    static constexpr double U2_hx = 36000;
    static constexpr double M_he2_s = 117;
    static constexpr double M_he2_ss = 672;
    static constexpr double c_p_sss = 1300;
    static constexpr double u2_L = 850;
    static constexpr double u2_H = 900;
    static constexpr double v2_L = 743;
    static constexpr double v2_H = 873;
    static constexpr double C1_2 = V_he2_s;
    static constexpr double C2_2 = -U2_hx / (M_he2_s * c_p_ss);
    static constexpr double C3_2 = V_he2_ss;
    static constexpr double C4_2 = U2_hx / (M_he2_ss * c_p_sss);

    Eigen::VectorXd u2_init = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd v2_init = Eigen::VectorXd::Zero(N);

    // Reactivity
    static constexpr double alpha_f = -5.904e-5;
    static constexpr double alpha_g = -6.624e-5;
    static constexpr double tau_l = 16.73;
    static constexpr double tau_c = 8.46;
    static constexpr double max_rho_change = 1e-4;
    double rho_0_value = 3.3e-5;

    // Transport Delays
    static constexpr double tau_hx_c = 0;
    static constexpr double tau_c_hx = 0;
    static constexpr double tau_hx_r = 0;
    static constexpr double tau_r_hx = 0;
    static constexpr double tau_r_pp = 0;
    static constexpr double tau_pp_r = 0;

    // Initial Conditions
    static constexpr double Ts_in = bc_s0;
    static constexpr double Ts_out = bc_sL;
    static constexpr double Tss_in = v_L;
    static constexpr double Tss_out = v_H;
    static constexpr double Tsss_in = v2_L;
    static constexpr double Tsss_out = v2_H;
};

void initialize_parameters(Parameters& params);

#endif // PARAMETERS_HPP
