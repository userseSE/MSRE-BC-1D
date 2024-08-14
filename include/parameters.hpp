#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <array>

constexpr int N = 200;  // spatial descretization
constexpr int Nx = N;  // spatial descretization

// Neutronics
extern const double dt;
extern const double L;
extern const int N;
extern const double dz;
extern const double V;
extern const double D;
extern const double sigma_a;
extern const double nu_sigma_f;
extern const std::array<double, 6> beta;
extern const double Beta;
extern const double delta;
extern const std::array<double, 6> lambda_i;
extern double phi_0[N];
extern double c0[N];

// Thermal-Hydraulics
extern const double c_p_s;
extern const double Vc;
extern const double Ms;
extern const double Mg;
extern const double gamma;
extern const double U;
extern const double c_p_g;
extern const double bc_s0;
extern const double bc_sL;
extern const double bc_g0;
extern const double bc_gL;
extern double initialS[N];
extern double initialG[N];

// Heat Exchanger 1
extern const int Nx;
extern const double L_HX;
extern const double dx;
extern const double V_he_s;
extern const double V_he_ss;
extern const double U_hx;
extern const double M_he_s;
extern const double M_he_ss;
extern const double c_p_ss;
extern double u_init[Nx];
extern double v_init[Nx];

// Heat Exchanger 2
extern const double L_HX2;
extern const double V_he2_s;
extern const double V_he2_ss;
extern const double U2_hx;
extern const double M_he2_s;
extern const double M_he2_ss;
extern const double c_p_sss;
extern double u2_init[Nx];
extern double v2_init[Nx];

// Reactivity
extern const double alpha_f;
extern const double alpha_g;
extern const double tau_l;
extern const double tau_c;

// Transport Delays
extern const double tau_hx_c;
extern const double tau_c_hx;
extern const double tau_hx_r;
extern const double tau_r_hx;
extern const double tau_r_pp;
extern const double tau_pp_r;

// Initial Conditions
extern const double Ts_in;
extern const double Ts_out;
extern const double Tss_in;
extern const double Tss_out;
extern const double Tsss_in;
extern const double Tsss_out;

#endif // PARAMETERS_HPP
