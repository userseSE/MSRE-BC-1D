#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <Eigen/Dense>
#include <array>
#include <numeric>

const int N = 200;  // spatial discretization
const int Nx = N;   // spatial discretization

// Neutronics
extern double dt;
extern const double L;
extern const double A;
extern const double flux_to_power;
extern const double dz;
extern const double V;
extern const double D;
extern const double sigma_a;
extern double nu_sigma_f;  // Remove const to make it modifiable
extern const double sigma_f;
extern const std::array<double, 6> beta;
extern const double Beta;
extern const std::array<double, 6> lambda_i;
extern Eigen::VectorXd phi_0;  // Use Eigen::VectorXd
extern Eigen::VectorXd c1;
extern Eigen::VectorXd c2;
extern Eigen::VectorXd c3;
extern Eigen::VectorXd c4;
extern Eigen::VectorXd c5;
extern Eigen::VectorXd c6;
extern double t0;
extern double t1;

// Thermal-Hydraulics
extern const double c_p_s;
extern const double Vc;
extern const double Ms;
extern const double Mg;
extern const double gamma;
extern const double U_sg;
extern const double U_gs;
extern const double c_p_g;
extern const double bc_s0;
extern const double bc_sL;
extern const double bc_g0;
extern const double bc_gL;
extern Eigen::VectorXd initialS;  // Use Eigen::VectorXd
extern Eigen::VectorXd initialG;  // Use Eigen::VectorXd

// Heat Exchanger 1
extern const double L_HX;
extern const double dx;
extern const double V_he_s;
extern const double V_he_ss;
extern const double U_hx;
extern const double M_he_s;
extern const double M_he_ss;
extern const double c_p_ss;
extern const double u_L;
extern const double u_H;
extern const double v_L;
extern const double v_H;
extern Eigen::VectorXd u_init;  // Use Eigen::VectorXd
extern Eigen::VectorXd v_init;  // Use Eigen::VectorXd

// Heat Exchanger 2
extern const double L_HX2;
extern const double V_he2_s;
extern const double V_he2_ss;
extern const double U2_hx;
extern const double M_he2_s;
extern const double M_he2_ss;
extern const double c_p_sss;
extern const double u2_L;
extern const double u2_H;
extern const double v2_L;
extern const double v2_H;
extern Eigen::VectorXd u2_init;  // Use Eigen::VectorXd
extern Eigen::VectorXd v2_init;  // Use Eigen::VectorXd

// Reactivity
extern const double alpha_f;
extern const double alpha_g;
extern const double tau_l;
extern const double tau_c;
extern const double max_rho_change;
extern double rho_0_value;  // Remove const to make it modifiable

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

void initialize_parameters();

#endif // PARAMETERS_HPP
