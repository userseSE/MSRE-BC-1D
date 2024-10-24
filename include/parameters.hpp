#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

// Spatial discretization constants
constexpr int N = 200;  // spatial discretization
constexpr int Nx = N;   // spatial discretization
constexpr int time_span = 1000;
constexpr double rho_insertion = 0.0;  // pcm, 50 * N
constexpr int length_th = 2 * N;
constexpr int length_neutr = 8 * N;
constexpr int length_hx = 2 * Nx;

struct Parameters {
    // seperate submodules into different structs
    // Neutronics
    static constexpr double dt = 0.01;
    static constexpr double L = 172;
    static constexpr double A = 4094;
    static constexpr double flux_to_power = 3.12e10;
    static constexpr double dz = L / (N - 1);
    static constexpr double V1 = 1.103497 * 1e7;   //1.103497 * 1e7;
    static constexpr double V2 = 2.2e5;
    static constexpr double D1 = 0.96343 * 7;  //0.96343 * 7;
    static constexpr double D2 = 0.282;
    static constexpr double sigma_a1 = 0.002161939172413793;  //0.002161939172413793;  0.001382
    static constexpr double sigma_a2 = 0.030;       //0.0054869
    static constexpr double nu_sigma_f1 = 0.00220605; //0.0044120763; 0.00440605
    static constexpr double nu_sigma_f2 = 0.00220;
    static constexpr double sigma_s12 = 0.013;      //0.0023
    static constexpr double sigma_f = (nu_sigma_f1 + nu_sigma_f2)/2.41;//0.004411764705882353 / 2.41;
    static constexpr double beta[] = {0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088};
    static constexpr double Beta = 0.000228+ 0.000788+ 0.000664+ 0.000736+ 0.000136+ 0.000088;
    // static constexpr double Beta = 0.000228 + 0.000788 + 0.000664 + 0.000736 + 0.000136 + 0.000088;
    static constexpr double lambda_i[] = {0.0126, 0.0337, 0.139, 0.325, 1.13, 2.5};

    static constexpr double t0 = 0;
    static constexpr double t1 = 1;

    // TODO: init put into another struct, initial conditions set as const
    double phi1_0[N] = {0};
    double phi2_0[N] = {0};
    double c1[N] = {0};
    double c2[N] = {0};
    double c3[N] = {0};
    double c4[N] = {0};
    double c5[N] = {0};
    double c6[N] = {0};

    double A1[N][N] = {0};
    double A2[N][N] = {0};
    double B1[N][N] = {0};
    double B2[N][N] = {0};

    // Thermal-Hydraulics
    static constexpr double c_p_s = 2090;
    static constexpr double Vc = 20;
    static constexpr double Ms = 1448;
    static constexpr double Mg = 3687;
    static constexpr double gamma = 0.93;
    static constexpr double U_sg = 36000;
    static constexpr double U_gs = 18000;
    static constexpr double c_p_g = 1757;
    static constexpr double bc_s0 = 900;
    static constexpr double bc_sL = 1000;
    static constexpr double bc_g0 = 950;
    static constexpr double bc_gL = 1050;
    static constexpr double fixed_boundary_sL = 950.0;
    static constexpr double fixed_boundary_gL = 970.0;
    // static constexpr double bc_s0 = 500;
    // static constexpr double bc_sL = 600;
    // static constexpr double bc_g0 = 550;
    // static constexpr double bc_gL = 650;
    static constexpr double a_th = Vc;
    static constexpr double b_th = U_gs / (Ms * c_p_s);
    static constexpr double c_th = U_sg / (Mg * c_p_g);
    static constexpr double d_th = L * gamma / (Ms * c_p_s);
    static constexpr double e_th = L * (1 - gamma) / (Mg * c_p_g);

    double initialS[N] = {0};
    double initialG[N] = {0};
    double AT[N][N] = {0};

    // Heat Exchanger 1
    static constexpr double L_HX = 200;
    static constexpr double dx = L_HX / (N - 1);
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

    double u_init[N] = {0};
    double v_init[N] = {0};
    double A_HX[Nx][Nx] = {0};

    // TODO: cosimulation - generate and use 3 sets of temperatures from brayton cycle, min, max, nominal in simulink
    // Heat Exchanger 2
    static constexpr double L_HX2 = 200;
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

    double u2_init[N] = {0};
    double v2_init[N] = {0};
    double A_HX2[Nx][Nx] = {0};

    // Reactivity
    static constexpr double alpha_f = -5.904e-5;
    static constexpr double alpha_g = -6.624e-5;
    static constexpr double tau_l = 16.73;
    static constexpr double tau_c = 8.46;
    static constexpr double max_rho_change = 1e-3;
    static constexpr double rho_0_value = 0.0; //3.3e-5;

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

    //ode_solver

};

void initialize_parameters(Parameters& params);

#endif // PARAMETERS_HPP
