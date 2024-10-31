#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

// Spatial discretization constants
constexpr int N = 200;  // spatial discretization
constexpr int Nx = N;   // spatial discretization
constexpr int time_span = 250;
constexpr float rho_insertion = 0.0;  // pcm, 50 * N
constexpr int length_th = 2 * N;
constexpr int length_neutr = 8 * N;
constexpr int length_hx = 2 * Nx;

// for CSR
const int MAX_NON_ZERO = 3 * N - 4; // Maximum number of non-zero entries in CSR
// Define CSR structure using arrays
struct CSRMatrix {
    float values[MAX_NON_ZERO];    // Non-zero values
    int col_indices[MAX_NON_ZERO];  // Column indices of non-zero values
    int row_pointers[N + 1];        // Row pointers (N + 1 for boundary condition)
};

struct Parameters {
    // seperate submodules into different structs
    // Neutronics
    static constexpr float dt = 0.01f;
    static constexpr float L = 172.0f;
    static constexpr float A = 4094.0f;
    static constexpr float flux_to_power = 3.12e10f;
    static constexpr float dz = L / (N - 1);
    static constexpr float V1 = 1.103497e7f;   //1.103497 * 1e7;
    static constexpr float V2 = 2.2e5f;
    static constexpr float D1 = 6.74401f;  //0.96343 * 7;
    static constexpr float D2 = 0.282f;
    static constexpr float sigma_a1 = 0.002161939172413793f;  //0.002161939172413793;  0.001382
    static constexpr float sigma_a2 = 0.030f;       //0.0054869
    static constexpr float nu_sigma_f1 = 0.00220605f; //0.0044120763; 0.00440605
    static constexpr float nu_sigma_f2 = 0.00220f;
    static constexpr float sigma_s12 = 0.013f;      //0.0023
    static constexpr float sigma_f = (nu_sigma_f1 + nu_sigma_f2)/2.41f;//0.004411764705882353 / 2.41;
    static constexpr float beta[] = {0.000228f, 0.000788f, 0.000664f, 0.000736f, 0.000136f, 0.000088f};
    static constexpr float Beta = 0.000228f+ 0.000788f+ 0.000664f+ 0.000736f+ 0.000136f+ 0.000088f;
    static constexpr float lambda_i[] = {0.0126f, 0.0337f, 0.139f, 0.325f, 1.13f, 2.5f};

    static constexpr float t0 = 0;
    static constexpr float t1 = 1;

    // TODO: init put into another struct, initial conditions set as const
    float phi1_0[N] = {0};
    float phi2_0[N] = {0};
    float c1[N] = {0};
    float c2[N] = {0};
    float c3[N] = {0};
    float c4[N] = {0};
    float c5[N] = {0};
    float c6[N] = {0};

    float A1[N][N] = {0};
    float A2[N][N] = {0};
    float B1[N][N] = {0};
    float B2[N][N] = {0};
    CSRMatrix B1_csr, B2_csr;

    // Thermal-Hydraulics
    static constexpr float c_p_s = 2090;
    static constexpr float Vc = 20;
    static constexpr float Ms = 1448;
    static constexpr float Mg = 3687;
    static constexpr float gamma = 0.93;
    static constexpr float U_sg = 36000;
    static constexpr float U_gs = 18000;
    static constexpr float c_p_g = 1757;
    static constexpr float bc_s0 = 900;
    static constexpr float bc_sL = 1000;
    static constexpr float bc_g0 = 950;
    static constexpr float bc_gL = 1050;
    static constexpr float fixed_boundary_sL = 950.0;
    static constexpr float fixed_boundary_gL = 970.0;
    // static constexpr float bc_s0 = 500;
    // static constexpr float bc_sL = 600;
    // static constexpr float bc_g0 = 550;
    // static constexpr float bc_gL = 650;
    static constexpr float a_th = Vc;
    static constexpr float b_th = U_gs / (Ms * c_p_s);
    static constexpr float c_th = U_sg / (Mg * c_p_g);
    static constexpr float d_th = L * gamma / (Ms * c_p_s);
    static constexpr float e_th = L * (1 - gamma) / (Mg * c_p_g);

    float initialS[N] = {0};
    float initialG[N] = {0};
    float AT[N][N] = {0};
    CSRMatrix AT_csr;

    // Heat Exchanger 1
    static constexpr float L_HX = 200;
    static constexpr float dx = L_HX / (N - 1);
    static constexpr float V_he_s = 171.2;
    static constexpr float V_he_ss = 105.7;
    static constexpr float U_hx = 82800;
    static constexpr float M_he_s = 342;
    static constexpr float M_he_ss = 117;
    static constexpr float c_p_ss = 2416;
    static constexpr float u_L = 900;
    static constexpr float u_H = 950;
    static constexpr float v_L = 850;
    static constexpr float v_H = 900;
    static constexpr float C1_1 = V_he_s;
    static constexpr float C2_1 = -U_hx / (M_he_s * c_p_s);
    static constexpr float C3_1 = V_he_ss;
    static constexpr float C4_1 = U_hx / (M_he_ss * c_p_ss);

    float u_init[N] = {0};
    float v_init[N] = {0};
    float A_HX[Nx][Nx] = {0};
    CSRMatrix A_HX_csr;

    // TODO: cosimulation - generate and use 3 sets of temperatures from brayton cycle, min, max, nominal in simulink
    // Heat Exchanger 2
    static constexpr float L_HX2 = 200;
    static constexpr float V_he2_s = 105.7;
    static constexpr float V_he2_ss = 45.3;
    static constexpr float U2_hx = 36000;
    static constexpr float M_he2_s = 117;
    static constexpr float M_he2_ss = 672;
    static constexpr float c_p_sss = 1300;
    static constexpr float u2_L = 850;
    static constexpr float u2_H = 900;
    static constexpr float v2_L = 743;
    static constexpr float v2_H = 873;
    static constexpr float C1_2 = V_he2_s;
    static constexpr float C2_2 = -U2_hx / (M_he2_s * c_p_ss);
    static constexpr float C3_2 = V_he2_ss;
    static constexpr float C4_2 = U2_hx / (M_he2_ss * c_p_sss);

    float u2_init[N] = {0};
    float v2_init[N] = {0};
    float A_HX2[Nx][Nx] = {0};
    CSRMatrix A_HX2_csr;

    // Reactivity
    static constexpr float alpha_f = -5.904e-5;
    static constexpr float alpha_g = -6.624e-5;
    static constexpr float tau_l = 16.73;
    static constexpr float tau_c = 8.46;
    static constexpr float max_rho_change = 1e-3;
    static constexpr float rho_0_value = 0.0; //3.3e-5;

    // Transport Delays
    static constexpr float tau_hx_c = 0;
    static constexpr float tau_c_hx = 0;
    static constexpr float tau_hx_r = 0;
    static constexpr float tau_r_hx = 0;
    static constexpr float tau_r_pp = 0;
    static constexpr float tau_pp_r = 0;

    // Initial Conditions
    static constexpr float Ts_in = bc_s0;
    static constexpr float Ts_out = bc_sL;
    static constexpr float Tss_in = v_L;
    static constexpr float Tss_out = v_H;
    static constexpr float Tsss_in = v2_L;
    static constexpr float Tsss_out = v2_H;

    //ode_solver

};

void initialize_parameters(Parameters& params);

#endif // PARAMETERS_HPP
