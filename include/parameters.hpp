#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

// Spatial discretization constants
#define N 200  // spatial discretization
#define Nx 200   // spatial discretization
#define time_span 250
#define rho_insertion 0.0  // pcm, 50 * N
#define length_th 400  //2 * N
#define length_neutr 1600  //8 * N;
#define length_hx 400 // 2 * Nx;

// for CSR
const int MAX_NON_ZERO = 3 * N - 4; // Maximum number of non-zero entries in CSR
// Define CSR structure using arrays
struct CSRMatrix {
    float values[MAX_NON_ZERO];    // Non-zero values
    int col_indices[MAX_NON_ZERO];  // Column indices of non-zero values
    int row_pointers[N + 1];        // Row pointers (N + 1 for boundary condition)
};

struct Param_Neutronics {
    // seperate submodules into different structs
    // Neutronics
    static constexpr float dt = 0.01f;
    static constexpr float L = 172.0f;
    static constexpr float A = 4094.0f;
    static constexpr float flux_to_power = 3.12e10f;
    static constexpr float dz = 0.864321608040201f; // L / (N - 1);
    static constexpr float V1 = 1.103497e7f;   //1.103497 * 1e7;
    static constexpr float V2 = 2.2e5f;
    static constexpr float D1 = 6.74401f;  //0.96343 * 7;
    static constexpr float D2 = 0.282f;
    static constexpr float sigma_a1 = 0.002161939172413793f;  //0.002161939172413793;  0.001382
    static constexpr float sigma_a2 = 0.030f;       //0.0054869
    static constexpr float nu_sigma_f1 = 0.00220605f; //0.0044120763; 0.00440605
    static constexpr float nu_sigma_f2 = 0.00220f;
    static constexpr float sigma_s12 = 0.013f;      //0.0023
    static constexpr float sigma_f = 0.0018282365f;//(nu_sigma_f1 + nu_sigma_f2)/2.41f
    static constexpr float beta[] = {0.000228f, 0.000788f, 0.000664f, 0.000736f, 0.000136f, 0.000088f};
    static constexpr float Beta = 0.00264f; //0.000228f+ 0.000788f+ 0.000664f+ 0.000736f+ 0.000136f+ 0.000088f;
    static constexpr float lambda_i[] = {0.0126f, 0.0337f, 0.139f, 0.325f, 1.13f, 2.5f};

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
};
struct Param_Thermal {
    // Thermal-Hydraulics
    static constexpr float L = 172.0f;
    static constexpr float dz = 0.864321608040201f;
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

    static constexpr float a_th = 20.0; //Vc
    static constexpr float b_th = 0.0059478178f; //U_gs / (Ms * c_p_s);
    static constexpr float c_th = 0.00555722014f; //U_sg / (Mg * c_p_g);
    static constexpr float d_th = 0.000052856274f; //L * gamma / (Ms * c_p_s);
    static constexpr float e_th = 1.8585814053252681e-6f; //L * (1 - gamma) / (Mg * c_p_g);

    float initialS[N] = {0};
    float initialG[N] = {0};
    float AT[N][N] = {0};
    CSRMatrix AT_csr;
};
struct Param_HX1 {
    // Heat Exchanger 1
    static constexpr float L_HX = 200;
    static constexpr float dx = 1.0050251256281406; // L_HX / (Nx - 1)
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
    
    static constexpr float C1_1 = 171.2; // V_he_s;
    static constexpr float C2_1 = -0.11583983883152858f;// -U_hx / (M_he_s * c_p_s);
    static constexpr float C3_1 = 105.7f;// V_he_ss;
    static constexpr float C4_1 = 0.29291900152827305f; // U_hx / (M_he_ss * c_p_ss);

    float u_init[N] = {0};
    float v_init[N] = {0};
    float A_HX[Nx][Nx] = {0};
    CSRMatrix A_HX_csr;
};
struct Param_HX2 {
    // TODO: cosimulation - generate and use 3 sets of temperatures from brayton cycle, min, max, nominal in simulink
    // Heat Exchanger 2
    static constexpr float L_HX2 = 200.0;
    static constexpr float dx = 1.0050251256281406; // L_HX / (Nx - 1)
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

    static constexpr float C1_2 = 105.7; // V_he2_s;
    static constexpr float C2_2 = -0.1273560876209883f; // -U2_hx / (M_he2_s * c_p_ss);
    static constexpr float C3_2 = 45.3; // V_he2_ss;
    static constexpr float C4_2 = 0.04120879120879121f; // U2_hx / (M_he2_ss * c_p_sss);

    float u2_init[N] = {0};
    float v2_init[N] = {0};
    float A_HX2[Nx][Nx] = {0};
    CSRMatrix A_HX2_csr;
};
struct Param_React {
    // Reactivity
    static constexpr float alpha_f = -5.904e-5;
    static constexpr float alpha_g = -6.624e-5;
    static constexpr float tau_l = 16.73;
    static constexpr float tau_c = 8.46;
    static constexpr float max_rho_change = 1e-3;
    static constexpr float rho_0_value = 0.0; //3.3e-5;

    static constexpr float bc_s0 = 900;
    static constexpr float bc_sL = 1000;
    static constexpr float bc_g0 = 950;
    static constexpr float bc_gL = 1050;
    static constexpr float L = 172.0f;

    float initialS[N] = {0};
    float initialG[N] = {0};
};

    // Initial Conditions
    // static constexpr float Ts_in = bc_s0;
    // static constexpr float Ts_out = bc_sL;
    // static constexpr float Tss_in = v_L;
    // static constexpr float Tss_out = v_H;
    // static constexpr float Tsss_in = v2_L;
    // static constexpr float Tsss_out = v2_H;

    //ode_solver

void initialize_neutronics(Param_Neutronics& params);
void initialize_thermal_hydraulics(Param_Thermal& params);
void initialize_heat_exchanger_1(Param_HX1& params);
void initialize_heat_exchanger_2(Param_HX2& params);
void initialize_reactivity(Param_React& params);

#endif // PARAMETERS_HPP
