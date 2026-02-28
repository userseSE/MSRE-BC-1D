#include "HX2.hpp"
#include "parameters.hpp"
#include "ode_solver_HX2.hpp"

//void pde_to_ode_hx2(float t, const float y[length_hx], float dydt[length_hx], Param_HX2 &params) {
//    // Split the input state vector y into u and v
//    const float* u = &y[0];    // u is the first half of y
//    const float* v = &y[Nx];   // v is the second half of y
//
//    float du_dt[Nx] = {0};    // Array to store the derivative of u
//    float dv_dt[Nx] = {0};    // Array to store the derivative of v
//
//    // Compute du_dt and dv_dt using the provided constants
//    for (int i = 0; i < Nx; ++i) {
//#pragma HLS UNROLL
//        for (int idx = params.A_HX2_csr.row_pointers[i]; idx < params.A_HX2_csr.row_pointers[i + 1]; ++idx) {
//#pragma HLS UNROLL
//            int j = params.A_HX2_csr.col_indices[idx];
//            du_dt[i] += params.A_HX2_csr.values[idx] * u[j];
//            dv_dt[i] += params.A_HX2_csr.values[idx] * v[j];
//        }
//        du_dt[i] += params.C2_2 * (u[i] - v[i]);
//        dv_dt[i] += params.C4_2 * (u[i] - v[i]);
//    }
//
//    // Apply time-varying boundary conditions
//    du_dt[0] = params.u2_L - u[0];
//    du_dt[Nx - 1] = params.u2_H - u[Nx - 1];
//    dv_dt[0] = 0; // Fixed condition: no change at v[0]
//    dv_dt[Nx - 1] = 0; // Fixed condition: no change at v[Nx-1]
//
//    // Populate dydt with the derivatives
//    for (int i = 0; i < Nx; ++i) {
//#pragma HLS UNROLL
//        dydt[i] = du_dt[i];        // First half for du_dt
//        dydt[Nx + i] = dv_dt[i];   // Second half for dv_dt
//    }
//}

void pde_to_ode_hx2(float t, const float y[length_hx], float dydt[length_hx], Param_HX2 &params) {
    // Split the input state vector y into u and v
    const float* u = &y[0];    // u is the first half of y
    const float* v = &y[Nx];   // v is the second half of y

    float du_dt[Nx] = {0};    // Array to store the derivative of u
    float dv_dt[Nx] = {0};    // Array to store the derivative of v
    float temp_dydt[length_hx] = {0}; // Local cache for dydt

    // Compute du_dt and dv_dt using the provided constants
    for (int i = 0; i < Nx; ++i) {
#pragma HLS UNROLL
        for (int idx = params.A_HX2_csr.row_pointers[i]; idx < params.A_HX2_csr.row_pointers[i + 1]; ++idx) {
#pragma HLS UNROLL
            int j = params.A_HX2_csr.col_indices[idx];
            du_dt[i] += params.A_HX2_csr.values[idx] * u[j];
            dv_dt[i] += params.A_HX2_csr.values[idx] * v[j];
        }
        du_dt[i] += params.C2_2 * (u[i] - v[i]);
        dv_dt[i] += params.C4_2 * (u[i] - v[i]);

        // Cache results in temp_dydt
        temp_dydt[i] = du_dt[i];          // First half for du_dt
        temp_dydt[Nx + i] = dv_dt[i];    // Second half for dv_dt
    }

    // Apply time-varying boundary conditions on temp_dydt
    temp_dydt[0] = params.u2_L - u[0];
    temp_dydt[Nx - 1] = params.u2_H - u[Nx - 1];
    temp_dydt[Nx] = 0; // Fixed condition: no change at v[0]
    temp_dydt[2 * Nx - 1] = 0; // Fixed condition: no change at v[Nx-1]

    // Write results back to dydt
    for (int i = 0; i < length_hx; ++i) {
#pragma HLS UNROLL
        dydt[i] = temp_dydt[i];
    }
}

void HX2(float y_hx2[length_hx], float Ts_HX2_L, int step) {
    Param_HX2 params;
    // Set boundary conditions
    y_hx2[Nx + 1] = params.v2_L;    // Fixed boundary condition at v[0]
    y_hx2[length_hx - 1] = params.v2_H; // Fixed boundary condition at v[L]

    // Initial condition vector y0
    float y0[length_hx];
    if (step == 0) {
        for (int i = 0; i < Nx; ++i) {
#pragma HLS UNROLL
            y0[i] = params.u2_init[i];
            y0[Nx + i] = params.v2_init[i];
        }
    } else {
        for (int i = 0; i < length_hx; ++i) {
#pragma HLS UNROLL
            y0[i] = y_hx2[i];
        }
    }
    y0[Nx - 1] = Ts_HX2_L;

    // Solve the system of ODEs
    ode_solver_hx2(y0, pde_to_ode_hx2, step, params);

    // Update y_hx2 with the solution from the ODE solver
    for (int i = 0; i < length_hx; ++i) {
#pragma HLS UNROLL
        y_hx2[i] = y0[i];  // Update the original array with the new values
    }
}
