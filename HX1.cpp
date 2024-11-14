#include "HX1.hpp"
#include "parameters.hpp"
#include "ode_solver_HX1.hpp"

void pde_to_ode_hx1(float t, const float y[length_hx], float dydt[length_hx], Param_HX1 &params) {
    // Split the input state vector y into u and v
    const float* u = &y[0];    // u is the first half of y
    const float* v = &y[Nx];   // v is the second half of y

    float du_dt[Nx] = {0};    // Array to store the derivative of u
    float dv_dt[Nx] = {0};    // Array to store the derivative of v

    // Compute du_dt and dv_dt using the provided constants
    for (int i = 0; i < Nx; ++i) {
#pragma HLS UNROLL
        for (int idx = params.A_HX_csr.row_pointers[i]; idx < params.A_HX_csr.row_pointers[i + 1]; ++idx) {
#pragma HLS UNROLL
            int j = params.A_HX_csr.col_indices[idx];
            du_dt[i] += params.A_HX_csr.values[idx] * u[j];
            dv_dt[i] += params.A_HX_csr.values[idx] * v[j];
        }
        du_dt[i] += params.C2_1 * (u[i] - v[i]);
        dv_dt[i] += params.C4_1 * (u[i] - v[i]);
        dydt[i] = du_dt[i];        // First half for du_dt
        dydt[Nx + i] = dv_dt[i];   // Second half for dv_dt
    }

    // Apply time-varying boundary conditions
    dydt[0] = params.u_L - dydt[0];
    dydt[Nx - 1] = params.u_H - dydt[Nx - 1];
    dydt[Nx] = params.v_L - dydt[Nx];
    dydt[2*Nx - 1] = params.v_H - dydt[2*Nx - 1];
}

void HX1(float y_hx1[length_hx], float Ts_HX1_L, int step) {
    Param_HX1 params;

    // Initial condition vector y0
    float y0[length_hx];
    if (step == 0) {
        for (int i = 0; i < Nx; ++i) {
#pragma HLS UNROLL
            y0[i] = params.u_init[i];
            y0[Nx + i] = params.v_init[i];
        }
    } else {
        for (int i = 0; i < length_hx; ++i) {
#pragma HLS UNROLL
            y0[i] = y_hx1[i];
        }
    }
    y0[Nx - 1] = Ts_HX1_L;

    // Solve the system of ODEs
    ode_solver_hx1(y0, pde_to_ode_hx1, step, params);

    // Update y_hx1 with the solution from the ODE solver
    for (int i = 0; i < length_hx; ++i) {
#pragma HLS UNROLL
        y_hx1[i] = y0[i];  // Update the original array with the new values
    }
}

