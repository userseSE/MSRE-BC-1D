#include "HX2.hpp"
#include "parameters.hpp"
#include "ode_solver_HX.hpp"

void pde_to_ode_hx2(double t, const double y[length_hx], double dydt[length_hx], Parameters &params) {
    // Split the input state vector y into u and v
    const double* u = &y[0];    // u is the first half of y
    const double* v = &y[Nx];   // v is the second half of y

    double du_dt[Nx] = {0};    // Array to store the derivative of u
    double dv_dt[Nx] = {0};    // Array to store the derivative of v

    // Compute du_dt and dv_dt using the provided constants
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nx; ++j) {
            du_dt[i] += params.C1_2 * params.A_HX2[i][j] * u[j];
            dv_dt[i] += params.C3_2 * params.A_HX2[i][j] * v[j];
        }
        du_dt[i] += params.C2_2 * (u[i] - v[i]);
        dv_dt[i] += params.C4_2 * (u[i] - v[i]);
    }

    // Apply time-varying boundary conditions
    du_dt[0] = params.u2_L - u[0];
    du_dt[Nx - 1] = params.u2_H - u[Nx - 1];
    dv_dt[0] = 0; // Fixed condition: no change at v[0]
    dv_dt[Nx - 1] = 0; // Fixed condition: no change at v[Nx-1]

    // Populate dydt with the derivatives
    for (int i = 0; i < Nx; ++i) {
        dydt[i] = du_dt[i];        // First half for du_dt
        dydt[Nx + i] = dv_dt[i];   // Second half for dv_dt
    }
}

void HX2(double y_hx2[length_hx], double Ts_HX2_L, int step, Parameters &params) {
    // Set boundary conditions
    y_hx2[Nx + 1] = params.v2_L;    // Fixed boundary condition at v[0]
    y_hx2[length_hx - 1] = params.v2_H; // Fixed boundary condition at v[L]

    double u[Nx], v[Nx];
    if (step == 0) {
        for (int i = 0; i < Nx; ++i) {
            u[i] = params.u2_init[i];
            v[i] = params.v2_init[i];
        }
    } else {
        for (int i = 0; i < Nx; ++i) {
            u[i] = y_hx2[i];
            v[i] = y_hx2[Nx + i];
        }
    }
    u[Nx - 1] = Ts_HX2_L;

    // Initial condition vector y0
    double y0[length_hx];
    if (step == 0) {
        for (int i = 0; i < Nx; ++i) {
            y0[i] = params.u2_init[i];
            y0[Nx + i] = params.v2_init[i];
        }
    } else {
        for (int i = 0; i < length_hx; ++i) {
            y0[i] = y_hx2[i];
        }
    }

    // Solve the system of ODEs
    ode_solver_hx(y0, pde_to_ode_hx2, step, params);

    // Update y_hx2 with the solution from the ODE solver
    for (int i = 0; i < length_hx; ++i) {
        y_hx2[i] = y0[i];  // Update the original array with the new values
    }
}
