#include <iostream>
#include "HX1.hpp"
#include "parameters.hpp"
#include "ode_solver_HX.hpp"

void pde_to_ode_hx1(double t, const double y[length_hx], double dydt[length_hx], Parameters &params) {
    // Split the input state vector y into u and v
    const double* u = &y[0];    // u is the first half of y
    const double* v = &y[Nx];   // v is the second half of y

    double du_dt[Nx] = {0};    // Array to store the derivative of u
    double dv_dt[Nx] = {0};    // Array to store the derivative of v

    // Compute du_dt and dv_dt using the provided constants
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nx; ++j) {
            du_dt[i] += params.C1_1 * params.A_HX[i][j] * u[j];
            dv_dt[i] += params.C3_1 * params.A_HX[i][j] * v[j];
        }
        du_dt[i] += params.C2_1 * (u[i] - v[i]);
        dv_dt[i] += params.C4_1 * (u[i] - v[i]);
    }

    // Apply time-varying boundary conditions
    du_dt[0] = params.u_L - u[0];
    du_dt[Nx - 1] = params.u_H - u[Nx - 1];
    dv_dt[0] = params.v_L - v[0];
    dv_dt[Nx - 1] = params.v_H - v[Nx - 1];

    // Populate dydt with the derivatives
    for (int i = 0; i < Nx; ++i) {
        dydt[i] = du_dt[i];        // First half for du_dt
        dydt[Nx + i] = dv_dt[i];   // Second half for dv_dt
    }
}

void HX1(double y_hx1[length_hx], double Ts_HX1_L, double Tss_HX1_0, int step, Parameters &params) {
    std::cout<<"HX1 is called"<<std::endl;
    // Set boundary conditions
    double u[Nx], v[Nx];
    if (step == 0) {
        for (int i = 0; i < Nx; ++i) {
            u[i] = params.u_init[i];
            v[i] = params.v_init[i];
        }
    } else {
        for (int i = 0; i < Nx; ++i) {
            u[i] = y_hx1[i];
            v[i] = y_hx1[Nx + i];
        }
    }
    u[Nx - 1] = Ts_HX1_L;
    v[0] = Tss_HX1_0;

    // Initial condition vector y0
    double y0[length_hx];
    if (step == 0) {
        for (int i = 0; i < Nx; ++i) {
            y0[i] = params.u_init[i];
            y0[Nx + i] = params.v_init[i];
        }
    } else {
        for (int i = 0; i < length_hx; ++i) {
            y0[i] = y_hx1[i];
        }
    }

    // Solve the system of ODEs
    ode_solver_HX(y0, pde_to_ode_hx1, step, params);

    // Update y_hx1 with the solution from the ODE solver
    for (int i = 0; i < length_hx; ++i) {
        y_hx1[i] = y0[i];  // Update the original array with the new values
    }
}

