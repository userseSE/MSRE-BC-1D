#include "HX2.hpp"
#include "parameters.hpp"
#include "ode_solver.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// Include necessary linear algebra headers
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <boost/numeric/odeint.hpp>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace boost::numeric::odeint;



// Discretize the spatial domain
MatrixXd A_HX2(Nx, Nx);

void initialize_A_HX2(Parameters& params) {
    double dx = params.dx;
    A_HX2 = MatrixXd::Zero(Nx, Nx);
    for (int i = 0; i < Nx; ++i) {
        A_HX2(i, i) = -2.0 / (dx * dx);
        if (i < Nx - 1) {
            A_HX2(i, i + 1) = 1.0 / (dx * dx);
        }
        if (i > 0) {
            A_HX2(i, i - 1) = 1.0 / (dx * dx);
        }
    }

    A_HX2(0, 0) = -1.0 / (dx * dx);
    A_HX2(Nx - 1, Nx - 1) = -1.0 / (dx * dx);
}

VectorXd HX2(VectorXd& y_hx2, double Ts_HX2_L, int step, Parameters& params) {
    // Ensure A_HX2 is initialized
    static bool initialized = false;
    if (!initialized) {
        initialize_A_HX2(params);
        initialized = true;
    }

    // Set boundary conditions
    y_hx2[Nx + 1] = params.v2_L;    // Fixed boundary condition at v[0]
    y_hx2[2*Nx - 1] = params.v2_H; // Fixed boundary condition at v[L]

    VectorXd u(Nx), v(Nx);
    if (step == 0) {
        u = params.u2_init;  // Use Eigen::VectorXd directly
        v = params.v2_init;
    } else {
        u = y_hx2.head(Nx);
        v = y_hx2.tail(Nx);
    }
    u[Nx - 1] = Ts_HX2_L;
    // v[0] = Tss_HX2_0;

    // Define the ODE system for the second heat exchanger model compatible with odeint
    std::function<void(double, const VectorXd &, VectorXd &)> pde_to_ode_hx2 = 
    [&](double t, const VectorXd &y, VectorXd &dydt) {
        // Split the input state vector y into u and v
        VectorXd u = y.head(Nx);
        VectorXd v = y.tail(Nx);

        // Compute du_dt and dv_dt using the provided constants
        VectorXd du_dt = params.C1_2 * (A_HX2 * u) + params.C2_2 * (u - v);
        VectorXd dv_dt = params.C3_2 * (A_HX2 * v) + params.C4_2 * (u - v);

        // // Apply time-varying boundary conditions
        du_dt[0] = params.u2_L - u[0];
        du_dt[Nx - 1] = params.u2_H - u[Nx - 1];
        // dv_dt[0] = v2_L - v[0];
        // dv_dt[Nx - 1] = v2_H - v[Nx - 1];
        dv_dt[0] = 0; // Fixed condition: no change at v[0]
        dv_dt[Nx - 1] = 0; // Fixed condition: no change at v[Nx-1]

        // double k = 0.05;  // Heat transfer coefficient to ambient (example value)
        // double ambient_temp = 300.0;  // Ambient temperature in Kelvin
        // du_dt[0] = -k * (u[0] - ambient_temp);
        // du_dt[N - 1] = -k * (u[Nx - 1] - ambient_temp);
        // // dv_dt[0] = -k * (v[0] - ambient_temp);
        // // dv_dt[N - 1] = -k * (v[Nx - 1] - ambient_temp);
        // dv_dt[0] = 0; // Fixed condition: no change at v[0]
        // dv_dt[Nx - 1] = 0; // Fixed condition: no change at v[Nx-1]

        // Populate dydt with the derivatives
        dydt.head(Nx) = du_dt;
        dydt.tail(Nx) = dv_dt;
    };

    // Initial condition vector
    VectorXd y0(2 * Nx);
    if (step == 0) {
        y0.head(Nx) = params.u2_init;
        y0.tail(Nx) = params.v2_init;
    } else {
        y0 = y_hx2;
    }

    // Solve the ODE system
    VectorXd solution_y_hx2 = ode_solver(y0, pde_to_ode_hx2, step);

    // Return the solution at the last time step
    return solution_y_hx2;
}
