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

const double C1 = V_he2_s;
const double C2 = -U2_hx / (M_he2_s * c_p_ss);
const double C3 = V_he2_ss;
const double C4 = U2_hx / (M_he2_ss * c_p_sss);

// Discretize the spatial domain
MatrixXd A_HX2(Nx, Nx);

void initialize_A_HX2() {
    A_HX2 = MatrixXd::Zero(Nx, Nx);
    for (int i = 0; i < Nx; ++i) {
        A_HX2(i, i) = -1.0;
        if (i < Nx - 1) {
            A_HX2(i, i + 1) = 1.0 / dx;
        }
    }
    A_HX2(0, 0) = 1.0 / dx;
    A_HX2(Nx - 1, Nx - 1) = 1.0 / dx;
}

VectorXd HX2(VectorXd& y_hx2, double Ts_HX2_L, double Tss_HX2_0, int step) {
    // Ensure A_HX2 is initialized
    static bool initialized = false;
    if (!initialized) {
        initialize_A_HX2();
        initialized = true;
    }

    // Set boundary conditions
    VectorXd u(Nx), v(Nx);
    if (step == 0) {
        u = u2_init;  // Use Eigen::VectorXd directly
        v = v2_init;
    } else {
        u = y_hx2.head(Nx);
        v = y_hx2.tail(Nx);
    }
    u[Nx - 1] = Ts_HX2_L;
    v[0] = Tss_HX2_0;

    // Define the ODE system for the second heat exchanger model compatible with odeint
    std::function<void(double, const VectorXd &, VectorXd &)> pde_to_ode_hx2 = 
    [&](double t, const VectorXd &y, VectorXd &dydt) {
        // Split the input state vector y into u and v
        VectorXd u = y.head(Nx);
        VectorXd v = y.tail(Nx);

        // Compute du_dt and dv_dt using the provided constants
        VectorXd du_dt = C1 * (A_HX2 * u) + C2 * (u - v);
        VectorXd dv_dt = C3 * (A_HX2 * v) + C4 * (u - v);

        // Apply time-varying boundary conditions
        du_dt[0] = u2_L - u[0];
        du_dt[Nx - 1] = u2_H - u[Nx - 1];
        dv_dt[0] = v2_L - v[0];
        dv_dt[Nx - 1] = v2_H - v[Nx - 1];

        // Populate dydt with the derivatives
        dydt.head(Nx) = du_dt;
        dydt.tail(Nx) = dv_dt;
    };

    // Initial condition vector
    VectorXd y0(2 * Nx);
    if (step == 0) {
        y0.head(Nx) = u2_init;
        y0.tail(Nx) = v2_init;
    } else {
        y0 = y_hx2;
    }

    // Solve the ODE system
    VectorXd solution_y_hx2 = ode_solver(y0, pde_to_ode_hx2);

    // Return the solution at the last time step
    return solution_y_hx2;
}
