#include "HX1.hpp"
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

const double C1 = V_he_s;
const double C2 = -U_hx / (M_he_s * c_p_s);
const double C3 = V_he_ss;
const double C4 = U_hx / (M_he_ss * c_p_ss);

// Discretize the spatial domain
MatrixXd A_HX(Nx, Nx);

void initialize_A_HX() {
    A_HX = MatrixXd::Zero(Nx, Nx);
    for (int i = 0; i < Nx; ++i) {
        A_HX(i, i) = -2.0 / (dx * dx);
        if (i < Nx - 1) {
            A_HX(i, i + 1) = 1.0 / (dx * dx);
        }
        if (i > 0) {
            A_HX(i, i - 1) = 1.0 / (dx * dx);
        }
    }

    A_HX(0, 0) = -1.0 / (dx * dx);
    A_HX(Nx - 1, Nx - 1) = -1.0 / (dx * dx);
}

VectorXd HX1(VectorXd& y_hx1, double Ts_HX1_L, double Tss_HX1_0, int step) {
    // Ensure A_HX is initialized
    static bool initialized = false;
    if (!initialized) {
        initialize_A_HX();
        initialized = true;
    }

    // Set boundary conditions
    VectorXd u(Nx), v(Nx);
    if (step == 0) {
        u = u_init;  // Use Eigen::VectorXd directly
        v = v_init;
    } else {
        u = y_hx1.head(Nx);
        v = y_hx1.tail(Nx);
    }
    u[Nx - 1] = Ts_HX1_L;
    v[0] = Tss_HX1_0;

    // Define the ODE system for the heat exchanger compatible with odeint
    std::function<void(double, const VectorXd &, VectorXd &)> pde_to_ode_hx1 = 
    [&](double t, const VectorXd &y, VectorXd &dydt) {
        // Split the input state vector y into u and v
        VectorXd u = y.head(Nx);
        VectorXd v = y.tail(Nx);

        // Compute du_dt and dv_dt using the provided constants
        VectorXd du_dt = C1 * (A_HX * u) + C2 * (u - v);
        VectorXd dv_dt = C3 * (A_HX * v) + C4 * (u - v);

        // Apply time-varying boundary conditions
        du_dt[0] = u_L - u[0];
        du_dt[Nx - 1] = u_H - u[Nx - 1];
        dv_dt[0] = v_L - v[0];
        dv_dt[Nx - 1] = v_H - v[Nx - 1];

        // Populate dydt with the derivatives
        dydt.head(Nx) = du_dt;
        dydt.tail(Nx) = dv_dt;
    };

    // Initial condition vector
    VectorXd y0(2 * Nx);
    if (step == 0) {
        y0.head(Nx) = u_init;
        y0.tail(Nx) = v_init;
    } else {
        y0 = y_hx1;
    }

    // Solve the system of ODEs
    VectorXd solution_y_hx1 = ode_solver(y0, pde_to_ode_hx1);

    // Return the solution at the last time step
    return solution_y_hx1;
}
