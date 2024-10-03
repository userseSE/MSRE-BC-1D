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
// using Eigen::SparseMatrix;

const double C1 = V_he_s;
const double C2 = -U_hx / (M_he_s * c_p_s);
const double C3 = V_he_ss;
const double C4 = U_hx / (M_he_ss * c_p_ss);

// Discretize the spatial domain
MatrixXd A_HX(Nx, Nx);

void initialize_A_HX() {
    A_HX = MatrixXd::Zero(Nx, Nx);
    for (int i = 0; i < Nx; ++i) {
        A_HX(i, i) = -2.0 / (dx*dx);
        if (i < Nx - 1) {
            A_HX(i, i + 1) = 1.0 / (dx*dx);
        }
        if (i > 0) {
            A_HX(i, i - 1) = 1.0 / (dx*dx);
        }
    }

    A_HX(0, 0) = -1.0 / (dx*dx);
    A_HX(Nx - 1, Nx - 1) = -1.0 / (dx*dx);

}

std::vector<std::vector<double>> HX1(std::vector<std::vector<double>> &y_hx1, 
                                     double Ts_HX1_L, double Tss_HX1_0, int step) {
    // Ensure A_HX is initialized
    static bool initialized = false;
    if (!initialized) {
        initialize_A_HX();
        initialized = true;
    }

    // Set boundary conditions
    std::vector<double> u(Nx), v(Nx);
    if (step == 0) {
        std::copy(hx1_initial_conditions[0].begin(), hx1_initial_conditions[0].end(), u.begin());
        std::copy(hx1_initial_conditions[1].begin(), hx1_initial_conditions[1].end(), v.begin());
    } else {
        std::copy(y_hx1[0].begin(), y_hx1[0].end(), u.begin());
        std::copy(y_hx1[1].begin(), y_hx1[1].end(), v.begin());
    }
    u.back() = Ts_HX1_L;
    v.front() = Tss_HX1_0;

    // Define the ODE system for the heat exchanger compatible with odeint
    auto pde_to_ode_hx1 = [&](const std::vector<std::vector<double>> &y, 
                              std::vector<std::vector<double>> &dydt, double t) {
        // Split the input state vector y into u and v
        Eigen::Map<const VectorXd> u(y[0].data(), Nx);
        Eigen::Map<const VectorXd> v(y[1].data(), Nx);

        // Compute du_dt and dv_dt using the provided constants
        VectorXd du_dt = C1 * (A_HX * u) + C2 * (u - v);
        VectorXd dv_dt = C3 * (A_HX * v) + C4 * (u - v);

        // Apply time-varying boundary conditions
        du_dt[0] = u_L - u[0];
        du_dt[Nx - 1] = u_H - u[Nx - 1];
        dv_dt[0] = v_L - v[0];
        dv_dt[Nx - 1] = v_H - v[Nx - 1];

        // Copy the derivatives to the dydt vector
        std::copy(du_dt.data(), du_dt.data() + Nx, dydt[0].begin());
        std::copy(dv_dt.data(), dv_dt.data() + Nx, dydt[1].begin());
    };

    // Initial condition vector (2D)
    std::vector<std::vector<double>> y0(2, std::vector<double>(Nx));
    if (step == 0) {
        std::copy(hx1_initial_conditions[0].begin(), hx1_initial_conditions[0].end(), y0[0].begin());
        std::copy(hx1_initial_conditions[1].begin(), hx1_initial_conditions[1].end(), y0[1].begin());
    } else {
        y0 = y_hx1;
    }

    // Solve the ODE system using 2D vectors
    std::vector<std::vector<double>> solution_y_hx1 = ode_solver(y0, pde_to_ode_hx1);

    // Return the solution at the last time step
    return solution_y_hx1;
}
