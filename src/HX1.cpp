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

using Eigen::VectorXd;
using Eigen::MatrixXd;
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
        A_HX(i, i) = -1.0;
        if (i < Nx - 1) {
            A_HX(i, i + 1) = 1.0 / dx;
        }
    }
    A_HX(0, 0) = 1.0 / dx;
    A_HX(Nx - 1, Nx - 1) = 1.0 / dx;
}

std::vector<double> HX1(std::vector<double>& y_hx1, double Ts_HX1_L, double Tss_HX1_0, int step) {
    // Ensure A_HX is initialized
    static bool initialized = false;
    if (!initialized) {
        initialize_A_HX();
        initialized = true;
    }

    // Set boundary conditions
    std::vector<double> u(Nx), v(Nx);
    if (step == 0) {
        std::copy(u_init, u_init + Nx, u.begin());
        std::copy(v_init, v_init + Nx, v.begin());
    } else {
        std::copy(y_hx1.begin(), y_hx1.begin() + Nx, u.begin());
        std::copy(y_hx1.begin() + Nx, y_hx1.end(), v.begin());
    }
    u.back() = Ts_HX1_L;
    v.front() = Tss_HX1_0;

    // Define the ODE system for the heat exchanger
    auto pde_to_ode_hx1 = [&](double t, const std::vector<double>& y) {
        VectorXd u = Eigen::Map<const VectorXd>(y.data(), Nx);
        VectorXd v = Eigen::Map<const VectorXd>(y.data() + Nx, Nx);

        VectorXd du_dt = C1 * (A_HX * u) + C2 * (u - v);
        VectorXd dv_dt = C3 * (A_HX * v) + C4 * (u - v);

        // Apply time-varying boundary conditions
        du_dt[0] = u_L - u[0];
        du_dt[Nx - 1] = u_H - u[Nx - 1];
        dv_dt[0] = v_L - v[0];
        dv_dt[Nx - 1] = v_H - v[Nx - 1];

        std::vector<double> dydt(2 * Nx);
        std::copy(du_dt.data(), du_dt.data() + Nx, dydt.begin());
        std::copy(dv_dt.data(), dv_dt.data() + Nx, dydt.begin() + Nx);

        return dydt;
    };

    // auto pde_to_ode_hx1 = [&](double t, const std::vector<double>& y, std::vector<double>& dydt) {
    // // Create Eigen vectors from the input state vector
    //     Eigen::Map<const Eigen::VectorXd> u(y.data(), Nx);
    //     Eigen::Map<const Eigen::VectorXd> v(y.data() + Nx, Nx);

    //     // Calculate du/dt and dv/dt
    //     Eigen::VectorXd du_dt = C1 * (A_HX * u) + C2 * (u - v);
    //     Eigen::VectorXd dv_dt = C3 * (A_HX * v) + C4 * (u - v);

    //     // Apply time-varying boundary conditions
    //     du_dt[0] = u_L - u[0];
    //     du_dt[Nx - 1] = u_H - u[Nx - 1];
    //     dv_dt[0] = v_L - v[0];
    //     dv_dt[Nx - 1] = v_H - v[Nx - 1];

    //     // Copy the results into the output derivative vector
    //     std::copy(du_dt.data(), du_dt.data() + Nx, dydt.begin());
    //     std::copy(dv_dt.data(), dv_dt.data() + Nx, dydt.begin() + Nx);
    // };

    // Initial condition vector
    std::vector<double> y0(2 * Nx);
    if (step == 0) {
        std::copy(u_init, u_init + Nx, y0.begin());
        std::copy(v_init, v_init + Nx, y0.begin() + Nx);
    } else {
        y0 = y_hx1;
    }

    // Solve the ODE system
    // std::vector<double> solution_y_hx1 = ode_solver(y0, {}, pde_to_ode_hx1, dt);

    // std::vector<double> y0_hx1 = { /* initial conditions */ };
    // double dt_hx1 = 0.01;
    std::vector<double> solution_y_hx1 = ode_solver(y0, {}, pde_to_ode_hx1, dt);

    // Return the solution at the last time step
    return solution_y_hx1;
}
