#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>

#include "radau.hpp"

// Define the type for the ODE function
using ODEFunction = std::function<std::vector<double>(double t, const std::vector<double>& y)>;

// Function to compute the Jacobian matrix numerically
Eigen::MatrixXd computeJacobian(const ODEFunction& f, double t, const std::vector<double>& y) {
    double epsilon = 1e-8;
    int n = y.size();
    Eigen::MatrixXd J(n, n);
    std::vector<double> f0 = f(t, y);

    for (int j = 0; j < n; ++j) {
        std::vector<double> y_eps = y;
        y_eps[j] += epsilon;
        std::vector<double> f_eps = f(t, y_eps);
        for (int i = 0; i < n; ++i) {
            J(i, j) = (f_eps[i] - f0[i]) / epsilon;
        }
    }

    return J;
}

// Newton-Raphson method for solving nonlinear equations
std::vector<double> newtonSolve(const Eigen::MatrixXd& J, const std::vector<double>& F, const std::vector<double>& y_guess, double tol, int max_iter) {

    Eigen::VectorXd y = Eigen::Map<const Eigen::VectorXd>(y_guess.data(), y_guess.size());
    Eigen::VectorXd F_vec = Eigen::Map<const Eigen::VectorXd>(F.data(), F.size());
    Eigen::VectorXd delta_y = J.fullPivLu().solve(-F_vec);

    int iter = 0;
    while (delta_y.norm() > tol && iter < max_iter) {
        y += delta_y;
        delta_y = J.fullPivLu().solve(-F_vec);
        iter++;
    }

    return std::vector<double>(y.data(), y.data() + y.size());
}

// The Radau solver function
std::vector<double> radau_solver(ODEFunction f, double t0, double tf, std::vector<double> y0) {
    int n = y0.size();
    double dt = (tf - t0) / 5;  // Initial time step, adjust as needed
    std::vector<double> y = y0;
    double t = t0;

    double tol = 1e-5;
    int max_iter = 30;
    
    while (t < tf) {
        // Predict the next step (initial guess for Newton)
        std::vector<double> y_guess = y;
        std::vector<double> f_y = f(t, y);

        // Compute the Jacobian matrix
        Eigen::MatrixXd J = computeJacobian(f, t, y);

        // Solve the nonlinear system using Newton-Raphson
        std::vector<double> y_next = newtonSolve(J, f_y, y_guess, tol, max_iter);

        // Update the solution
        y = y_next;
        t += dt;
    }

    return y;
}

