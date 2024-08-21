#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <functional>

#include "radau.hpp"

// Define the type for the ODE function
using ODEFunction = std::function<void(double, const std::vector<double>&, std::vector<double>&)>;

std::vector<std::vector<double>> radau_solver(ODEFunction f, double t0, double tf, std::vector<double> y0, double dt) {
    int N = y0.size();
    int steps = static_cast<int>((tf - t0) / dt);
    std::vector<std::vector<double>> solution(steps + 1, std::vector<double>(N));

    // Radau coefficients (Radau IIA)
    double c2 = (4.0 - sqrt(6.0)) / 10.0;
    double c3 = (4.0 + sqrt(6.0)) / 10.0;
    Eigen::Matrix3d A;
    A << 5.0/12, -1.0/12, 0,
         3.0/4,  1.0/4,  0,
         1.0/4,  5.0/12, 1.0/12;
    Eigen::Vector3d b(5.0/12, 3.0/4, 1.0/4);

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd J(N, N); // Jacobian matrix

    std::vector<double> y = y0;
    solution[0] = y0;

    for (int i = 0; i < steps; ++i) {
        double t = t0 + i * dt;

        // Nonlinear system for Radau (requires solving a system of equations)
        Eigen::VectorXd k1(N), k2(N), k3(N);
        Eigen::VectorXd y_new(N), F(N);
        
        // Initial guess (Euler step)
        std::vector<double> y_temp(N);
        f(t, y, y_temp);
        Eigen::VectorXd y_dot = Eigen::Map<Eigen::VectorXd>(y_temp.data(), N);
        k1 = y_dot;

        // Newton-Raphson iteration for solving implicit equations
        for (int j = 0; j < 10; ++j) { // Max 10 iterations
            Eigen::VectorXd g1 = y_dot - (1.0 / dt) * (I - A(0, 0) * J) * k1 - A(0, 1) * k2 - A(0, 2) * k3;
            Eigen::VectorXd g2 = y_dot - (1.0 / dt) * (I - A(1, 0) * J) * k1 - A(1, 1) * k2 - A(1, 2) * k3;
            Eigen::VectorXd g3 = y_dot - (1.0 / dt) * (I - A(2, 0) * J) * k1 - A(2, 1) * k2 - A(2, 2) * k3;

            // Solve the linear system using LU decomposition
            Eigen::MatrixXd G(N, 3 * N);
            G << J * k1, J * k2, J * k3;

            Eigen::VectorXd g(N * 3);
            g << g1, g2, g3;

            Eigen::VectorXd dk = G.colPivHouseholderQr().solve(g);

            k1 -= dk.segment(0, N);
            k2 -= dk.segment(N, N);
            k3 -= dk.segment(2 * N, N);

            if (dk.norm() < 1e-8) break; // Convergence criteria
        }

        // Update solution
        // y_new = y + dt * (b(0) * k1 + b(1) * k2 + b(2) * k3);
        // y_new = y + dt * (b(0) * k1 + b(1) * k2 + b(2) * k3);
        // Eigen::VectorXd y_new_vec = y + dt * (b(0) * k1 + b(1) * k2 + b(2) * k3);
        // Ensure y, k1, k2, and k3 are Eigen::VectorXd
        // Eigen::VectorXd y_new_vec = y + dt * (b(0) * k1 + b(1) * k2 + b(2) * k3);
        Eigen::VectorXd y_new_vec = Eigen::Map<Eigen::VectorXd>(y.data(), y.size()) + dt * (b[0] * k1 + b[1] * k2 + b[2] * k3);
        
        std::vector<double> y_new_vector(y_new_vec.data(), y_new_vec.data() + y_new_vec.size());
        
        y = std::vector<double>(y_new.data(), y_new.data() + N);
        solution[i + 1] = y;
    }

    return solution;
}
