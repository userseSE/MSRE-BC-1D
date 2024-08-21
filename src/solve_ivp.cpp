#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

#include "solve_ivp.hpp"

// Define a type for the ODE function
using ODEFunction = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;
// using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

// An example Runge-Kutta solver (e.g., RK45)
void runge_kutta_45(ODEFunction f, double t0, double tf, std::vector<double>& y0, double dt) {
    double t = t0;
    std::vector<double> y = y0;
    size_t n = y.size();
    std::vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), y_temp(n);

    while (t < tf) {
        f(t, y, k1);
        for (size_t i = 0; i < n; ++i) {
            y_temp[i] = y[i] + dt * k1[i] / 4.0;
        }
        
        f(t + dt / 4.0, y_temp, k2);
        for (size_t i = 0; i < n; ++i) {
            y_temp[i] = y[i] + dt * (k1[i] / 32.0 + k2[i] * 9.0 / 32.0);
        }

        f(t + 3.0 * dt / 8.0, y_temp, k3);
        for (size_t i = 0; i < n; ++i) {
            y_temp[i] = y[i] + dt * (1932.0 * k1[i] / 2197.0 - 7200.0 * k2[i] / 2197.0 + 7296.0 * k3[i] / 2197.0);
        }

        f(t + 12.0 * dt / 13.0, y_temp, k4);
        for (size_t i = 0; i < n; ++i) {
            y_temp[i] = y[i] + dt * (439.0 * k1[i] / 216.0 - 8.0 * k2[i] + 3680.0 * k3[i] / 513.0 - 845.0 * k4[i] / 4104.0);
        }

        f(t + dt, y_temp, k5);
        for (size_t i = 0; i < n; ++i) {
            y_temp[i] = y[i] - 8.0 * k1[i] / 27.0 + 2.0 * k2[i] - 3544.0 * k3[i] / 2565.0 + 1859.0 * k4[i] / 4104.0 - 11.0 * k5[i] / 40.0;
        }

        f(t + dt / 2.0, y_temp, k6);
        for (size_t i = 0; i < n; ++i) {
            y[i] += dt * (16.0 * k1[i] / 135.0 + 6656.0 * k3[i] / 12825.0 + 28561.0 * k4[i] / 56430.0 - 9.0 * k5[i] / 50.0 + 2.0 * k6[i] / 55.0);
        }

        t += dt;
    }

    y0 = y;
}

// Function to simulate `solve_ivp`
std::vector<std::vector<double>> solve_ivp(ODEFunction f, std::pair<double, double> t_span, const std::vector<double>& y0, double dt) {
    double t0 = t_span.first;
    double tf = t_span.second;
    std::vector<double> y = y0;
    std::vector<std::vector<double>> solution;

    runge_kutta_45(f, t0, tf, y, dt);

    solution.push_back(y);
    return solution;
}

// Example usage
// void example_ode(double t, const std::vector<double>& y, std::vector<double>& dydt) {
//     dydt[0] = -0.5 * y[0];
// }

// int main() {
//     double t0 = 0.0;
//     double tf = 10.0;
//     double dt = 0.1;
//     std::vector<double> y0 = {2.0};
//     std::pair<double, double> t_span = {t0, tf};

//     auto solution = solve_ivp(example_ode, t_span, y0, dt);

//     for (const auto& state : solution) {
//         std::cout << state[0] << std::endl;
//     }

//     return 0;
// }
