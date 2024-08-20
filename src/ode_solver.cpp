#include "ode_solver.hpp"
#include "parameters.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <functional>
#include <cmath>


// Function to implement the Radau method
// Note: Implementing Radau or any advanced ODE solver from scratch can be very complex.
// Here we'll assume a simple Euler method for demonstration. You would need to implement
// or use a library for Radau if required.
// using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                 const std::vector<double>& bc, 
//                                 ODEFunction vector_to_be_solved) {
//     // Initial condition vector
//     std::vector<double> y = ic;
//     std::vector<std::vector<double>> solution;
//     std::vector<double> results;

//     // Time vector
//     double t = 0.0;
//     double dt = dt;
//     double t_end = 1.0;  // Define the time until which you want to solve the ODE
//     double tolerance = 1e-1;  // Tolerance for the RK45 method

//     solution.push_back(y);  // Initial state

//     while (t < t_end) {
//         std::vector<double> k1 = vector_to_be_solved(t, y);
//         std::vector<double> y_temp(y.size());

//         for (size_t i = 0; i < y.size(); ++i) {
//             y_temp[i] = y[i] + 0.25 * dt * k1[i];
//         }
//         std::vector<double> k2 = vector_to_be_solved(t + 0.25 * dt, y_temp);

//         for (size_t i = 0; i < y.size(); ++i) {
//             y_temp[i] = y[i] + (3.0 / 32.0) * dt * k1[i] + (9.0 / 32.0) * dt * k2[i];
//         }
//         std::vector<double> k3 = vector_to_be_solved(t + (3.0 / 8.0) * dt, y_temp);

//         for (size_t i = 0; i < y.size(); ++i) {
//             y_temp[i] = y[i] + (1932.0 / 2197.0) * dt * k1[i] - (7200.0 / 2197.0) * dt * k2[i] + (7296.0 / 2197.0) * dt * k3[i];
//         }
//         std::vector<double> k4 = vector_to_be_solved(t + (12.0 / 13.0) * dt, y_temp);

//         for (size_t i = 0; i < y.size(); ++i) {
//             y_temp[i] = y[i] + (439.0 / 216.0) * dt * k1[i] - 8.0 * dt * k2[i] + (3680.0 / 513.0) * dt * k3[i] - (845.0 / 4104.0) * dt * k4[i];
//         }
//         std::vector<double> k5 = vector_to_be_solved(t + dt, y_temp);

//         for (size_t i = 0; i < y.size(); ++i) {
//             y_temp[i] = y[i] - (8.0 / 27.0) * dt * k1[i] + 2.0 * dt * k2[i] - (3544.0 / 2565.0) * dt * k3[i] + (1859.0 / 4104.0) * dt * k4[i] - (11.0 / 40.0) * dt * k5[i];
//         }
//         std::vector<double> k6 = vector_to_be_solved(t + 0.5 * dt, y_temp);

//         // Compute the next value of y using RK45
//         std::vector<double> y_next(y.size());
//         for (size_t i = 0; i < y.size(); ++i) {
//             y_next[i] = y[i] + dt * (16.0 / 135.0 * k1[i] + 6656.0 / 12825.0 * k3[i] + 28561.0 / 56430.0 * k4[i] - 9.0 / 50.0 * k5[i] + 2.0 / 55.0 * k6[i]);
//         }

//         // Estimate the error using the difference between the 4th and 5th order approximations
//         std::vector<double> y_err(y.size());
//         for (size_t i = 0; i < y.size(); ++i) {
//             y_err[i] = dt * (1.0 / 360.0 * k1[i] - 128.0 / 4275.0 * k3[i] - 2197.0 / 75240.0 * k4[i] + 1.0 / 50.0 * k5[i] + 2.0 / 55.0 * k6[i]);
//         }

//         // Calculate the maximum error and determine if the step size should be adjusted
//         double max_err = *std::max_element(y_err.begin(), y_err.end());
//         double safety_factor = 0.9;

//         if (max_err > tolerance) {
//             // Error is too large, reduce time step
//             dt *= safety_factor * std::pow(tolerance / max_err, 0.25);
//         } else {
//             // Accept the step and move forward
//             y = y_next;
//             t += dt;
//             solution.push_back(y);

//             if (max_err < tolerance / 10.0) {
//                 // Error is very small, increase time step
//                 dt *= safety_factor * std::pow(tolerance / max_err, 0.2);
//             }
//         }
//     }

//     std::cout << "Number of steps: " << solution.size() << std::endl;

//     results = solution.back();
//     return results;
// }
std::vector<double> ode_solver(const std::vector<double>& ic, 
                                            const std::vector<double>& bc, 
                                            ODEFunction vector_to_be_solved) {
    // Initial condition vector
    std::vector<double> y = ic;
    std::vector<std::vector<double>> solution;
    std::vector<double> results;
    // Time vector
    int time_steps = static_cast<int>(1.0 / dt);
    double t = 0.0;
    
    solution.push_back(y);  // Initial state

    // Simple Euler integration loop for demonstration
    for (int step = 0; step < time_steps; ++step) {
        std::vector<double> dy = vector_to_be_solved(t, y);
        
        for (size_t i = 0; i < y.size(); ++i) {
            y[i] += dt * dy[i];  // Euler step
        }
        
        t += dt;
        solution.push_back(y);
    }
    // std::cout<<solution.size()<<std::endl;

    results=solution.back();
    // return results;
    return ic;
}
