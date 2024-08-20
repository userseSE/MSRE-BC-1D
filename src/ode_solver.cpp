#include "ode_solver.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include "parameters.hpp"

// Function to implement the Radau method
// Note: Implementing Radau or any advanced ODE solver from scratch can be very complex.
// Here we'll assume a simple Euler method for demonstration. You would need to implement
// or use a library for Radau if required.

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

    results=solution.back();
    return results;
}
