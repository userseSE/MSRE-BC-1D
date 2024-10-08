#include "reactivity.hpp"
#include "parameters.hpp"
#include "src/Core/Matrix.h"
#include <Eigen/Dense>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>

using Eigen::VectorXd;

VectorXd reactivity(const VectorXd& temperature_fuel_r, const VectorXd& temperature_graphite_r, 
                    const VectorXd& temperature_fuel, const VectorXd& temperature_graphite, 
                    int step, int time_span, double rho_insertion) {
    
    // Calculate the initial reactivity (rho_0_value)
    double rho_0_value = std::accumulate(beta.begin(), beta.end(), 0.0);
    for (int i = 0; i < 6; ++i) {
        rho_0_value -= beta[i] / (1.0 + (1.0 / (lambda_i[i] * tau_c)) * (1.0 - std::exp(-lambda_i[i] * tau_l)));
    }

    // Initialize reactivity vectors
    VectorXd rho_0 = VectorXd::Constant(N, rho_0_value);

    // Determine the reactivity insertion based on the step
    double reactdata = (step < time_span / 2) ? 0.0 : rho_insertion * 1e-5;
    VectorXd react = VectorXd::Constant(N, reactdata);

    // Calculate reactivity feedback
    VectorXd rho_feedback = VectorXd::Zero(N);
    for (int i = 0; i < N; ++i) {
        double fuel_r = temperature_fuel_r[i];
        double fuel = temperature_fuel[i];
        double graphite_r = temperature_graphite_r[i];
        double graphite = temperature_graphite[i];

        // Calculate fuel and graphite differences, then sum
        double fuel_diff = (fuel_r - fuel) * alpha_f;
        double graphite_diff = (graphite_r - graphite) * alpha_g;
        double result = fuel_diff + graphite_diff;

        // Clip (clamp) the result between -max_rho_change and max_rho_change
        rho_feedback[i] = std::clamp(result, -max_rho_change, max_rho_change);
    }

    std::cout << "rho_feedback: " << rho_feedback[0] << std::endl;

    // Calculate total reactivity
    VectorXd rho = rho_0 + rho_feedback;
    // VectorXd rho = VectorXd::Constant(N,0.0);
    // std::cout << "rho: " << rho[0] << std::endl;

    return rho;
}
