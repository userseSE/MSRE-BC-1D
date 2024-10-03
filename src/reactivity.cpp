#include "reactivity.hpp"
#include "parameters.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iostream>

std::vector<double> reactivity(const std::vector<double>& temperature_fuel_r, const std::vector<double>& temperature_graphite_r, 
                const std::vector<double>& temperature_fuel, const std::vector<double>& temperature_graphite, 
                  int step, int time_span, double rho_insertion) {
    
    // rho_0, rho_fb, rho_insertion 3 vectors
    // std::cout<<"test rho"<<std::endl;
    double rho_0_value = std::accumulate(beta.begin(), beta.end(), 0.0);
    for (int i = 0; i < 6; ++i) {
        rho_0_value -= beta[i] / (1.0 + (1.0 / (lambda_i[i] * tau_c)) * (1.0 - std::exp(-lambda_i[i] * tau_l)));
    }

    // std::vector<double> rho_0(N, rho_0_value);
    std::vector<double> rho_0 (N, 0.0);
    // Determine the reactivity insertion based on the step
    std::vector<double> reactdata = {0, rho_insertion * 1e-5};
    std::vector<double> react;
    if (step < time_span / 2) {
        react.assign(N, reactdata[0]);
    } else {
        react.assign(N, reactdata[1]);
    }

    // Calculate reactivity feedback
    std::vector<double> rho_feedback(N, 0.0);

    for (size_t i = 0; i < N; ++i) {
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

    std::cout<<"rho_feedback: "<< rho_feedback.front()<<std::endl;

    // Calculate total reactivity
    std::vector<double> rho(N, 0.0);
    std::transform(rho_0.begin(), rho_0.end(), rho_feedback.begin(), rho.begin(), std::plus<double>());
    std::cout<<"rho: "<< rho.front()<<std::endl;
    return rho;
}
