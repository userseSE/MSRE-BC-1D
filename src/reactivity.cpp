#include "reactivity.hpp"
#include "parameters.hpp"
#include <vector>
#include <numeric>

double reactivity(const std::vector<double>& temperature_fuel, 
                  const std::vector<double>& temperature_graphite, 
                  int step, int time_span, double rho_insertion) {
    
    double rho_0 = std::accumulate(beta.begin(), beta.end(), 0.0);
    for (int i = 0; i < 6; ++i) {
        rho_0 -= beta[i] / (1.0 + (1.0 / (lambda_i[i] * tau_c)) * (1.0 - std::exp(-lambda_i[i] * tau_l)));
    }
    
    // Determine the reactivity insertion based on the step
    std::vector<double> reactdata = {50, rho_insertion * 1e-5};
    double react;
    if (step < time_span / 2) {
        react = reactdata[0];
    } else {
        react = reactdata[1];
    }

    // Calculate reactivity feedback
    double rho_feedback_fuel = std::inner_product(initialS, initialS + N, temperature_fuel.begin(), 0.0, std::plus<double>(), 
                                                  [](double init, double temp) { return (init - temp) * alpha_f; });
    double rho_feedback_graphite = std::inner_product(initialG, initialG + N, temperature_graphite.begin(), 0.0, std::plus<double>(), 
                                                      [](double init, double temp) { return (init - temp) * alpha_g; });
    double rho_feedback = (rho_feedback_fuel + rho_feedback_graphite) / N;

    // Calculate total reactivity
    double rho = rho_0 + rho_feedback + react;

    return rho;
}
