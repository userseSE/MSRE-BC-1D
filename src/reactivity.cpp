#include "reactivity.hpp"
#include "parameters.hpp"
// #include "src/Core/Matrix.h"
#include <Eigen/Dense>
// #include <numeric>
// #include <algorithm>
// #include <cmath>
#include <iostream>

using Eigen::VectorXd;

VectorXd reactivity(const VectorXd& temperature_fuel, const VectorXd& temperature_graphite, 
                    int step, int time_span, double rho_insertion, Parameters& params) {
    
    // // Calculate the initial reactivity (rho_0_value)
    // double rho_0_value = std::accumulate(beta.begin(), beta.end(), 0.0);
    // for (int i = 0; i < 6; ++i) {
    //     rho_0_value -= beta[i] / (1.0 + (1.0 / (lambda_i[i] * tau_c)) * (1.0 - std::exp(-lambda_i[i] * tau_l)));
    // }
    VectorXd temperature_fuel_r = params.initialS;
    VectorXd temperature_graphite_r = params.initialG;
    double max_rho_change = params.max_rho_change;

    // if (step > time_span / 2 + 500) {
    //     params.rho_0_value = -(rho_insertion)*1e-5;
    // }
    // Initialize reactivity vectors
    VectorXd rho_0 = VectorXd::Constant(N, params.rho_0_value);
    // std::cout << rho_0_value << std::endl;

    // Determine the reactivity insertion based on the step
    double reactdata = (step < time_span / 2) ? 0.0 : (rho_insertion / N) * 1e-5;
    VectorXd react = VectorXd::Constant(N, reactdata);

    // Calculate reactivity feedback
    VectorXd rho_feedback = VectorXd::Zero(N);
    for (int i = 0; i < N; ++i) {
        double fuel_r = temperature_fuel_r[i];
        double fuel = temperature_fuel[i];
        double graphite_r = temperature_graphite_r[i];
        double graphite = temperature_graphite[i];

        // Calculate fuel and graphite differences, then sum
        double fuel_diff = (fuel - fuel_r) * params.alpha_f;
        double graphite_diff = (graphite - graphite_r) * params.alpha_g;
        double result = fuel_diff + graphite_diff;

        // if (step > time_span / 2 && step < time_span / 2 + 1000) {
        //     // rho_feedback[i] = result;
        //     rho_feedback[i] = std::clamp(result, -1e-3, 1e-3);
        // } else {
        //     // Clip (clamp) the result between -max_rho_change and max_rho_change
        //     rho_feedback[i] = std::clamp(result, -max_rho_change, max_rho_change);
        // }

        // Clip (clamp) the result between -max_rho_change and max_rho_change
        // rho_feedback[i] = std::clamp(result, -max_rho_change, max_rho_change);
        rho_feedback[i] = result;
    }

    // std::cout << "rho_feedback: " << rho_feedback[0] << std::endl;

    // Calculate total reactivity
    VectorXd rho = rho_0 + rho_feedback + react;
    // if (step == time_span/2 || step == time_span/2+1){
    //     std::cout << "rho_fb: " << step<<": "<<rho_feedback.sum()/N << std::endl;
    //     std::cout << "rho_0: " << step<<": "<<rho_0.sum()/N << std::endl;
    //     std::cout << "react: " << step<<": "<<react.sum()/N << std::endl;
    //     std::cout << "rho: " << step<<": "<<rho.sum()/N << std::endl;
    //     std::cout << "rho_fb[N/2]: " << step<<": "<<rho_feedback[N/2]<< std::endl;
    //     std::cout << "rho_0: " << step<<": "<<rho_0[N/2]<< std::endl;
    //     std::cout << "react: " << step<<": "<<react[N/2] << std::endl;
    //     std::cout << "rho: " << step<<": "<<rho[N/2] << std::endl;
    // }
    // VectorXd rho = VectorXd::Constant(N,0.0);
    // std::cout << "rho: " << rho[0] << std::endl;
    // rho = VectorXd::Constant(N, 0.0);
    return rho;
}
