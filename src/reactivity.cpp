#include "reactivity.hpp"
#include "parameters.hpp"

void reactivity(const float temperature_fuel[N], const float temperature_graphite[N], 
                int step, float rho[N]) {
    Param_React params;
    float rho_feedback[N]; // Array for reactivity feedback
    float react[N];   // Array for reactivity insertion
    float max_rho_change = params.max_rho_change;

    // Determine the reactivity insertion based on the step
    float reactdata = (step < time_span / 2) ? 0.0 : (rho_insertion / N) * 1e-5;

    // Initialize reactivity insertion and reactivity feedback arrays
    for (int i = 0; i < N; ++i) { 
        react[i] = reactdata;
        rho_feedback[i] = 0.0;
    }

    // Calculate reactivity feedback based on temperature differences
    for (int i = 0; i < N; ++i) {
        // Calculate fuel and graphite temperature differences
        float fuel_diff = (temperature_fuel[i] - params.initialS[i]) * params.alpha_f;
        float graphite_diff = (temperature_graphite[i] - params.initialG[i]) * params.alpha_g;

        // Sum the fuel and graphite contributions to reactivity feedback
        rho_feedback[i] = fuel_diff + graphite_diff;
        // std::cout<<graphite_diff<<std::endl;
    }

    // Calculate total reactivity: rho_0 + rho_feedback + react
    for (int i = 0; i < N; ++i) {
        rho[i] = params.rho_0_value + rho_feedback[i] + react[i];
    }
}

