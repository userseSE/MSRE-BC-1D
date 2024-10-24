#include "parameters.hpp"
#ifndef REACTIVITY_HPP
#define REACTIVITY_HPP

void reactivity(const double temperature_fuel[N], const double temperature_graphite[N], 
                int step, int time_span, Parameters& params, double rho[N]);

#endif // REACTIVITY_HPP
