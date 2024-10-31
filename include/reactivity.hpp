#include "parameters.hpp"
#ifndef REACTIVITY_HPP
#define REACTIVITY_HPP

void reactivity(const float temperature_fuel[N], const float temperature_graphite[N], 
                int step, Param_React& params, float rho[N]);

#endif // REACTIVITY_HPP
