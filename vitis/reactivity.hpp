#ifndef REACTIVITY_HPP
#define REACTIVITY_HPP

#include "parameters.hpp"

void reactivity(const float temperature_fuel[N], const float temperature_graphite[N],
                int step, float rho[N], float rho_insertion);

#endif // REACTIVITY_HPP
