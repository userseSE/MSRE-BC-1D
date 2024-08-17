#ifndef REACTIVITY_HPP
#include <vector>
#define REACTIVITY_HPP

double reactivity(const std::vector<double>& temperature_fuel, 
                  const std::vector<double>& temperature_graphite, 
                  int step, int time_span, double rho_insertion);

#endif // REACTIVITY_HPP
