#ifndef REACTIVITY_HPP
#include <vector>
#define REACTIVITY_HPP

std::vector<double> reactivity(const std::vector<double>& temperature_fuel_r, const std::vector<double>& temperature_hraphite_r, const std::vector<double>& temperature_fuel, 
                  const std::vector<double>& temperature_graphite, 
                  int step, int time_span, double rho_insertion);

#endif // REACTIVITY_HPP
