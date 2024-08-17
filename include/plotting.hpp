#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <vector>

void plot_results(const std::vector<double>& rho_matrix, 
                  const std::vector<double>& phi_middle_matrix,
                  const std::vector<std::vector<double>>& ci_middle_matrix, 
                  const std::vector<double>& temperature_fuel_middle_matrix);

#endif // PLOTTING_HPP
