#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <vector>

void save_results(const std::vector<double>& rho_matrix, 
                  const std::vector<double>& phi_middle_matrix,
                  const std::vector<std::vector<double> >& ci_middle_matrix, 
                  const std::vector<double>& temperature_fuel_middle_matrix);
void save_spacial_results(const std::vector<double>& phi, 
                          const std::vector<double>& ci, 
                          const std::vector<double>& temperature_fuel, 
                          const std::vector<double>& temperature_graphite, 
                          const std::vector<double>& Ts_HX1, 
                          const std::vector<double>& Tss_HX1, 
                          const std::vector<double>& Tss_HX2, 
                          const std::vector<double>& Tsss_HX2);
#endif // PLOTTING_HPP
