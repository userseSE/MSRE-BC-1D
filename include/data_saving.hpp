#ifndef DATA_SAVING_HPP
#define DATA_SAVING_HPP

#include <Eigen/Core>
#include <string>
#include "parameters.hpp"

void save_results(const double rho_matrix[time_span], 
                  const double phi1_middle_matrix[time_span],
                  const double phi2_middle_matrix[time_span],
                  const double ci_middle_matrix[time_span][6], 
                  const double temperature_fuel_middle_matrix[time_span],
                  const std::string& folder);

void save_spacial_results(const double rho[N], 
                          const double y_n[length_neutr], 
                          const double y_th[length_th], 
                          const double y_hx1[length_hx], 
                          const double y_hx2[length_hx], 
                          const std::string& folder);

#endif // DATA_SAVING_HPP