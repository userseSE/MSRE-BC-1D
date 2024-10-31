#ifndef DATA_SAVING_HPP
#define DATA_SAVING_HPP

#include <Eigen/Core>
#include <string>
#include "parameters.hpp"

void save_results(const float rho_matrix[time_span], 
                  const float phi_middle_matrix[time_span * 2],
                  const float ci_middle_matrix[time_span * 6], 
                  const float temperature_fuel_middle_matrix[time_span],
                  const std::string& folder);

void save_spacial_results(const float rho[N], 
                          const float y_n[length_neutr], 
                          const float y_th[length_th], 
                          const float y_hx1[length_hx], 
                          const float y_hx2[length_hx], 
                          const std::string& folder);

#endif // DATA_SAVING_HPP