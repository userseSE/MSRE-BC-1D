#ifndef THERMAL_HYDRAULICS_HPP
#define THERMAL_HYDRAULICS_HPP

#include <vector>

// Function to solve the thermal hydraulics problem
std::vector<double> thermal_hydraulics(std::vector<double>& y_th, 
                                       const std::vector<double>& q_prime, 
                                       double Ts_core_0, 
                                       int step);

#endif // THERMAL_HYDRAULICS_HPP
