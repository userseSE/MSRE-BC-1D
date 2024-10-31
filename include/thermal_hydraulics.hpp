#ifndef THERMAL_HYDRAULICS_HPP
#define THERMAL_HYDRAULICS_HPP

#include "parameters.hpp"

// Function to solve the thermal hydraulics problem
void thermal_hydraulics(float y_th[length_th], const float q_prime[N], float Ts_core_0, int step, Param_Thermal &params);

#endif // THERMAL_HYDRAULICS_HPP
