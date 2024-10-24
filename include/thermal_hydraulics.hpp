#ifndef THERMAL_HYDRAULICS_HPP
#define THERMAL_HYDRAULICS_HPP

#include "parameters.hpp"

// Function to solve the thermal hydraulics problem
void thermal_hydraulics(double y_th[length_th], const double q_prime[N], double Ts_core_0, int step, Parameters &params);

#endif // THERMAL_HYDRAULICS_HPP
