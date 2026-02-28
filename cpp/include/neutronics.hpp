#ifndef NEUTRONICS_HPP
#define NEUTRONICS_HPP

#include "parameters.hpp"

// Function to solve the neutronics problem
void neutronics(float y_n[length_neutr], const float rho[N], int step, Param_Neutronics& params, float min_step_neutr);

#endif // NEUTRONICS_HPP
