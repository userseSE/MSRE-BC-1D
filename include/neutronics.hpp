#ifndef NEUTRONICS_HPP
#define NEUTRONICS_HPP

#include "parameters.hpp"

// Function to solve the neutronics problem
void neutronics(double y_n[length_neutr], const double rho[N], int step, Parameters &params);

#endif // NEUTRONICS_HPP
