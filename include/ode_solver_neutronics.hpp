#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer)(float, const float[length_neutr], float[length_neutr], Param_Neutronics &params, const float Keff[N]);

// Update the signature to use VectorXd instead of state_type
void ode_solver_neutr(float y[length_neutr], OdeFuncPointer ode_func, int step, Param_Neutronics &params, const float Keff[N]);

#endif // ODE_SOLVER_HPP
