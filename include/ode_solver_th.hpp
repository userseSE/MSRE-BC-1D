#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer_th)(float, const float[length_neutr], float[length_neutr], Param_Thermal &params, const float q_prime[N]);

// Update the signature to use VectorXd instead of state_type
void ode_solver_th(float y[length_th], OdeFuncPointer_th ode_func, int step, Param_Thermal &params, const float q_prime[N]);

#endif // ODE_SOLVER_HPP
