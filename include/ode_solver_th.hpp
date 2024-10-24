#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer)(double, const double[length_neutr], double[length_neutr], Parameters &params, const double q_prime[N]);

// Update the signature to use VectorXd instead of state_type
void ode_solver_th(double y[length_th], OdeFuncPointer ode_func, int step, Parameters &params, const double q_prime[N]);

#endif // ODE_SOLVER_HPP
