#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer)(double, const double[length_hx], double[length_hx], Parameters &params);

// Update the signature to use VectorXd instead of state_type
void ode_solver_hx(double y[length_hx], OdeFuncPointer ode_func, int step, Parameters &params);

#endif // ODE_SOLVER_HPP
