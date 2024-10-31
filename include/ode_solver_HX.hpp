#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer_hx)(float, const float[length_hx], float[length_hx], Parameters &params);

// Update the signature to use VectorXd instead of state_type
void ode_solver_hx(float y[length_hx], OdeFuncPointer_hx ode_func, int step, Parameters &params);

#endif // ODE_SOLVER_HPP
