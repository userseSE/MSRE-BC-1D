#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer_hx1)(float, const float[length_hx], float[length_hx], Param_HX1 &params);

// Update the signature to use VectorXd instead of state_type
void ode_solver_hx1(float y[length_hx], OdeFuncPointer_hx1 ode_func, int step, Param_HX1 &params, float min_step);

#endif // ODE_SOLVER_HPP
