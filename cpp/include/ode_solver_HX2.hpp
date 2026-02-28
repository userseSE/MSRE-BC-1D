#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include "parameters.hpp"

typedef void (*OdeFuncPointer_hx2)(float, const float[length_hx], float[length_hx], Param_HX2 &params);

// Update the signature to use VectorXd instead of state_type
void ode_solver_hx2(float y[length_hx], OdeFuncPointer_hx2 ode_func, int step, Param_HX2 &params);

#endif // ODE_SOLVER_HPP
