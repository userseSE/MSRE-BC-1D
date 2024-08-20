#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <vector>
#include <functional>
#include <functional>

// Function signature for the ODE system to be solved
using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

// Solver function declaration
std::vector<double> ode_solver(const std::vector<double>& ic, 
                                            const std::vector<double>& bc, 
                                            ODEFunction vector_to_be_solved);

#endif // ODE_SOLVER_HPP
