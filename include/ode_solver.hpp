#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <vector>
#include <functional>
#include <functional>

// Function signature for the ODE system to be solved
// using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;
// using ODEFunction = std::function<void(double, const std::vector<double>&, std::vector<double>&)>;

// using ODEFunction = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;
// using ODEFunction = std::function<void(double, const std::vector<double>&, std::vector<double>&)>;

// void ode_function_wrapper(double t, double y[], double yp[], ODEFunction& func);
using ODEFunction = std::function<std::vector<double>(double t, const std::vector<double>& y)>;

// Solver function declaration
std::vector<double> ode_solver(const std::vector<double>& y0, 
                                            const std::vector<double>& bc, 
                                            ODEFunction vector_to_be_solved, double dt);

#endif // ODE_SOLVER_HPP
