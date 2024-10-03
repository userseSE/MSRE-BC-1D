#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <vector>
#include <functional>
#include <functional>
#include <boost/numeric/odeint.hpp>

// using ODEFunction = std::function<std::vector<double>(double t, const std::vector<double>& y)>;
using ODEFunction = std::function<void(const std::vector<double>& y, std::vector<double>& dydt, double t)>;

// Solver function declaration
std::vector<double> ode_solver(const std::vector<double> &y0,
                               ODEFunction vector_to_be_solved);

#endif // ODE_SOLVER_HPP
