#ifndef ODE_SOLVER_HPP
#define ODE_SOLVER_HPP

#include <Eigen/Dense>
#include <functional>

using Eigen::VectorXd;

// Update the signature to use VectorXd instead of state_type
VectorXd ode_solver(const VectorXd &y0, const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func, int step);

#endif // ODE_SOLVER_HPP
