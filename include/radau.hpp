#ifndef RADAU_HPP
# include <vector>
# include <functional>
# include <Eigen/Dense>
# define RADAU_HPP

using ODEFunction = std::function<std::vector<double>(double t, const std::vector<double>& y)>;
Eigen::MatrixXd computeJacobian(const ODEFunction& f, double t, const std::vector<double>& y);
std::vector<double> newtonSolve(const Eigen::MatrixXd& J, const std::vector<double>& F, const std::vector<double>& y_guess, double tol, int max_iter);
std::vector<double> radau_solver(ODEFunction f, double t0, double tf, std::vector<double> y0);
#endif


