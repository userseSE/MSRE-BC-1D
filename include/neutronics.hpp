#ifndef NEUTRONICS_HPP
#define NEUTRONICS_HPP

#include "parameters.hpp"
#include <vector>
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::VectorXd state_type;

// Function to solve the neutronics problem
Eigen::VectorXd neutronics(const state_type &y_n, const std::vector<double> &rho, int step, Parameters& params);

#endif // NEUTRONICS_HPP
