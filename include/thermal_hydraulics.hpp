#ifndef THERMAL_HYDRAULICS_HPP
#define THERMAL_HYDRAULICS_HPP

#include "parameters.hpp"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::VectorXd;

// Function to solve the thermal hydraulics problem
VectorXd thermal_hydraulics(VectorXd &y_th, const VectorXd &q_prime, double Ts_core_0, int step, Parameters& params);

#endif // THERMAL_HYDRAULICS_HPP
