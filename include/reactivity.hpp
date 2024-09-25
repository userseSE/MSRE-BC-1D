#include "parameters.hpp"
#ifndef REACTIVITY_HPP
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define REACTIVITY_HPP

using Eigen::VectorXd;
VectorXd reactivity(const VectorXd& temperature_fuel, const VectorXd& temperature_graphite, 
                    int step, int time_span, double rho_insertion, Parameters& params);

#endif // REACTIVITY_HPP
