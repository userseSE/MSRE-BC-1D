#ifndef HX1_HPP
#define HX1_HPP

#include "parameters.hpp"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// Function to solve the heat exchanger 1 problem
VectorXd HX1(VectorXd& y_hx1, double Ts_HX1_L, double Tss_HX1_0, int step, Parameters& params); ;

#endif // HX1_HPP
