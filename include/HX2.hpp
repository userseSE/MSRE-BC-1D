#ifndef HX2_HPP
#define HX2_HPP

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::VectorXd;
using Eigen::MatrixXd;

// Function to solve the heat exchanger 2 problem
VectorXd HX2(VectorXd& y_hx2, double Ts_HX2_L, double Tss_HX2_0, int step) ;

#endif // HX2_HPP
