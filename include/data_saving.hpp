#ifndef PLOTTING_HPP
#define PLOTTING_HPP

#include <vector>
#include <Eigen/Core>

void save_results(const Eigen::VectorXd& rho_matrix, 
                  const Eigen::VectorXd& phi_middle_matrix,
                  const Eigen::MatrixXd& ci_middle_matrix, 
                  const Eigen::VectorXd& temperature_fuel_middle_matrix);
void save_spacial_results(const Eigen::VectorXd& phi, 
                          const Eigen::VectorXd& ci, 
                          const Eigen::VectorXd& temperature_fuel, 
                          const Eigen::VectorXd& temperature_graphite, 
                          const Eigen::VectorXd& Ts_HX1, 
                          const Eigen::VectorXd& Tss_HX1, 
                          const Eigen::VectorXd& Tss_HX2, 
                          const Eigen::VectorXd& Tsss_HX2);
#endif // PLOTTING_HPP
