#ifndef DATA_SAVING_HPP
#define DATA_SAVING_HPP

#include <Eigen/Core>
#include <string>

void save_results(const Eigen::VectorXd& rho_matrix, 
                  const Eigen::VectorXd& phi1_middle_matrix,
                  const Eigen::VectorXd& phi2_middle_matrix,
                  const Eigen::MatrixXd& ci_middle_matrix, 
                  const Eigen::VectorXd& temperature_fuel_middle_matrix,
                  const std::string& folder);

void save_spacial_results(const Eigen::VectorXd& phi1, 
                          const Eigen::VectorXd& phi2, 
                          const Eigen::VectorXd& ci, 
                          const Eigen::VectorXd& rho,
                          const Eigen::VectorXd& temperature_fuel, 
                          const Eigen::VectorXd& temperature_graphite, 
                          const Eigen::VectorXd& Ts_HX1, 
                          const Eigen::VectorXd& Tss_HX1, 
                          const Eigen::VectorXd& Tss_HX2, 
                          const Eigen::VectorXd& Tsss_HX2,
                          const std::string& folder);

#endif // DATA_SAVING_HPP