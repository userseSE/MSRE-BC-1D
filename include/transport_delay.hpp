#ifndef TRANSPORT_DELAY_HPP
#define TRANSPORT_DELAY_HPP

#include <vector>
#include <Eigen/Core>

double transport_delay(double T0, int time_delay, double initial_output, 
                       Eigen::VectorXd& buffer, int step);

#endif // TRANSPORT_DELAY_HPP
