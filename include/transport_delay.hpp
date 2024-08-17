#ifndef TRANSPORT_DELAY_HPP
#define TRANSPORT_DELAY_HPP

#include <vector>

double transport_delay(double T0, int time_delay, double initial_output, 
                       std::vector<double>& buffer, int step);

#endif // TRANSPORT_DELAY_HPP
