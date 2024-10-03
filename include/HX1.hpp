#ifndef HX1_HPP
#define HX1_HPP

#include <vector>

// Function to solve the heat exchanger 1 problem
std::vector<std::vector<double>> HX1(std::vector<std::vector<double>> &y_hx1, 
                                     double Ts_HX1_L, double Tss_HX1_0, int step);

#endif // HX1_HPP
