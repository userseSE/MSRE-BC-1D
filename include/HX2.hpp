#ifndef HX2_HPP
#define HX2_HPP

#include <vector>

// Function to solve the heat exchanger 2 problem
std::vector<std::vector<double>> HX2(std::vector<std::vector<double>> &y_hx2, 
                                     double Ts_HX2_L, double Tss_HX2_0, int step);

#endif // HX2_HPP
