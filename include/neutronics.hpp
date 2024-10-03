#ifndef NEUTRONICS_HPP
#define NEUTRONICS_HPP

#include <vector>
#include <functional>

// Function to solve the neutronics problem
std::pair<std::vector<std::vector<double>>, std::vector<double>> neutronics(const std::vector<std::vector<double>> &y_n, const std::vector<double> &rho,
           int step);

#endif // NEUTRONICS_HPP
