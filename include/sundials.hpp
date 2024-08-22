#include<iostream>
#include<functional>

using ODEFunction = std::function<std::vector<double>(double t, const std::vector<double>& y)>;

std::vector<double> sundials(const std::vector<double>& y0, 
                               const std::vector<double>& bc, 
                               ODEFunction vector_to_be_solved, 
                               double dt);