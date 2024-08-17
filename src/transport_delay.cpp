#include "transport_delay.hpp"
#include <vector>

double transport_delay(double T0, int time_delay, double initial_output, 
                       std::vector<double>& buffer, int step) {
    
    double T1;
    
    if (step < time_delay) {
        T1 = initial_output;
        buffer.push_back(T0);
    } else {
        T1 = buffer.front();
        buffer.erase(buffer.begin());  // Remove the first element
        buffer.push_back(T0);
    }
    
    return T1;
}
