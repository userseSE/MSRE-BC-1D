#include <Eigen/Dense>

double transport_delay(double T0, int time_delay, double initial_output, 
                       Eigen::VectorXd& buffer, int step) {
    
    double T1;

    if (step < time_delay) {
        T1 = initial_output;
        Eigen::VectorXd new_buffer(buffer.size() + 1);
        new_buffer.head(buffer.size()) = buffer;
        new_buffer[buffer.size()] = T0;
        buffer = new_buffer;
    } else {
        T1 = buffer[0];  // Front of the buffer
        buffer.segment(0, buffer.size() - 1) = buffer.tail(buffer.size() - 1);
        buffer[buffer.size() - 1] = T0;
    }
    
    return T1;
}

