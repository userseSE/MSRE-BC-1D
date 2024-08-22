#include "ode_solver.hpp"
#include "parameters.hpp"
// #include "rkf45.hpp"
// #include "solve_ivp.hpp"
#include "radau.hpp"
#include "sundials.hpp"

#include <vector>
// #include <cmath>
// #include <algorithm>
#include <iostream>
#include <functional>
// #include <boost/numeric/odeint.hpp>


// using namespace boost::numeric::odeint;

// // Define the type for the state vector
// typedef std::vector<double> state_type;

// // Define the ODE system (dy/dt = f(t, y))
// void ode_system(const state_type &y, state_type &dy, double t, ODEFunction vector_to_be_solved) {
//     dy = vector_to_be_solved(t, y);
// }

// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                 const std::vector<double>& bc, 
//                                 ODEFunction vector_to_be_solved, 
//                                 double dt) {

//     double t_end = 1.0;
//     state_type y = ic;  // Initial conditions
//     std::vector<double> results;

//     // Integrate the ODEs using odeint with the specified time step
//     runge_kutta4<state_type> stepper;
//     integrate_const(stepper, 
//                     std::bind(ode_system, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, vector_to_be_solved),
//                     y, 
//                     0.0, t_end, dt);

//     results = y;  // The final state after integration
//     return results;
// }

// using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                const std::vector<double>& bc, 
//                                ODEFunction vector_to_be_solved, 
//                                double dt) {

//     double tolerance = 1e-6;
//     int max_iterations = 100;
//     // Initial condition vector
//     std::vector<double> y = ic;
//     std::vector<std::vector<double>> solution;
//     std::vector<double> results;
//     int time_steps = static_cast<int>(1.0 / dt);
//     double t = 0.0;

//     solution.push_back(y);  // Initial state

//     // Backward Euler integration loop
//     for (int step = 0; step < time_steps; ++step) {
//         std::vector<double> y_next = y;  // Initial guess for y_next

//         for (int iter = 0; iter < max_iterations; ++iter) {
//             std::vector<double> dy = vector_to_be_solved(t + dt, y_next);
//             std::vector<double> y_new(y.size());

//             bool converged = true;
//             for (size_t i = 0; i < y.size(); ++i) {
//                 y_new[i] = y[i] + dt * dy[i];
//                 if (std::fabs(y_new[i] - y_next[i]) > tolerance) {
//                     converged = false;
//                 }
//             }

//             y_next = y_new;

//             if (converged) {
//                 break;
//             }
//         }

//         // Update the time and state
//         t += dt;
//         y = y_next;
//         solution.push_back(y);
//     }

//     // std::cout << "Number of steps: " << solution.size() << std::endl;

//     results = solution.back();
//     return results;
// }
// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                             const std::vector<double>& bc, 
//                                             ODEFunction vector_to_be_solved, double dt) {
//     // Initial condition vector
//     std::vector<double> y = ic;
//     std::vector<std::vector<double>> solution;
//     std::vector<double> results;
//     // Time vector
//     int time_steps = static_cast<int>(1.0 / dt);
//     double t = 0.0;
    
//     solution.push_back(y);  // Initial state

//     // Simple Euler integration loop for demonstration
//     for (int step = 0; step < time_steps; ++step) {
//         std::vector<double> dy = vector_to_be_solved(t, y);
        
//         for (size_t i = 0; i < y.size(); ++i) {
//             y[i] += dt * dy[i];  // Euler step
//         }
        
//         t += dt;
//         solution.push_back(y);
//     }
//     // std::cout<<solution.size()<<std::endl;

//     results=solution.back();
//     // return results;
//     return results;
// }
// Define the ODE function type
// typedef void (* ODEFunction)(double t, const double y[], double yp[]);
// Define the ODE function type as a standard function
// using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

// // Global variable to hold the current ODE function
// ODEFunction* current_ode_function = nullptr;

// // Wrapper to convert std::function to a function pointer compatible with RKF45
// void ode_function_wrapper(double t, double y[], double yp[]) {
//     std::vector<double> y_vec(y, y + 2);  // Assuming y has 2 elements (adjust based on your system)
//     std::vector<double> yp_vec = (*current_ode_function)(t, y_vec);
//     for (size_t i = 0; i < yp_vec.size(); ++i) {
//         yp[i] = yp_vec[i];
//     }
// }

// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                 const std::vector<double>& bc, 
//                                 ODEFunction vector_to_be_solved, double dt) {
//     // Set the global pointer to the current ODE function
//     current_ode_function = &vector_to_be_solved;

//     // Initial condition vector
//     std::vector<double> y = ic;
//     std::vector<double> yp(y.size(), 0.0);  // Derivatives vector
//     std::vector<std::vector<double>> solution;
//     std::vector<double> results;

//     // Time vector
//     double t = 0.0;
//     double tout = t + dt;
//     double relerr = 1.0e-1;
//     double abserr = 1.0e-1;
//     int flag = 1;  // Initial call to RKF45

//     solution.push_back(y);  // Initial state

//     // RKF45 integration loop
//     while (t < 1.0) {
//         int status = r8_rkf45(ode_function_wrapper, y.size(), y.data(), yp.data(), &t, tout, &relerr, abserr, flag);
        
//         if (status != 2 && status != -2) {
//             std::cerr << "RKF45 failed with flag = " << status << std::endl;
//             break;
//         }
        
//         solution.push_back(y);
//         tout = t + dt;  // Update target time
//     }

//     results = solution.back();
//     return results;
// }

// #include "solve_ivp.hpp" // Include the header for solve_ivp

// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                const std::vector<double>& bc, 
//                                std::function<std::vector<double>(double, const std::vector<double>&)> vector_to_be_solved, 
//                                double dt) {
//     // Initial condition vector
//     std::vector<double> y = ic;
//     std::pair<double, double> t_span(0.0, dt);
//     std::vector<std::vector<double>> solution;

//     // // Wrap the ODE function to match the signature expected by solve_ivp
//     // auto wrapper_function = [&](double t, const std::vector<double>& y, std::vector<double>& yp) {
//     //     yp = vector_to_be_solved(t, y);
//     // };

//     // Call solve_ivp with the wrapped function
//     solution = solve_ivp(vector_to_be_solved, t_span, y, dt);

//     return solution.back();
// }

// std::vector<double> ode_solver(const std::vector<double>& ic, 
//                                const std::vector<double>& bc, 
//                                ODEFunction vector_to_be_solved, double dt) {
//     // Time span (assuming you want to integrate from 0 to 1)
//     std::pair<double, double> t_span = {0.0, 1.0};

//     // Call solve_ivp with the provided ODE function and initial conditions
//     auto solution = solve_ivp(vector_to_be_solved, t_span, ic, dt);

//     // Return the last computed value (at the final time step)
//     return solution.back();
// }
std::vector<double> ode_solver(const std::vector<double>& y0, 
                               const std::vector<double>& bc, 
                               ODEFunction vector_to_be_solved, 
                               double dt) {
    std::vector<double> result;
    // std::vector<std::vector<double>> solution;

    result = radau_solver(vector_to_be_solved, 0.0, 1.0, y0);
    // result= sundials(y0, {}, vector_to_be_solved, dt);
    return result;
}


