#ifndef SOLVE_IVP_HPP
#define SOLVE_IVP_HPP

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

// Type alias for the ODE function signature
using ODEFunction = std::function<void(double t, const std::vector<double>& y, std::vector<double>& dydt)>;
// using ODEFunction = std::function<void(double, const std::vector<double>&, std::vector<double>&)>;
// using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

// Function declarations

/**
 * @brief Perform a single integration step using the Runge-Kutta-Fehlberg method (RK45).
 * 
 * @param f The ODE system to solve.
 * @param t0 The initial time.
 * @param tf The final time.
 * @param y0 The initial state vector.
 * @param dt The time step size.
 */
void runge_kutta_45(ODEFunction f, double t0, double tf, std::vector<double>& y0, double dt);

/**
 * @brief Solve an initial value problem for a system of ODEs.
 * 
 * @param f The ODE system to solve.
 * @param t_span A pair representing the start and end times.
 * @param y0 The initial state vector.
 * @param dt The time step size.
 * @return std::vector<std::vector<double>> The solution vector containing the state at each time step.
 */
std::vector<std::vector<double>> solve_ivp(ODEFunction f, std::pair<double, double> t_span, const std::vector<double>& y0, double dt);

#endif // SOLVE_IVP_HPP
