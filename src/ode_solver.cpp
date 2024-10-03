#include "ode_solver.hpp"
#include "parameters.hpp"
// #include "rkf45.hpp"
// // #include "solve_ivp.hpp"
// #include "radau.hpp"

#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>
// #include "ode/ode_rkf_32.h"

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

typedef std::vector<std::vector<double>> state_type;  // 2D vector

std::vector<std::vector<double>> ode_solver(const std::vector<std::vector<double>> &y0,
                                            ODEFunction vector_to_be_solved) {
  // Copy initial state y0 to result
  std::vector<std::vector<double>> result = y0;

  std::cout << "y0: " << y0[0][0] << std::endl;  // Adjusted to access 2D vector

  // Instantiate the Rosenbrock4 stepper (for stiff problems)
  rosenbrock4<state_type> stepper;

  // Integrate the system using adaptive time-stepping
  integrate_adaptive(stepper, vector_to_be_solved, result, t0, t1, dt);

  return result;
}

