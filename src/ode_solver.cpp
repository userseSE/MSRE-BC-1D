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

typedef std::vector<double> state_type;

std::vector<double> ode_solver(const std::vector<double> &y0,
                               ODEFunction vector_to_be_solved) {
  // Copy initial state y0 to result
  std::vector<double> result = y0;

  std::cout << "y0: " << y0.front() << std::endl;

  typedef runge_kutta_dopri5<state_type> stepper_type;

  // Use Runge-Kutta-Dopri5 stepper for adaptive time-stepping
  runge_kutta_dopri5<state_type> stepper;
  
  // Integrate the system with adaptive time-stepping
  integrate_adaptive(stepper, vector_to_be_solved, result, t0, t1, dt);

  // TODO: TRY Implicit method

  return result;
}
