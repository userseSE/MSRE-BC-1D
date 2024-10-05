#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <functional>
#include <iostream>

// #include <cvode/cvode.h>             // Main CVODE integrator header
// #include <nvector/nvector_serial.h>   // N_Vector for serial operations
// #include <cvode/cvode_dense.h>        // Dense linear solver
// #include <sundials/sundials_types.h>  // Definitions of realtype

#include "ode_solver.hpp"
#include "parameters.hpp"

using state_type = VectorXd;
using namespace boost::numeric::odeint;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// // Custom Rosenbrock Stepper
// class Rosenbrock4Eigen {
// public:
//     using state_type = VectorXd;

//     // Constructor for tolerances
//     Rosenbrock4Eigen(double atol = 1.0e-6, double rtol = 1.0e-9)
//         : atol_(atol), rtol_(rtol) {}

//     // Perform a single step of the Rosenbrock method
//     void do_step(std::function<void(const state_type &, state_type &,
//     double)> &ode_func,
//                  state_type &y, double t, double dt) {
//         // Define necessary variables for the Rosenbrock method
//         state_type k1 = state_type::Zero(y.size());
//         state_type k2 = state_type::Zero(y.size());

//         // Compute the first stage of Rosenbrock (using ode_func)
//         ode_func(y, k1, t);  // First evaluation of the ODE function

//         // Perform a single Rosenbrock step (simplified)
//         state_type y_next = y + dt * k1;  // Euler's method approximation for
//         illustration y = y_next;
//     }

// private:
//     double atol_, rtol_;
// };

// // ODE solver function using Rosenbrock method for Eigen::VectorXd
// VectorXd ode_solver(const VectorXd &y0, std::function<void(const VectorXd &,
// VectorXd &, double)> &ode_func) {
//     double t0 = 0.0;
//     double t_end = 1.0;

//     // Create a copy of the initial condition
//     VectorXd y = y0;

//     // Define the custom Rosenbrock stepper
//     Rosenbrock4Eigen stepper;

//     // Manual time-stepping loop
//     double t = t0;
//     while (t < t_end) {
//         stepper.do_step(ode_func, y, t, dt);  // One step of the custom
//         Rosenbrock method t += dt;
//     }

//     // Return the solution
//     return y;
// }
// class RK23Eigen {
// public:
//     using state_type = VectorXd;

//     // Constructor for tolerances
//     RK23Eigen(double atol = 1.0e-3, double rtol = 1.0e-6)
//         : atol_(atol), rtol_(rtol) {}

//     // Perform a single step of the RK23 method
//     void do_step(std::function<void(const state_type &, state_type &,
//     double)> &ode_func,
//                  state_type &y, double &t, double &dt) {
//         // Define necessary variables for RK23 method
//         state_type k1 = state_type::Zero(y.size());
//         state_type k2 = state_type::Zero(y.size());
//         state_type k3 = state_type::Zero(y.size());
//         state_type k4 = state_type::Zero(y.size());

//         // Compute the stages of RK23
//         ode_func(y, k1, t);                       // k1 = f(t, y)
//         ode_func(y + dt * 0.5 * k1, k2, t + dt / 2.0);  // k2 = f(t + dt/2, y
//         + dt/2 * k1) ode_func(y + dt * (0.75 * k1), k3, t + dt * 3.0 / 4.0);
//         // k3 = f(t + 3dt/4, y + 3dt/4 * k1)

//         // Compute RK2 solution (low-order)
//         state_type y_rk2 = y + dt * (k1 * 2.0 / 9.0 + k2 * 1.0 / 3.0 + k3
//         * 4.0 / 9.0);

//         // Compute RK3 solution (higher-order)
//         ode_func(y_rk2, k4, t + dt);  // k4 = f(t + dt, y_rk2)
//         state_type y_rk3 = y + dt * (k1 * 7.0 / 24.0 + k2 * 1.0 / 4.0 + k3
//         * 1.0 / 3.0 + k4 / 8.0);

//         // Error estimate (difference between RK2 and RK3 solutions)
//         state_type error_estimate = y_rk3 - y_rk2;

//         // Compute the maximum error tolerance (based on relative and
//         absolute tolerances) double error_norm =
//         error_estimate.lpNorm<Eigen::Infinity>(); double tolerance = rtol_ *
//         y.lpNorm<Eigen::Infinity>() + atol_;

//         // Adjust the step size based on the error
//         if (error_norm > tolerance) {
//             dt *= std::max(0.5, 0.9 * std::pow(tolerance / error_norm, 1.0
//             / 3.0));  // Decrease step size
//         } else {
//             y = y_rk3;  // Accept the higher-order solution
//             t += dt;    // Update time
//             dt *= std::min(2.0, 0.9 * std::pow(tolerance / error_norm, 1.0
//             / 3.0));  // Increase step size
//         }
//     }

// private:
//     double atol_, rtol_;
// };

// VectorXd ode_solver(const VectorXd &y0, std::function<void(const VectorXd &,
// VectorXd &, double)> &ode_func) {
//     double t0 = 0.0;
//     double t_end = 1.0;

//     // Create a copy of the initial condition
//     VectorXd y = y0;

//     // Define the RK23 stepper
//     RK23Eigen stepper;

//     // Manual time-stepping loop
//     double t = t0;
//     double dt = 0.01;  // Initial guess for time step size
//     while (t < t_end) {
//         if (t + dt > t_end) dt = t_end - t;  // Adjust final step size
//         stepper.do_step(ode_func, y, t, dt);  // One step of RK23 method
//     }

//     // Return the solution
//     return y;
// }

class RungeKutta {
public:
  RungeKutta(double rtol = 1e-3, double atol = 1e-6,
             double max_step = std::numeric_limits<double>::infinity())
      : rtol_(rtol), atol_(atol), max_step_(max_step), min_factor_(0.2),
        max_factor_(10.0), safety_(0.9) {
    initialize_coefficients();
  }

  // Function to perform a single RK step
  void rk_step(const std::function<void(double, const state_type &,
                                        state_type &)> &ode_func,
               double t, state_type &y, state_type &f, double h, MatrixXd &K) {

                // std::cout << "rk_step: t = " << t << std::endl;
// std::cout << "y before step: " << y.transpose() << std::endl;

    // Stages (s=0 is the initial stage)
    K.col(0) = f; // Set K[0] = f

    // Loop over each stage to compute K values
    for (int s = 1; s < A.rows(); ++s) {
      state_type dy = state_type::Zero(y.size());
      for (int j = 0; j < s; ++j) {
        dy += A(s, j) * K.col(j);
      }
      dy *= h;
      state_type temp_y = y + dy; // Ensure we pass a VectorXd, not a block

      // Fix here: create a temporary VectorXd from K.col(s)
      VectorXd temp_Ks = K.col(s); // Convert Eigen block to a VectorXd
      ode_func(t + C(s) * h, temp_y,
               temp_Ks); // Pass temp_Ks as the third argument

      // Assign the result back to K.col(s) after the function call
      K.col(s) = temp_Ks;
    }

    // Compute y_new using the B coefficients (Runge-Kutta solution)
    state_type y_new = y;
    for (int s = 0; s < B.size(); ++s) {
      y_new += h * B(s) * K.col(s);
    }

    // Update y to y_new
    y = y_new;

    // Compute new f for the next step
    VectorXd temp_f = K.col(K.cols() - 1);
    ode_func(t + h, y, temp_f);
    K.col(K.cols() - 1) = temp_f;
  }

  // Function to compute the error norm using the E vector
  double estimate_error_norm(const MatrixXd &K, double h, const state_type &y,
                             const state_type &y_new) {
    state_type scale = state_type::Constant(y.size(), atol_) +
                       rtol_ * y.cwiseAbs().cwiseMax(y_new.cwiseAbs());
    state_type error = state_type::Zero(y.size());
    for (int s = 0; s < E.size(); ++s) {
      error += E(s) * K.col(s);
    }
    error *= h;
    return error.cwiseAbs().cwiseQuotient(scale).maxCoeff();
  }

  // Main ODE solving function
  void solve(const std::function<void(double, const state_type &, state_type &)> &ode_func,
             state_type &y) {
    double t0=0.0; 
    double t1=1.0;
    double t = t0;
    double h = max_step_; // Initial step size
    state_type f(y.size());
    ode_func(t, y, f); // Compute initial f

    MatrixXd K(y.size(), A.rows() + 1); // Runge-Kutta K stages

    while (t < t1) {
      if (t + h > t1) {
        h = t1 - t; // Adjust for final step
      }

      state_type y_old = y;
      rk_step(ode_func, t, y, f, h, K);

      double error_norm = estimate_error_norm(K, h, y_old, y);

      if (error_norm <= 1.0) {
        // Step accepted, move to the next time step
        t += h;
        h *= std::min(
            max_factor_,
            safety_ * std::pow(error_norm, -1.0 / 4.0)); // Increase step size
      } else {
        // Step rejected, reduce step size
        h *= std::max(min_factor_, safety_ * std::pow(error_norm, -1.0 / 4.0));
        y = y_old; // Restore the previous y
      }
    }
  }

private:
  double rtol_, atol_, max_step_, min_factor_, max_factor_, safety_;

  // Coefficients for the Runge-Kutta method
  MatrixXd A = MatrixXd::Zero(16, 16); // A matrix (16x16)
  VectorXd B = VectorXd::Zero(12);     // B vector
  VectorXd C = VectorXd::Zero(16);     // C vector
  VectorXd E = VectorXd::Zero(13);     // E vector for error estimation

  // Initialize coefficients for RK method
  void initialize_coefficients() {
    // Fill in the coefficients here, using your Python-provided coefficients
    // C << 0.0, 0.0526001519587677318785587544488,
    // 0.0789002279381515978178381316732, /* continue with C values */
    // ... similarly initialize B, E, and the A matrix
  }
};

VectorXd ode_solver(const VectorXd &y0,
    const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func) {
  
  double t0 = 0.0;
  double t_end = 1.0;

  // Create a copy of the initial condition
  VectorXd y = y0;

  // Define the RK23 stepper
  // Instantiate the solver
  RungeKutta solver;

  // Solve from t=0 to t=1
//   std::cout << "Before solve, y: " << y.transpose() << std::endl;
  solver.solve(ode_func, y);
// std::cout << "After solve, y: " << y.transpose() << std::endl;

  return y;
}