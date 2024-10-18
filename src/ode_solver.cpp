#include <Eigen/Dense>
#include <fstream>
#include <functional>
#include <iostream>

#include "ode_solver.hpp"
#include "parameters.hpp"

using state_type = VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Define a class for the Runge-Kutta-Fehlberg (RKF45) method
class RungeKuttaFehlberg45 {
public:
  double tol;          // Tolerance for adaptive step-size control
  double h_min, h_max; // Min and max step sizes

  RungeKuttaFehlberg45(double tolerance = 1e-8, double min_step = 1e-10,
                       double max_step = 0.1)
      : tol(tolerance), h_min(min_step), h_max(max_step) {}

  void solve(
      const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func,
      VectorXd &y, int step, double t0 = 0.0, double t1 = 1.0) {
    double t = t0;
    double h = (t1 - t0) / 100; // Initial step size (can be adjusted)

    while (t < t1) {
      if (t + h > t1)
        h = t1 - t; // Adjust final step size to reach t1

      // Perform a single RKF45 step
      VectorXd y_new, error_estimate;
      rkf45_step(ode_func, t, y, h, y_new, error_estimate);

      // Estimate the error and adjust step size
      double error_norm = error_estimate.norm() / y.norm();
      double safety_factor = 0.9; // Safety factor for step-size control
      double scale = safety_factor * std::pow(tol / error_norm, 0.25);

      if (error_norm < tol) {
        // Accept the step
        t += h;
        y = y_new;
        // std::cout << "t = " << t << ", h = " << h << ", y = " <<
        // y.transpose() << std::endl;
      }

      // Update step size for the next iteration
      h *= std::clamp(scale, 0.5,2.0); // Clamp the scaling to avoid drastic step-size changes
      h = std::clamp(h, h_min, h_max); // Ensure step size stays within bounds

      // if (step >= 1000 && step < 1001) {
      //     std::cout << "error_norm: "<< error_norm <<std::endl;
      //     std::cout <<"h: "<< h <<std::endl;
      // }
      // if (step >= 9999 && step < 10001) {
      //   // Open a file in append mode
      //   std::ofstream outputFile("output.txt", std::ios::app);

      //   if (outputFile.is_open()) {
      //     outputFile << "error_norm: " << error_norm << std::endl;
      //     outputFile << "h: " << h << std::endl;
      //     outputFile.close(); // Close the file
      //   } else {
      //     std::cerr << "Unable to open file";
      //   }
      // }
    }
  }

private:
  void rkf45_step(
      const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func,
      double t, const VectorXd &y, double h, VectorXd &y_new,
      VectorXd &error_estimate) {
    // RKF45 coefficients
    static const double a2 = 0.25, a3 = 3.0 / 8.0, a4 = 12.0 / 13.0, a5 = 1.0,
                        a6 = 0.5;
    static const double b2 = 0.25;
    static const double b3[] = {3.0 / 32.0, 9.0 / 32.0};
    static const double b4[] = {1932.0 / 2197.0, -7200.0 / 2197.0,
                                7296.0 / 2197.0};
    static const double b5[] = {439.0 / 216.0, -8.0, 3680.0 / 513.0,
                                -845.0 / 4104.0};
    static const double b6[] = {-8.0 / 27.0, 2.0, -3544.0 / 2565.0,
                                1859.0 / 4104.0, -11.0 / 40.0};

    static const double c[] = {
        16.0 / 135.0,      0.0,         6656.0 / 12825.0,
        28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0}; // 5th order solution
    static const double c_star[] = {25.0 / 216.0,    0.0,
                                    1408.0 / 2565.0, 2197.0 / 4104.0,
                                    -1.0 / 5.0,      0.0}; // 4th order solution

    // Compute the stages
    VectorXd k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()),
        k5(y.size()), k6(y.size());

    ode_func(t, y, k1);
    ode_func(t + a2 * h, y + b2 * k1 * h, k2);
    ode_func(t + a3 * h, y + b3[0] * k1 * h + b3[1] * k2 * h, k3);
    ode_func(t + a4 * h, y + b4[0] * k1 * h + b4[1] * k2 * h + b4[2] * k3 * h,
             k4);
    ode_func(t + a5 * h,
             y + b5[0] * k1 * h + b5[1] * k2 * h + b5[2] * k3 * h +
                 b5[3] * k4 * h,
             k5);
    ode_func(t + a6 * h,
             y + b6[0] * k1 * h + b6[1] * k2 * h + b6[2] * k3 * h +
                 b6[3] * k4 * h + b6[4] * k5 * h,
             k6);

    // Compute the 4th and 5th order estimates
    y_new = y + h * (c[0] * k1 + c[2] * k3 + c[3] * k4 + c[4] * k5 + c[5] * k6);
    VectorXd y_star = y + h * (c_star[0] * k1 + c_star[2] * k3 +
                               c_star[3] * k4 + c_star[4] * k5);

    // Estimate the local error
    error_estimate = y_new - y_star;
  }
};

VectorXd ode_solver(
    const VectorXd &y0,
    const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func,
    int step) {

  double t0 = 0.0;
  double t_end = 1.0;
  // std::cout << "t0: " << t0 << std::endl;

  // Create a copy of the initial condition
  VectorXd y = y0;

  // Define the RK23 stepper
  RungeKuttaFehlberg45 solver;

  // std::cout << "Before solve, y: " << std::endl;
  // Solve from t=0 to t=1
  solver.solve(ode_func, y, step);

  return y;
}
