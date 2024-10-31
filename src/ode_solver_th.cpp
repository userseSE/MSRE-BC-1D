#include "ode_solver_th.hpp"
#include "parameters.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <Eigen>

// Define a class for the Runge-Kutta-Fehlberg (RKF45) method
class RungeKuttaFehlberg45_th {
public:
  double tol;          // Tolerance for adaptive step-size control
  double h_min, h_max; // Min and max step sizes

  RungeKuttaFehlberg45_th(double tolerance = 1e-8, double min_step = 1e-10,
                       double max_step = 0.1) : tol(tolerance), h_min(min_step), h_max(max_step) {}

  void solve(const OdeFuncPointer_th ode_func, double y[length_th], int step,
             Parameters &params, const double q_prime[N], double t0 = 0.0,
             double t1 = 1.0) {
    // std::cout << "Neutronics ODE solver called" << std::endl;
    double t = t0;
    double h = (t1 - t0) / 100; // Initial step size (can be adjusted)

    while (t < t1) {
      if (t + h > t1)
        h = t1 - t; // Adjust final step size to reach t1

      // Perform a single RKF45 step
      double y_new[length_th], error_estimate[length_th];
      rkf45_step_th(ode_func, t, y, params, q_prime, h, y_new, error_estimate);

      // Estimate the error and adjust step size
      double error_estimate_norm = 0.0;
      double y_norm = 0.0;
      error_estimate_norm= stableNorm(error_estimate, length_th);
      y_norm=stableNorm(y, length_th);
      double error_norm = error_estimate_norm / y_norm;

      double safety_factor = 0.9; // Safety factor for step-size control
      double scale = safety_factor * std::pow(tol / error_norm, 0.25);

      if (error_norm < tol) {
        t += h;
        for (int i = 0; i < length_th; ++i) {
          y[i] = y_new[i];
        }
      }
      // Update step size for the next iteration
      h *= std::clamp(scale, 0.5,2.0);
      h = std::clamp(h, h_min, h_max);
    }
  }

private:
  double stableNorm(const double* arr, int N) {
    Eigen::Map<const Eigen::VectorXd> eigen_vec(arr, N); // Map array to Eigen vector
    return eigen_vec.stableNorm();
}

  void clamp(double &value, double low, double high) {
    if (value < low) {
      value = low;
    } else if (value > high) {
      value = high;
    }
  }

  void rkf45_step_th(const OdeFuncPointer_th ode_func, double t,
                  const double y[length_th], Parameters &params,
                  const double q_prime[N], double h, double y_new[length_th],
                  double error_estimate[length_th]) {
    // std::cout << "Neutronics RKF45 step called" << std::endl;

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
    double k1[length_th], k2[length_th], k3[length_th],
        k4[length_th], k5[length_th], k6[length_th];
    double result[length_th];

    // Ensure k1, k2, etc. are properly initialized to zero
    for (int i = 0; i < length_th; ++i) {
      k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = k6[i] = 0.0;
    }

    // First stage (k1)
    ode_func(t, y, k1, params, q_prime);
    // Compute intermediate results for k2
    for (int i = 0; i < length_th; ++i) {
      result[i] = y[i] + b2 * h * k1[i];
    }
    ode_func(t + a2 * h, result, k2, params, q_prime);
    // Compute intermediate results for k3
    for (int i = 0; i < length_th; ++i) {
      result[i] = y[i] + b3[0] * h * k1[i] + b3[1] * h * k2[i];
    }
    ode_func(t + a3 * h, result, k3, params, q_prime);

    // Compute intermediate results for k4
    for (int i = 0; i < length_th; ++i) {
      result[i] =
          y[i] + b4[0] * h * k1[i] + b4[1] * h * k2[i] + b4[2] * h * k3[i];
    }
    ode_func(t + a4 * h, result, k4, params, q_prime);

    // Compute intermediate results for k5
    for (int i = 0; i < length_th; ++i) {
      result[i] = y[i] + b5[0] * h * k1[i] + b5[1] * h * k2[i] + b5[2] * h * k3[i] + b5[3] * h * k4[i];
    }
    ode_func(t + a5 * h, result, k5, params, q_prime);

    // Compute intermediate results for k6
    for (int i = 0; i < length_th; ++i) {
      result[i] = y[i] + b6[0] * h * k1[i] + b6[1] * h * k2[i] + b6[2] * h * k3[i] + b6[3] * h * k4[i] + b6[4] * h * k5[i];
    }
    ode_func(t + a6 * h, result, k6, params, q_prime);

    // Compute the 4th and 5th order estimates (y_new and y_star)
    double y_star[length_th];
    for (int i = 0; i < length_th; ++i) {
      y_new[i] = y[i] + h * (c[0] * k1[i] + c[2] * k3[i] + c[3] * k4[i] + c[4] * k5[i] + c[5] * k6[i]);
      y_star[i] = y[i] + h * (c_star[0] * k1[i] + c_star[2] * k3[i] + c_star[3] * k4[i] + c_star[4] * k5[i]);
    }

    // Estimate the local error
    for (int i = 0; i < length_th; ++i) {
      error_estimate[i] = y_new[i] - y_star[i];
    }
  }
};
void ode_solver_th(double y[length_th], OdeFuncPointer_th ode_func, int step,
                      Parameters &params, const double q_prime[N]) {

  RungeKuttaFehlberg45_th solver;

  solver.solve(ode_func, y, step, params, q_prime);
}