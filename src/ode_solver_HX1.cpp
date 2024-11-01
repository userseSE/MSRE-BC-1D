#include "ode_solver_HX1.hpp"
#include "parameters.hpp"
#include <cmath>
// Define a class for the Runge-Kutta-Fehlberg (RKF45) method
class RungeKuttaFehlberg45_hx {
public:
  float tol;          // Tolerance for adaptive step-size control
  float h_min, h_max; // Min and max step sizes

  RungeKuttaFehlberg45_hx(float tolerance = 1e-8, float min_step = 1e-10,
                       float max_step = 0.1) : tol(tolerance), h_min(min_step), h_max(max_step) {}

  void solve(const OdeFuncPointer_hx1 ode_func, float y[length_hx], int step,
             Param_HX1 &params, float t0 = 0.0,
             float t1 = 1.0) {
    // std::cout << "Neutronics ODE solver called" << std::endl;
    float t = t0;
    float h = (t1 - t0) / 100; // Initial step size (can be adjusted)

    while (t < t1) {
      if (t + h > t1)
        h = t1 - t; // Adjust final step size to reach t1

      // Perform a single RKF45 step
      float y_new[length_hx], error_estimate[length_hx];
      rkf45_step_hx(ode_func, t, y, params , h, y_new, error_estimate);

      // Estimate the error and adjust step size
      float error_estimate_norm = 0.0;
      float y_norm = 0.0;
      // error_estimate_norm= stableNorm(error_estimate, length_hx);
      // y_norm=stableNorm(y, length_hx);
      calculateNorm(error_estimate, error_estimate_norm);
      calculateNorm(y, y_norm);
      float error_norm = error_estimate_norm / y_norm;

      float safety_factor = 0.9; // Safety factor for step-size control
      float scale = safety_factor * std::pow(tol / error_norm, 0.25);

      if (error_norm < tol) {
        t += h;
        for (int i = 0; i < length_hx; ++i) {
          y[i] = y_new[i];
        }
      }
      float h_temp=0.0;
      clamp(h_temp, scale, 0.5,2.0);
      h *= h_temp;
      clamp(h,h, h_min, h_max);
    }
  }

private:
void clamp(float& result, const float& value, const float& min_val, const float& max_val) {
    if (value < min_val) result= min_val;
    if (value > max_val) result= max_val;
    if (value >=min_val && value<=max_val) result= value;
}
void calculateNorm(const float arr[length_hx], float &norm) {
    norm = 0.0;
    for (int i = 0; i < length_hx; ++i) {
        norm += arr[i] * arr[i];
    }
    norm = std::sqrt(norm);
}
  // void clamp(float &value, float low, float high) {
  //   if (value < low) {
  //     value = low;
  //   } else if (value > high) {
  //     value = high;
  //   }
  // }

  void rkf45_step_hx(const OdeFuncPointer_hx1 ode_func, float t,
                  const float y[length_hx], Param_HX1 &params, float h, float y_new[length_hx],
                  float error_estimate[length_hx]) {
    // std::cout << "Neutronics RKF45 step called" << std::endl;

    // RKF45 coefficients
    static const float a2 = 0.25, a3 = 3.0 / 8.0, a4 = 12.0 / 13.0, a5 = 1.0,
                        a6 = 0.5;
    static const float b2 = 0.25;
    static const float b3[] = {3.0 / 32.0, 9.0 / 32.0};
    static const float b4[] = {1932.0 / 2197.0, -7200.0 / 2197.0,
                                7296.0 / 2197.0};
    static const float b5[] = {439.0 / 216.0, -8.0, 3680.0 / 513.0,
                                -845.0 / 4104.0};
    static const float b6[] = {-8.0 / 27.0, 2.0, -3544.0 / 2565.0,
                                1859.0 / 4104.0, -11.0 / 40.0};

    static const float c[] = {
        16.0 / 135.0,      0.0,         6656.0 / 12825.0,
        28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0}; // 5th order solution
    static const float c_star[] = {25.0 / 216.0,    0.0,
                                    1408.0 / 2565.0, 2197.0 / 4104.0,
                                    -1.0 / 5.0,      0.0}; // 4th order solution

    // Compute the stages
    float k1[length_hx], k2[length_hx], k3[length_hx],
        k4[length_hx], k5[length_hx], k6[length_hx];
    float result[length_hx];

    // Ensure k1, k2, etc. are properly initialized to zero
    for (int i = 0; i < length_hx; ++i) {
      k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = k6[i] = 0.0;
    }

    // First stage (k1)
    ode_func(t, y, k1, params);
    // for(int i=0; i<length_hx; ++i){
    //   std::cout << "k1[" << i << "]: " << k1[i] << std::endl;
    // }
    // Compute intermediate results for k2
    for (int i = 0; i < length_hx; ++i) {
      result[i] = y[i] + b2 * h * k1[i];
    }
    ode_func(t + a2 * h, result, k2, params);
    // for(int i=0; i<length_hx; ++i){
    //   std::cout << "k2[" << i << "]: " << k2[i] << std::endl;
    // }

    // Compute intermediate results for k3
    for (int i = 0; i < length_hx; ++i) {
      result[i] = y[i] + b3[0] * h * k1[i] + b3[1] * h * k2[i];
    }
    ode_func(t + a3 * h, result, k3, params);

    // Compute intermediate results for k4
    for (int i = 0; i < length_hx; ++i) {
      result[i] = y[i] + b4[0] * h * k1[i] + b4[1] * h * k2[i] + b4[2] * h * k3[i];
    }
    ode_func(t + a4 * h, result, k4, params);

    // Compute intermediate results for k5
    for (int i = 0; i < length_hx; ++i) {
      result[i] = y[i] + b5[0] * h * k1[i] + b5[1] * h * k2[i] + b5[2] * h * k3[i] + b5[3] * h * k4[i];
    }
    ode_func(t + a5 * h, result, k5, params);

    // Compute intermediate results for k6
    for (int i = 0; i < length_hx; ++i) {
      result[i] = y[i] + b6[0] * h * k1[i] + b6[1] * h * k2[i] + b6[2] * h * k3[i] + b6[3] * h * k4[i] + b6[4] * h * k5[i];
    }
    ode_func(t + a6 * h, result, k6, params);

    // Compute the 4th and 5th order estimates (y_new and y_star)
    float y_star[length_hx];
    for (int i = 0; i < length_hx; ++i) {
      y_new[i] = y[i] + h * (c[0] * k1[i] + c[2] * k3[i] + c[3] * k4[i] + c[4] * k5[i] + c[5] * k6[i]);
      y_star[i] = y[i] + h * (c_star[0] * k1[i] + c_star[2] * k3[i] + c_star[3] * k4[i] + c_star[4] * k5[i]);
    }
    // Estimate the local error
    for (int i = 0; i < length_hx; ++i) {
      error_estimate[i] = y_new[i] - y_star[i];
    }
  }
};
void ode_solver_hx1(float y[length_hx], OdeFuncPointer_hx1 ode_func, int step,
                      Param_HX1 &params) {

  RungeKuttaFehlberg45_hx solver;

  solver.solve(ode_func, y, step, params);
}