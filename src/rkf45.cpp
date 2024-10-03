// #include <vector>
// #include <iostream>
// #include <cmath>
// #include <algorithm> // for std::max_element

// #include "rkf45.hpp"

// template <typename Func>
// void rk45_step(Func&& f, std::vector<double>& y, double& t, double& dt) {

//     double tol = 1e-6; // 定义容忍度
//     // 定义Runge-Kutta-Fehlberg 4(5)的系数
//     static constexpr double a2 = 0.25, a3 = 3.0 / 8.0, a4 = 12.0 / 13.0, a5 = 1.0, a6 = 0.5;
//     static constexpr double b21 = 0.25;
//     static constexpr double b31 = 3.0 / 32.0, b32 = 9.0 / 32.0;
//     static constexpr double b41 = 1932.0 / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0;
//     static constexpr double b51 = 439.0 / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0;
//     static constexpr double b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;

//     static constexpr double c1 = 16.0 / 135.0, c3 = 6656.0 / 12825.0, c4 = 28561.0 / 56430.0;
//     static constexpr double c5 = -9.0 / 50.0, c6 = 2.0 / 55.0;
    
//     static constexpr double d1 = 1.0 / 360.0, d3 = -128.0 / 4275.0, d4 = -2197.0 / 75240.0;
//     static constexpr double d5 = 1.0 / 50.0, d6 = 2.0 / 55.0;

//     size_t N = y.size(); // 获取向量大小

//     // 工作向量
//     std::vector<double> k1(N), k2(N), k3(N), k4(N), k5(N), k6(N), y_temp(N), error(N);

//     // 计算k1
//     f(t, y, k1);
//     for (size_t i = 0; i < N; ++i) {
//         y_temp[i] = y[i] + dt * b21 * k1[i];
//     }

//     // 计算k2
//     f(t + a2 * dt, y_temp, k2);
//     for (size_t i = 0; i < N; ++i) {
//         y_temp[i] = y[i] + dt * (b31 * k1[i] + b32 * k2[i]);
//     }

//     // 计算k3
//     f(t + a3 * dt, y_temp, k3);
//     for (size_t i = 0; i < N; ++i) {
//         y_temp[i] = y[i] + dt * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
//     }

//     // 计算k4
//     f(t + a4 * dt, y_temp, k4);
//     for (size_t i = 0; i < N; ++i) {
//         y_temp[i] = y[i] + dt * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
//     }

//     // 计算k5
//     f(t + a5 * dt, y_temp, k5);
//     for (size_t i = 0; i < N; ++i) {
//         y_temp[i] = y[i] + dt * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
//     }

//     // 计算k6
//     f(t + a6 * dt, y_temp, k6);

//     // 计算y和误差
//     for (size_t i = 0; i < N; ++i) {
//         double y4 = y[i] + dt * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i]);
//         double y5 = y[i] + dt * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i]);
//         error[i] = std::abs(y5 - y4);
//         y[i] = y5;  // 使用高阶解更新y
//     }

//     // 调整步长
//     double max_error = *std::max_element(error.begin(), error.end());
//     if (max_error > tol) {
//         // 如果误差过大，减小步长
//         dt *= 0.9 * std::pow(tol / max_error, 0.25);
//     } else {
//         // 如果误差在容许范围内，前进时间
//         t += dt;
//         if (max_error < tol / 10) {
//             // 如果误差非常小，增大步长
//             dt *= 0.9 * std::pow(tol / max_error, 0.2);
//         }
//     }
// }