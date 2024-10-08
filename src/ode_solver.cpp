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
//         state_type y_next = y + dt * k1;  // Euler's method approximation for illustration 
//         y = y_next;
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
//         stepper.do_step(ode_func, y, t, dt);  // One step of the custom Rosenbrock method 
//         t += dt;
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
//         ode_func(y + dt * 0.5 * k1, k2, t + dt / 2.0);  // k2 = f(t + dt/2, y + dt/2 * k1) 
//         ode_func(y + dt * (0.75 * k1), k3, t + dt * 3.0 / 4.0); // k3 = f(t + 3dt/4, y + 3dt/4 * k1)

//         // Compute RK2 solution (low-order)
//         state_type y_rk2 = y + dt * (k1 * 2.0 / 9.0 + k2 * 1.0 / 3.0 + k3
//         * 4.0 / 9.0);

//         // Compute RK3 solution (higher-order)
//         ode_func(y_rk2, k4, t + dt);  // k4 = f(t + dt, y_rk2)
//         state_type y_rk3 = y + dt * (k1 * 7.0 / 24.0 + k2 * 1.0 / 4.0 + k3
//         * 1.0 / 3.0 + k4 / 8.0);

//         // Error estimate (difference between RK2 and RK3 solutions)
//         state_type error_estimate = y_rk3 - y_rk2;

//         // Compute the maximum error tolerance (based on relative and absolute tolerances) 
//         double error_norm =
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



// class RungeKutta {
// public:
//   RungeKutta(double rtol = 1e-3, double atol = 1e-6,
//              double max_step = std::numeric_limits<double>::infinity())
//       : rtol_(rtol), atol_(atol), max_step_(max_step), min_factor_(0.2),
//         max_factor_(10.0), safety_(0.9) {
//     initialize_coefficients();
//   }

//   // Function to perform a single RK step
//   void rk_step(const std::function<void(double, const state_type &,
//                                         state_type &)> &ode_func,
//                double t, state_type &y, state_type &f, double h, MatrixXd &K) {
//     K.col(0) = f; // Set K[0] = f

//     // Loop over each stage to compute K values
//     for (int s = 1; s < A.rows(); ++s) {
//       state_type dy = state_type::Zero(y.size());
//       for (int j = 0; j < s; ++j) {
//         dy += A(s, j) * K.col(j);
//       }
//       dy *= h;
//       state_type temp_y = y + dy; // Ensure we pass a VectorXd, not a block

//       // Fix here: create a temporary VectorXd from K.col(s)
//       VectorXd temp_Ks = K.col(s); // Convert Eigen block to a VectorXd
//       ode_func(t + C(s) * h, temp_y,
//                temp_Ks); // Pass temp_Ks as the third argument

//       // Assign the result back to K.col(s) after the function call
//       K.col(s) = temp_Ks;
//     }

//     // Compute y_new using the B coefficients (Runge-Kutta solution)
//     state_type y_new = y;
//     for (int s = 0; s < B.size(); ++s) {
//       y_new += h * B(s) * K.col(s);
//     }

//     // Update y to y_new
//     y = y_new;

//     // Compute new f for the next step
//     VectorXd temp_f = K.col(K.cols() - 1);
//     ode_func(t + h, y, temp_f);
//     K.col(K.cols() - 1) = temp_f;
//   }

//   // Function to compute the error norm using the E vector
//   double estimate_error_norm(const MatrixXd &K, double h, const state_type &y,
//                              const state_type &y_new) {
//     // Compute the scaling factor
//     state_type scale = state_type::Constant(y.size(), atol_) +
//                        rtol_ * y.cwiseAbs().cwiseMax(y_new.cwiseAbs());

//     // Estimate error using E3 and E5
//     state_type err5 = state_type::Zero(y.size());
//     state_type err3 = state_type::Zero(y.size());

//     for (int s = 0; s < E5.size(); ++s) {
//       err5 += E5(s) * K.col(s);
//     }
//     for (int s = 0; s < E3.size(); ++s) {
//       err3 += E3(s) * K.col(s);
//     }

//     // Calculate correction factor
//     double err5_norm = err5.cwiseAbs().cwiseQuotient(scale).norm();
//     double err3_norm = err3.cwiseAbs().cwiseQuotient(scale).norm();
//     double denom = std::hypot(err5_norm, 0.1 * err3_norm);

//     if (denom > 0.0) {
//       double correction_factor = err5_norm / denom;
//       return h * err5_norm * correction_factor;
//     }

//     return 0.0;
//   }

//   // Main ODE solving function
//   void solve(const std::function<void(double, const state_type &, state_type &)>
//                  &ode_func,
//              state_type &y) {
//     double t0 = 0.0;
//     double t1 = 1.0;
//     double t = t0;
//     double h = max_step_; // Initial step size
//     state_type f(y.size());
//     ode_func(t, y, f); // Compute initial f

//     MatrixXd K(y.size(), A.rows() + 1); // Runge-Kutta K stages

//     while (t < t1) {
//       if (t + h > t1) {
//         h = t1 - t; // Adjust for final step
//       }

//       state_type y_old = y;
//       rk_step(ode_func, t, y, f, h, K);

//       double error_norm = estimate_error_norm(K, h, y_old, y);

//       if (error_norm <= 1.0) {
//         // Step accepted, move to the next time step
//         t += h;
//         h *= std::min(
//             max_factor_,
//             safety_ * std::pow(error_norm, -1.0 / 4.0)); // Increase step size
//       } else {
//         // Step rejected, reduce step size
//         h *= std::max(min_factor_, safety_ * std::pow(error_norm, -1.0 / 4.0));
//         y = y_old; // Restore the previous y
//       }
//     }
//   }

// private:
//   double rtol_, atol_, max_step_, min_factor_, max_factor_, safety_;

//   // Coefficients for the Runge-Kutta method
//   MatrixXd A = MatrixXd::Zero(16, 16); // A matrix (16x16)
//   VectorXd B = VectorXd::Zero(12);     // B vector
//   VectorXd C = VectorXd::Zero(16);     // C vector
//   VectorXd E3 =
//       VectorXd::Zero(13); // E3 vector for third-order error estimation
//   VectorXd E5 =
//       VectorXd::Zero(13); // E5 vector for fifth-order error estimation

//   // Initialize coefficients for RK method
//   void initialize_coefficients() {
//     std::cout << "Test C" << std::endl;
//     // Initialize C vector (16 elements)
//     C << 0.0, 0.0526001519587677318785587544488,
//         0.0789002279381515978178381316732, 0.118350341907227396726757197510,
//         0.281649658092772603273242802490, 0.333333333333333333333333333333,
//         0.25, 0.307692307692307692307692307692,
//         0.651282051282051282051282051282, 0.6, 0.857142857142857142857142857142,
//         1.0, 1.0, 0.1, 0.2, 0.777777777777777777777777777778;
//     std::cout << "Test A" << std::endl;
//     // Initialize A matrix (16x16 elements, N_STAGES_EXTENDED)
//     A.setZero(); // Initialize the entire matrix to zero
//     A(1, 0) = 5.26001519587677318785587544488e-2;
//     A(2, 0) = 1.97250569845378994544595329183e-2;
//     A(2, 1) = 5.91751709536136983633785987549e-2;
//     A(3, 0) = 2.95875854768068491816892993775e-2;
//     A(3, 2) = 8.87627564304205475450678981324e-2;
//     A(4, 0) = 2.41365134159266685502369798665e-1;
//     A(4, 2) = -8.84549479328286085344864962717e-1;
//     A(4, 3) = 9.24834003261792003115737966543e-1;
//     A(5, 0) = 3.7037037037037037037037037037e-2;
//     A(5, 3) = 1.70828608729473871279604482173e-1;
//     A(5, 4) = 1.25467687566822425016691814123e-1;
//     A(6, 0) = 3.7109375e-2;
//     A(6, 3) = 1.70252211019544039314978060272e-1;
//     A(6, 4) = 6.02165389804559606850219397283e-2;
//     A(6, 5) = -1.7578125e-2;
//     A(7, 0) = 3.70920001185047927108779319836e-2;
//     A(7, 3) = 1.70383925712239993810214054705e-1;
//     A(7, 4) = 1.07262030446373284651809199168e-1;
//     A(7, 5) = -1.53194377486244017527936158236e-2;
//     A(7, 6) = 8.27378916381402288758473766002e-3;
//     A(8, 0) = 6.24110958716075717114429577812e-1;
//     A(8, 3) = -3.36089262944694129406857109825;
//     A(8, 4) = -8.68219346841726006818189891453e-1;
//     A(8, 5) = 2.75920996994467083049415600797e1;
//     A(8, 6) = 2.01540675504778934086186788979e1;
//     A(8, 7) = -4.34898841810699588477366255144e1;
//     A(9, 0) = 4.77662536438264365890433908527e-1;
//     A(9, 3) = -2.48811461997166764192642586468;
//     A(9, 4) = -5.90290826836842996371446475743e-1;
//     A(9, 5) = 2.12300514481811942347288949897e1;
//     A(9, 6) = 1.52792336328824235832596922938e1;
//     A(9, 7) = -3.32882109689848629194453265587e1;
//     A(9, 8) = -2.03312017085086261358222928593e-2;
//     A(10, 0) = -9.3714243008598732571704021658e-1;
//     A(10, 3) = 5.18637242884406370830023853209;
//     A(10, 4) = 1.09143734899672957818500254654;
//     A(10, 5) = -8.14978701074692612513997267357;
//     A(10, 6) = -1.85200656599969598641566180701e1;
//     A(10, 7) = 2.27394870993505042818970056734e1;
//     A(10, 8) = 2.49360555267965238987089396762;
//     A(10, 9) = -3.0467644718982195003823669022;
//     A(11, 0) = 2.27331014751653820792359768449;
//     A(11, 3) = -1.05344954667372501984066689879e1;
//     A(11, 4) = -2.00087205822486249909675718444;
//     A(11, 5) = -1.79589318631187989172765950534e1;
//     A(11, 6) = 2.79488845294199600508499808837e1;
//     A(11, 7) = -2.85899827713502369474065508674;
//     A(11, 8) = -8.87285693353062954433549289258;
//     A(11, 9) = 1.23605671757943030647266201528e1;
//     A(11, 10) = 6.43392746015763530355970484046e-1;
//     A(12, 0) = 5.42937341165687622380535766363e-2;
//     A(12, 5) = 4.45031289275240888144113950566;
//     A(12, 6) = 1.89151789931450038304281599044;
//     A(12, 7) = -5.8012039600105847814672114227;
//     A(12, 8) = 3.1116436695781989440891606237e-1;
//     A(12, 9) = -1.52160949662516078556178806805e-1;
//     A(12, 10) = 2.01365400804030348374776537501e-1;
//     A(12, 11) = 4.47106157277725905176885569043e-2;
//     A(13, 0) = 5.61675022830479523392909219681e-2;
//     A(13, 6) = 2.53500210216624811088794765333e-1;
//     A(13, 7) = -2.46239037470802489917441475441e-1;
//     A(13, 8) = -1.24191423263816360469010140626e-1;
//     A(13, 9) = 1.5329179827876569731206322685e-1;
//     A(13, 10) = 8.20105229563468988491666602057e-3;
//     A(13, 11) = 7.56789766054569976138603589584e-3;
//     A(13, 12) = -8.298e-3;
//     A(14, 0) = 3.18346481635021405060768473261e-2;
//     A(14, 5) = 2.83009096723667755288322961402e-2;
//     A(14, 6) = 5.35419883074385676223797384372e-2;
//     A(14, 7) = -5.49237485713909884646569340306e-2;
//     A(14, 10) = -1.08347328697249322858509316994e-4;
//     A(14, 11) = 3.82571090835658412954920192323e-4;
//     A(14, 12) = -3.40465008687404560802977114492e-4;
//     A(14, 13) = 1.41312443674632500278074618366e-1;
//     A(15, 0) = -4.28896301583791923408573538692e-1;
//     A(15, 5) = -4.69762141536116384314449447206;
//     A(15, 6) = 7.68342119606259904184240953878;
//     A(15, 7) = 4.06898981839711007970213554331;
//     A(15, 8) = 3.56727187455281109270669543021e-1;
//     A(15, 12) = -1.39902416515901462129418009734e-3;
//     A(15, 13) = 2.9475147891527723389556272149;
//     A(15, 14) = -9.15095847217987001081870187138;
//     std::cout << "Test B" << std::endl;
//     // Initialize B vector (12 elements)
//     B = A.row(12).head(12);

//     // E3.resize(13);
//     std::cout << "Test E3" << std::endl;
//     // Initialize E vector for error estimation (13 elements, interpolator
//     // power)
//     for (int i = 0; i < B.size(); ++i) {
//       E3[i] = B[i];
//     }
//     // Modify specific elements of E3 based on Python code
//     E3[0] -= 0.244094488188976377952755905512;
//     E3[8] -= 0.733846688281611857341361741547;
//     E3[11] -= 0.220588235294117647058823529412e-1;

//     // Initialize E5 elements
//     E5[0] = 0.1312004499419488073250102996e-1;
//     E5[5] = -0.1225156446376204440720569753e+1;
//     E5[6] = -0.4957589496572501915214079952;
//     E5[7] = 0.1664377182454986536961530415e+1;
//     E5[8] = -0.3503288487499736816886487290;
//     E5[9] = 0.3341791187130174790297318841;
//     E5[10] = 0.8192320648511571246570742613e-1;
//     E5[11] = -0.2235530786388629525884427845e-1;

//     // std::cout << "Test Done" << std::endl;
//   }
// };

// const double SAFETY = 0.9;
// const double MIN_FACTOR = 0.2;
// const double MAX_FACTOR = 10.0;

// // Function to take a single RK23 step
// void rk_step(const std::function<void(double, const VectorXd &, VectorXd &)> &fun, 
//              double t, const VectorXd &y, const VectorXd &f, double h, 
//              const Eigen::MatrixXd &A, const VectorXd &B, const VectorXd &C, 
//              std::vector<VectorXd> &K, VectorXd &y_new, VectorXd &f_new) {
    
//     K[0] = f;  // Store f(t, y) as the first stage
//     for (int s = 1; s < A.rows(); ++s) {
//         VectorXd dy = VectorXd::Zero(y.size());
//         for (int j = 0; j < s; ++j) {
//             dy += K[j] * A(s, j);
//         }
//         dy *= h;
//         fun(t + C[s] * h, y + dy, K[s]);
//     }

//     // Compute the new solution
//     y_new = y;
//     for (int i = 0; i < B.size(); ++i) {
//         y_new += h * B[i] * K[i];
//     }

//     // Evaluate f at the new time step
//     fun(t + h, y_new, f_new);
//     K.back() = f_new;  // Store the new derivative
// }
// VectorXd ode_solver(
//     const VectorXd &y0,
//     const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func) {

//   double t0 = 0.0;
//   double t_end = 1.0;
//   std::cout << "t0: " << t0 << std::endl;
//   // Create a copy of the initial condition
//   VectorXd y = y0;
//   //   std::cout << "y0: " << y.transpose() << std::endl;
//   // Define the RK23 stepper
//   // Instantiate the solver
//   RungeKutta solver;
//   std::cout << "Before solve, y: " << std::endl;
//   // Solve from t=0 to t=1
//   // std::cout << "Before solve, y: " << y.transpose() << std::endl;
//   solver.solve(ode_func, y);
//   // std::cout << "After solve, y: " << y.transpose() << std::endl;

//   return y;
// }




// // RK23 class
// class RungeKutta23 {
// public:
//     int order = 3;
//     int error_estimator_order = 2;
//     int n_stages = 3;

//     // Coefficients for the RK23 method
//     VectorXd C = (VectorXd(3) << 0, 0.5, 0.75).finished();
//     Eigen::MatrixXd A = (Eigen::MatrixXd(3, 3) << 0, 0, 0,
//                                                   0.5, 0, 0,
//                                                   0, 0.75, 0).finished();
//     VectorXd B = (VectorXd(3) << 2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0).finished();
//     VectorXd E = (VectorXd(4) << 5.0 / 72.0, -1.0 / 12.0, -1.0 / 9.0, 1.0 / 8.0).finished();

//     double h_abs = 1e-2;  // Initial step size
//     double rtol = 1e-3;
//     double atol = 1e-6;

//     RungeKutta23() {}

//     bool solve(const std::function<void(double, const VectorXd &, VectorXd &)> &fun, 
//                VectorXd &y) {
//         // Time variables
//         double t = 0.0;
//         double t_end = 1.0;

//         // Error scaling
//         double error_exponent = -1.0 / (error_estimator_order + 1);

//         // Derivative at the initial time
//         VectorXd f(y.size());
//         fun(t, y, f);

//         // Storage for RK stages
//         std::vector<VectorXd> K(n_stages + 1, VectorXd::Zero(y.size()));

//         bool step_accepted = false;
//         while (t < t_end) {
//             // Make sure step does not exceed the final time
//             double h = std::min(h_abs, t_end - t);
            
//             VectorXd y_new(y.size());
//             VectorXd f_new(y.size());

//              rk_step(fun, t, y, f, h, A, B, C, K, y_new, f_new);

//             // Error estimate and scaling
//             VectorXd scale = VectorXd::Constant(y.size(), atol) + VectorXd::Constant(y.size(), rtol).cwiseProduct(y_new.cwiseAbs().cwiseMax(y.cwiseAbs()));
//             VectorXd error = (y_new - y).cwiseQuotient(scale);
//             double error_norm = error.norm();

//             if (error_norm < 1.0) {
//                 // Step accepted, update variables
//                 t += h;
//                 y = y_new;
//                 f = f_new;

//                 // Step size adjustment
//                 double factor = std::min(MAX_FACTOR, SAFETY * std::pow(error_norm, error_exponent));
//                 h_abs *= factor;
//                 step_accepted = true;
//             } else {
//                 // Step rejected, reduce step size
//                 h_abs *= std::max(MIN_FACTOR, SAFETY * std::pow(error_norm, error_exponent));
//                 step_accepted = false;
//             }
//         }

//         return step_accepted;
//     }
// };

// Define a class for the Runge-Kutta-Fehlberg (RKF45) method
class RungeKuttaFehlberg45 {
public:
    double tol; // Tolerance for adaptive step-size control
    double h_min, h_max; // Min and max step sizes

    RungeKuttaFehlberg45(double tolerance = 1e-7, double min_step = 1e-10, double max_step = 0.1) 
        : tol(tolerance), h_min(min_step), h_max(max_step) {}

    void solve(const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func, 
               VectorXd &y, double t0 = 0.0, double t1 = 1.0) {
        double t = t0;
        double h = (t1 - t0) / 100; // Initial step size (can be adjusted)

        while (t < t1) {
            if (t + h > t1) h = t1 - t; // Adjust final step size to reach t1
            
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
                // std::cout << "t = " << t << ", h = " << h << ", y = " << y.transpose() << std::endl;
            }

            // Update step size for the next iteration
            h *= std::clamp(scale, 0.5, 2.0); // Clamp the scaling to avoid drastic step-size changes
            h = std::clamp(h, h_min, h_max); // Ensure step size stays within bounds
        }
    }

private:
    void rkf45_step(const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func, 
                    double t, const VectorXd& y, double h, 
                    VectorXd &y_new, VectorXd &error_estimate) {
        // RKF45 coefficients
        static const double a2 = 0.25, a3 = 3.0/8.0, a4 = 12.0/13.0, a5 = 1.0, a6 = 0.5;
        static const double b2 = 0.25;
        static const double b3[] = { 3.0/32.0, 9.0/32.0 };
        static const double b4[] = { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0 };
        static const double b5[] = { 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0 };
        static const double b6[] = { -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0 };

        static const double c[] = { 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 }; // 5th order solution
        static const double c_star[] = { 25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0 }; // 4th order solution

        // Compute the stages
        VectorXd k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()), k5(y.size()), k6(y.size());
        
        ode_func(t, y, k1);
        ode_func(t + a2 * h, y + b2 * k1 * h, k2);
        ode_func(t + a3 * h, y + b3[0] * k1 * h + b3[1] * k2 * h, k3);
        ode_func(t + a4 * h, y + b4[0] * k1 * h + b4[1] * k2 * h + b4[2] * k3 * h, k4);
        ode_func(t + a5 * h, y + b5[0] * k1 * h + b5[1] * k2 * h + b5[2] * k3 * h + b5[3] * k4 * h, k5);
        ode_func(t + a6 * h, y + b6[0] * k1 * h + b6[1] * k2 * h + b6[2] * k3 * h + b6[3] * k4 * h + b6[4] * k5 * h, k6);

        // Compute the 4th and 5th order estimates
        y_new = y + h * (c[0] * k1 + c[2] * k3 + c[3] * k4 + c[4] * k5 + c[5] * k6);
        VectorXd y_star = y + h * (c_star[0] * k1 + c_star[2] * k3 + c_star[3] * k4 + c_star[4] * k5);

        // Estimate the local error
        error_estimate = y_new - y_star;
    }
};

VectorXd ode_solver(const VectorXd &y0, const std::function<void(double, const VectorXd &, VectorXd &)> &ode_func) {

    double t0 = 0.0;
    double t_end = 1.0;
    // std::cout << "t0: " << t0 << std::endl;

    // Create a copy of the initial condition
    VectorXd y = y0;

    // Define the RK23 stepper
    RungeKuttaFehlberg45 solver;

    // std::cout << "Before solve, y: " << std::endl;
    // Solve from t=0 to t=1
    solver.solve(ode_func, y);

    return y;
}

