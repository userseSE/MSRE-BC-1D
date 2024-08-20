#include "thermal_hydraulics.hpp"
#include "parameters.hpp"
#include "ode_solver.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// Include necessary linear algebra headers
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;

const double a_th = Vc;
const double b_th = U / (Ms * c_p_s);
const double c_th = U / (Mg * c_p_g);
const double d_th = L * gamma / (Ms * c_p_s);
const double e_th = L * (1 - gamma) / (Mg * c_p_g);

// Discretize the spatial domain
SparseMatrix<double> AT(N, N);

void initialize_AT() {
    AT = SparseMatrix<double>(N, N);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(3 * N - 2);
    
    for (int i = 0; i < N; ++i) {
        tripletList.emplace_back(i, i, -1.0);
        if (i < N - 1) {
            tripletList.emplace_back(i, i + 1, 1.0 / dz);
        }
    }
    
    AT.setFromTriplets(tripletList.begin(), tripletList.end());
    
    AT.coeffRef(0, 0) = 1.0 / dz;
    AT.coeffRef(N - 1, N - 1) = -1.0 / dz;
}

std::vector<double> thermal_hydraulics(std::vector<double>& y_th, 
                                       const std::vector<double>& q_prime, 
                                       double Ts_core_0, 
                                       int step) {
    // Ensure AT is initialized
    static bool initialized = false;
    if (!initialized) {
        initialize_AT();
        initialized = true;
    }

    // Set boundary conditions
    if (step == 0) {
        y_th[0] = Ts_core_0;
    } else {
        y_th[0] = Ts_core_0;
    }

    // Define the ODE system for the thermal hydraulics
    auto pde_to_ode_th = [&](double t, const std::vector<double>& y) {
       
        VectorXd temperature_fuel = Eigen::Map<const VectorXd>(y.data(), N);
        VectorXd temperature_graphite = Eigen::Map<const VectorXd>(y.data() + N, N);
        
        // std::cout<<(AT * temperature_fuel).size()<<std::endl;
        // std::cout<<VectorXd::Map(q_prime.data(), N).size()<<std::endl;
        // std::cout<<"test pde_to_ode_th"<<std::endl;
        VectorXd temperature_fuel_dt = a_th * (AT * temperature_fuel) + b_th * (temperature_graphite - temperature_fuel) + d_th * VectorXd::Map(q_prime.data(), N);
        // std::cout << "Matrix AT size: " << AT.rows() << "x" << AT.cols() << std::endl;
        // std::cout << "Vector temperature_fuel size: " << temperature_fuel.size() << std::endl;
        VectorXd temperature_graphite_dt = c_th * (temperature_fuel - temperature_graphite) + e_th * VectorXd::Map(q_prime.data(), N);
        // std::cout<<"test pde_to_ode_th"<<std::endl;
        // Apply time-varying boundary conditions
        temperature_fuel_dt[0] = bc_s0 - temperature_fuel[0];
        temperature_fuel_dt[N - 1] = bc_sL - temperature_fuel[N - 1];
        temperature_graphite_dt[0] = bc_g0 - temperature_graphite[0];
        temperature_graphite_dt[N - 1] = bc_gL - temperature_graphite[N - 1];

        std::vector<double> dydt(2 * N);
        std::copy(temperature_fuel_dt.data(), temperature_fuel_dt.data() + N, dydt.begin());
        std::copy(temperature_graphite_dt.data(), temperature_graphite_dt.data() + N, dydt.begin() + N);
        
        return dydt;
    };

    // Initial condition vector
    std::vector<double> y0(2 * N);
    if (step == 0) {
        std::copy(initialS, initialS + N, y0.begin());
        std::copy(initialG, initialG + N, y0.begin() + N);
    } else {
        y0 = y_th;
    }
    // std::cout<<"size of y0 in th"<<y0.size()<<std::endl;
    // Solve the ODE system
    std::vector<double> solution_y_th = ode_solver(y0, {}, pde_to_ode_th);
    // std::cout<<"size of y in th"<<solution_y_th.size()<<std::endl;
    // Return the solution at the last time step
    return solution_y_th;
}
