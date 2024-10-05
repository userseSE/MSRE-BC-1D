#include "data_saving.hpp"
#include "parameters.hpp"

#include <filesystem>  // C++17
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <string>

namespace fs = std::filesystem;

void export_to_csv(const Eigen::VectorXd& data, const std::string& filename) {
    std::cout << "Exporting to CSV file: " << filename << std::endl;
    fs::create_directories("dataset");  // Create "dataset" folder if it doesn't exist
    std::ofstream file("../dataset/" + filename);
    for (int i = 0; i < data.size(); ++i) {
        file << data[i] << "\n";
    }
    file.close();
}

void export_matrix_to_csv(const Eigen::MatrixXd& matrix, const std::string& filename) {
    fs::create_directories("dataset");  // Create "dataset" folder if it doesn't exist
    std::ofstream file("../dataset/" + filename);
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            file << matrix(i, j);
            if (j < matrix.cols() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
}

void save_results(const Eigen::VectorXd& rho_matrix, 
                  const Eigen::VectorXd& phi_middle_matrix,
                  const Eigen::MatrixXd& ci_middle_matrix, 
                  const Eigen::VectorXd& temperature_fuel_middle_matrix) {
    // Export each data vector to a CSV file
    export_to_csv(rho_matrix, "rho_matrix.csv");
    export_to_csv(phi_middle_matrix, "phi_middle_matrix.csv");
    export_matrix_to_csv(ci_middle_matrix, "ci_middle_matrix.csv");
    export_to_csv(temperature_fuel_middle_matrix, "temperature_fuel_middle_matrix.csv");
    std::cout << "Results saved to CSV files." << std::endl;
}

void save_spacial_results(const Eigen::VectorXd& phi, 
                          const Eigen::VectorXd& ci, 
                          const Eigen::VectorXd& temperature_fuel, 
                          const Eigen::VectorXd& temperature_graphite, 
                          const Eigen::VectorXd& Ts_HX1, 
                          const Eigen::VectorXd& Tss_HX1, 
                          const Eigen::VectorXd& Tss_HX2, 
                          const Eigen::VectorXd& Tsss_HX2) {
    // Export each spatial data vector to a CSV file
    export_to_csv(phi, "phi.csv");

    // Reshape and transpose ci for export
    Eigen::MatrixXd segmented_ci = Eigen::Map<const Eigen::MatrixXd>(ci.data(), N, 6);
    // Eigen::MatrixXd transposed_ci = segmented_ci.transpose();
    export_matrix_to_csv(segmented_ci, "ci.csv");

    // Combine temperature_fuel and temperature_graphite into a matrix for export
    Eigen::MatrixXd temperature_core(N, 2);
    temperature_core.col(0) = temperature_fuel;
    temperature_core.col(1) = temperature_graphite;
    export_matrix_to_csv(temperature_core, "temperature_core.csv");

    // Combine Ts_HX1 and Tss_HX1 into a matrix for export
    Eigen::MatrixXd temperature_HX1(Nx, 2);
    temperature_HX1.col(0) = Ts_HX1;
    temperature_HX1.col(1) = Tss_HX1;
    export_matrix_to_csv(temperature_HX1, "temperature_HX1.csv");

    // Combine Tss_HX2 and Tsss_HX2 into a matrix for export
    Eigen::MatrixXd temperature_HX2(Nx, 2);
    temperature_HX2.col(0) = Tss_HX2;
    temperature_HX2.col(1) = Tsss_HX2;
    export_matrix_to_csv(temperature_HX2, "temperature_HX2.csv");
}
