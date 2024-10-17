#include "data_saving.hpp"
#include "parameters.hpp"

#include <filesystem>  // C++17
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <string>

namespace fs = std::filesystem;

void export_to_csv(const Eigen::VectorXd& data, const std::string& folder, const std::string& filename) {
    std::cout << "Exporting to CSV file: " << folder << "/" << filename << std::endl;
    fs::create_directories(folder);  // Create folder if it doesn't exist
    std::ofstream file(folder + "/" + filename);
    for (int i = 0; i < data.size(); ++i) {
        file << data[i] << "\n";
    }
    file.close();
}

void export_matrix_to_csv(const Eigen::MatrixXd& matrix, const std::string& folder, const std::string& filename) {
    fs::create_directories(folder);  // Create folder if it doesn't exist
    std::ofstream file(folder + "/" + filename);
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
                  const Eigen::VectorXd& phi1_middle_matrix,
                  const Eigen::VectorXd& phi2_middle_matrix,
                  const Eigen::MatrixXd& ci_middle_matrix, 
                  const Eigen::VectorXd& temperature_fuel_middle_matrix,
                  const std::string& folder) {
    // Export each data vector to a CSV file in the specific folder
    export_to_csv(rho_matrix, folder, "rho_matrix.csv");
    Eigen::MatrixXd phi_middle(phi1_middle_matrix.size(), 2);
    phi_middle.col(0) = phi1_middle_matrix;
    phi_middle.col(1) = phi2_middle_matrix;
    export_matrix_to_csv(phi_middle, folder, "phi_middle_matrix.csv");
    export_matrix_to_csv(ci_middle_matrix, folder, "ci_middle_matrix.csv");
    export_to_csv(temperature_fuel_middle_matrix, folder, "temperature_fuel_middle_matrix.csv");
    std::cout << "Results saved to CSV files in " << folder << "." << std::endl;
}

void save_spacial_results(const Eigen::VectorXd& phi1, 
                          const Eigen::VectorXd& phi2, 
                          const Eigen::VectorXd& ci, 
                          const Eigen::VectorXd& rho,
                          const Eigen::VectorXd& temperature_fuel, 
                          const Eigen::VectorXd& temperature_graphite, 
                          const Eigen::VectorXd& Ts_HX1, 
                          const Eigen::VectorXd& Tss_HX1, 
                          const Eigen::VectorXd& Tss_HX2, 
                          const Eigen::VectorXd& Tsss_HX2,
                          const std::string& folder) {
    // Export each spatial data vector to a CSV file in the specific folder
    Eigen::MatrixXd phi(phi1.size(), 2);
    phi.col(0) = phi1;
    phi.col(1) = phi2;
    export_matrix_to_csv(phi, folder, "phi.csv");

    export_to_csv(rho, folder, "rho.csv");

    // Reshape and transpose ci for export
    Eigen::MatrixXd segmented_ci = Eigen::Map<const Eigen::MatrixXd>(ci.data(), ci.size() / 6, 6);
    export_matrix_to_csv(segmented_ci, folder, "ci.csv");

    // Combine temperature_fuel and temperature_graphite into a matrix for export
    Eigen::MatrixXd temperature_core(temperature_fuel.size(), 2);
    temperature_core.col(0) = temperature_fuel;
    temperature_core.col(1) = temperature_graphite;
    export_matrix_to_csv(temperature_core, folder, "temperature_core.csv");

    // Combine Ts_HX1 and Tss_HX1 into a matrix for export
    Eigen::MatrixXd temperature_HX1(Ts_HX1.size(), 2);
    temperature_HX1.col(0) = Ts_HX1;
    temperature_HX1.col(1) = Tss_HX1;
    export_matrix_to_csv(temperature_HX1, folder, "temperature_HX1.csv");

    // Combine Tss_HX2 and Tsss_HX2 into a matrix for export
    Eigen::MatrixXd temperature_HX2(Tss_HX2.size(), 2);
    temperature_HX2.col(0) = Tss_HX2;
    temperature_HX2.col(1) = Tsss_HX2;
    export_matrix_to_csv(temperature_HX2, folder, "temperature_HX2.csv");
    std::cout << "Spatial results saved to CSV files in " << folder << "." << std::endl;
}
