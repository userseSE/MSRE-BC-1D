#include "data_saving.hpp"
#include "parameters.hpp"

#include <filesystem>  // C++17
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <string>

namespace fs = std::filesystem;

// Helper function to export a 1D array to CSV
void export_to_csv(const double* data, int size, const std::string& folder, const std::string& filename) {
    std::cout << "Exporting to CSV file: " << folder << "/" << filename << std::endl;
    fs::create_directories(folder);  // Create folder if it doesn't exist
    std::ofstream file(folder + "/" + filename);
    for (int i = 0; i < size; ++i) {
        file << data[i] << "\n";
    }
    file.close();
}

// Helper function to export a 2D array to CSV
void export_matrix_to_csv(const double* data, int rows, int cols, const std::string& folder, const std::string& filename) {
    fs::create_directories(folder);  // Create folder if it doesn't exist
    std::ofstream file(folder + "/" + filename);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file << data[i * cols + j];  // Access 2D array stored in 1D
            if (j < cols - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
}

void save_results(const double rho_matrix[time_span], 
                  const double phi1_middle_matrix[time_span],
                  const double phi2_middle_matrix[time_span],
                  const double ci_middle_matrix[time_span][6], 
                  const double temperature_fuel_middle_matrix[time_span],
                  const std::string& folder) {
    // Export each data vector to a CSV file in the specific folder
    export_to_csv(rho_matrix, time_span, folder, "rho_matrix.csv");

    // Combine phi1_middle_matrix and phi2_middle_matrix into one matrix for export
    double phi_middle[time_span][2];
    for (int i = 0; i < time_span; ++i) {
        phi_middle[i][0] = phi1_middle_matrix[i];
        phi_middle[i][1] = phi2_middle_matrix[i];
    }
    export_matrix_to_csv(&phi_middle[0][0], time_span, 2, folder, "phi_middle_matrix.csv");

    // Export ci_middle_matrix
    export_matrix_to_csv(&ci_middle_matrix[0][0], time_span, 6, folder, "ci_middle_matrix.csv");

    // Export temperature_fuel_middle_matrix
    export_to_csv(temperature_fuel_middle_matrix, time_span, folder, "temperature_fuel_middle_matrix.csv");

    std::cout << "Results saved to CSV files in " << folder << "." << std::endl;
}

void save_spacial_results(const double rho[N], 
                          const double y_n[length_neutr], 
                          const double y_th[length_th], 
                          const double y_hx1[length_hx], 
                          const double y_hx2[length_hx], 
                          const std::string& folder) {
    // Export the 1D arrays to CSV
    export_to_csv(rho, N, folder, "rho.csv");
    export_to_csv(y_n, length_neutr, folder, "y_n.csv");
    export_to_csv(y_th, length_th, folder, "y_th.csv");
    export_to_csv(y_hx1, length_hx, folder, "y_hx1.csv");
    export_to_csv(y_hx2, length_hx, folder, "y_hx2.csv");

    std::cout << "Spatial results saved to CSV files in " << folder << "." << std::endl;
}
