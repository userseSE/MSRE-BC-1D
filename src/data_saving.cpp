#include "data_saving.hpp"
#include "parameters.hpp"

#include <filesystem> // C++17
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

namespace fs = std::filesystem;

void export_to_csv(const std::vector<double>& data, const std::string& filename) {
    // std::cout << "Exporting to CSV file: " << filename << std::endl;
    fs::create_directories("dataset");  // Create "dataset" folder if it doesn't exist
    std::ofstream file("../dataset/" + filename);
    for (const auto& value : data) {
        file << value << "\n";
    }
    file.close();
}

void export_matrix_to_csv(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    fs::create_directories("dataset");  // Create "dataset" folder if it doesn't exist
    std::ofstream file("../dataset/" + filename);
    for (const auto& row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
    // std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;
}

void save_results(const std::vector<double>& rho_matrix, 
                  const std::vector<double>& phi_middle_matrix,
                  const std::vector<std::vector<double>>& ci_middle_matrix, 
                  const std::vector<double>& temperature_fuel_middle_matrix) {
    // Export each data vector to a CSV file
    export_to_csv(rho_matrix, "rho_matrix.csv");
    export_to_csv(phi_middle_matrix, "phi_middle_matrix.csv");
    export_matrix_to_csv(ci_middle_matrix, "ci_middle_matrix.csv");
    export_to_csv(temperature_fuel_middle_matrix, "temperature_fuel_middle_matrix.csv");
    // std::cout << "Results saved to CSV files." << std::endl;
}

void save_spacial_results(const std::vector<double>& phi, 
                          const std::vector<double>& ci, 
                          const std::vector<double>& temperature_fuel, 
                          const std::vector<double>& temperature_graphite, 
                          const std::vector<double>& Ts_HX1, 
                          const std::vector<double>& Tss_HX1, 
                          const std::vector<double>& Tss_HX2, 
                          const std::vector<double>& Tsss_HX2) {
    // Export each spatial data vector to a CSV file
    export_to_csv(phi, "phi.csv");
    std::vector<std::vector<double>> segmented_ci(6, std::vector<double>(200));
    for (int i = 0; i < 6; ++i) {
        std::copy(ci.begin() + i * 200, ci.begin() + (i + 1) * 200, segmented_ci[i].begin());
    }
    std::vector<std::vector<double>> transposed_ci(200, std::vector<double>(6));
    for (size_t i = 0; i < segmented_ci.size(); ++i) {
        for (size_t j = 0; j < segmented_ci[0].size(); ++j) {
            transposed_ci[j][i] = segmented_ci[i][j];
        }
    }
    export_matrix_to_csv(transposed_ci, "ci.csv");

    std::vector<std::vector<double>> temperature_core = {temperature_fuel, temperature_graphite};
    std::vector<std::vector<double>> transposed_temperature_core(200, std::vector<double>(2));
    for (size_t i = 0; i < temperature_core.size(); ++i) {
        for (size_t j = 0; j < temperature_core[0].size(); ++j) {
            transposed_temperature_core[j][i] = temperature_core[i][j];
        }
    }
    export_matrix_to_csv(transposed_temperature_core, "temperature_core.csv");

    std::vector<std::vector<double>> temperature_HX1 = {Ts_HX1, Tss_HX1};
    std::vector<std::vector<double>> transposed_temperature_HX1(200, std::vector<double>(2));
    for (size_t i = 0; i < temperature_HX1.size(); ++i) {
        for (size_t j = 0; j < temperature_HX1[0].size(); ++j) {
            transposed_temperature_HX1[j][i] = temperature_HX1[i][j];
        }
    }
    export_matrix_to_csv(transposed_temperature_HX1, "temperature_HX1.csv");

    std::vector<std::vector<double>> temperature_HX2 = {Tss_HX2, Tsss_HX2};
    std::vector<std::vector<double>> transposed_temperature_HX2(200, std::vector<double>(2));
    for (size_t i = 0; i < temperature_HX2.size(); ++i) {
        for (size_t j = 0; j < temperature_HX2[0].size(); ++j) {
            transposed_temperature_HX2[j][i] = temperature_HX2[i][j];
        }
    }
    export_matrix_to_csv(transposed_temperature_HX2, "temperature_HX2.csv");
}
