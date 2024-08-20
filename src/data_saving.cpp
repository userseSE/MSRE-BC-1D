#include "data_saving.hpp"
#include "parameters.hpp"

#include <filesystem> // C++17
#include <fstream>
#include <vector>
#include <string>

namespace fs = std::filesystem;

void export_to_csv(const std::vector<double>& data, const std::string& filename) {
    fs::create_directories("dataset");  // Create "dataset" folder if it doesn't exist
    std::ofstream file("dataset/" + filename);
    for (const auto& value : data) {
        file << value << "\n";
    }
    file.close();
}

void export_matrix_to_csv(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    fs::create_directories("dataset");  // Create "dataset" folder if it doesn't exist
    std::ofstream file("dataset/" + filename);
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
    export_to_csv(ci, "ci.csv");
    export_to_csv(temperature_fuel, "temperature_fuel.csv");
    export_to_csv(temperature_graphite, "temperature_graphite.csv");
    export_to_csv(Ts_HX1, "Ts_HX1.csv");
    export_to_csv(Tss_HX1, "Tss_HX1.csv");
    export_to_csv(Tss_HX2, "Tss_HX2.csv");
    export_to_csv(Tsss_HX2, "Tsss_HX2.csv");
}
