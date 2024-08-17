#include "plotting.hpp"
#include <matplotlibcpp.h> // Assuming matplotlibcpp is used for plotting
#include "parameters.hpp"

namespace plt = matplotlibcpp;

void plot_results(const std::vector<double>& rho_matrix, 
                  const std::vector<double>& phi_middle_matrix,
                  const std::vector<std::vector<double>>& ci_middle_matrix, 
                  const std::vector<double>& temperature_fuel_middle_matrix) {
    
    std::vector<double> z(N);
    for (int i = 0; i < N; ++i) {
        z[i] = static_cast<double>(i) * L / (N - 1);
    }

    plt::figure_size(1400, 600);

    // Plot 1: Neutron Flux
    plt::subplot(2, 2, 1);
    plt::plot(z, phi_middle_matrix);
    plt::title("Neutron Flux");

    // Plot 2: Delayed Neutron Precursors
    plt::subplot(2, 2, 2);
    for (int i = 0; i < 6; ++i) {
        std::vector<double> ci(N);
        for (int j = 0; j < N; ++j) {
            ci[j] = ci_middle_matrix[j][i];
        }
        plt::plot(z, ci);
    }
    plt::title("Delayed Neutron Precursors");

    // Plot 3: Reactivity
    plt::subplot(2, 2, 3);
    plt::plot(rho_matrix);
    plt::title("Reactivity");

    // Plot 4: Temperature in the Core
    plt::subplot(2, 2, 4);
    plt::plot(z, temperature_fuel_middle_matrix);
    plt::title("Temperature in the Core");

    plt::tight_layout();
    plt::save("neutron_flux_and_precursors.png");

    // Additional plots can be created in the same manner...
}
