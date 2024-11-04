#include "HX1.hpp"
#include "HX2.hpp"
#include "neutronics.hpp"
#include "parameters.hpp"
#include "reactivity.hpp"
#include "thermal_hydraulics.hpp"

#include <iostream>

int main() {
  Param_Neutronics params_neutr;
  // Initialization
  float rho[N];
  float y_n[length_neutr];
  float y_th[length_th];
  float y_hx1[length_hx];
  float y_hx2[length_hx];

  for (int step = 0; step < time_span; ++step) {
    std::cout << "Time step: " << step << std::endl;
    neutronics(y_n, rho, step, params_neutr);
    thermal_hydraulics(y_th, y_n,y_hx1[0], step);
    HX1(y_hx1, y_th[N - 1], step);
    HX2(y_hx2, y_hx1[Nx - 1], step);
    reactivity(y_th, y_th + N, step, rho);
  }
  return 0;
}
