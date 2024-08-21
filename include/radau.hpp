# include <vector>
# include <functional>

using ODEFunction = std::function<void(double, const std::vector<double>&, std::vector<double>&)>;

std::vector<std::vector<double>> radau_solver(ODEFunction f, double t0, double tf, std::vector<double> y0, double dt);
