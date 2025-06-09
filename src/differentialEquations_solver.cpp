#include "NumLibCpp/differentialEquations_solver.hpp"

namespace NumLibCpp {

double rk4_solve(std::function<double(double, double)> f, double x0, double y0, double x_target, int num_steps) {
    if (num_steps <= 0) {
        throw std::invalid_argument("Liczba krokow musi byc dodatnia.");
    }

    double h = (x_target - x0) / num_steps;
    double x = x0;
    double y = y0;

    for (int i = 0; i < num_steps; ++i) {
        double k1 = h * f(x, y);
        double k2 = h * f(x + 0.5 * h, y + 0.5 * k1);
        double k3 = h * f(x + 0.5 * h, y + 0.5 * k2);
        double k4 = h * f(x + h, y + k3);

        y = y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        x = x + h; // lub x = x0 + (i+1)*h dla większej precyzji
    }
    return y;
}

} // namespace NumLibCpp