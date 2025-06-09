#ifndef NUMLIBCPP_ODE_SOLVER_H
#define NUMLIBCPP_ODE_SOLVER_H

#include <functional>
#include <vector>
#include <stdexcept> // Dla std::invalid_argument

namespace NumLibCpp {

/**
 * @brief Rozwiązuje równanie różniczkowe zwyczajne pierwszego rzędu y' = f(x, y)
 *        metodą Rungego-Kutty 4. rzędu (RK4).
 *
 * @param f Funkcja f(x, y) definiująca równanie różniczkowe.
 * @param x0 Początkowa wartość x.
 * @param y0 Początkowa wartość y (warunek początkowy y(x0) = y0).
 * @param x_target Wartość x, dla której szukane jest rozwiązanie y.
 * @param num_steps Liczba kroków całkowania (musi być dodatnia).
 * @return double Wartość y w punkcie x_target.
 * @throws std::invalid_argument Jeśli `num_steps <= 0`.
 *
 * @example
 * @code
 * #include <NumLibCpp/ode_solver.h>
 * #include <iostream>
 * #include <functional>
 * #include <cmath> // Dla std::exp
 *
 * // Równanie y' = y, y(0) = 1. Rozwiązanie analityczne: y(x) = exp(x)
 * double f_exp(double x, double y) {
 *     (void)x; // x nie jest używane w tym konkretnym f
 *     return y;
 * }
 *
 * int main() {
 *     try {
 *         double y_at_1 = NumLibCpp::rk4_solve(f_exp, 0.0, 1.0, 1.0, 100);
 *         std::cout << "RK4 dla y'=y, y(0)=1: y(1) = " << y_at_1 << std::endl; // Oczekiwane: ~exp(1) = 2.718
 *     } catch (const std::exception& e) {
 *         std::cerr << "Blad: " << e.what() << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
double rk4_solve(std::function<double(double, double)> f, double x0, double y0, double x_target, int num_steps);

} // namespace NumLibCpp

#endif //NUMLIBCPP_ODE_SOLVER_H