#ifndef NUMLIBCPP_NONLINEAR_SOLVER_H
#define NUMLIBCPP_NONLINEAR_SOLVER_H

#include <functional>
#include <stdexcept> // Dla std::runtime_error, std::invalid_argument
#include <limits>    // Dla std::numeric_limits

namespace NumLibCpp {

/**
 * @brief Znajduje pierwiastek równania nieliniowego f(x) = 0 metodą siecznych.
 *
 * Metoda siecznych jest iteracyjna i wymaga dwóch początkowych przybliżeń pierwiastka.
 * "Metoda krzywej linii" jest interpretowana jako metoda siecznych (Secant Method).
 *
 * @param func Funkcja f(x), której pierwiastka szukamy.
 * @param x0 Pierwsze przybliżenie początkowe.
 * @param x1 Drugie przybliżenie początkowe (różne od x0).
 * @param tol Tolerancja błędu (kryterium zbieżności, |x_new - x_old| < tol lub |f(x_new)| < tol).
 * @param max_iter Maksymalna liczba iteracji.
 * @return double Przybliżona wartość pierwiastka.
 * @throws std::invalid_argument Jeśli `tol <= 0` lub `max_iter <= 0`.
 * @throws std::runtime_error Jeśli metoda nie zbiegnie w `max_iter` iteracjach, lub jeśli `f(x1) - f(x0)` jest zbyt bliskie zeru (dzielenie przez zero).
 *
 * @example
 * @code
 * #include <NumLibCpp/nonlinear_solver.h>
 * #include <iostream>
 * #include <functional>
 * #include <cmath> // Dla std::cos
 *
 * // f(x) = x^2 - 2. Pierwiastki to sqrt(2) i -sqrt(2)
 * double func_quad(double x) { return x * x - 2.0; }
 *
 * int main() {
 *     try {
 *         double root = NumLibCpp::secant_method(func_quad, 1.0, 2.0, 1e-6, 100);
 *         std::cout << "Pierwiastek x^2 - 2 = 0 (metoda siecznych): " << root << std::endl; // Oczekiwane: ~1.41421
 *
 *         // f(x) = cos(x) - x. Pierwiastek ~0.739
 *         root = NumLibCpp::secant_method([](double x){ return std::cos(x) - x; }, 0.0, 1.0, 1e-7, 50);
 *         std::cout << "Pierwiastek cos(x) - x = 0: " << root << std::endl; // Oczekiwane: ~0.739085
 *     } catch (const std::exception& e) {
 *         std::cerr << "Blad: " << e.what() << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
double secant_method(std::function<double(double)> func, double x0, double x1, double tol, int max_iter);

} // namespace NumLibCpp

#endif //NUMLIBCPP_NONLINEAR_SOLVER_H