#ifndef NUMLIBCPP_DIFFERENTIATION_H
#define NUMLIBCPP_DIFFERENTIATION_H

#include <functional>
#include <stdexcept> // Dla std::invalid_argument

namespace NumLibCpp {

/**
 * @brief Oblicza numerycznie pochodną funkcji `func` w punkcie `x` używając różnicy centralnej.
 *
 * Formuła: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
 *
 * @param func Funkcja do zróżniczkowania, przyjmująca `double` i zwracająca `double`.
 * @param x Punkt, w którym obliczana jest pochodna.
 * @param h Mały krok (szerokość przedziału różniczkowania). Powinien być mały, ale nie za mały, aby uniknąć błędów numerycznych.
 * @return double Przybliżona wartość pochodnej.
 * @throws std::invalid_argument Jeśli `h <= 0`.
 *
 * @example
 * @code
 * #include <NumLibCpp/differentiation.h>
 * #include <iostream>
 * #include <functional>
 * #include <cmath> // Dla std::sin, std::cos
 *
 * // f(x) = x^3, f'(x) = 3x^2
 * double cubic_func(double x) { return x * x * x; }
 *
 * int main() {
 *     try {
 *         double deriv_at_2 = NumLibCpp::central_difference(cubic_func, 2.0, 1e-5);
 *         std::cout << "Pochodna x^3 w x=2: " << deriv_at_2 << std::endl; // Oczekiwane: ~12.0
 *
 *         // f(x) = sin(x), f'(x) = cos(x). f'(0) = cos(0) = 1
 *         deriv_at_2 = NumLibCpp::central_difference([](double x){ return std::sin(x);}, 0.0, 1e-5);
 *         std::cout << "Pochodna sin(x) w x=0: " << deriv_at_2 << std::endl; // Oczekiwane: ~1.0
 *     } catch (const std::exception& e) {
 *         std::cerr << "Blad: " << e.what() << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
double central_difference(std::function<double(double)> func, double x, double h);

} // namespace NumLibCpp

#endif //NUMLIBCPP_DIFFERENTIATION_H