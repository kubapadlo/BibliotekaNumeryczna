#ifndef NUMLIBCPP_INTEGRATION_H
#define NUMLIBCPP_INTEGRATION_H

#include <functional>
#include <stdexcept> // Dla std::invalid_argument

namespace NumLibCpp {

/**
 * @brief Oblicza całkę oznaczoną funkcji `func` na przedziale `[a, b]` używając złożonej metody Simpsona.
 *
 * @param func Funkcja do całkowania, przyjmująca `double` i zwracająca `double`.
 * @param a Dolna granica całkowania.
 * @param b Górna granica całkowania.
 * @param n Liczba podprzedziałów (musi być parzysta i dodatnia).
 * @return double Przybliżona wartość całki.
 * @throws std::invalid_argument Jeśli `n` nie jest parzyste, `n <= 0`, lub `a > b`.
 *
 * @example
 * @code
 * #include <NumLibCpp/integration.h>
 * #include <iostream>
 * #include <functional>
 * #include <cmath> // Dla M_PI i std::sin
 *
 * double my_function(double x) { return x * x; } // f(x) = x^2
 *
 * int main() {
 *     try {
 *         // Całka z x^2 od 0 do 1 (wynik analityczny: 1/3)
 *         double integral_val = NumLibCpp::simpson_integrate(my_function, 0.0, 1.0, 100);
 *         std::cout << "Calka z x^2 od 0 do 1: " << integral_val << std::endl; // Oczekiwane: ~0.33333
 *
 *         // Całka z sin(x) od 0 do PI (wynik analityczny: 2)
 *         integral_val = NumLibCpp::simpson_integrate([](double x){ return std::sin(x); }, 0.0, M_PI, 100);
 *         std::cout << "Calka z sin(x) od 0 do PI: " << integral_val << std::endl; // Oczekiwane: ~2.0
 *     } catch (const std::exception& e) {
 *         std::cerr << "Blad: " << e.what() << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
double simpson_integrate(std::function<double(double)> func, double a, double b, int n);

} // namespace NumLibCpp

#endif //NUMLIBCPP_INTEGRATION_H