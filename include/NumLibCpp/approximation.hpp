#ifndef NUMLIBCPP_APPROXIMATION_H
#define NUMLIBCPP_APPROXIMATION_H

#include <vector>
#include <utility>   // Dla std::pair
#include <functional> // Dla std::function
#include <stdexcept> // Dla std::invalid_argument, std::runtime_error

namespace NumLibCpp {
    
/**
 * @brief Wykonuje aproksymację funkcji f(x) wielomianem stopnia 'degree' na przedziale [a, b]
 *        metodą najmniejszych kwadratów w sensie ciągłym.
 *
 * Funkcja znajduje współczynniki c_0, c_1, ..., c_degree wielomianu
 * P(x) = c_0 + c_1*x + ... + c_degree*x^degree, który minimalizuje całkę
 * ∫[a,b] (f(x) - P(x))^2 dx.
 * Współczynniki są rozwiązaniem układu równań Ac = d, gdzie:
 * A_ij = ∫[a,b] x^(i+j) dx
 * d_i  = ∫[a,b] f(x)*x^i dx
 * Całki są obliczane numerycznie metodą Simpsona.
 * Układ równań jest rozwiązywany metodą eliminacji Gaussa.
 *
 * @param func_to_approx Funkcja f(x) do aproksymowania.
 * @param a Dolna granica przedziału aproksymacji.
 * @param b Górna granica przedziału aproksymacji.
 * @param degree Stopień wielomianu aproksymującego (musi być >= 0).
 * @param num_simpson_intervals Liczba podprzedziałów dla metody Simpsona (musi być parzysta i dodatnia).
 * @return std::vector<double> Wektor współczynników wielomianu [c_0, c_1, ..., c_degree].
 * @throws std::invalid_argument Jeśli `degree < 0`, `a >= b`, lub `num_simpson_intervals` jest niepoprawne.
 * @throws std::runtime_error Jeśli rozwiązanie układu równań napotka problemy (np. macierz osobliwa).
 *
 * @example
 * @code
 * #include <NumLibCpp/approximation.h>
 * #include <NumLibCpp/integration.h> // Potrzebne dla Simpsona, ale tu jest używany wewnętrznie
 * #include <NumLibCpp/linear_algebra.h> // Potrzebne dla Gaussa, ale tu jest używany wewnętrznie
 * #include <iostream>
 * #include <vector>
 * #include <cmath>
 * #include <iomanip>
 *
 * double my_func(double x) { return std::exp(x) * std::cos(5 * x) - std::pow(x, 3); }
 * double poly_eval(double x, const std::vector<double>& coeffs) {
 *     double res = 0.0;
 *     for (size_t i = 0; i < coeffs.size(); ++i) res += coeffs[i] * std::pow(x, i);
 *     return res;
 * }
 *
 * int main() {
 *     double a = -1.0, b_val = 2.0; // b_val zamiast b zeby nie kolidowac z para (a,b) z linear_least_squares
 *     int degree = 3;
 *     int simpson_intervals = 100;
 *     std::cout << std::fixed << std::setprecision(6);
 *
 *     try {
 *         std::vector<double> coeffs = NumLibCpp::polynomial_approximation(my_func, a, b_val, degree, simpson_intervals);
 *         std::cout << "Wspolczynniki wielomianu aproksymujacego stopnia " << degree << ":" << std::endl;
 *         for (size_t i = 0; i < coeffs.size(); ++i) {
 *             std::cout << "c" << i << " = " << coeffs[i] << std::endl;
 *         }
 *
 *         // Test w kilku punktach
 *         std::cout << "\nPorownanie w kilku punktach (x, f(x), P(x)):" << std::endl;
 *         for (double x_test = a; x_test <= b_val; x_test += (b_val - a) / 4.0) {
 *             std::cout << x_test << "\t" << my_func(x_test) << "\t" << poly_eval(x_test, coeffs) << std::endl;
 *         }
 *     } catch (const std::exception& e) {
 *         std::cerr << "Blad aproksymacji wielomianowej: " << e.what() << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
std::vector<double> polynomial_approximation(
    std::function<double(double)> func_to_approx,
    double a,
    double b,
    int degree,
    int num_simpson_intervals
);

} // namespace NumLibCpp

#endif //NUMLIBCPP_APPROXIMATION_H
