#ifndef NUMLIBCPP_INTERPOLATION_H
#define NUMLIBCPP_INTERPOLATION_H

#include <vector>
#include <stdexcept> // Dla std::invalid_argument

namespace NumLibCpp {

/**
 * @brief Oblicza wartość interpolowaną w punkcie x_interp używając wielomianu Lagrange'a.
 *
 * @param x_nodes Wektor współrzędnych x znanych punktów (węzłów). Muszą być unikalne.
 * @param y_nodes Wektor współrzędnych y znanych punktów (węzłów).
 * @param x_interp Punkt, w którym ma być obliczona wartość interpolowana.
 * @return double Wartość interpolowana y w punkcie x_interp.
 * @throws std::invalid_argument Jeśli `x_nodes` i `y_nodes` mają różne rozmiary, są puste, lub `x_nodes` zawierają powtarzające się wartości.
 *
 * @example
 * @code
 * #include <NumLibCpp/interpolation.h>
 * #include <iostream>
 * #include <vector>
 *
 * int main() {
 *     std::vector<double> x_nodes = {0.0, 1.0, 2.0};
 *     std::vector<double> y_nodes = {0.0, 1.0, 4.0}; // y = x^2
 *     try {
 *         double y_interp = NumLibCpp::lagrange_interpolate(x_nodes, y_nodes, 1.5);
 *         std::cout << "Interpolowana wartosc w x=1.5: " << y_interp << std::endl; // Oczekiwane: 2.25
 *     } catch (const std::exception& e) {
 *         std::cerr << "Blad: " << e.what() << std::endl;
 *     }
 *     return 0;
 * }
 * @endcode
 */
double lagrange_interpolate(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double x_interp);

} // namespace NumLibCpp

#endif //NUMLIBCPP_INTERPOLATION_H