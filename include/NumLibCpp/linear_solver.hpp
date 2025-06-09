#ifndef NUMLIBCPP_LINEAR_SOLVER_HPP
#define NUMLIBCPP_LINEAR_SOLVER_HPP

#include <vector>
#include <stdexcept> // Dla std::runtime_error, std::invalid_argument

namespace NumLibCpp {

/**
 * @brief Rozwiązuje układ równań liniowych Ax = b metodą eliminacji Gaussa-Jordana.
 *
 * Funkcja modyfikuje kopie macierzy A i wektora b, aby doprowadzić macierz A do postaci jednostkowej
 * (pełna eliminacja Gaussa-Jordana), a wynikowy wektor b zawiera rozwiązanie x.
 *
 * @param A Kwadratowa macierz współczynników (NxN).
 * @param b Wektor wyrazów wolnych (N).
 * @return Wektor x będący rozwiązaniem układu równań.
 * @throws std::invalid_argument Jeśli macierz A nie jest kwadratowa lub rozmiary są niezgodne.
 * @throws std::runtime_error Jeśli macierz A jest osobliwa (nieodwracalna).
 */
std::vector<double> gauss_elimination(std::vector<std::vector<double>> A, std::vector<double> b);

} // namespace NumLibCpp

#endif // NUMLIBCPP_LINEAR_SOLVER_HPP
