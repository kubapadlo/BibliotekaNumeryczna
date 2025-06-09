#include "NumLibCpp/approximation.hpp"
#include "NumLibCpp/integration.hpp"    // Dla simpson_integrate
#include "NumLibCpp/linear_solver.hpp" // Dla gauss_elimination
#include <numeric> // Dla std::accumulate
#include <cmath>   // Dla std::pow
#include <limits>  // Dla std::numeric_limits

namespace NumLibCpp {

std::vector<double> polynomial_approximation(
    std::function<double(double)> func_to_approx,
    double a,
    double b,
    int degree,
    int num_simpson_intervals) {

    if (degree < 0) {
        throw std::invalid_argument("Stopien wielomianu (degree) musi byc nieujemny.");
    }
    if (a >= b) { // simpson_integrate obsłuży a > b, ale a == b da 0, co może być problemem dla macierzy
        throw std::invalid_argument("Dolna granica calkowania 'a' musi byc mniejsza niz gorna granica 'b'.");
    }
    // Sprawdzenia dla num_simpson_intervals są wewnątrz simpson_integrate

    int matrix_size = degree + 1;
    std::vector<std::vector<double>> A_matrix(matrix_size, std::vector<double>(matrix_size));
    std::vector<double> d_vector(matrix_size);

    // Budowanie macierzy A_matrix
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            int power = i + j;
            auto integrand_A = [power](double x_val) {
                if (power == 0) return 1.0; // x^0 = 1
                return std::pow(x_val, static_cast<double>(power));
            };
            // simpson_integrate rzuci wyjątkiem, jeśli num_simpson_intervals jest niepoprawne
            A_matrix[i][j] = simpson_integrate(integrand_A, a, b, num_simpson_intervals);
        }
    }

    // Budowanie wektora d_vector
    for (int i = 0; i <= degree; ++i) {
        int power = i;
        auto integrand_d = [&func_to_approx, power](double x_val) {
            if (power == 0) return func_to_approx(x_val); // f(x) * x^0
            return func_to_approx(x_val) * std::pow(x_val, static_cast<double>(power));
        };
        d_vector[i] = simpson_integrate(integrand_d, a, b, num_simpson_intervals);
    }

    // Rozwiązanie układu Ac = d
    // Funkcja gauss_elimination przyjmuje przez wartość, więc tworzy kopie
    try {
        return gauss_elimination(A_matrix, d_vector);
    } catch (const std::runtime_error& e) {
        // Przechwycenie błędu z Gaussa (np. macierz osobliwa) i rzucenie dalej
        // z bardziej kontekstowym komunikatem lub po prostu rzucenie dalej
        throw std::runtime_error(std::string("Blad podczas rozwiazywania ukladu rownan dla aproksymacji: ") + e.what());
    }
}

} // namespace NumLibCpp
