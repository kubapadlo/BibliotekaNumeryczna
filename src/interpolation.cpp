#include "NumLibCpp/interpolation.hpp" // Dołączenie odpowiedniego nagłówka
#include <vector>
#include <stdexcept>
#include <cmath> // Dla std::abs, przydatne do sprawdzania duplikatów

// Otwieramy przestrzeń nazw, aby definicja pasowała do deklaracji
namespace NumLibCpp {

/**
 * @brief Implementacja funkcji obliczającej wartość interpolowaną metodą Lagrange'a.
 */
double lagrange_interpolate(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double x_interp) {
    
    // --- Walidacja danych wejściowych ---

    if (x_nodes.size() != y_nodes.size()) {
        throw std::invalid_argument("Vectors x_nodes and y_nodes must have the same size.");
    }

    if (x_nodes.empty()) {
        throw std::invalid_argument("Input vectors cannot be empty.");
    }

    // Sprawdzenie, czy węzły x są unikalne (opcjonalne, ale bardzo zalecane)
    for (size_t i = 0; i < x_nodes.size(); ++i) {
        for (size_t j = i + 1; j < x_nodes.size(); ++j) {
            // Używamy małej tolerancji na wypadek problemów z precyzją zmiennoprzecinkową
            if (std::abs(x_nodes[i] - x_nodes[j]) < 1e-9) {
                throw std::invalid_argument("x_nodes must contain unique values.");
            }
        }
    }
    
    // --- Logika interpolacji Lagrange'a ---

    double interpolated_value = 0.0;
    size_t n = x_nodes.size();

    for (size_t i = 0; i < n; ++i) {
        double basis_polynomial = 1.0; // l_i(x)
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                basis_polynomial *= (x_interp - x_nodes[j]) / (x_nodes[i] - x_nodes[j]);
            }
        }
        interpolated_value += y_nodes[i] * basis_polynomial;
    }

    return interpolated_value;
}

} // namespace NumLibCpp