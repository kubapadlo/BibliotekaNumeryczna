#include "NumLibCpp/differentiation.hpp"
#include <cmath> // Dla std::abs

namespace NumLibCpp {

double central_difference(std::function<double(double)> func, double x, double h) {
    if (h <= 0.0) {
        throw std::invalid_argument("Krok h musi byc dodatni.");
    }
    // Dla bardzo małych h, (x+h) może być równe x numerycznie.
    // Zabezpieczenie przed tym jest skomplikowane i zależy od precyzji double.
    // Tutaj zakładamy, że h jest "rozsądnie" małe.
    return (func(x + h) - func(x - h)) / (2.0 * h);
}

} // namespace NumLibCpp