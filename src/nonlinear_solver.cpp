#include "NumLibCpp/nonlinear_solver.hpp"
#include <cmath> // Dla std::abs, std::fabs

namespace NumLibCpp {

double secant_method(std::function<double(double)> func, double x0, double x1, double tol, int max_iter) {
    if (tol <= 0.0) {
        throw std::invalid_argument("Tolerancja musi byc dodatnia.");
    }
    if (max_iter <= 0) {
        throw std::invalid_argument("Maksymalna liczba iteracji musi byc dodatnia.");
    }
    if (std::abs(x0 - x1) < std::numeric_limits<double>::epsilon()) {
         throw std::invalid_argument("Poczatkowe przyblizenia x0 i x1 musza byc rozne.");
    }


    double fx0 = func(x0);
    double fx1 = func(x1);

    for (int i = 0; i < max_iter; ++i) {
        if (std::abs(fx1 - fx0) < std::numeric_limits<double>::epsilon() * 100) { // Mnożnik dla bezpieczeństwa
            // To może oznaczać, że f(x1) i f(x0) są bardzo blisko,
            // co może prowadzić do dzielenia przez bardzo małą liczbę, lub
            // że znaleźliśmy płaski region funkcji.
            // Jeśli fx1 jest bliskie 0, to x1 jest prawdopodobnie pierwiastkiem.
            if (std::abs(fx1) < tol) return x1;
            throw std::runtime_error("Dzielenie przez wartosc bliska zeru (fx1 - fx0). Funkcja moze byc plaska w poblizu przyblizen.");
        }

        double x_next = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        if (std::abs(x_next - x1) < tol || std::abs(func(x_next)) < tol) {
            return x_next;
        }

        x0 = x1;
        fx0 = fx1;
        x1 = x_next;
        fx1 = func(x1);
    }

    throw std::runtime_error("Metoda siecznych nie zbiegla w maksymalnej liczbie iteracji.");
}

} // namespace NumLibCpp