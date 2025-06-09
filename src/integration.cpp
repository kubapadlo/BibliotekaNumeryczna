#include "NumLibCpp/integration.hpp"
#include <cmath> // Dla std::abs

namespace NumLibCpp {

double simpson_integrate(std::function<double(double)> func, double a, double b, int n) {
    if (n <= 0) {
        throw std::invalid_argument("Liczba podprzedzialow n musi byc dodatnia.");
    }
    if (n % 2 != 0) {
        throw std::invalid_argument("Liczba podprzedzialow n musi byc parzysta dla metody Simpsona.");
    }
    if (a == b) { // Całka po punkcie jest 0
        return 0.0;
    }
    if (a > b) { // Odwracamy granice i znak
        return -simpson_integrate(func, b, a, n);
    }


    double h = (b - a) / n;
    double sum = func(a) + func(b); // f(x_0) + f(x_n)

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2 == 1) { // Nieparzyste indeksy (x_1, x_3, ...)
            sum += 4 * func(x);
        } else { // Parzyste indeksy (x_2, x_4, ...)
            sum += 2 * func(x);
        }
    }

    return sum * h / 3.0;
}

} // namespace NumLibCpp