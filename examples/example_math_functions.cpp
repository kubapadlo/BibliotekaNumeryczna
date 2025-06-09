#include <NumLibCpp/linear_solver.hpp>
#include <NumLibCpp/interpolation.hpp>
#include <NumLibCpp/approximation.hpp> // Dodane dla polynomial_approximation
#include <NumLibCpp/integration.hpp>
#include <NumLibCpp/differentiation.hpp>
#include <NumLibCpp/nonlinear_solver.hpp>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Funkcja do ewaluacji wielomianu dla przykładu aproksymacji
double eval_polynomial(double x, const std::vector<double>& coeffs) {
    double result = 0.0;
    double term = 1.0; // x^0
    for (double coeff : coeffs) {
        result += coeff * term;
        term *= x;
    }
    return result;
}


int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "--- NumLibCpp: Przyklady funkcji matematycznych ---" << std::endl << std::endl;

    // 1. Rozwiązywanie układu równań liniowych (Gauss) - bez zmian
    std::cout << "1. Rozwiazywanie ukladu rownan Ax=b metoda Gaussa:" << std::endl;
    std::vector<std::vector<double>> A_gauss = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    std::vector<double> b_gauss = {8, -11, -3};
    std::cout << "   Macierz A:" << std::endl;
    for(const auto& row : A_gauss) {
        std::cout << "     | ";
        for(double val : row) std::cout << std::setw(5) << val << " ";
        std::cout << "|" << std::endl;
    }
    std::cout << "   Wektor b: { ";
    for(size_t i=0; i<b_gauss.size(); ++i) std::cout << b_gauss[i] << (i==b_gauss.size()-1 ? "" : ", ");
    std::cout << "}" << std::endl;
    try {
        std::vector<double> x_gauss = NumLibCpp::gauss_elimination(A_gauss, b_gauss);
        std::cout << "   Rozwiazanie x = { ";
        for (size_t i = 0; i < x_gauss.size(); ++i) {
            std::cout << x_gauss[i] << (i == x_gauss.size() - 1 ? "" : ", ");
        }
        std::cout << "}" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "   Blad Gaussa: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    // 2. Interpolacja Lagrange'a - bez zmian
    std::cout << "2. Interpolacja Lagrange'a dla f(x) = x^2:" << std::endl;
    std::vector<double> x_nodes_interp = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y_nodes_interp;
    for(double x_val : x_nodes_interp) y_nodes_interp.push_back(x_val * x_val);
    double x_interp_val = 1.5;
    std::cout << "   Wezly x: { "; for(size_t i=0; i<x_nodes_interp.size(); ++i) std::cout << x_nodes_interp[i] << (i==x_nodes_interp.size()-1 ? "" : ", "); std::cout << "}" << std::endl;
    std::cout << "   Wezly y: { "; for(size_t i=0; i<y_nodes_interp.size(); ++i) std::cout << y_nodes_interp[i] << (i==y_nodes_interp.size()-1 ? "" : ", "); std::cout << "}" << std::endl;
    try {
        double y_interp = NumLibCpp::lagrange_interpolate(x_nodes_interp, y_nodes_interp, x_interp_val);
        std::cout << "   Wartosc interpolowana w x = " << x_interp_val << ": " << y_interp << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "   Blad interpolacji: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    // 3. Aproksymacja wielomianowa
    std::cout << "3. Aproksymacja wielomianowa funkcji f(x) = exp(x)*cos(5x) - x^3:" << std::endl;
    auto func_to_approx = [](double x) { return std::exp(x) * std::cos(5 * x) - std::pow(x, 3); };
    double approx_a = -1.0;
    double approx_b = 2.0;
    int approx_degree = 4; // Stopień wielomianu aproksymującego
    int simpson_intervals_approx = 200; // Liczba podprzedziałów dla Simpsona
    
    std::cout << "   Aproksymacja na przedziale [" << approx_a << ", " << approx_b << "] wielomianem stopnia " << approx_degree << "." << std::endl;
    try {
        std::vector<double> coeffs = NumLibCpp::polynomial_approximation(func_to_approx, approx_a, approx_b, approx_degree, simpson_intervals_approx);
        std::cout << "   Wspolczynniki wielomianu aproksymujacego P(x):" << std::endl;
        for (size_t i = 0; i < coeffs.size(); ++i) {
            std::cout << "     c" << i << " = " << coeffs[i] << (coeffs[i] >= 0 ? " " : "") << (i == coeffs.size() - 1 ? "" : "\n");
        }
        std::cout << std::endl;

        std::cout << "   Porownanie f(x) z P(x) w kilku punktach:" << std::endl;
        std::cout << "   x      | f(x)     | P(x)     | Blad_abs" << std::endl;
        std::cout << "   -------|----------|----------|----------" << std::endl;
        for (double x_test = approx_a; x_test <= approx_b + 1e-9; x_test += (approx_b - approx_a) / 5.0) {
            double fx = func_to_approx(x_test);
            double px = eval_polynomial(x_test, coeffs);
            std::cout << "   " << std::setw(6) << x_test << " | "
                      << std::setw(8) << fx << " | "
                      << std::setw(8) << px << " | "
                      << std::setw(8) << std::abs(fx - px) << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "   Blad aproksymacji wielomianowej: " << e.what() << std::endl;
    }
    std::cout << std::endl;


    // 4. Całkowanie numeryczne (Simpson) - było 3, teraz 4
    std::cout << "4. Calkowanie numeryczne f(x) = sin(x) od 0 do PI metoda Simpsona:" << std::endl;
    auto func_sin = [](double val) { return std::sin(val); };
    try {
        double integral_sin = NumLibCpp::simpson_integrate(func_sin, 0.0, M_PI, 100);
        std::cout << "   Calka = " << integral_sin << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "   Blad calkowania: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    // 5. Różniczkowanie numeryczne - było 4, teraz 5
    std::cout << "5. Rozniczkowanie numeryczne f(x) = x^3 w x=2:" << std::endl;
    auto func_cubic = [](double x) { return x * x * x; };
    try {
        double derivative = NumLibCpp::central_difference(func_cubic, 2.0, 1e-5);
        std::cout << "   Pochodna f'(2) = " << derivative << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "   Blad rozniczkowania: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    // 6. Rozwiązywanie równań nieliniowych (Metoda Siecznych) - było 5, teraz 6
    std::cout << "6. Rozwiazywanie rownania nieliniowego cos(x) - x = 0 metoda siecznych:" << std::endl;
    auto func_cos_x_eq_x = [](double x) { return std::cos(x) - x; };
    try {
        double root = NumLibCpp::secant_method(func_cos_x_eq_x, 0.0, 1.0, 1e-7, 100);
        std::cout << "   Pierwiastek x = " << root << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "   Blad metody siecznych: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    return 0;
}