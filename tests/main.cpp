#include "test_helpers.hpp" // Nasz plik z makrami i podstawowymi nagłówkami

// Dołączamy wszystkie nagłówki z biblioteki, które będziemy testować
#include "NumLibCpp/linear_solver.hpp"
#include "NumLibCpp/interpolation.hpp"
#include "NumLibCpp/approximation.hpp"
#include "NumLibCpp/integration.hpp"
#include "NumLibCpp/differentialEquations_solver.hpp"
#include "NumLibCpp/nonlinear_solver.hpp"
#include "NumLibCpp/differentiation.hpp"

// --- 1. Testy algebry liniowej ---
void test_linear_algebra() {
    // Poprawny przypadek
    std::vector<std::vector<double>> A1 = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    std::vector<double> b1 = {8, -11, -3};
    std::vector<double> x1 = NumLibCpp::gauss_elimination(A1, b1);
    ASSERT_NEAR(x1[0], 2.0, 1e-9);
    ASSERT_NEAR(x1[1], 3.0, 1e-9);
    ASSERT_NEAR(x1[2], -1.0, 1e-9);
    std::cout << "  gauss_elimination (correct): PASSED" << std::endl;

    // Błędny przypadek (macierz osobliwa)
    std::vector<std::vector<double>> A2 = {{1, 1}, {1, 1}};
    std::vector<double> b2 = {2, 3};
    ASSERT_THROW(NumLibCpp::gauss_elimination(A2, b2), std::runtime_error);
    std::cout << "  gauss_elimination (singular matrix): PASSED" << std::endl;
}

// --- 2. Testy interpolacji ---
void test_interpolation() {
    // Poprawny przypadek (interpolacja kwadratowa)
    std::vector<double> x1 = {0.0, 1.0, 2.0};
    std::vector<double> y1 = {0.0, 1.0, 4.0}; // y = x^2
    double y_interp = NumLibCpp::lagrange_interpolate(x1, y1, 1.5);
    ASSERT_NEAR(y_interp, 2.25, 1e-9);
    std::cout << "  lagrange_interpolate (correct): PASSED" << std::endl;

    // Błędny przypadek (różne rozmiary wektorów)
    std::vector<double> x2 = {0.0, 1.0};
    std::vector<double> y2 = {0.0, 1.0, 2.0};
    ASSERT_THROW(NumLibCpp::lagrange_interpolate(x2, y2, 0.5), std::invalid_argument);
    std::cout << "  lagrange_interpolate (mismatched sizes): PASSED" << std::endl;
}

// --- 3. Testy aproksymacji ---
void test_approximation() {
    // Poprawny przypadek (aproksymacja wielomianowa)
    auto func = [](double x) { return x * x; }; // f(x) = x^2
    std::vector<double> coeffs = NumLibCpp::polynomial_approximation(func, -1.0, 1.0, 2, 100);
    ASSERT_NEAR(coeffs[2], 1.0, 1e-6); // Współczynnik przy x^2
    std::cout << "  polynomial_approximation (correct): PASSED" << std::endl;

    // Błędny przypadek (niepoprawne argumenty)
    ASSERT_THROW(NumLibCpp::polynomial_approximation(func, 0.0, 1.0, -1, 100), std::invalid_argument);
    std::cout << "  polynomial_approximation (invalid args): PASSED" << std::endl;
}

// --- 4. Testy całkowania ---
void test_integration() {
    // Poprawny przypadek (całka z x^2)
    auto func = [](double x) { return x * x; };
    double result = NumLibCpp::simpson_integrate(func, 0.0, 1.0, 100);
    ASSERT_NEAR(result, 1.0/3.0, 1e-7);
    std::cout << "  simpson_integrate (correct): PASSED" << std::endl;

    // Błędny przypadek (nieparzysta liczba przedziałów)
    ASSERT_THROW(NumLibCpp::simpson_integrate(func, 0.0, 1.0, 99), std::invalid_argument);
    std::cout << "  simpson_integrate (odd intervals): PASSED" << std::endl;
}

// --- 5. Testy równań różniczkowych ---
void test_ode_solver() {
    // Poprawny przypadek (y' = y, y(0) = 1 => y(1) = e)
    auto f_exp = [](double x, double y) { (void)x; return y; };
    double result = NumLibCpp::rk4_solve(f_exp, 0.0, 1.0, 1.0, 100);
    ASSERT_NEAR(result, std::exp(1.0), 1e-6);
    std::cout << "  rk4_solve (correct): PASSED" << std::endl;

    // Błędny przypadek (zerowa liczba kroków)
    ASSERT_THROW(NumLibCpp::rk4_solve(f_exp, 0.0, 1.0, 1.0, 0), std::invalid_argument);
    std::cout << "  rk4_solve (zero steps): PASSED" << std::endl;
}

// --- 6. Testy rozwiązywania równań nieliniowych ---
void test_nonlinear_solver() {
    // Poprawny przypadek (pierwiastek z x^2 - 2)
    auto func_sqrt2 = [](double x) { return x * x - 2.0; };
    double root = NumLibCpp::secant_method(func_sqrt2, 1.0, 2.0, 1e-7, 100);
    ASSERT_NEAR(root, std::sqrt(2.0), 1e-7);
    std::cout << "  secant_method (correct): PASSED" << std::endl;

    // Błędny przypadek (maksymalna liczba iteracji osiągnięta)
    auto func_no_real_root = [](double x){ return x*x + 1.0; };
    ASSERT_THROW(NumLibCpp::secant_method(func_no_real_root, -10.0, 10.0, 1e-5, 5), std::runtime_error);
    std::cout << "  secant_method (max iterations): PASSED" << std::endl;
}

// --- 7. Testy różniczkowania ---
void test_differentiation() {
    // Poprawny przypadek (pochodna z x^3)
    auto cubic_func = [](double x) { return x * x * x; };
    double derivative = NumLibCpp::central_difference(cubic_func, 2.0, 1e-5);
    ASSERT_NEAR(derivative, 12.0, 1e-6); // 3*x^2 dla x=2 to 12
    std::cout << "  central_difference (correct): PASSED" << std::endl;

    // Błędny przypadek (niepoprawny krok h)
    ASSERT_THROW(NumLibCpp::central_difference(cubic_func, 2.0, 0.0), std::invalid_argument);
    std::cout << "  central_difference (invalid h): PASSED" << std::endl;
}

// --- Główna funkcja uruchamiająca testy ---
int main() {
    struct TestCase {
        std::string name;
        std::function<void()> function;
    };
    std::vector<TestCase> tests;
    int tests_passed = 0;
    int tests_failed = 0;

    auto ADD_TEST = [&](const std::string& name, std::function<void()> func) {
        tests.push_back({name, func});
    };
    
    std::cout << std::fixed << std::setprecision(8);

    ADD_TEST("LinearAlgebra", test_linear_algebra);
    ADD_TEST("Interpolation", test_interpolation);
    ADD_TEST("Approximation", test_approximation);
    ADD_TEST("Integration", test_integration);
    ADD_TEST("ODESolver", test_ode_solver);
    ADD_TEST("NonlinearSolver", test_nonlinear_solver);
    ADD_TEST("Differentiation", test_differentiation);

    std::cout << "Running " << tests.size() << " test suites." << std::endl;

    for (const auto& test_case : tests) {
        std::cout << "==============================================" << std::endl;
        std::cout << "[ RUN      ] " << test_case.name << std::endl;
        try {
            test_case.function();
            std::cout << "[       OK ] " << test_case.name << " suite finished." << std::endl;
            tests_passed++;
        } catch (const std::exception& e) {
            std::cout << "[  FAILED  ] " << test_case.name << " suite failed! Exception: " << e.what() << std::endl;
            tests_failed++;
        } catch (...) {
            std::cout << "[  FAILED  ] " << test_case.name << " suite failed! Unknown exception." << std::endl;
            tests_failed++;
        }
    }

    std::cout << "\n==============================================" << std::endl;
    std::cout << "                 TEST SUMMARY" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Total test suites: " << tests.size() << std::endl;
    std::cout << "Passed: " << tests_passed << std::endl;
    std::cout << "Failed: " << tests_failed << std::endl;
    std::cout << "==============================================" << std::endl;

    return (tests_failed == 0) ? 0 : 1;
}