#include "NumLibCpp/differentialEquations_solver.hpp"
#include <NumLibCpp/approximation.hpp>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>

// Problem inżynierski 1: Chłodzenie obiektu (Prawo Newtona)
// dT/dt = -k * (T - T_env)
// Gdzie:
// T - temperatura obiektu
// T_env - temperatura otoczenia
// k - stała chłodzenia

// Funkcja dla ODE solvera
double cooling_law(double t, double T, double k_val, double T_env_val) {
    (void)t; // Czas nie jest jawnie używany w tej postaci równania
    return -k_val * (T - T_env_val);
}

// Problem inżynierski 2: Aproksymacja danych pomiarowych
// Załóżmy, że mamy pomiary naprężenia (sigma) w funkcji odkształcenia (epsilon)
// i chcemy znaleźć moduł Younga (E) z prawa Hooke'a (sigma = E * epsilon),
// co jest formą y = ax (gdzie b=0). Nasza funkcja aproksymacji jest y = ax + b,
// więc jeśli b wyjdzie bliskie 0, to 'a' będzie naszym E.

int main() {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "--- NumLibCpp: Przykłady problemów inżynierskich ---" << std::endl << std::endl;

    // 1. Symulacja chłodzenia obiektu
    std::cout << "1. Symulacja chlodzenia obiektu (prawo Newtona):" << std::endl;
    double T_initial = 100.0;  // [C] Początkowa temperatura obiektu
    double T_environment = 20.0; // [C] Temperatura otoczenia
    double k_cooling = 0.05;   // [1/min] Stała chłodzenia
    double time_target = 10.0; // [min] Czas, dla którego chcemy znać temperaturę
    int num_steps = 100;

    // Użycie std::bind do przekazania dodatkowych parametrów k i T_env
    auto cooling_ode_func = std::bind(cooling_law,
                                      std::placeholders::_1, // dla t
                                      std::placeholders::_2, // dla T
                                      k_cooling, T_environment);
    try {
        double T_final = NumLibCpp::rk4_solve(cooling_ode_func, 0.0, T_initial, time_target, num_steps);
        std::cout << "   Temperatura poczatkowa: " << T_initial << " C" << std::endl;
        std::cout << "   Temperatura otoczenia: " << T_environment << " C" << std::endl;
        std::cout << "   Stala chlodzenia k: " << k_cooling << " 1/min" << std::endl;
        std::cout << "   Temperatura po " << time_target << " minutach: " << T_final << " C" << std::endl;
        // Analityczne rozwiązanie: T(t) = T_env + (T_initial - T_env) * exp(-k*t)
        double T_analytical = T_environment + (T_initial - T_environment) * std::exp(-k_cooling * time_target);
        std::cout << "   Temperatura analityczna: " << T_analytical << " C" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "   Blad symulacji chlodzenia: " << e.what() << std::endl;
    }
    std::cout << std::endl;

    return 0;
}
