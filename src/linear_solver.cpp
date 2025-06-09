#include "NumLibCpp/linear_solver.hpp"
#include <cmath> // Dla std::abs, std::fabs
#include <algorithm> // Dla std::swap
#include <limits>    // Dla std::numeric_limits

namespace NumLibCpp {

std::vector<double> gauss_elimination(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();
    if (n == 0) {
        throw std::invalid_argument("Macierz A nie moze byc pusta.");
    }
    if (A[0].size() != static_cast<size_t>(n)) {
        throw std::invalid_argument("Macierz A musi byc kwadratowa.");
    }
    if (b.size() != static_cast<size_t>(n)) {
        throw std::invalid_argument("Rozmiar wektora b musi byc zgodny z rozmiarem macierzy A.");
    }

    for (int i = 0; i < n; ++i) {
        // Częściowy wybór elementu głównego (pivot)
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                max_row = k;
            }
        }
        std::swap(A[i], A[max_row]);
        std::swap(b[i], b[max_row]);

        // Sprawdzenie osobliwości
        if (std::abs(A[i][i]) < std::numeric_limits<double>::epsilon()) {
            throw std::runtime_error("Macierz jest osobliwa lub bliska osobliwej.");
        }

        // Normalizacja wiersza pivotowego
        for (int k = i + 1; k < n; ++k) {
            A[i][k] /= A[i][i];
        }
        b[i] /= A[i][i];
        A[i][i] = 1.0; // Ustawienie elementu diagonalnego na 1 (opcjonalne, ale czytelne)


        // Eliminacja dla pozostałych wierszy
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = i; j < n; ++j) { // Zaczynamy od kolumny i, bo A[k][<i] już są 0
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }

    // W tym momencie macierz A powinna być jednostkowa, a wektor b zawierać rozwiązanie
    // Jeśli nie normalizowaliśmy wiersza pivotowego (A[i][i] nie było dzielone),
    // to potrzebne jest podstawienie wsteczne. Przy pełnej eliminacji Gaussa-Jordana
    // jak tutaj (gdzie dążymy do macierzy jednostkowej), b jest rozwiązaniem.
    // Jednak dla pewności, wykonajmy podstawienie wsteczne na wypadek, gdyby A nie była idealnie jednostkowa
    // z powodu błędów numerycznych lub jeśli nie zrobiliśmy pełnej eliminacji Gaussa-Jordana.
    // W naszym przypadku, b już powinno zawierać rozwiązanie, ale poniższy kod jest bardziej ogólny dla eliminacji Gaussa.
    // UWAGA: Powyższa implementacja wykonuje eliminację Gaussa-Jordana, więc `b` jest rozwiązaniem.
    // Poniższe podstawienie wsteczne jest dla standardowej eliminacji Gaussa (górnotrójkątna macierz).
    // Dla Gaussa-Jordana `b` jest już gotowe.

    // W obecnej implementacji (Gauss-Jordan), b jest rozwiązaniem.
    // Jeśli chcemy klasycznego Gaussa (górnotrójkątna) i potem podstawienie wsteczne:
    // To pętla eliminacji wyglądałaby tak:
    /*
    for (int i = 0; i < n; ++i) {
        // ... (pivotowanie jak wyżej) ...
        // ... (sprawdzenie osobliwości jak wyżej) ...
        for (int k = i + 1; k < n; ++k) { // Eliminacja tylko pod diagonalą
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    // Podstawienie wsteczne
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        if (std::abs(A[i][i]) < std::numeric_limits<double>::epsilon()) { // Ponowne sprawdzenie, choć nie powinno tu dojść
             throw std::runtime_error("Macierz osobliwa podczas podstawiania wstecznego.");
        }
        x[i] /= A[i][i];
    }
    return x;
    */
    // Ponieważ zaimplementowałem Gauss-Jordana, gdzie celem jest macierz jednostkowa,
    // wektor 'b' po transformacjach jest bezpośrednio wektorem rozwiązania 'x'.
    return b;
}

} // namespace NumLibCpp