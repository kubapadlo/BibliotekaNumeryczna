# BibliotekaNumeryczna

# NumLibCpp - Biblioteka Numeryczna w C++

NumLibCpp to biblioteka C++ dostarczająca implementacje popularnych metod obliczeniowych.

## Funkcjonalności

Biblioteka implementuje następujące metody:

*   **Rozwiązywanie układów równań liniowych:** Metoda eliminacji Gaussa.
*   **Interpolacja:** Interpolacja Lagrange'a.
*   **Aproksymacja:** Metoda najmniejszych kwadratów dla aproksymacji liniowej.
*   **Całkowanie numeryczne:** Złożona metoda Simpsona.
*   **Rozwiązywanie równań różniczkowych zwyczajnych:** Metoda Rungego-Kutty 4. rzędu (RK4).
*   **Rozwiązywanie równań nieliniowych:** Metoda Krzywej Lini
*   **Różniczkowanie numeryczne:** Metoda różnic centralnych.

## Struktura Projektu

NumLibCpp/

├── CMakeLists.txt       # Główny plik budujący CMake

├── README.md            # Ten plik

├── include/NumLibCpp/   # Pliki nagłówkowe (.h)

├── src/                 # Pliki źródłowe (.cpp)

├── tests/               # Testy jednostkowe

└── examples/            # Przykłady użycia

## Wymagania

*   Kompilator C++17
*   CMake

## Budowanie Biblioteki

1.  **Sklonuj repozytorium (jeśli jest w repozytorium git):**
    ```bash
    git clone <URL_REPOZYTORIUM>
    cd NumLibCpp
    ```

2.  **Utwórz katalog budowania i skonfiguruj CMake:**
    ```bash
    mkdir build
    cd build
    cmake ..
    ```

3.  **Skompiluj projekt:**
    ```bash
    cmake --build .
    ```

4.  **(Opcjonalnie) Zainstaluj bibliotekę:**
    ```bash
    cmake --install . --prefix /sciezka/do/instalacji
    ```

## Uruchamianie Testów

Po zbudowaniu projektu, testy można uruchomić z katalogu `build`:

```bash
./tests/Debug/run_all_tests.exe 
```

## Uruchamianie przykładów użycia
```
./examples/Debug/math_functions.exe
./examples/Debug/engineering_problem.exe
```
