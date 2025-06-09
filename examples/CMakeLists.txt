# Ten plik buduje dwa oddzielne programy demonstracyjne.

# --- Przykład 1: Problem Inżynierski ---
# Definiujemy cel wykonywalny o nazwie "engineering_problem"
# który zostanie zbudowany z pliku "example_engineering_problem.cpp"
add_executable(engineering_problem example_engineering_problem.cpp)

# Linkujemy ten program z naszą biblioteką NumLibCpp, aby
# mógł używać jej funkcji (rk4_solve, linear_least_squares, etc.).
target_link_libraries(engineering_problem PRIVATE NumLibCpp)


# --- Przykład 2: Funkcje Matematyczne ---
# Definiujemy drugi cel wykonywalny o nazwie "math_functions"
# zbudowany z pliku "example_math_functions.cpp"
add_executable(math_functions example_math_functions.cpp)

# Ten program również musi być zlinkowany z naszą biblioteką.
target_link_libraries(math_functions PRIVATE NumLibCpp)


# (Opcjonalnie, ale dobra praktyka) Definiowanie miejsca instalacji
# Jeśli uruchomisz "cmake --install .", te pliki zostaną skopiowane
# do katalogu bin/ w miejscu instalacji.
include(GNUInstallDirs)
install(TARGETS engineering_problem math_functions
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)