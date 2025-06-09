# Plik CMake dla programu testowego

add_executable(run_all_tests main.cpp)

# Linkuj nasz program z biblioteką NumLibCpp
# Dzięki temu w main.cpp będziemy mogli używać funkcji z biblioteki
target_link_libraries(run_all_tests PRIVATE NumLibCpp)