cmake_minimum_required(VERSION 3.15)
project(NumLibCpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Biblioteka ---
add_library(NumLibCpp
    src/linear_solver.cpp
    src/interpolation.cpp
    src/approximation.cpp
    src/integration.cpp
    src/differentialEquations_solver.cpp
    src/nonlinear_solver.cpp
    src/differentiation.cpp
)

target_include_directories(NumLibCpp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include> 
)

# --- Testy ---
add_subdirectory(tests)

# --- Przykłady ---
add_subdirectory(examples)

# --- Instalacja (opcjonalnie, ale dobra praktyka) ---
include(GNUInstallDirs)
install(TARGETS NumLibCpp
    EXPORT NumLibCppTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY include/NumLibCpp/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/NumLibCpp)

install(EXPORT NumLibCppTargets
    FILE NumLibCppTargets.cmake
    NAMESPACE NumLibCpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NumLibCpp
)

# Aby inne projekty CMake mogły łatwo znaleźć bibliotekę
export(PACKAGE NumLibCpp)