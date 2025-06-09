#pragma once // Zabezpieczenie przed wielokrotnym dołączeniem

#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <typeinfo> // dla typeid

// Definicja M_PI dla MSVC i innych kompilatorów
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- JEDNA, SPÓJNA DEFINICJA MAKR TESTOWYCH ---

#define ASSERT_TRUE(condition) \
    if (!(condition)) { \
        std::cerr << "\nASSERT_TRUE FAILED: (" #condition ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        throw std::runtime_error("Assertion failed"); \
    }

#define ASSERT_THROW(expression, ExceptionType) \
    try { \
        (expression); \
        std::cerr << "\nASSERT_THROW FAILED: Expected " #ExceptionType " from (" #expression ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        throw std::runtime_error("Assertion failed: Expected exception"); \
    } catch (const ExceptionType&) { \
        /* Test passed */ \
    } catch (const std::exception& e) { \
        std::cerr << "\nASSERT_THROW FAILED: Expected " #ExceptionType " but got " << typeid(e).name() << " from (" #expression ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        throw std::runtime_error("Assertion failed: Wrong exception type"); \
    }

#define ASSERT_NEAR(val1, val2, tolerance) \
    if (std::abs((val1) - (val2)) > (tolerance)) { \
        std::cerr << std::fixed << std::setprecision(10) \
                  << "\nASSERT_NEAR FAILED: |" #val1 " - " #val2 "| > " #tolerance \
                  << " (values: |" << (val1) << " - " << (val2) << "| = " << std::abs((val1)-(val2)) << ")" \
                  << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        throw std::runtime_error("Assertion failed"); \
    }