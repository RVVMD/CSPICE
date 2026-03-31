#ifndef CSPICE_TEST_FRAMEWORK_H
#define CSPICE_TEST_FRAMEWORK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

/* Test result tracking */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

/* Assertion macros */
#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        printf("  FAIL: %s:%d - Condition failed: %s\n", __FILE__, __LINE__, #cond); \
        return 0; \
    } \
} while(0)

#define ASSERT_FALSE(cond) do { \
    if ((cond)) { \
        printf("  FAIL: %s:%d - Expected false but was true: %s\n", __FILE__, __LINE__, #cond); \
        return 0; \
    } \
} while(0)

#define ASSERT_EQ(expected, actual) do { \
    if ((expected) != (actual)) { \
        printf("  FAIL: %s:%d - Expected %d but got %d\n", __FILE__, __LINE__, (int)(expected), (int)(actual)); \
        return 0; \
    } \
} while(0)

#define ASSERT_DOUBLE_EQ(expected, actual, tolerance) do { \
    double diff = fabs((expected) - (actual)); \
    if (diff > (tolerance)) { \
        printf("  FAIL: %s:%d - Expected %.10f but got %.10f (diff: %.10e, tol: %.10e)\n", \
               __FILE__, __LINE__, (expected), (actual), diff, (tolerance)); \
        return 0; \
    } \
} while(0)

#define ASSERT_NULL(ptr) do { \
    if ((ptr) != NULL) { \
        printf("  FAIL: %s:%d - Expected NULL but got %p\n", __FILE__, __LINE__, (void*)(ptr)); \
        return 0; \
    } \
} while(0)

#define ASSERT_NOT_NULL(ptr) do { \
    if ((ptr) == NULL) { \
        printf("  FAIL: %s:%d - Expected non-NULL but got NULL\n", __FILE__, __LINE__); \
        return 0; \
    } \
} while(0)

/* Test runner macros */
#define RUN_TEST(test_func) do { \
    tests_run++; \
    printf("Running %s... ", #test_func); \
    if (test_func()) { \
        tests_passed++; \
        printf("PASSED\n"); \
    } else { \
        tests_failed++; \
        printf("FAILED\n"); \
    } \
} while(0)

/* Test function prototype */
typedef int (*test_func_t)(void);

/* Print test summary */
#define PRINT_TEST_SUMMARY() do { \
    printf("\n=== Test Summary ===\n"); \
    printf("Total:  %d\n", tests_run); \
    printf("Passed: %d\n", tests_passed); \
    printf("Failed: %d\n", tests_failed); \
    if (tests_failed == 0) { \
        printf("All tests PASSED!\n"); \
    } else { \
        printf("Some tests FAILED!\n"); \
    } \
} while(0)

/* Helper function to compare doubles with tolerance */
static inline bool doubles_equal(double a, double b, double tolerance) {
    return fabs(a - b) <= tolerance;
}

/* Calculate relative error */
static inline double relative_error(double expected, double actual) {
    if (expected == 0.0) return fabs(actual);
    return fabs((actual - expected) / expected);
}

#endif /* CSPICE_TEST_FRAMEWORK_H */
