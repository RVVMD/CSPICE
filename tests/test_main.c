/*
 * Main Test Runner for CSPICE
 * Runs all unit tests and validation tests
 */

#include "../mna_solver.h"
#include "test_framework.h"
#include <stdio.h>
#include <stdlib.h>

/* External test runners */
extern int run_passive_tests(void);
extern int run_transformer_tests(void);
extern int run_source_tests(void);

/* Forward declarations */
int run_validation_tests(void);
int run_netlist_tests(void);
extern int run_transformer_comprehensive_tests(void);
extern int run_saturated_transformer_realistic_tests(void);

int main(int argc, char** argv) {
    printf("===========================================\n");
    printf("   CSPICE Circuit Simulation Test Suite\n");
    printf("===========================================\n\n");

    /* Run all test categories */
    run_passive_tests();
    run_transformer_tests();
    run_source_tests();
    run_transformer_comprehensive_tests();
    run_saturated_transformer_realistic_tests();
    run_validation_tests();
    run_netlist_tests();

    /* Print final summary */
    PRINT_TEST_SUMMARY();

    return tests_failed > 0 ? 1 : 0;
}
