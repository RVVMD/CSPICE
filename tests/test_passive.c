/*
 * Unit Tests for Passive Elements (R, L, C)
 * Tests basic functionality and accuracy of resistor, capacitor, and inductor models
 */

#include "../mna_solver.h"
#include "test_framework.h"
#include "../src/solver/dc.h"
#include "../src/solver/transient.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Test 1: Simple resistor voltage divider */
static int test_resistor_voltage_divider(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Create nodes: V1(1-0) -> R1(1-2) -> R2(2-0) */
    int node0 = 0;  /* Ground */
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* 10V source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 10.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* R1 = 2kΩ */
    ComponentHandle r1;
    status = mna_add_resistor(&solver, node1, node2, 2000.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* R2 = 3kΩ */
    ComponentHandle r2;
    status = mna_add_resistor(&solver, node2, node0, 3000.0, &r2);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run DC analysis */
    status = mna_solve_dc(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* V(node2) should be 10V * (3k / (2k + 3k)) = 6V */
    double v_node2 = mna_get_node_voltage(&solver, node2);
    ASSERT_DOUBLE_EQ(6.0, v_node2, 1e-6);
    
    /* V(node1) should be 10V */
    double v_node1 = mna_get_node_voltage(&solver, node1);
    ASSERT_DOUBLE_EQ(10.0, v_node1, 1e-6);
    
    /* Current through series circuit: I = 10V / 5kΩ = 2mA (may be negative due to convention) */
    double i_v1 = mna_get_component_current(&solver, v1);
    ASSERT_TRUE(fabs(fabs(i_v1) - 0.002) < 1e-5);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 2: Resistor validation - reject negative/zero values */
static int test_resistor_validation(void) {
    MNASolver solver;
    MNAStatus status;
    ComponentHandle handle;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    
    /* Test negative resistance */
    status = mna_add_resistor(&solver, node1, node0, -100.0, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test zero resistance */
    status = mna_add_resistor(&solver, node1, node0, 0.0, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test very small but positive resistance (should succeed) */
    status = mna_add_resistor(&solver, node1, node0, 1e-6, &handle);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 3: Capacitor validation */
static int test_capacitor_validation(void) {
    MNASolver solver;
    MNAStatus status;
    ComponentHandle handle;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    
    /* Test negative capacitance */
    status = mna_add_capacitor(&solver, node1, node0, -1e-6, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test zero capacitance */
    status = mna_add_capacitor(&solver, node1, node0, 0.0, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 4: Inductor validation */
static int test_inductor_validation(void) {
    MNASolver solver;
    MNAStatus status;
    ComponentHandle handle;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    
    /* Test negative inductance */
    status = mna_add_inductor(&solver, node1, node0, -1e-3, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    /* Test zero inductance */
    status = mna_add_inductor(&solver, node1, node0, 0.0, &handle);
    ASSERT_TRUE(status == MNA_INVALID_PARAMETER);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 5: Parallel resistors */
static int test_parallel_resistors(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    
    /* 12V source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 12.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Three parallel resistors: 1kΩ, 2kΩ, 3kΩ */
    ComponentHandle r1, r2, r3;
    status = mna_add_resistor(&solver, node1, node0, 1000.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node1, node0, 2000.0, &r2);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node1, node0, 3000.0, &r3);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run DC analysis */
    status = mna_solve_dc(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Total current: I = V/R_eq where 1/R_eq = 1/1k + 1/2k + 1/3k = 11/6k */
    /* R_eq = 6k/11 = 545.45Ω */
    /* I_total = 12V / 545.45Ω = 22mA */
    double i_total = mna_get_component_current(&solver, v1);
    double expected_current = 12.0 / (1000.0 * 2000.0 * 3000.0 / 
        (2000.0 * 3000.0 + 1000.0 * 3000.0 + 1000.0 * 2000.0));
    ASSERT_DOUBLE_EQ(expected_current, fabs(i_total), 1e-6);
    
    /* Individual currents: I1 = 12mA, I2 = 6mA, I3 = 4mA */
    double i_r1 = mna_get_component_current(&solver, r1);
    ASSERT_DOUBLE_EQ(0.012, fabs(i_r1), 1e-6);
    
    double i_r2 = mna_get_component_current(&solver, r2);
    ASSERT_DOUBLE_EQ(0.006, fabs(i_r2), 1e-6);
    
    double i_r3 = mna_get_component_current(&solver, r3);
    ASSERT_DOUBLE_EQ(0.004, fabs(i_r3), 1e-6);
    
    mna_destroy(&solver);
    return 1;
}

/* Run all passive element tests */
int run_passive_tests(void) {
    printf("\n=== Passive Element Tests ===\n");
    
    RUN_TEST(test_resistor_voltage_divider);
    RUN_TEST(test_resistor_validation);
    RUN_TEST(test_capacitor_validation);
    RUN_TEST(test_inductor_validation);
    RUN_TEST(test_parallel_resistors);
    
    return tests_failed == 0;
}
