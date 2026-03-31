/*
 * Netlist-based Tests
 * Tests using circuit construction to verify end-to-end functionality
 */

#include "../mna_solver.h"
#include "test_framework.h"
#include "../src/solver/dc.h"
#include "../src/solver/transient.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Test 1: Simple RC circuit */
static int test_rc_circuit(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* V1 1 0 DC 10V */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 10.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* S1 1 2 0.001 (switch) */
    ComponentHandle s1;
    status = mna_add_switch(&solver, node1, node2, 0.001, &s1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* R1 2 0 1k */
    ComponentHandle r1;
    status = mna_add_resistor(&solver, node2, node0, 1000.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* C1 2 0 1u */
    ComponentHandle c1;
    status = mna_add_capacitor(&solver, node2, node0, 1e-6, &c1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run transient */
    mna_init_transient(&solver);
    
    /* Simulate for 5ms */
    for (int i = 0; i < 500; i++) {
        status = mna_solve_transient_step(&solver, 10e-6);
        ASSERT_TRUE(status == MNA_SUCCESS);
    }
    
    /* Capacitor should be charged to ~10V */
    double v_final = mna_get_node_voltage(&solver, node2);
    ASSERT_DOUBLE_EQ(10.0, v_final, 0.5);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 2: RLC series circuit */
static int test_rlc_series_circuit(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    int node3 = mna_create_node(&solver);
    
    /* V1 1 0 DC 5V */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 5.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* S1 1 2 */
    ComponentHandle s1;
    status = mna_add_switch(&solver, node1, node2, 0.001, &s1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* R1 2 3 100 */
    ComponentHandle r1;
    status = mna_add_resistor(&solver, node2, node3, 100.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* L1 3 0 100mH */
    ComponentHandle l1;
    status = mna_add_inductor(&solver, node3, node0, 0.1, &l1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* C1 2 0 10uF (parallel to R+L for stability) */
    ComponentHandle c1;
    status = mna_add_capacitor(&solver, node2, node0, 10e-6, &c1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run transient */
    mna_init_transient(&solver);
    
    /* Simulate */
    for (int i = 0; i < 200; i++) {
        status = mna_solve_transient_step(&solver, 10e-6);
        ASSERT_TRUE(status == MNA_SUCCESS);
    }
    
    /* Verify simulation completes without divergence */
    double v_node2 = mna_get_node_voltage(&solver, node2);
    ASSERT_TRUE(fabs(v_node2) < 100.0);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 3: Transformer with load */
static int test_transformer_load(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* V1 1 0 DC 24V */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 24.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* XFMR 1 0 2 0 2:1 */
    ComponentHandle xfmr;
    status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 2.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Rload 2 0 50 */
    ComponentHandle r_load;
    status = mna_add_resistor(&solver, node2, node0, 50.0, &r_load);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run DC */
    status = mna_solve_dc(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Secondary voltage should be 24V/2 = 12V */
    double v_sec = mna_get_node_voltage(&solver, node2);
    ASSERT_DOUBLE_EQ(12.0, v_sec, 1.0);
    
    /* Secondary current: 12V/50Ω = 240mA */
    double i_sec = mna_get_component_current(&solver, r_load);
    ASSERT_DOUBLE_EQ(0.24, fabs(i_sec), 0.02);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 4: Bridge circuit */
static int test_bridge_circuit(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    int node3 = mna_create_node(&solver);
    int node4 = mna_create_node(&solver);
    
    /* 10V source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 10.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Bridge resistors */
    ComponentHandle r1, r2, r3, r4;
    status = mna_add_resistor(&solver, node1, node2, 1000.0, &r1);  /* R1 */
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node2, node0, 2000.0, &r2);  /* R2 */
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node1, node3, 2000.0, &r3);  /* R3 */
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node3, node0, 4000.0, &r4);  /* R4 */
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run DC */
    status = mna_solve_dc(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Balanced bridge: R1/R2 = R3/R4 -> 1k/2k = 2k/4k = 0.5 */
    /* V_node2 = 10V * 2k/(1k+2k) = 6.67V */
    /* V_node3 = 10V * 4k/(2k+4k) = 6.67V */
    /* V_diff should be 0 for balanced bridge */
    double v2 = mna_get_node_voltage(&solver, node2);
    double v3 = mna_get_node_voltage(&solver, node3);
    
    ASSERT_TRUE(fabs(fabs(v2) - 6.667) < 0.5);
    ASSERT_TRUE(fabs(fabs(v3) - 6.667) < 0.5);
    ASSERT_TRUE(fabs(v2 - v3) < 0.1);  /* Balanced */
    
    mna_destroy(&solver);
    return 1;
}

/* Test 5: Current divider */
static int test_current_divider(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    
    /* 100mA current source */
    ComponentHandle i1;
    status = mna_add_current_source(&solver, node1, node0, 0.1, &i1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Three parallel resistors */
    ComponentHandle r1, r2, r3;
    status = mna_add_resistor(&solver, node1, node0, 300.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node1, node0, 600.0, &r2);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_resistor(&solver, node1, node0, 600.0, &r3);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Run DC */
    status = mna_solve_dc(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Total conductance: 1/300 + 1/600 + 1/600 = 1/150 */
    /* R_eq = 150Ω */
    /* V = I*R = 0.1A * 150Ω = 15V (may be negative due to convention) */
    double v_node1 = mna_get_node_voltage(&solver, node1);
    ASSERT_TRUE(fabs(fabs(v_node1) - 15.0) < 0.5);
    
    /* Currents: I1 = 15V/300Ω = 50mA, I2 = I3 = 15V/600Ω = 25mA */
    double i_r1 = mna_get_component_current(&solver, r1);
    ASSERT_DOUBLE_EQ(0.05, fabs(i_r1), 0.001);
    
    double i_r2 = mna_get_component_current(&solver, r2);
    ASSERT_DOUBLE_EQ(0.025, fabs(i_r2), 0.001);
    
    mna_destroy(&solver);
    return 1;
}

/* Run all netlist tests */
int run_netlist_tests(void) {
    printf("\n=== Netlist-based Tests ===\n");
    
    RUN_TEST(test_rc_circuit);
    RUN_TEST(test_rlc_series_circuit);
    RUN_TEST(test_transformer_load);
    RUN_TEST(test_current_divider);
    
    return tests_failed == 0;
}
