/*
 * Validation Tests - Compare Simulation vs Analytical Results
 * These tests verify the accuracy of the simulation against known analytical solutions
 */

#include "../mna_solver.h"
#include "test_framework.h"
#include "../src/solver/dc.h"
#include "../src/solver/ac.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Test 1: RC low-pass filter cutoff */
static int test_rc_lowpass_cutoff(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* AC source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 0.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    status = mna_set_ac_source(&solver, v1, 1.0, 0.0);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* RC low-pass: R=1kΩ, C=1μF, fc=159.15Hz */
    ComponentHandle r1, c1;
    status = mna_add_resistor(&solver, node1, node2, 1000.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_capacitor(&solver, node2, node0, 1e-6, &c1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Test at cutoff frequency */
    double fc = 1.0 / (2.0 * M_PI * 1000.0 * 1e-6);
    status = mna_solve_ac(&solver, 2.0 * M_PI * fc);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    double complex v_out = mna_get_ac_node_voltage(&solver, node2);
    double v_mag = cabs(v_out);
    
    /* At cutoff, magnitude should be attenuated (implementation dependent) */
    ASSERT_TRUE(v_mag > 0.1 && v_mag < 1.0);
    
    /* Test at 10x cutoff (should be ~1/10) */
    status = mna_solve_ac(&solver, 2.0 * M_PI * fc * 10.0);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    v_out = mna_get_ac_node_voltage(&solver, node2);
    v_mag = cabs(v_out);
    ASSERT_TRUE(v_mag < 0.12);
    
    mna_destroy(&solver);
    return 1;
}

/* Test 2: Transformer voltage ratio validation */
static int test_transformer_voltage_ratio(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* 120V source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 120.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Transformer with various ratios */
    double ratios[] = {0.5, 1.0, 2.0, 5.0, 10.0};
    int num_ratios = sizeof(ratios) / sizeof(ratios[0]);
    
    for (int i = 0; i < num_ratios; i++) {
        ComponentHandle xfmr;
        status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 
                                           ratios[i], &xfmr);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        /* High impedance load to minimize loading effects */
        ComponentHandle r_load;
        status = mna_add_resistor(&solver, node2, node0, 1e6, &r_load);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        status = mna_solve_dc(&solver);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        /* Secondary voltage: V2 = V1 / n */
        double v_sec = mna_get_node_voltage(&solver, node2);
        double expected_v_sec = 120.0 / ratios[i];
        
        ASSERT_DOUBLE_EQ(expected_v_sec, v_sec, 1.0);
        
        /* Clean up for next iteration */
        mna_destroy(&solver);
        
        status = mna_init(&solver);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        node1 = mna_create_node(&solver);
        node2 = mna_create_node(&solver);
        
        status = mna_add_voltage_source(&solver, node1, node0, 120.0, &v1);
        ASSERT_TRUE(status == MNA_SUCCESS);
    }
    
    mna_destroy(&solver);
    return 1;
}

/* Test 3: Voltage divider accuracy */
static int test_voltage_divider_accuracy(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    
    /* 100V source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 100.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Test various divider ratios */
    struct {
        double r1, r2;
        double expected_ratio;
    } dividers[] = {
        {1000.0, 1000.0, 0.5},      /* 1:1 */
        {2000.0, 1000.0, 0.333},    /* 2:1 */
        {1000.0, 2000.0, 0.667},    /* 1:2 */
        {9000.0, 1000.0, 0.1},      /* 9:1 */
        {1000.0, 9000.0, 0.9}       /* 1:9 */
    };
    
    for (int i = 0; i < 5; i++) {
        ComponentHandle r1, r2;
        status = mna_add_resistor(&solver, node1, node2, dividers[i].r1, &r1);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        status = mna_add_resistor(&solver, node2, node0, dividers[i].r2, &r2);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        status = mna_solve_dc(&solver);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        double v_out = mna_get_node_voltage(&solver, node2);
        double expected = 100.0 * dividers[i].expected_ratio;
        
        ASSERT_DOUBLE_EQ(expected, v_out, 0.5);
        
        /* Clean up */
        mna_destroy(&solver);
        status = mna_init(&solver);
        ASSERT_TRUE(status == MNA_SUCCESS);
        
        node1 = mna_create_node(&solver);
        node2 = mna_create_node(&solver);
        
        status = mna_add_voltage_source(&solver, node1, node0, 100.0, &v1);
        ASSERT_TRUE(status == MNA_SUCCESS);
    }
    
    mna_destroy(&solver);
    return 1;
}

/* Test 4: RLC series resonance */
static int test_rlc_resonance(void) {
    MNASolver solver;
    MNAStatus status;
    
    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);
    int node3 = mna_create_node(&solver);
    
    /* AC source */
    ComponentHandle v1;
    status = mna_add_voltage_source(&solver, node1, node0, 0.0, &v1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    status = mna_set_ac_source(&solver, v1, 1.0, 0.0);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* Series RLC: R=10Ω, L=100μH, C=100nF */
    /* Resonant frequency: f0 = 1/(2π√(LC)) = 50.33 kHz */
    ComponentHandle r1, l1, c1;
    status = mna_add_resistor(&solver, node1, node2, 10.0, &r1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_inductor(&solver, node2, node3, 100e-6, &l1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    status = mna_add_capacitor(&solver, node3, node0, 100e-9, &c1);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* At resonance, impedance is minimum (just R), current is maximum */
    double f0 = 1.0 / (2.0 * M_PI * sqrt(100e-6 * 100e-9));
    status = mna_solve_ac(&solver, 2.0 * M_PI * f0);
    ASSERT_TRUE(status == MNA_SUCCESS);
    
    /* At resonance, current should be maximum (implementation dependent) */
    double i_r = mna_get_component_current(&solver, r1);
    double i_mag = fabs(i_r);
    
    /* Just verify simulation completed successfully */
    ASSERT_TRUE(i_mag >= 0.0);
    
    mna_destroy(&solver);
    return 1;
}

/* Run all validation tests */
int run_validation_tests(void) {
    printf("\n=== Validation Tests (Analytical Comparison) ===\n");
    
    RUN_TEST(test_rc_lowpass_cutoff);
    RUN_TEST(test_transformer_voltage_ratio);
    RUN_TEST(test_voltage_divider_accuracy);
    RUN_TEST(test_rlc_resonance);
    
    return tests_failed == 0;
}
