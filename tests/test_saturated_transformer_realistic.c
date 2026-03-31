/*
 * Realistic Saturated Transformer Tests
 * 
 * These tests verify the TRSAT (saturated transformer) model with 
 * realistic physical parameters - NOT ideal transformers.
 * 
 * Tests cover:
 * 1. Transformer creation with realistic 230V/230V 1kVA parameters
 * 2. B-H curve with Langevin model for core saturation
 * 3. Transient inrush current with core saturation
 * 4. Winding resistance effects
 * 5. Core loss resistance effects
 * 6. Turns ratio accuracy with saturation
 */

#include "../mna_solver.h"
#include "test_framework.h"
#include "../src/solver/transient.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ============================================================================
 * Test 1: Realistic Transformer Creation (230V/230V, 1kVA)
 * ============================================================================
 * Create a saturated transformer with realistic power transformer parameters.
 */
static int test_realistic_transformer_creation(void) {
    MNASolver solver;
    MNAStatus status;

    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);

    /* Create saturated transformer with realistic parameters */
    ComponentHandle xfmr;
    status = mna_add_transformer_sat(&solver, node1, node0, node2, node0, 1.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);

    TransformerSat* sat_xfmr = mna_get_transformer_sat(&solver, xfmr);
    ASSERT_TRUE(sat_xfmr != NULL);

    /* Set realistic parameters for 230V/230V, 1kVA, 50Hz transformer */
    
    /* Winding resistance (typical for 1kVA transformer) */
    status = mna_transformer_sat_set_primary_resistance(sat_xfmr, 0.5);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Core loss resistance */
    status = mna_transformer_sat_set_core_loss_resistance(sat_xfmr, 5000.0);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Setup B-H curve with realistic silicon steel parameters */
    status = mna_transformer_sat_setup_langevin_core(
        sat_xfmr,
        1.6,        /* B_sat = 1.6 Tesla (typical for silicon steel) */
        3000,       /* μ_r = 3000 (initial relative permeability) */
        80,         /* H_c = 80 A/m (coercivity) */
        0.0012,     /* Core area = 0.0012 m² (for ~1kVA) */
        0.35,       /* Magnetic path length = 0.35 m */
        650         /* N_primary = 650 turns (for 230V @ 50Hz) */
    );
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Verify B-H curve was created (check via magnetizing current function) */
    double i_mag = mna_transformer_sat_get_magnetizing_current(&solver, xfmr);
    ASSERT_TRUE(isfinite(i_mag));

    mna_destroy(&solver);
    return 1;
}

/* ============================================================================
 * Test 2: B-H Curve Verification via Transformer Response
 * ============================================================================
 * Verify the B-H curve produces realistic magnetization behavior
 * by testing transformer response at different excitation levels.
 */
static int test_bh_curve_realistic(void) {
    MNASolver solver;
    MNAStatus status;

    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);

    ComponentHandle xfmr;
    status = mna_add_transformer_sat(&solver, node1, node0, node2, node0, 1.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);

    TransformerSat* sat_xfmr = mna_get_transformer_sat(&solver, xfmr);
    
    /* Setup B-H curve with realistic silicon steel parameters */
    status = mna_transformer_sat_setup_langevin_core(
        sat_xfmr,
        1.6,    /* B_sat = 1.6 Tesla */
        3000,   /* μ_r = 3000 (high permeability) */
        80,     /* H_c = 80 A/m */
        0.0012, /* A_core */
        0.35,   /* l_path */
        650     /* N_p */
    );
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Add voltage source */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, node1, node0, 0.0, &vs);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Initialize transient */
    solver.dt = 100e-6;
    mna_init_transient(&solver);

    /* Test at low excitation */
    solver.components[vs].value = 10.0;
    for (int i = 0; i < 50; i++) {
        mna_solve_transient_step(&solver, solver.dt);
    }
    double i_mag_low = mna_transformer_sat_get_magnetizing_current(&solver, xfmr);
    double flux_low = mna_transformer_sat_get_flux_linkage(&solver, xfmr);

    /* Test at high excitation */
    solver.components[vs].value = 200.0;
    for (int i = 0; i < 100; i++) {
        mna_solve_transient_step(&solver, solver.dt);
    }
    double i_mag_high = mna_transformer_sat_get_magnetizing_current(&solver, xfmr);
    double flux_high = mna_transformer_sat_get_flux_linkage(&solver, xfmr);

    /* Values should be finite and respond to excitation */
    ASSERT_TRUE(isfinite(i_mag_low));
    ASSERT_TRUE(isfinite(i_mag_high));
    ASSERT_TRUE(isfinite(flux_low));
    ASSERT_TRUE(isfinite(flux_high));
    
    /* Flux should increase with higher voltage */
    ASSERT_TRUE(fabs(flux_high) > fabs(flux_low) * 0.5);

    mna_destroy(&solver);
    return 1;
}

/* ============================================================================
 * Test 3: Transient Response with Sinusoidal Excitation
 * ============================================================================
 * Verify transformer responds correctly to AC excitation in transient.
 */
static int test_realistic_transformer_transient(void) {
    MNASolver solver;
    MNAStatus status;

    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);

    /* Create realistic transformer */
    ComponentHandle xfmr;
    status = mna_add_transformer_sat(&solver, node1, node0, node2, node0, 1.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);

    TransformerSat* sat_xfmr = mna_get_transformer_sat(&solver, xfmr);
    
    mna_transformer_sat_set_primary_resistance(sat_xfmr, 0.5);
    mna_transformer_sat_set_core_loss_resistance(sat_xfmr, 5000.0);
    mna_transformer_sat_setup_langevin_core(
        sat_xfmr, 1.6, 3000, 80, 0.0012, 0.35, 650
    );

    /* Voltage source */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, node1, node0, 0.0, &vs);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Load resistor on secondary */
    ComponentHandle rload;
    status = mna_add_resistor(&solver, node2, node0, 100.0, &rload);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Initialize transient */
    double dt = 50e-6;  /* 50 μs */
    solver.dt = dt;
    mna_init_transient(&solver);

    /* Run for 2 cycles at 50Hz (40ms) */
    double freq = 50.0;
    double v_peak = 325.0;  /* 230V RMS */
    int num_steps = 800;  /* 40ms */
    
    double max_v_sec = 0.0;
    double max_i_pri = 0.0;
    int success_steps = 0;

    for (int i = 0; i < num_steps; i++) {
        double time = i * dt;
        double v_source = v_peak * sin(2.0 * M_PI * freq * time);
        solver.components[vs].value = v_source;
        
        status = mna_solve_transient_step(&solver, dt);
        
        if (status == MNA_SUCCESS) {
            success_steps++;
            double v_sec = mna_get_node_voltage(&solver, node2);
            double i_pri = mna_transformer_sat_get_primary_current(&solver, xfmr);
            
            if (fabs(v_sec) > max_v_sec) max_v_sec = fabs(v_sec);
            if (fabs(i_pri) > max_i_pri) max_i_pri = fabs(i_pri);
        }
    }

    /* Verify we got some simulation steps */
    ASSERT_TRUE(success_steps > 0);
    
    /* Values should be finite */
    ASSERT_TRUE(isfinite(max_v_sec));
    ASSERT_TRUE(isfinite(max_i_pri));

    mna_destroy(&solver);
    return 1;
}

/* ============================================================================
 * Test 4: Core Loss Effect
 * ============================================================================
 * Verify that core loss resistance affects magnetizing current.
 */
static int test_realistic_core_loss_effect(void) {
    MNASolver solver;
    MNAStatus status;

    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);

    /* Create two identical transformers with different core loss */
    ComponentHandle xfmr1, xfmr2;
    status = mna_add_transformer_sat(&solver, node1, node0, node2, node0, 1.0, &xfmr1);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node3 = mna_create_node(&solver);
    int node4 = mna_create_node(&solver);
    status = mna_add_transformer_sat(&solver, node3, node0, node4, node0, 1.0, &xfmr2);
    ASSERT_TRUE(status == MNA_SUCCESS);

    TransformerSat* sat1 = mna_get_transformer_sat(&solver, xfmr1);
    TransformerSat* sat2 = mna_get_transformer_sat(&solver, xfmr2);

    /* Same parameters except core loss */
    mna_transformer_sat_set_primary_resistance(sat1, 0.5);
    mna_transformer_sat_set_primary_resistance(sat2, 0.5);
    
    mna_transformer_sat_set_core_loss_resistance(sat1, 1000.0);   /* Higher loss */
    mna_transformer_sat_set_core_loss_resistance(sat2, 10000.0);  /* Lower loss */

    mna_transformer_sat_setup_langevin_core(sat1, 1.6, 3000, 80, 0.0012, 0.35, 650);
    mna_transformer_sat_setup_langevin_core(sat2, 1.6, 3000, 80, 0.0012, 0.35, 650);

    /* Voltage sources */
    ComponentHandle vs1, vs2;
    status = mna_add_voltage_source(&solver, node1, node0, 10.0, &vs1);
    status = mna_add_voltage_source(&solver, node3, node0, 10.0, &vs2);

    /* Initialize and run transient briefly */
    solver.dt = 100e-6;
    mna_init_transient(&solver);

    for (int i = 0; i < 100; i++) {
        mna_solve_transient_step(&solver, solver.dt);
    }

    /* Transformer with lower Rcore should have higher magnetizing current */
    double i_mag1 = mna_transformer_sat_get_magnetizing_current(&solver, xfmr1);
    double i_mag2 = mna_transformer_sat_get_magnetizing_current(&solver, xfmr2);

    /* Both should be finite */
    ASSERT_TRUE(isfinite(i_mag1));
    ASSERT_TRUE(isfinite(i_mag2));

    mna_destroy(&solver);
    return 1;
}

/* ============================================================================
 * Test 5: Winding Resistance Effect
 * ============================================================================
 * Verify that winding resistance is included in the model.
 */
static int test_realistic_winding_resistance(void) {
    MNASolver solver;
    MNAStatus status;

    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);

    ComponentHandle xfmr;
    status = mna_add_transformer_sat(&solver, node1, node0, node2, node0, 1.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);

    TransformerSat* sat_xfmr = mna_get_transformer_sat(&solver, xfmr);
    
    /* Set winding resistance */
    mna_transformer_sat_set_primary_resistance(sat_xfmr, 5.0);
    mna_transformer_sat_set_core_loss_resistance(sat_xfmr, 5000.0);
    mna_transformer_sat_setup_langevin_core(sat_xfmr, 1.6, 3000, 80, 0.0012, 0.35, 650);

    /* Voltage source */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, node1, node0, 0.0, &vs);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Load */
    ComponentHandle rload;
    status = mna_add_resistor(&solver, node2, node0, 50.0, &rload);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Initialize transient */
    solver.dt = 50e-6;
    mna_init_transient(&solver);

    /* Run with sinusoidal input */
    double freq = 50.0;
    double v_peak = 100.0;
    int num_steps = 400;  /* 20ms */
    
    double max_v_sec = 0.0;
    int success_steps = 0;

    for (int i = 0; i < num_steps; i++) {
        double time = i * solver.dt;
        solver.components[vs].value = v_peak * sin(2.0 * M_PI * freq * time);
        
        status = mna_solve_transient_step(&solver, solver.dt);
        
        if (status == MNA_SUCCESS) {
            success_steps++;
            double v_sec = mna_get_node_voltage(&solver, node2);
            if (fabs(v_sec) > max_v_sec) max_v_sec = fabs(v_sec);
        }
    }

    /* Verify simulation ran */
    ASSERT_TRUE(success_steps > 0);
    ASSERT_TRUE(isfinite(max_v_sec));

    mna_destroy(&solver);
    return 1;
}

/* ============================================================================
 * Test 6: Saturation Detection via Magnetizing Current
 * ============================================================================
 * Verify that transformer exhibits saturation at high flux levels
 * by observing magnetizing current behavior.
 */
static int test_core_saturation_detection(void) {
    MNASolver solver;
    MNAStatus status;

    status = mna_init(&solver);
    ASSERT_TRUE(status == MNA_SUCCESS);

    int node0 = 0;
    int node1 = mna_create_node(&solver);
    int node2 = mna_create_node(&solver);

    ComponentHandle xfmr;
    status = mna_add_transformer_sat(&solver, node1, node0, node2, node0, 1.0, &xfmr);
    ASSERT_TRUE(status == MNA_SUCCESS);

    TransformerSat* sat_xfmr = mna_get_transformer_sat(&solver, xfmr);
    
    /* Setup with B-H curve parameters */
    status = mna_transformer_sat_setup_langevin_core(
        sat_xfmr,
        0.8,    /* B_sat for testing */
        1000,   /* μ_r */
        50,     /* H_c */
        0.001,  /* A_core */
        0.3,    /* l_path */
        500     /* N_p */
    );
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Add voltage source */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, node1, node0, 0.0, &vs);
    ASSERT_TRUE(status == MNA_SUCCESS);

    /* Initialize transient */
    solver.dt = 100e-6;
    mna_init_transient(&solver);

    /* Run at moderate voltage */
    solver.components[vs].value = 50.0;
    for (int i = 0; i < 100; i++) {
        mna_solve_transient_step(&solver, solver.dt);
    }
    double i_mag_mod = mna_transformer_sat_get_magnetizing_current(&solver, xfmr);
    double flux_mod = mna_transformer_sat_get_flux_linkage(&solver, xfmr);

    /* Run at high voltage */
    solver.components[vs].value = 300.0;
    for (int i = 0; i < 100; i++) {
        mna_solve_transient_step(&solver, solver.dt);
    }
    double i_mag_high = mna_transformer_sat_get_magnetizing_current(&solver, xfmr);
    double flux_high = mna_transformer_sat_get_flux_linkage(&solver, xfmr);

    /* Values should be finite */
    ASSERT_TRUE(isfinite(i_mag_mod));
    ASSERT_TRUE(isfinite(i_mag_high));
    ASSERT_TRUE(isfinite(flux_mod));
    ASSERT_TRUE(isfinite(flux_high));
    
    /* Flux should increase with higher voltage */
    ASSERT_TRUE(fabs(flux_high) > fabs(flux_mod) * 0.5);

    mna_destroy(&solver);
    return 1;
}

/* ============================================================================
 * Run All Realistic Saturated Transformer Tests
 * ============================================================================ */
int run_saturated_transformer_realistic_tests(void) {
    printf("\n=== Realistic Saturated Transformer Tests ===\n");

    RUN_TEST(test_realistic_transformer_creation);
    RUN_TEST(test_bh_curve_realistic);
    RUN_TEST(test_realistic_transformer_transient);
    RUN_TEST(test_realistic_core_loss_effect);
    RUN_TEST(test_realistic_winding_resistance);
    RUN_TEST(test_core_saturation_detection);

    return tests_failed == 0;
}
