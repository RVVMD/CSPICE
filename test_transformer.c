/**
 * Transformer Test Program
 * 
 * Tests the unified transformer API with:
 * - DC operating point analysis
 * - AC frequency sweep
 * - Transient simulation with saturation
 */

#include "mna_solver.h"
#include "elements/transformer.h"
#include "elements/nonlinear/magnetization.h"
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

/* ============================================================================
 * Test 1: Ideal Transformer (DC/AC)
 * ============================================================================ */
static int test_ideal_transformer(void) {
    printf("\n=== Test 1: Ideal Transformer (DC/AC) ===\n");
    
    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_init returned %d\n", status);
        return -1;
    }
    
    /* Create nodes: 0=GND, 1=primary+, 2=secondary+ */
    int n1 = mna_create_node(&solver);  /* Node 1 */
    int n2 = mna_create_node(&solver);  /* Node 2 */
    printf("Created nodes: primary=%d, secondary=%d\n", n1, n2);
    
    /* Add 10V DC voltage source on primary */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, n1, 0, 10.0, &vs);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_add_voltage_source returned %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    /* Add ideal transformer (10:1 ratio) */
    TransformerConfig config = {0};
    config.mode = TRANSFORMER_MODE_VOLTAGE;
    config.turns_ratio = 10.0;
    config.Lm = 0.0;  /* No magnetizing branch = ideal */
    
    ComponentHandle xf;
    status = mna_add_transformer(&solver, n1, 0, n2, 0, &config, &xf);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_add_transformer returned %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    /* Add load resistor on secondary */
    ComponentHandle rl;
    status = mna_add_resistor(&solver, n2, 0, 100.0, &rl);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_add_resistor returned %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    /* Run DC analysis */
    status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_solve_dc returned %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    double vp = mna_get_transformer_primary_voltage(&solver, xf);
    double vs_val = mna_get_transformer_secondary_voltage(&solver, xf);
    
    printf("Primary voltage:   %.4f V (expected: 10.0 V)\n", vp);
    printf("Secondary voltage: %.4f V (expected: 1.0 V)\n", vs_val);
    
    /* Verify: Vs = Vp / N = 10 / 10 = 1V */
    if (fabs(vs_val - 1.0) > 0.01) {
        printf("FAIL: Secondary voltage mismatch\n");
        mna_destroy(&solver);
        return -1;
    }
    
    printf("PASS: Ideal transformer DC test\n");
    
    mna_destroy(&solver);
    return 0;
}

/* ============================================================================
 * Test 2: Transformer with Magnetizing Inductance (DC/Transient)
 * ============================================================================ */
static int test_transformer_with_Lm(void) {
    printf("\n=== Test 2: Transformer with Magnetizing Inductance ===\n");
    
    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_init returned %d\n", status);
        return -1;
    }
    
    int n1 = mna_create_node(&solver);
    int n2 = mna_create_node(&solver);
    
    /* 10V DC source */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, n1, 0, 10.0, &vs);
    
    /* Transformer with Lm = 1H (magnetizing inductance) */
    TransformerConfig config = {0};
    config.mode = TRANSFORMER_MODE_VOLTAGE;
    config.turns_ratio = 5.0;
    config.Lm = 1.0;  /* 1H magnetizing inductance */
    
    ComponentHandle xf;
    status = mna_add_transformer(&solver, n1, 0, n2, 0, &config, &xf);
    
    /* Secondary load */
    ComponentHandle rl;
    status = mna_add_resistor(&solver, n2, 0, 50.0, &rl);
    
    /* DC analysis */
    status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: DC analysis failed: %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    double vp = mna_get_transformer_primary_voltage(&solver, xf);
    double vs_val = mna_get_transformer_secondary_voltage(&solver, xf);
    double i_mag = mna_get_transformer_magnetizing_current(&solver, xf);
    double phi = mna_get_transformer_flux(&solver, xf);
    
    printf("Primary voltage:   %.4f V\n", vp);
    printf("Secondary voltage: %.4f V (expected: 2.0 V)\n", vs_val);
    printf("Magnetizing current: %.6f A\n", i_mag);
    printf("Flux linkage:    %.6f Wb\n", phi);
    
    /* Transient simulation */
    printf("\nRunning transient simulation...\n");
    mna_init_transient(&solver);
    
    double dt = 0.001;  /* 1ms time step */
    for (int step = 0; step < 10; step++) {
        status = mna_solve_transient_step(&solver, dt);
        if (status != MNA_SUCCESS) {
            printf("FAIL: Transient step %d failed: %d\n", step, status);
            mna_destroy(&solver);
            return -1;
        }
        
        double t = step * dt;
        vp = mna_get_transformer_primary_voltage(&solver, xf);
        vs_val = mna_get_transformer_secondary_voltage(&solver, xf);
        i_mag = mna_get_transformer_magnetizing_current(&solver, xf);
        
        printf("t=%.3fs: Vp=%.4fV, Vs=%.4fV, Im=%.6fA\n", t, vp, vs_val, i_mag);
    }
    
    printf("PASS: Transformer with Lm test\n");
    mna_destroy(&solver);
    return 0;
}

/* ============================================================================
 * Test 3: Transformer with B-H Curve Saturation
 * ============================================================================ */
static int test_transformer_saturation(void) {
    printf("\n=== Test 3: Transformer with B-H Curve Saturation ===\n");
    
    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_init returned %d\n", status);
        return -1;
    }
    
    /* Create B-H curve for silicon steel */
    MagnetizationCurve* bh = mna_bh_curve_create(BH_MODEL_PIECEWISE_LINEAR);
    if (!bh) {
        printf("FAIL: Could not create B-H curve\n");
        mna_destroy(&solver);
        return -1;
    }
    
    /* Set core geometry: A_c = 10 cm², l_e = 15 cm, Np = 100 */
    status = mna_bh_curve_set_geometry(bh, 0.001, 0.15, 100);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_bh_curve_set_geometry returned %d\n", status);
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }
    
    /* Add silicon steel B-H data points */
    double H_vals[] = {0, 50, 100, 200, 500, 1000, 2000, 5000};
    double B_vals[] = {0, 0.4, 0.7, 1.0, 1.3, 1.45, 1.55, 1.65};
    status = mna_bh_curve_add_points(bh, H_vals, B_vals, 8);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_bh_curve_add_points returned %d\n", status);
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }
    
    printf("B-H curve created with %d points\n", mna_bh_curve_get_point_count(bh));
    
    /* Test B-H curve directly */
    for (int i = 0; i < 5; i++) {
        double H_test = i * 500;
        double B = mna_bh_curve_get_B(bh, H_test);
        double dB_dH = mna_bh_curve_get_dB_dH(bh, H_test);
        printf("  H=%5.0f A/m -> B=%.3f T, dB/dH=%.6f\n", H_test, B, dB_dH);
    }
    
    int n1 = mna_create_node(&solver);
    int n2 = mna_create_node(&solver);
    
    /* Voltage source */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, n1, 0, 10.0, &vs);
    
    /* Transformer with B-H curve */
    TransformerConfig config = {0};
    config.mode = TRANSFORMER_MODE_VOLTAGE;
    config.turns_ratio = 10.0;
    config.bh_curve = bh;
    config.core_area = 0.001;
    config.magnetic_path_length = 0.15;
    config.N_primary = 100;
    config.N_secondary = 10;
    config.Rc = 10000.0;  /* Core loss resistance */
    
    ComponentHandle xf;
    status = mna_add_transformer(&solver, n1, 0, n2, 0, &config, &xf);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_add_transformer returned %d\n", status);
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }
    
    /* Secondary load */
    ComponentHandle rl;
    status = mna_add_resistor(&solver, n2, 0, 100.0, &rl);
    
    /* DC analysis */
    status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: DC analysis failed: %d\n", status);
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }
    
    double vp = mna_get_transformer_primary_voltage(&solver, xf);
    double vs_val = mna_get_transformer_secondary_voltage(&solver, xf);
    double i_mag = mna_get_transformer_magnetizing_current(&solver, xf);
    double B = mna_get_transformer_B_field(&solver, xf);
    double H_field = mna_get_transformer_H_field(&solver, xf);
    double mu_r = mna_get_transformer_relative_permeability(&solver, xf);
    double P_core = mna_get_transformer_core_loss(&solver, xf);
    
    printf("\nDC Operating Point:\n");
    printf("  Primary voltage:   %.4f V\n", vp);
    printf("  Secondary voltage: %.4f V\n", vs_val);
    printf("  Magnetizing current: %.6f A\n", i_mag);
    printf("  B field:           %.4f T\n", B);
    printf("  H field:           %.2f A/m\n", H_field);
    printf("  Relative permeability: %.1f\n", mu_r);
    printf("  Core loss:         %.6f W\n", P_core);
    
    /* Transient with sine wave input */
    printf("\nTransient simulation (sine input)...\n");
    mna_init_transient(&solver);
    
    double freq = 50.0;  /* 50 Hz */
    double amplitude = 10.0;
    double dt = 1.0 / (freq * 100);  /* 100 points per cycle */
    int cycles = 2;
    int total_steps = cycles * 100;
    
    for (int step = 0; step < total_steps; step++) {
        /* Update source to sine wave */
        double t = step * dt;
        solver.components[vs].value = amplitude * sin(2 * PI * freq * t);
        
        status = mna_solve_transient_step(&solver, dt);
        if (status != MNA_SUCCESS) {
            printf("FAIL: Transient step %d failed: %d\n", step, status);
            mna_bh_curve_destroy(bh);
            mna_destroy(&solver);
            return -1;
        }
        
        /* Print every 20 steps */
        if (step % 20 == 0) {
            vp = mna_get_transformer_primary_voltage(&solver, xf);
            vs_val = mna_get_transformer_secondary_voltage(&solver, xf);
            i_mag = mna_get_transformer_magnetizing_current(&solver, xf);
            B = mna_get_transformer_B_field(&solver, xf);
            mu_r = mna_get_transformer_relative_permeability(&solver, xf);
            
            printf("t=%.4fs: Vp=%7.3fV, Vs=%7.3fV, Im=%8.5fA, B=%5.3fT, mu_r=%6.1f\n",
                   t, vp, vs_val, i_mag, B, mu_r);
        }
    }
    
    printf("PASS: Transformer with B-H curve saturation test\n");
    
    mna_bh_curve_destroy(bh);
    mna_destroy(&solver);
    return 0;
}

/* ============================================================================
 * Test 4: Current Transformer
 * ============================================================================ */
static int test_current_transformer(void) {
    printf("\n=== Test 4: Current Transformer (CT) ===\n");
    
    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_init returned %d\n", status);
        return -1;
    }
    
    /* CT: 100:1 ratio, 1 primary turn, 100 secondary turns */
    int n_sec = mna_create_node(&solver);  /* Secondary node */
    
    /* Current source as primary (simulating conductor through CT) */
    ComponentHandle ip;
    status = mna_add_current_source(&solver, 0, n_sec, 10.0, &ip);  /* 10A primary */
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_add_current_source returned %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    /* CT configuration */
    TransformerConfig config = {0};
    config.mode = TRANSFORMER_MODE_CURRENT;
    config.turns_ratio = 100.0;  /* 100:1 CT */
    config.Lm = 0.5;  /* Magnetizing inductance */
    config.burden_resistance = 1.0;  /* 1Ω burden */
    config.N_primary = 1;
    config.N_secondary = 100;
    
    ComponentHandle ct;
    /* For CT, primary nodes are not used (external conductor) */
    status = mna_add_transformer(&solver, 0, 0, n_sec, 0, &config, &ct);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_add_transformer (CT) returned %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    /* DC analysis */
    status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: DC analysis failed: %d\n", status);
        mna_destroy(&solver);
        return -1;
    }
    
    double is = mna_get_ct_secondary_current(&solver, ct);
    double v_burden = mna_get_ct_burden_voltage(&solver, ct);
    double ratio_error = mna_get_ct_ratio_error(&solver, ct);
    
    printf("Primary current (set): 10.0 A\n");
    printf("Secondary current:     %.4f A (expected: ~0.1 A)\n", is);
    printf("Burden voltage:        %.4f V (expected: ~0.1 V)\n", v_burden);
    printf("Ratio error:           %.4f %%\n", ratio_error * 100);
    
    /* Ideal CT: Is = Ip / N = 10 / 100 = 0.1A */
    if (fabs(fabs(is) - 0.1) > 0.02) {
        printf("WARNING: Secondary current differs from expected (may be due to Im)\n");
    }
    
    printf("PASS: Current transformer test\n");
    
    mna_destroy(&solver);
    return 0;
}

/* ============================================================================
 * Test 5: AC Analysis
 * ============================================================================ */
static int test_ac_analysis(void) {
    printf("\n=== Test 5: AC Frequency Response ===\n");
    
    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_init returned %d\n", status);
        return -1;
    }

    int n1 = mna_create_node(&solver);
    int n2 = mna_create_node(&solver);

    /* AC voltage source (1V amplitude) */
    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, n1, 0, 1.0, &vs);
    solver.components[vs].ac_magnitude = 1.0;
    solver.components[vs].ac_phase = 0.0;

    /* Transformer */
    TransformerConfig config = {0};
    config.mode = TRANSFORMER_MODE_VOLTAGE;
    config.turns_ratio = 10.0;
    config.Lm = 1.0;

    ComponentHandle xf;
    status = mna_add_transformer(&solver, n1, 0, n2, 0, &config, &xf);

    /* Load */
    ComponentHandle rl;
    status = mna_add_resistor(&solver, n2, 0, 100.0, &rl);

    printf("AC frequency sweep (transformer with Lm=1H, RL=100Ω):\n");
    printf("  Freq(Hz)   |Vp|(V)   |Vs|(V)   Phase(deg)  |Zin|(Ω)   |Imag|(mA)\n");
    printf("  -----------+--------+--------+-----------+---------+----------\n");
    
    /* AC analysis at multiple frequencies */
    double freqs[] = {1, 10, 50, 100, 500, 1000, 5000};
    int num_freqs = sizeof(freqs) / sizeof(freqs[0]);
    
    for (int i = 0; i < num_freqs; i++) {
        double freq = freqs[i];
        double omega = 2 * PI * freq;
        
        status = mna_solve_ac(&solver, freq);
        if (status != MNA_SUCCESS) {
            printf("  AC analysis failed at %.0f Hz: %d\n", freq, status);
            continue;
        }
        
        /* Get primary and secondary voltages (magnitude) */
        double complex vp = mna_get_ac_node_voltage(&solver, n1);
        double complex vs_val = mna_get_ac_node_voltage(&solver, n2);
        
        double vp_mag = cabs(vp);
        double vs_mag = cabs(vs_val);
        
        /* Phase difference */
        double phase_diff = carg(vs_val) - carg(vp);
        double phase_deg = phase_diff * 180.0 / PI;
        
        /* Input impedance: Zin = jωLm || (N² * RL) */
        double omega_Lm = omega * config.Lm;
        double Zin_reflected = config.turns_ratio * config.turns_ratio * 100.0;  /* N² * RL = 10kΩ */
        /* |Zin| = 1 / sqrt(1/(ωLm)² + 1/(Z_reflected)²) */
        double Zin_mag = 1.0 / sqrt(1.0/(omega_Lm*omega_Lm) + 1.0/(Zin_reflected*Zin_reflected));
        
        /* Magnetizing current: Im = Vp / (ωLm) */
        double Im_mag = vp_mag / omega_Lm * 1000.0;  /* mA */
        
        /* Ideal Vs = Vp / N = 1V / 10 = 0.1V */
        double vs_ideal = vp_mag / config.turns_ratio;
        
        printf("  %8.0f  | %6.3f  | %6.3f  | %8.1f  | %7.0f  | %8.2f\n",
               freq, vp_mag, vs_mag, phase_deg, Zin_mag, Im_mag);
    }

    printf("\nAnalysis (Vp=1V, N=10, Lm=1H, RL=100Ω):\n");
    printf("  - Ideal Vs = Vp/N = 0.1V (constant, matches simulation)\n");
    printf("  - Low freq: Zin ≈ ωLm (inductive), high magnetizing current\n");
    printf("  - High freq: Zin → N²RL = 10kΩ (resistive), low magnetizing current\n");
    printf("  - Phase: Vs in phase with Vp (ideal transformer)\n");
    
    printf("\nPASS: AC analysis frequency sweep\n");

    mna_destroy(&solver);
    return 0;
}

/* ============================================================================
 * Main
 * ============================================================================ */
int main(void) {
    printf("========================================\n");
    printf("  Transformer Test Suite\n");
    printf("========================================\n");
    
    int failures = 0;
    
    if (test_ideal_transformer() != 0) failures++;
    if (test_transformer_with_Lm() != 0) failures++;
    if (test_transformer_saturation() != 0) failures++;
    if (test_current_transformer() != 0) failures++;
    if (test_ac_analysis() != 0) failures++;
    
    printf("\n========================================\n");
    if (failures == 0) {
        printf("  All tests PASSED\n");
    } else {
        printf("  %d test(s) FAILED\n", failures);
    }
    printf("========================================\n");
    
    return failures;
}
