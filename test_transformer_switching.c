/**
 * Transformer Switching Transient Test
 *
 * Simulates a transformer being connected/disconnected via switches:
 * - t=0ms: Both switches open (transformer de-energized)
 * - t=10ms: Switch 1 closes (connects to source)
 * - t=50ms: Switch 2 closes (connects load)
 * - t=150ms: Switch 2 opens (disconnects load)
 * - t=200ms: Switch 1 opens (disconnects from source)
 *
 * Outputs all voltages and currents to CSV file.
 */

#include "mna_solver.h"
#include "elements/transformer.h"
#include "elements/nonlinear/magnetization.h"
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

int main(void) {
    printf("========================================\n");
    printf("  Transformer Switching Transient Test\n");
    printf("========================================\n\n");

    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: mna_init returned %d\n", status);
        return -1;
    }

    /* ==================== Circuit Topology ====================
     *
     *      S1          Transformer         S2
     *   o--/ /--+-----+      +-----+-----+--/ /--+-----o
     *   |       |     |      |     |     |       |     |
     *  Vs      ( )   === Lm  ( )   |    === RL    |    GND
     *   |       |     |      |     |     |       |
     *   +-------+-----+------+-   +-----+-------+
     *                           |
     *                          GND
     *
     * Nodes:
     *   0 = Ground
     *   1 = Source side of S1
     *   2 = Primary side of transformer
     *   3 = Secondary side of transformer / S2
     *   4 = Load side of S2
     */

    int n_source = mna_create_node(&solver);    /* Node 1: Source */
    int n_primary = mna_create_node(&solver);   /* Node 2: Primary */
    int n_secondary = mna_create_node(&solver); /* Node 3: Secondary */
    int n_load = mna_create_node(&solver);      /* Node 4: Load */

    printf("Created nodes:\n");
    printf("  Source node:  %d\n", n_source);
    printf("  Primary node: %d\n", n_primary);
    printf("  Secondary node: %d\n", n_secondary);
    printf("  Load node:    %d\n", n_load);

    /* ==================== Components ==================== */

    /* AC Voltage source: 230V RMS = 325V peak, 50Hz */
    double V_peak = 230.0 * sqrt(2.0);  /* 325V */
    double freq = 50.0;

    ComponentHandle vs;
    status = mna_add_voltage_source(&solver, n_source, 0, V_peak, &vs);
    if (status != MNA_SUCCESS) {
        printf("FAIL: Could not add voltage source\n");
        mna_destroy(&solver);
        return -1;
    }

    /* Source impedance: Rs (typical grid resistance) */
    /* Rs = 1.0Ω (limits inrush current) */
    /* Create intermediate node for source impedance */
    int n_rs = mna_create_node(&solver);  /* After voltage source + Rs */

    ComponentHandle rs;
    status = mna_add_resistor(&solver, n_source, n_rs, 1.0, &rs);

    /* Switch S1: Connects source (with impedance) to transformer primary */
    /* Initially OPEN (state=0), closes at t=7.5ms */
    ComponentHandle s1;
    status = mna_add_switch(&solver, n_rs, n_primary, 0.001, &s1);  /* 1mΩ when closed */
    if (status != MNA_SUCCESS) {
        printf("FAIL: Could not add switch S1\n");
        mna_destroy(&solver);
        return -1;
    }
    solver.components[s1].state = 0;  /* Initially open */

    /* Switch S2: Connects transformer secondary to load */
    /* Initially OPEN (state=0), closes at t=50ms */
    ComponentHandle s2;
    status = mna_add_switch(&solver, n_secondary, n_load, 0.001, &s2);
    if (status != MNA_SUCCESS) {
        printf("FAIL: Could not add switch S2\n");
        mna_destroy(&solver);
        return -1;
    }
    solver.components[s2].state = 0;  /* Initially open */

    /* Transformer: 230V/24V, 50Hz, with saturation */
    double turns_ratio = 230.0 / 24.0;  /* ~9.58:1 */

    /* Create B-H curve for silicon steel core */
    MagnetizationCurve* bh = mna_bh_curve_create(BH_MODEL_PIECEWISE_LINEAR);
    mna_bh_curve_set_geometry(bh, 0.002, 0.2, 500);  /* A_c=20cm², l_e=20cm, Np=500 */

    /* Silicon steel B-H data */
    double H_vals[] = {0, 50, 100, 200, 500, 1000, 2000, 5000};
    double B_vals[] = {0, 0.4, 0.7, 1.0, 1.3, 1.45, 1.55, 1.65};
    mna_bh_curve_add_points(bh, H_vals, B_vals, 8);

    TransformerConfig config = {0};
    config.mode = TRANSFORMER_MODE_VOLTAGE;
    config.turns_ratio = turns_ratio;
    config.bh_curve = bh;
    config.core_area = 0.002;
    config.magnetic_path_length = 0.2;
    config.N_primary = 500;
    config.N_secondary = 52;  /* 24V secondary */
    config.Rc = 50000.0;  /* Core loss resistance */
    config.initial_flux = 0.0;

    ComponentHandle xf;
    status = mna_add_transformer(&solver, n_primary, 0, n_secondary, 0, &config, &xf);
    if (status != MNA_SUCCESS) {
        printf("FAIL: Could not add transformer\n");
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }

    /* Load: 100W resistive load at 24V */
    /* R = V²/P = 24²/100 = 5.76Ω */
    ComponentHandle rl;
    status = mna_add_resistor(&solver, n_load, 0, 5.76, &rl);
    if (status != MNA_SUCCESS) {
        printf("FAIL: Could not add load resistor\n");
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }

    printf("\nCircuit parameters:\n");
    printf("  Source:         %.1f V peak, %.0f Hz\n", V_peak, freq);
    printf("  Source impedance: Rs=%.1fΩ (limits inrush current)\n", 1.0);
    printf("  Transformer:    230V/24V (N=%.2f)\n", turns_ratio);
    printf("  Load:           100W @ 24V (R=%.2f Ω)\n", 24.0*24.0/100.0);
    printf("  Core:           Silicon steel, A_c=20cm², l_e=20cm\n");

    /* ==================== Switching Schedule ==================== */
    printf("\nSwitching schedule (2 cycles over 500ms):\n");
    printf("  Cycle 1:\n");
    printf("    t=0ms:    S1=OPEN,  S2=OPEN   (transformer de-energized)\n");
    printf("    t=7.5ms:  S1=CLOSE, S2=OPEN   (energize at V=Vpeak!)\n");
    printf("    t=50ms:   S1=CLOSE, S2=CLOSE  (loaded operation)\n");
    printf("    t=150ms:  S1=CLOSE, S2=OPEN   (unload)\n");
    printf("    t=200ms:  S1=OPEN,  S2=OPEN   (de-energize)\n");
    printf("  Cycle 2:\n");
    printf("    t=207.5ms:S1=CLOSE, S2=OPEN   (re-energize at V=-Vpeak)\n");
    printf("    t=250ms:  S1=CLOSE, S2=CLOSE  (loaded operation)\n");
    printf("    t=350ms:  S1=CLOSE, S2=OPEN   (unload)\n");
    printf("    t=400ms:  S1=OPEN,  S2=OPEN   (de-energize)\n");
    printf("    t=500ms:  END\n");
    printf("\nNote: Source impedance Rs=1Ω limits inrush current\n");
    printf("      At V=Vpeak: flux starts at max, no DC offset, no saturation\n");
    printf("      At V=0: flux has DC offset, potential for saturation and inrush\n");

    /* ==================== Open CSV File ==================== */
    FILE* csv = fopen("transformer_transient.csv", "w");
    if (!csv) {
        printf("FAIL: Could not open output file\n");
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }

    /* CSV header */
    fprintf(csv, "time_ms,Vs,Vp,Vs_sec,Vload,I_primary,I_secondary,I_s2,Imag,Phi,B_core,P_core\n");
    fprintf(csv, "(ms),(V),(V),(V),(V),(A),(A),(A),(A),(Wb),(T),(W)\n");

    /* ==================== Transient Simulation ==================== */
    printf("\nRunning transient simulation...\n");

    double dt = 1e-5;  /* 10 μs time step (100kHz sampling) */
    double t_final = 0.5;  /* 500ms total (2 cycles) */
    int total_steps = (int)(t_final / dt);

    mna_init_transient(&solver);

    /* Set initial source value to 0 for DC initialization */
    solver.components[vs].value = 0.0;

    /* DC operating point */
    status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        printf("FAIL: DC analysis failed: %d\n", status);
        fclose(csv);
        mna_bh_curve_destroy(bh);
        mna_destroy(&solver);
        return -1;
    }

    printf("DC initialization complete.\n");
    printf("Starting transient...\n\n");

    /* Progress markers */
    int next_progress = 0;
    int progress_step = total_steps / 20;  /* 5% increments */

    for (int step = 0; step <= total_steps; step++) {
        double t = step * dt;
        double t_ms = t * 1000.0;

        /* Update voltage source: V = Vp * sin(2πft) */
        solver.components[vs].value = V_peak * sin(2 * PI * freq * t);

        /* Handle switching events (2 cycles) */
        int s1_state = 0, s2_state = 0;

        /* Cycle 1 - energize at voltage peak (t=7.5ms = 1/4 cycle) */
        if (t_ms >= 7.5 && t_ms < 200.0) s1_state = 1;   /* S1: 7.5-200ms */
        if (t_ms >= 50.0 && t_ms < 152.5) s2_state = 1;   /* S2: 50-150ms */

        /* Cycle 2 - energize at voltage negative peak (t=207.5ms = 1.25 cycles) */
        if (t_ms >= 207.5 && t_ms < 400.0) s1_state = 1;  /* S1: 207.5-400ms */
        if (t_ms >= 250.0 && t_ms < 350.0) s2_state = 1;  /* S2: 250-350ms */

        /* Update switch states (triggers state storage update) */
        if (solver.components[s1].state != s1_state) {
            mna_set_switch_state(&solver, s1, s1_state);
        }
        if (solver.components[s2].state != s2_state) {
            mna_set_switch_state(&solver, s2, s2_state);
        }

        /* Solve transient step */
        status = mna_solve_transient_step(&solver, dt);
        if (status != MNA_SUCCESS) {
            printf("FAIL: Transient step %d failed: %d\n", step, status);
            fclose(csv);
            mna_bh_curve_destroy(bh);
            mna_destroy(&solver);
            return -1;
        }

        /* Write data every 10 steps (100μs resolution in CSV) */
        if (step % 10 == 0) {
            /* Get voltages */
            double v_source = mna_get_node_voltage(&solver, n_source);
            double v_primary = mna_get_transformer_primary_voltage(&solver, xf);
            double v_secondary = mna_get_transformer_secondary_voltage(&solver, xf);
            double v_load = mna_get_node_voltage(&solver, n_load);

            /* Get currents using dedicated transformer API */
            double i_primary = mna_get_transformer_primary_current(&solver, xf);
            double i_secondary = mna_get_transformer_secondary_current(&solver, xf);
            /* Is2: current through S2 (load switch) - calculated from voltage drop */
            double i_s2 = (v_secondary - v_load) / (s2_state ? 0.001 : 1e9);
            double i_mag = mna_get_transformer_magnetizing_current(&solver, xf);

            /* Get core quantities */
            double phi = mna_get_transformer_flux(&solver, xf);
            double B = mna_get_transformer_B_field(&solver, xf);
            double P_core = mna_get_transformer_core_loss(&solver, xf);

            fprintf(csv, "%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.6f,%.8f,%.6f,%.6f\n",
                    t_ms, v_source, v_primary, v_secondary, v_load,
                    i_primary, i_secondary, i_s2, i_mag, phi, B, P_core);
        }

        /* Progress indicator */
        if (step >= next_progress) {
            printf("  Progress: %3d%% (t=%.1fms)\n",
                   (step * 100) / total_steps, t_ms);
            next_progress += progress_step;
        }
    }

    printf("  Progress: 100%% (t=%.1fms)\n", t_final * 1000.0);
    printf("\nSimulation complete.\n");

    fclose(csv);
    mna_bh_curve_destroy(bh);
    mna_destroy(&solver);

    printf("\nOutput written to: transformer_transient.csv\n");
    printf("========================================\n");
    printf("  Test PASSED\n");
    printf("========================================\n");

    return 0;
}
