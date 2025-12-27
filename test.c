#include "mna_solver_v2_1.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Forward declaration of the ideal op-amp stamping function
void ideal_opamp_stamp(MNASolver* solver,
                      const int* nodes,
                      int num_nodes,
                      void* user_data,
                      double time,
                      double dt);

// Simple voltage-controlled voltage source (VCVS) as a 4-terminal element
void vcvs_stamp(MNASolver* solver,
                const int* nodes,
                int num_nodes,
                void* user_data,
                double time,
                double dt) {
    if (num_nodes < 4) return;

    int in_plus = nodes[0];   // Input non-inverting
    int in_minus = nodes[1];  // Input inverting
    int out_plus = nodes[2];  // Output positive
    int out_minus = nodes[3]; // Output negative (usually ground)

    double gain = *((double*)user_data);

    // Add equation for Vout = gain * (Vin+ - Vin-)
    if (out_plus > 0) {
        // Output node equation
        MAT(solver, out_plus-1, out_plus-1) += 1.0;

        if (in_plus > 0) MAT(solver, out_plus-1, in_plus-1) -= gain;
        if (in_minus > 0) MAT(solver, out_plus-1, in_minus-1) += gain;

        // Add small conductance for numerical stability
        MAT(solver, out_plus-1, out_plus-1) += 1e-9;
    }

    // Add small conductances to input nodes to prevent floating
    if (in_plus > 0) MAT(solver, in_plus-1, in_plus-1) += 1e-12;
    if (in_minus > 0) MAT(solver, in_minus-1, in_minus-1) += 1e-12;
}

int main() {
    MNASolver solver;
    MNAStatus status;

    printf("== Multi-Domain Circuit Analysis Test ==\n");
    printf("This program tests both DC and AC solvers on the same circuit\n\n");

    printf("Initializing MNA solver...\n");
    status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("Failed to initialize solver: %d\n", status);
        return 1;
    }

    printf("Creating nodes...\n");
    int gnd = 0;            // Ground (node 0)
    int vin = mna_create_node(&solver);       // Input voltage
    int vout = mna_create_node(&solver);      // Output voltage
    int op_in_plus = mna_create_node(&solver);
    int op_in_minus = mna_create_node(&solver);

    printf("Adding components for an active low-pass filter...\n");
    // Input voltage source (1V for DC, 1V magnitude for AC)
    mna_add_voltage_source(&solver, vin, gnd, 1.0, NULL);
    mna_set_ac_source(&solver, solver.num_components-1, 1.0, 0.0); // 1V magnitude, 0 phase

    // Non-inverting amplifier with gain = 2
    mna_add_resistor(&solver, op_in_minus, gnd, 1000.0, NULL); // R1 = 1k
    mna_add_resistor(&solver, op_in_minus, vout, 1000.0, NULL); // R2 = 1k

    // RC low-pass filter at input
    mna_add_resistor(&solver, vin, op_in_plus, 1000.0, NULL); // R3 = 1k
    mna_add_capacitor(&solver, op_in_plus, gnd, 1e-6, NULL);   // C1 = 1uF

    // Add custom op-amp (ideal) as 3-terminal element
    double gain = 1e6;
    int opamp_nodes[3] = {op_in_plus, op_in_minus, vout};
    status = mna_add_custom_n_pole(&solver, opamp_nodes, 3, ideal_opamp_stamp, &gain, NULL);
    if (status != MNA_SUCCESS) {
        printf("Failed to add op-amp: %d\n", status);
        mna_destroy(&solver);
        return 1;
    }

    printf("\n=== DC Analysis ===\n");
    printf("Solving DC operating point...\n");
    status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        printf("DC solve failed: %d\n", status);
    } else {
        double v_in = mna_get_node_voltage(&solver, vin);
        double v_out = mna_get_node_voltage(&solver, vout);
        double v_plus = mna_get_node_voltage(&solver, op_in_plus);
        double v_minus = mna_get_node_voltage(&solver, op_in_minus);

        printf("DC Results:\n");
        printf("  Input voltage (Vin):  %.6f V\n", v_in);
        printf("  Output voltage (Vout): %.6f V\n", v_out);
        printf("  Op-amp non-inverting:  %.6f V\n", v_plus);
        printf("  Op-amp inverting:      %.6f V\n", v_minus);

        // For DC, capacitor is open circuit, so the gain should be 2 (non-inverting amp)
        double expected_output = 2.0;
        double error = fabs(v_out - expected_output) / expected_output * 100.0;
        printf("  Expected output: %.6f V, Error: %.2f%%\n", expected_output, error);

        if (error < 5.0) {
            printf("  DC analysis: PASSED\n");
        } else {
            printf("  DC analysis: FAILED\n");
        }
    }

    printf("\n=== AC Analysis ===\n");
    printf("Analyzing frequency response (10Hz to 100kHz)...\n");

    // Open file for AC results
    FILE* ac_file = fopen("ac_analysis_results.csv", "w");
    if (!ac_file) {
        printf("Failed to create AC analysis output file\n");
        mna_destroy(&solver);
        return 1;
    }

    fprintf(ac_file, "frequency,vout_magnitude,vout_phase\n");

    // Logarithmic frequency sweep
    double start_freq = 10.0;    // 10 Hz
    double end_freq = 100e3;    // 100 kHz
    int num_points = 50;

    printf("Frequency (Hz)\tMagnitude (V)\tPhase (deg)\n");
    printf("-------------------------------------------\n");

    for (int i = 0; i < num_points; i++) {
        // Logarithmic spacing
        double freq = start_freq * pow(end_freq/start_freq, (double)i/(num_points-1));

        status = mna_solve_ac(&solver, freq);
        if (status != MNA_SUCCESS) {
            printf("AC solve failed at %.1f Hz: %d\n", freq, status);
            continue;
        }

        double complex vout_ac = mna_get_ac_node_voltage(&solver, vout);
        double magnitude = cabs(vout_ac);
        double phase = carg(vout_ac) * 180.0 / M_PI;

        // Write to file
        fprintf(ac_file, "%.6e,%.6e,%.6e\n", freq, magnitude, phase);

        // Print key frequencies
        if (i == 0 || i == num_points/4 || i == num_points/2 || i == 3*num_points/4 || i == num_points-1) {
            printf("%9.1f\t%10.6f\t%8.2f\n", freq, magnitude, phase);
        }
    }

    fclose(ac_file);
    printf("\nAC analysis results saved to ac_analysis_results.csv\n");
    printf("Plot this file to see the low-pass filter response\n");

    printf("\n=== Additional Test: Second-Order Filter with Custom Element ===\n");
    printf("Creating a Sallen-Key low-pass filter circuit...\n");

    // Create new solver for second circuit
    MNASolver solver2;
    status = mna_init(&solver2);
    if (status != MNA_SUCCESS) {
        printf("Failed to initialize second solver: %d\n", status);
        mna_destroy(&solver);
        return 1;
    }

    // Nodes for Sallen-Key filter
    int gnd2 = 0;
    int vin2 = mna_create_node(&solver2);
    int vout2 = mna_create_node(&solver2);
    int node_a = mna_create_node(&solver2);
    int node_b = mna_create_node(&solver2);

    // Voltage source
    mna_add_voltage_source(&solver2, vin2, gnd2, 1.0, NULL);
    mna_set_ac_source(&solver2, solver2.num_components-1, 1.0, 0.0);

    // Sallen-Key components
    mna_add_resistor(&solver2, vin2, node_a, 10e3, NULL);  // R1 = 10k
    mna_add_capacitor(&solver2, node_a, node_b, 10e-9, NULL); // C1 = 10nF
    mna_add_resistor(&solver2, node_b, gnd2, 10e3, NULL);  // R2 = 10k
    mna_add_capacitor(&solver2, node_b, vout2, 10e-9, NULL); // C2 = 10nF

    // VCVS (voltage-controlled voltage source) for unity gain buffer
    double vcvs_gain = 1.0;
    int vcvs_nodes[4] = {node_b, gnd2, vout2, gnd2}; // in+, in-, out+, out-
    status = mna_add_custom_n_pole(&solver2, vcvs_nodes, 4, vcvs_stamp, &vcvs_gain, NULL);
    if (status != MNA_SUCCESS) {
        printf("Failed to add VCVS: %d\n", status);
        mna_destroy(&solver);
        mna_destroy(&solver2);
        return 1;
    }

    printf("\nSolving DC for Sallen-Key filter...\n");
    status = mna_solve_dc(&solver2);
    if (status != MNA_SUCCESS) {
        printf("DC solve failed for Sallen-Key: %d\n", status);
    } else {
        double vout_dc = mna_get_node_voltage(&solver2, vout2);
        printf("DC output voltage: %.6f V (should be 1V for unity gain)\n", vout_dc);
    }

    printf("\nAnalyzing AC frequency response...\n");
    FILE* ac_file2 = fopen("sallen_key_ac_results.csv", "w");
    if (!ac_file2) {
        printf("Failed to create Sallen-Key AC output file\n");
        mna_destroy(&solver);
        mna_destroy(&solver2);
        return 1;
    }

    fprintf(ac_file2, "frequency,vout_magnitude,vout_phase\n");

    // Calculate theoretical cutoff frequency for this Sallen-Key filter
    // fc = 1/(2*pi*R*C) for equal R and C values
    double r_val = 10e3;
    double c_val = 10e-9;
    double fc_theoretical = 1.0 / (2.0 * M_PI * r_val * c_val);
    printf("Theoretical cutoff frequency: %.2f Hz\n", fc_theoretical);

    printf("\nFrequency (Hz)\tMagnitude (V)\tPhase (deg)\n");
    printf("-------------------------------------------\n");

    for (int i = 0; i < num_points; i++) {
        double freq = start_freq * pow(end_freq/start_freq, (double)i/(num_points-1));

        status = mna_solve_ac(&solver2, freq);
        if (status != MNA_SUCCESS) {
            printf("AC solve failed at %.1f Hz for Sallen-Key: %d\n", freq, status);
            continue;
        }

        double complex vout_ac = mna_get_ac_node_voltage(&solver2, vout2);
        double magnitude = cabs(vout_ac);
        double phase = carg(vout_ac) * 180.0 / M_PI;

        fprintf(ac_file2, "%.6e,%.6e,%.6e\n", freq, magnitude, phase);

        if (i == 0 || i == num_points/4 || i == num_points/2 || i == 3*num_points/4 || i == num_points-1) {
            printf("%9.1f\t%10.6f\t%8.2f\n", freq, magnitude, phase);
        }
    }

    fclose(ac_file2);
    printf("\nSallen-Key AC analysis results saved to sallen_key_ac_results.csv\n");

    // Clean up
    printf("\nCleaning up...\n");
    mna_destroy(&solver);
    mna_destroy(&solver2);

    return 0;
}

// Simple ideal op-amp stamping function for testing
void ideal_opamp_stamp(MNASolver* solver,
                      const int* nodes,
                      int num_nodes,
                      void* user_data,
                      double time,
                      double dt) {
    if (num_nodes != 3) {
        fprintf(stderr, "Error: Op-amp requires exactly 3 terminals\n");
        return;
    }

    double gain = *((double*)user_data);
    int in_plus = nodes[0];   // Non-inverting input
    int in_minus = nodes[1];  // Inverting input
    int output = nodes[2];    // Output

    // Stamp output equation: Vout = gain * (V+ - V-)
    if (output > 0) {
        MAT(solver, output-1, output-1) += 1.0;

        if (in_plus > 0) MAT(solver, output-1, in_plus-1) -= gain;
        if (in_minus > 0) MAT(solver, output-1, in_minus-1) += gain;

        // Add small conductance for stability
        MAT(solver, output-1, output-1) += 1e-9;
    }

    // Add very small conductances to input nodes
    if (in_plus > 0) MAT(solver, in_plus-1, in_plus-1) += 1e-12;
    if (in_minus > 0) MAT(solver, in_minus-1, in_minus-1) += 1e-12;
}
