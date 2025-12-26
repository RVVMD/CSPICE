#include "mna_solver_v2_1.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Custom nonlinear comparator with hysteresis (Schmitt trigger)
void schmitt_trigger_stamp(MNASolver* solver,
                          const int* nodes,
                          int num_nodes,
                          void* user_data,
                          double time,
                          double dt) {
    if (num_nodes < 3) return; // Need at least in+, in-, out

    int in_plus = nodes[0];   // Non-inverting input
    int in_minus = nodes[1];  // Inverting input (threshold)
    int output = nodes[2];    // Output

    double* thresholds = (double*)user_data;
    double upper_thresh = thresholds[0];
    double lower_thresh = thresholds[1];
    double output_high = thresholds[2];
    double output_low = thresholds[3];

    // Get current output state (use last_value as state memory)
    double current_out = (output > 0 && solver->x) ? solver->x[output-1] : 0.0;
    int state = (current_out > (output_high + output_low) / 2) ? 1 : 0;

    // Get input voltages
    double v_plus = (in_plus > 0) ? solver->x[in_plus-1] : 0.0;
    double v_minus = (in_minus > 0) ? solver->x[in_minus-1] : 0.0;

    // Determine next state based on hysteresis
    // For the matrix stamping, we use a large gain to approximate the switching behavior
    const double gain = 1e6;

    if (output > 0) {
        if (state == 1) {
            // Currently high - switch when input goes below lower threshold
            if (v_plus < lower_thresh) {
                state = 0;
            }
            MAT(solver, output-1, output-1) = 1.0;
            if (in_plus > 0) MAT(solver, output-1, in_plus-1) = -gain;
            solver->b[output-1] = -gain * lower_thresh;
        } else {
            // Currently low - switch when input goes above upper threshold
            if (v_plus > upper_thresh) {
                state = 1;
            }
            MAT(solver, output-1, output-1) = 1.0;
            if (in_plus > 0) MAT(solver, output-1, in_plus-1) = -gain;
            solver->b[output-1] = -gain * upper_thresh;
        }
    }

    // Add small conductances to prevent floating nodes
    if (in_plus > 0) MAT(solver, in_plus-1, in_plus-1) += 1e-12;
    if (in_minus > 0) MAT(solver, in_minus-1, in_minus-1) += 1e-12;
    if (output > 0) MAT(solver, output-1, output-1) += 1e-9;
}

// 555 Timer Model as a custom n-pole element (8-pin package)
// Pins: 1-GND, 2-TRIG, 3-OUT, 4-RESET, 5-CTRL, 6-THRES, 7-DISCH, 8-VCC
void timer_555_stamp(MNASolver* solver,
                    const int* nodes,
                    int num_nodes,
                    void* user_data,
                    double time,
                    double dt) {
    if (num_nodes < 8) return;

    int gnd = nodes[0];      // Pin 1 - Ground
    int trig = nodes[1];     // Pin 2 - Trigger
    int out = nodes[2];      // Pin 3 - Output
    int reset = nodes[3];    // Pin 4 - Reset (active low)
    int ctrl = nodes[4];     // Pin 5 - Control voltage
    int thres = nodes[5];    // Pin 6 - Threshold
    int disch = nodes[6];    // Pin 7 - Discharge
    int vcc = nodes[7];      // Pin 8 - VCC

    double vcc_voltage = 5.0; // Standard 5V supply

    // Internal thresholds (2/3 Vcc and 1/3 Vcc, or based on ctrl pin)
    double threshold_upper = (ctrl > 0 && solver->x) ?
                            solver->x[ctrl-1] : (2.0/3.0) * vcc_voltage;
    double threshold_lower = threshold_upper / 2.0;

    // Get current voltages
    double v_trig = (trig > 0 && solver->x) ? solver->x[trig-1] : 0.0;
    double v_thres = (thres > 0 && solver->x) ? solver->x[thres-1] : 0.0;
    double v_reset = (reset > 0 && solver->x) ? solver->x[reset-1] : vcc_voltage; // Active low

    // Internal flip-flop state (approximation)
    static int ff_state = 1; // Start with output high

    // Flip-flop logic (simplified)
    if (v_reset < 0.7) { // Reset active (below ~0.7V)
        ff_state = 0;
    } else {
        if (v_trig < threshold_lower) {
            ff_state = 1; // Set
        }
        if (v_thres > threshold_upper) {
            ff_state = 0; // Reset
        }
    }

    // Output stage (push-pull)
    if (out > 0) {
        if (ff_state == 1) {
            // Output high (with small resistance for stability)
            MAT(solver, out-1, out-1) += 1e3; // 1k conductance
            solver->b[out-1] += 1e3 * vcc_voltage;
        } else {
            // Output low (connect to ground)
            MAT(solver, out-1, out-1) += 1e3; // 1k conductance
            solver->b[out-1] += 0.0;
        }
    }

    // Discharge transistor (open collector)
    if (disch > 0) {
        if (ff_state == 0) {
            // Discharge transistor ON (connect to ground)
            MAT(solver, disch-1, disch-1) += 1e2; // 100 S (10mΩ)
            solver->b[disch-1] += 0.0;
        } else {
            // Discharge transistor OFF (high impedance)
            MAT(solver, disch-1, disch-1) += 1e-9; // Very small conductance
        }
    }

    // Control pin input impedance
    if (ctrl > 0) {
        MAT(solver, ctrl-1, ctrl-1) += 1e-8; // 100MΩ input impedance
    }

    // Prevent floating nodes with small conductances
    if (trig > 0) MAT(solver, trig-1, trig-1) += 1e-12;
    if (thres > 0) MAT(solver, thres-1, thres-1) += 1e-12;
    if (reset > 0) MAT(solver, reset-1, reset-1) += 1e-12;
}

// Phase-Locked Loop (PLL) test circuit with voltage-controlled oscillator
int main() {
    MNASolver solver;
    MNAStatus status;

    printf("Initializing MNA solver for complex PLL circuit...\n");
    status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        printf("Failed to initialize solver: %d\n", status);
        return 1;
    }

    printf("Creating nodes for PLL circuit...\n");
    int gnd = 0;             // Ground (node 0)
    int vcc = mna_create_node(&solver);        // +5V power supply
    int vdd = mna_create_node(&solver);        // +12V for VCO

    // Reference signal path
    int ref_in = mna_create_node(&solver);     // Reference input
    int phase_detector_out = mna_create_node(&solver);

    // VCO control and output
    int vco_control = mna_create_node(&solver);
    int vco_out = mna_create_node(&solver);

    // Filter nodes
    int filter_node1 = mna_create_node(&solver);
    int filter_node2 = mna_create_node(&solver);

    // 555 Timer nodes (8 pins as described in timer_555_stamp)
    int timer_gnd = gnd;                     // Pin 1 - Ground
    int timer_trig = mna_create_node(&solver); // Pin 2 - Trigger
    int timer_out = mna_create_node(&solver);  // Pin 3 - Output
    int timer_reset = vcc;                   // Pin 4 - Reset (tied to VCC)
    int timer_ctrl = mna_create_node(&solver); // Pin 5 - Control
    int timer_thres = mna_create_node(&solver); // Pin 6 - Threshold
    int timer_disch = mna_create_node(&solver); // Pin 7 - Discharge
    int timer_vcc = vcc;                     // Pin 8 - VCC

    // Schmitt trigger nodes
    int schmitt_in_plus = mna_create_node(&solver);
    int schmitt_in_minus = mna_create_node(&solver);
    int schmitt_out = mna_create_node(&solver);

    printf("Adding power supplies...\n");
    mna_add_voltage_source(&solver, vcc, gnd, 5.0, NULL);   // 5V supply
    mna_add_voltage_source(&solver, vdd, gnd, 12.0, NULL);  // 12V supply

    printf("Adding reference signal source...\n");
    // 1kHz square wave reference signal (we'll change it during simulation)
    mna_add_voltage_source(&solver, ref_in, gnd, 0.0, NULL);

    printf("Adding passive components...\n");
    // Phase detector - simple XOR implementation with resistors
    mna_add_resistor(&solver, ref_in, phase_detector_out, 10e3, NULL); // 10k
    mna_add_resistor(&solver, vco_out, phase_detector_out, 10e3, NULL); // 10k

    // Loop filter (active low-pass)
    mna_add_resistor(&solver, phase_detector_out, filter_node1, 10e3, NULL); // 10k
    mna_add_capacitor(&solver, filter_node1, gnd, 100e-9, NULL); // 100nF
    mna_add_resistor(&solver, filter_node1, filter_node2, 1e3, NULL); // 1k
    mna_add_capacitor(&solver, filter_node2, gnd, 1e-6, NULL); // 1uF

    // VCO control components
    mna_add_resistor(&solver, filter_node2, vco_control, 10e3, NULL); // 10k

    // 555 Timer external components for VCO
    mna_add_resistor(&solver, vcc, timer_disch, 10e3, NULL); // R1 = 10k
    mna_add_resistor(&solver, timer_disch, timer_thres, 100e3, NULL); // R2 = 100k (variable)
    mna_add_capacitor(&solver, timer_thres, gnd, 10e-9, NULL); // C = 10nF

    // Connect 555 output to VCO output
    mna_add_resistor(&solver, timer_out, vco_out, 100, NULL); // Small resistor for isolation

    // Feedback divider (divide by 4)
    mna_add_resistor(&solver, vco_out, schmitt_in_plus, 1e6, NULL); // High impedance buffer

    // Configure Schmitt trigger thresholds (1/3 and 2/3 of VCC)
    double thresholds[4] = {3.33, 1.67, 5.0, 0.0}; // Upper, lower, high output, low output

    printf("Adding custom n-pole elements...\n");
    // Add 555 timer as 8-terminal element
    int timer_nodes[8] = {timer_gnd, timer_trig, timer_out, timer_reset,
                         timer_ctrl, timer_thres, timer_disch, timer_vcc};
    status = mna_add_custom_n_pole(&solver, timer_nodes, 8, timer_555_stamp, NULL, NULL);
    if (status != MNA_SUCCESS) {
        printf("Failed to add 555 timer: %d\n", status);
        mna_destroy(&solver);
        return 1;
    }

    // Add Schmitt trigger as comparator for frequency division
    int schmitt_nodes[3] = {schmitt_in_plus, schmitt_in_minus, schmitt_out};
    status = mna_add_custom_n_pole(&solver, schmitt_nodes, 3, schmitt_trigger_stamp, thresholds, NULL);
    if (status != MNA_SUCCESS) {
        printf("Failed to add Schmitt trigger: %d\n", status);
        mna_destroy(&solver);
        return 1;
    }

    // Connect Schmitt trigger reference voltage (midpoint)
    mna_add_voltage_source(&solver, schmitt_in_minus, gnd, 2.5, NULL); // 2.5V reference

    // Connect output to feedback path (divide by 2)
    mna_add_capacitor(&solver, schmitt_out, timer_trig, 100e-12, NULL); // Small coupling cap
    mna_add_resistor(&solver, timer_trig, gnd, 10e3, NULL); // Pull-down

    // Control voltage pin connection
    mna_add_resistor(&solver, vco_control, timer_ctrl, 1e6, NULL); // High impedance buffer
    mna_add_capacitor(&solver, timer_ctrl, gnd, 10e-9, NULL); // Filtering cap

    printf("Initializing transient analysis...\n");
    mna_init_transient(&solver);

    printf("Running PLL transient simulation...\n");
    double t = 0.0;
    double dt = 1e-7; // 100ns timestep
    double duration = 10e-3; // 10ms total duration (many cycles)
    int steps = (int)(duration / dt);
    int print_interval = (int)(0.001 / dt); // Print every 1ms

    // Open file for results
    FILE* output_file = fopen("pll_transient_results.csv", "w");
    if (!output_file) {
        printf("Failed to create output file\n");
        mna_destroy(&solver);
        return 1;
    }

    fprintf(output_file, "time,ref_input,vco_output,control_voltage,filter_voltage\n");

    // Find reference source handle
    ComponentHandle ref_source_handle = -1;
    for (int i = 0; i < solver.num_components; i++) {
        if (solver.components[i].node1 == ref_in && solver.components[i].node2 == gnd) {
            ref_source_handle = i;
            break;
        }
    }

    int success_count = 0;
    int fail_count = 0;

    // Simulation start time for performance measurement
    clock_t start_time = clock();

    for (int step = 0; step <= steps; step++) {
        // Generate 1kHz square wave reference (period = 1ms)
        double ref_period = 0.001; // 1ms = 1kHz
        double ref_voltage = ((fmod(t, ref_period) < ref_period/2) ? 5.0 : 0.0);
        solver.components[ref_source_handle].value = ref_voltage;

        // Solve this time step
        status = mna_solve_transient_step(&solver, dt);

        if (status == MNA_SUCCESS) {
            success_count++;
            double ref_val = mna_get_node_voltage(&solver, ref_in);
            double vco_val = mna_get_node_voltage(&solver, vco_out);
            double ctrl_val = mna_get_node_voltage(&solver, vco_control);
            double filter_val = mna_get_node_voltage(&solver, filter_node2);

            // Write to file (less frequently to keep file size manageable)
            if (step % 10 == 0) {
                fprintf(output_file, "%.6e,%.6e,%.6e,%.6e,%.6e\n",
                       t, ref_val, vco_val, ctrl_val, filter_val);
            }

            // Print progress periodically
            if (step % print_interval == 0) {
                printf("t = %.3f ms: Ref = %.3f V, VCO = %.3f V, Ctrl = %.3f V\n",
                       t * 1000.0, ref_val, vco_val, ctrl_val);

                // Check for convergence issues
                if (fabs(ctrl_val) > 100.0) {
                    printf("Warning: Control voltage is unusually high (%.3f V)\n", ctrl_val);
                }
            }
        } else {
            fail_count++;
            printf("Transient step failed at t = %.3f ms (step %d): status %d\n",
                   t * 1000.0, step, status);

            if (fail_count > 5) {
                printf("Too many failures, aborting simulation\n");
                break;
            }

            // Try with smaller timestep on failure
            dt /= 2.0;
            step--; // Retry this step
            continue;
        }

        t += dt;
    }

    clock_t end_time = clock();
    double sim_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    fclose(output_file);
    printf("\nPLL transient simulation completed:\n");
    printf("  Simulation time: %.2f seconds\n", sim_time);
    printf("  Circuit time simulated: %.3f ms\n", t * 1000.0);
    printf("  Successful steps: %d\n", success_count);
    printf("  Failed steps: %d\n", fail_count);
    printf("  Results saved to pll_transient_results.csv\n");

    // Print final voltages
    double final_ref = mna_get_node_voltage(&solver, ref_in);
    double final_vco = mna_get_node_voltage(&solver, vco_out);
    double final_ctrl = mna_get_node_voltage(&solver, vco_control);
    double final_filter = mna_get_node_voltage(&solver, filter_node2);

    printf("\nFinal steady-state voltages:\n");
    printf("  Reference input: %.4f V\n", final_ref);
    printf("  VCO output: %.4f V\n", final_vco);
    printf("  Control voltage: %.4f V\n", final_ctrl);
    printf("  Filter voltage: %.4f V\n", final_filter);

    printf("\nCleaning up...\n");
    mna_destroy(&solver);
    printf("Test completed!\n");

    return 0;
}
