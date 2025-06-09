#include <stdio.h>
#include <stdlib.h>
#include "cspice.h"

void print_results(double t, double *node_voltages) {
    // Print time, input AC voltage, and DC output voltage
    printf("%.6f, %.6f, %.6f\n", t, node_voltages[1], node_voltages[2]);
}

int main() {
    // Initialize circuit with proper parameters
    Circuit circuit;
    const double rel_tol = 1e-6;       // Relative tolerance
    const double abs_tol = 1e-9;       // Absolute tolerance
    const int max_iter = 100;          // Maximum Newton iterations
    const int init_comp_cap = 16;      // Initial component capacity
    const int init_node_cap = 8;       // Initial node capacity

    init_circuit(&circuit, rel_tol, abs_tol, max_iter, init_comp_cap, init_node_cap);

    // Set simulation parameters
    const double time_step = 10e-6;    // 10μs time step for 60Hz signal
    set_fixed_time_step(&circuit, time_step);

    // Circuit parameters
    const double ac_amplitude = 10.0;  // 10V peak (20V peak-to-peak)
    const double ac_frequency = 60.0;  // 60Hz
    const double cap_value = 100e-6;   // 100μF smoothing capacitor
    const double load_resistance = 1000.0; // 1kΩ load

    // Diode parameters (1N4007-like)
    const double diode_Is = 1e-14;     // Saturation current
    const double diode_N = 1.6;        // Emission coefficient

    // Create full-wave bridge rectifier circuit
    // Nodes: 0=GND, 1=AC input, 2=Output, 3-6=bridge internal nodes

    // AC voltage source (10V peak, 60Hz, 0V DC offset)
    double ac_params[] = {0.0, ac_amplitude, ac_frequency, 0.0}; // DC, Ampl, Freq, Phase
    if (!add_voltage_source(&circuit, "VAC", 1, 0, SINE_SOURCE, ac_params)) {
        fprintf(stderr, "Failed to add AC voltage source\n");
        free_circuit(&circuit);
        return EXIT_FAILURE;
    }

    // Diode bridge (D1-D4)
    if (!add_diode(&circuit, "D1", 1, 3, diode_Is, diode_N) ||  // AC+ to bridge top
        !add_diode(&circuit, "D2", 4, 1, diode_Is, diode_N) ||  // AC- to bridge bottom
        !add_diode(&circuit, "D3", 2, 4, diode_Is, diode_N) ||  // Output to bridge bottom
        !add_diode(&circuit, "D4", 3, 2, diode_Is, diode_N)) {  // Bridge top to output
        fprintf(stderr, "Failed to add diodes\n");
        free_circuit(&circuit);
        return EXIT_FAILURE;
    }

    // Load resistor
    if (!add_resistor(&circuit, "Rload", 2, 0, load_resistance)) {
        fprintf(stderr, "Failed to add load resistor\n");
        free_circuit(&circuit);
        return EXIT_FAILURE;
    }

    // Smoothing capacitor (initial voltage = 0V)
    if (!add_capacitor(&circuit, "Cfilter", 2, 0, cap_value, 0.0)) {
        fprintf(stderr, "Failed to add filter capacitor\n");
        free_circuit(&circuit);
        return EXIT_FAILURE;
    }

    // Run transient analysis for 5 AC cycles (~83.33ms)
    const double simulation_time = 5.0 / ac_frequency; // 5 cycles
    printf("Time (s), AC Input (V), DC Output (V)\n");
    transient_analysis(&circuit, simulation_time, print_results);

    // Clean up
    free_circuit(&circuit);
    return EXIT_SUCCESS;
}
