#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mna_solver.h"

// Diode parameters
#define DIODE_IS 1e-12     // Saturation current (A)
#define DIODE_N 1.0        // Ideality factor
#define DIODE_MAX_EXP 700.0 // Max exponent to avoid overflow
#define DIODE_MIN_EXP -300.0 // Min exponent for reverse bias

// Circuit parameters
#define V_PEAK 10.0        // AC source peak voltage (V)
#define FREQ 60.0          // AC source frequency (Hz)
#define R_LOAD 1000.0      // Load resistance (Ω)
#define C_FILTER 1e-3      // Filter capacitance (F)
#define TOTAL_CYCLES 5     // Number of AC cycles to simulate
#define NUM_POINTS 10000   // Number of simulation points

// Diode nonlinear function
void diode_func(MNASolver* solver, int comp_index, double voltage,
                double* current, double* conductance) {
    (void)solver;         // Unused
    (void)comp_index;     // Unused

    double Vt = MNA_VT;   // Thermal voltage
    double x = voltage / (DIODE_N * Vt);

    // Clamp exponent to avoid numerical issues
    if (x > DIODE_MAX_EXP) {
        double exp_x = exp(DIODE_MAX_EXP);
        *current = DIODE_IS * (exp_x - 1.0);
        *conductance = DIODE_IS / (DIODE_N * Vt) * exp_x;
    } else if (x < DIODE_MIN_EXP) {
        *current = -DIODE_IS;
        *conductance = DIODE_IS / (DIODE_N * Vt) * exp(DIODE_MIN_EXP);
    } else {
        double exp_x = exp(x);
        *current = DIODE_IS * (exp_x - 1.0);
        *conductance = DIODE_IS / (DIODE_N * Vt) * exp_x;
    }
}

int main() {
    MNASolver solver;
    mna_init(&solver);

    // Add components to the circuit
    // AC Voltage Source (between nodes 1 and 2)
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, 1, 2, 0.0);

    // Diode Bridge (D1-D4)
    mna_add_custom_nonlinear(&solver, 1, 3, diode_func, NULL, 0.0); // D1: anode=1, cathode=3
    mna_add_custom_nonlinear(&solver, 2, 3, diode_func, NULL, 0.0); // D2: anode=2, cathode=3
    mna_add_custom_nonlinear(&solver, 0, 1, diode_func, NULL, 0.0); // D3: anode=0 (GND), cathode=1
    mna_add_custom_nonlinear(&solver, 0, 2, diode_func, NULL, 0.0); // D4: anode=0 (GND), cathode=2

    // Load Resistor (between node 3 and GND)
    mna_add_component(&solver, MNA_RESISTOR, 3, 0, R_LOAD);

    // Filter Capacitor (between node 3 and GND)
    mna_add_component(&solver, MNA_CAPACITOR, 3, 0, C_FILTER);

    // Initialize transient analysis
    mna_init_transient(&solver);

    // Open CSV file for writing results
    FILE *fp = fopen("rectifier.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }
    fprintf(fp, "time,V1,V2,V3,I_R,I_C,I_D1,I_D2,I_D3,I_D4\n");

    // Simulation parameters
    double total_time = TOTAL_CYCLES / FREQ;
    double dt = total_time / NUM_POINTS;
    double t = 0.0;
    double V3_prev = 0.0; // Previous capacitor voltage (for current calculation)

    // Run transient simulation
    for (int step = 0; step < NUM_POINTS; step++) {
        // Update AC source: V = V_peak * sin(2*pi*f*t)
        solver.components[0].value = V_PEAK * sin(2 * M_PI * FREQ * t);

        // Solve one time step
        if (!mna_solve_transient_step(&solver, dt)) {
            fprintf(stderr, "Solver failed at step %d\n", step);
            break;
        }

        // Get node voltages
        double V1 = mna_get_node_voltage(&solver, 1);
        double V2 = mna_get_node_voltage(&solver, 2);
        double V3 = mna_get_node_voltage(&solver, 3);

        // Calculate load resistor current (Ohm's Law)
        double I_R = V3 / R_LOAD;

        // Calculate capacitor current: I_C = C * dV/dt ≈ C * (V3 - V3_prev)/dt
        double I_C = C_FILTER * (V3 - V3_prev) / dt;
        V3_prev = V3; // Update for next step

        // Calculate diode currents using diode function
        double I_D1, I_D2, I_D3, I_D4, g_dummy;
        diode_func(&solver, 0, V1 - V3, &I_D1, &g_dummy); // D1: V1-V3
        diode_func(&solver, 0, V2 - V3, &I_D2, &g_dummy); // D2: V2-V3
        diode_func(&solver, 0, -V1, &I_D3, &g_dummy);     // D3: -V1
        diode_func(&solver, 0, -V2, &I_D4, &g_dummy);     // D4: -V2

        // Write results to CSV
        fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                t, V1, V2, V3, I_R, I_C, I_D1, I_D2, I_D3, I_D4);

        t += dt; // Advance time
    }

    fclose(fp);
    printf("Simulation completed. Results saved to rectifier.csv\n");
    return 0;
}
