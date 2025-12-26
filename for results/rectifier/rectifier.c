#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mna_solver.h"

// Diode parameters for 1N456
#define DIODE_IS 2.5e-8      // Saturation current (A) - calculated from leakage current (25nA)
#define DIODE_N 1.0          // Ideality factor (unchanged, typical for Si diodes)
#define DIODE_MAX_EXP 700.0  // Unchanged (prevents overflow)
#define DIODE_MIN_EXP -300.0 // Unchanged (reverse bias)
#define DIODE_VJ 0.7         // Junction potential (V) - standard for Si
#define DIODE_M 0.5          // Grading coefficient (unchanged)
#define DIODE_CJ0 1e-12      // Zero-bias junction capacitance (F) - reasonable estimate

// Circuit parameters
#define V_PEAK 10.0        // AC source peak voltage (V)
#define FREQ 60.0          // AC source frequency (Hz)
#define R_LOAD 1000.0      // Load resistance (Ω)
#define C_FILTER 1e-4      // Filter capacitance (F)
#define TOTAL_CYCLES 5     // Number of AC cycles to simulate
#define NUM_POINTS 10000  // Number of simulation points

// Diode nonlinear function (updated for capacitive model)
void diode_func(MNASolver* solver, int comp_index, double voltage, double current,
                double* value1, double* value2, NonlinearType nl_type) {
    (void)solver;         // Unused
    (void)comp_index;     // Unused
    (void)current;        // Only used for inductor

    double Vt = MNA_VT;   // Thermal voltage

    if (nl_type == NONLINEAR_RESISTOR) {
        // Resistive part (standard diode equation)
        double x = voltage / (DIODE_N * Vt);

        // Clamp exponent to avoid numerical issues
        if (x > DIODE_MAX_EXP) {
            double exp_x = exp(DIODE_MAX_EXP);
            *value1 = DIODE_IS * (exp_x - 1.0);  // current
            *value2 = DIODE_IS / (DIODE_N * Vt) * exp_x; // conductance
        } else if (x < DIODE_MIN_EXP) {
            *value1 = -DIODE_IS; // current
            *value2 = DIODE_IS / (DIODE_N * Vt) * exp(DIODE_MIN_EXP); // conductance
        } else {
            double exp_x = exp(x);
            *value1 = DIODE_IS * (exp_x - 1.0); // current
            *value2 = DIODE_IS / (DIODE_N * Vt) * exp_x; // conductance
        }
    } else if (nl_type == NONLINEAR_CAPACITOR) {
        // Capacitive part (junction capacitance)
        // Clamp voltage to avoid singularity
        double v_clamped = (voltage < DIODE_VJ - 0.001) ? voltage : DIODE_VJ - 0.001;

        // Compute junction capacitance
        double Cj = DIODE_CJ0 / pow(1 - v_clamped/DIODE_VJ, DIODE_M);

        // Compute charge at this voltage
        *value1 = (DIODE_CJ0 * DIODE_VJ) / (1-DIODE_M) *
                  (1 - pow(1 - v_clamped/DIODE_VJ, 1-DIODE_M));
        *value2 = Cj;  // Capacitance
    }
}

int main() {
    MNASolver solver;
    mna_init(&solver);

    // Arrays to store component indices for diodes
    int diodeR_indices[4];   // Resistor part indices
    int diodeC_indices[4];   // Capacitor part indices

    // Add components to the circuit
    // AC Voltage Source (between nodes 1 and 2)
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, 1, 2, 0.0);

    // Diode Bridge (D1-D4) - now with capacitive effects
    // D1: between 1 (anode) and 3 (cathode)
    diodeR_indices[0] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 1, 3, NONLINEAR_RESISTOR,
                            diode_func, NULL, 0.0, 0.0);
    diodeC_indices[0] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 1, 3, NONLINEAR_CAPACITOR,
                            diode_func, NULL, 0.0, 0.0);

    // D2: between 2 (anode) and 3 (cathode)
    diodeR_indices[1] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 2, 3, NONLINEAR_RESISTOR,
                            diode_func, NULL, 0.0, 0.0);
    diodeC_indices[1] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 2, 3, NONLINEAR_CAPACITOR,
                            diode_func, NULL, 0.0, 0.0);

    // D3: between 0 (GND, anode) and 1 (cathode)
    diodeR_indices[2] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 0, 1, NONLINEAR_RESISTOR,
                            diode_func, NULL, 0.0, 0.0);
    diodeC_indices[2] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 0, 1, NONLINEAR_CAPACITOR,
                            diode_func, NULL, 0.0, 0.0);

    // D4: between 0 (GND, anode) and 2 (cathode)
    diodeR_indices[3] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 0, 2, NONLINEAR_RESISTOR,
                            diode_func, NULL, 0.0, 0.0);
    diodeC_indices[3] = solver.num_components;
    mna_add_custom_nonlinear(&solver, 0, 2, NONLINEAR_CAPACITOR,
                            diode_func, NULL, 0.0, 0.0);

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
    fprintf(fp, "time,V1,V2,V3,I_R,I_C,I_D1_res,I_D1_cap,I_D2_res,I_D2_cap,I_D3_res,I_D3_cap,I_D4_res,I_D4_cap\n");

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

        // Calculate diode currents
        double I_D1_res, I_D1_cap, I_D2_res, I_D2_cap;
        double I_D3_res, I_D3_cap, I_D4_res, I_D4_cap;
        double g_dummy;

        // Get resistive currents
        diode_func(&solver, 0, V1 - V3, 0, &I_D1_res, &g_dummy, NONLINEAR_RESISTOR);
        diode_func(&solver, 0, V2 - V3, 0, &I_D2_res, &g_dummy, NONLINEAR_RESISTOR);
        diode_func(&solver, 0, -V1, 0, &I_D3_res, &g_dummy, NONLINEAR_RESISTOR);
        diode_func(&solver, 0, -V2, 0, &I_D4_res, &g_dummy, NONLINEAR_RESISTOR);

        // Get capacitive currents from companion models
        Component* cap_comp;
        cap_comp = &solver.components[diodeC_indices[0]];
        I_D1_cap = cap_comp->trans_G_eq * (V1 - V3) + cap_comp->trans_I_eq;

        cap_comp = &solver.components[diodeC_indices[1]];
        I_D2_cap = cap_comp->trans_G_eq * (V2 - V3) + cap_comp->trans_I_eq;

        cap_comp = &solver.components[diodeC_indices[2]];
        I_D3_cap = cap_comp->trans_G_eq * (-V1) + cap_comp->trans_I_eq;

        cap_comp = &solver.components[diodeC_indices[3]];
        I_D4_cap = cap_comp->trans_G_eq * (-V2) + cap_comp->trans_I_eq;

        // Write results to CSV
        fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                t, V1, V2, V3, I_R, I_C,
                I_D1_res, I_D1_cap, I_D2_res, I_D2_cap,
                I_D3_res, I_D3_cap, I_D4_res, I_D4_cap);

        t += dt; // Advance time
    }

    fclose(fp);
    printf("Simulation completed. Results saved to rectifier.csv\n");
    return 0;
}
