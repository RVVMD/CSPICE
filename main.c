#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "mna_solver.h"

// Diode parameters
#define DIODE_IS 1e-12  // Saturation current (A)
#define DIODE_N 1.0     // Ideality factor
#define DIODE_VT 0.02585 // Thermal voltage

void diode_func(MNASolver* solver, int comp_index,
                double vd, double* current, double* conductance) {
    // Improved diode model with reverse bias protection
    const double I_s = DIODE_IS;
    const double n = DIODE_N;
    const double V_T = DIODE_VT;

    if (vd < -3 * n * V_T) {
        // Reverse bias approximation
        *current = -I_s;
        *conductance = MNA_MIN_CONDUCTANCE;
    } else if (vd > 0.7) {
        // Forward bias approximation
        *current = I_s * exp(vd/(n*V_T));
        *conductance = *current / (n*V_T);
    } else {
        // Full exponential model
        *current = I_s * (exp(vd/(n*V_T)) - 1);
        *conductance = I_s/(n*V_T) * exp(vd/(n*V_T));
    }
}

int main() {
    MNASolver solver;
    mna_init(&solver);

    // Create circuit:
    // DC Source (1V) + AC Source (1V) -> Resistor (1k) -> Diode -> Ground
    // Nodes: 0=Ground, 1=Source, 2=DiodeAnode

    // Add components
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, 0, 1, 1.0);  // DC=1V
    mna_add_component(&solver, MNA_RESISTOR, 1, 2, 1000);       // 1k resistor
    int diode_idx = solver.num_components;  // Store diode index
    mna_add_custom_nonlinear(&solver, 2, 0, diode_func, NULL, 0.0);  // Better initial guess

    // Configure AC source (1V magnitude, 0° phase)
    mna_set_ac_source(&solver, 0, 1.0, 0.0);

    // Solve DC operating point
    printf("Solving DC operating point...\n");
    int dc_success = mna_solve_dc(&solver);

    // Print DC results
    printf("\nDC Solution %s:\n", dc_success ? "SUCCESS" : "FAILED");
    printf("Node 1 (Source) Voltage: %.3f V\n", mna_get_node_voltage(&solver, 1));
    printf("Node 2 (Diode Anode) Voltage: %.3f V\n", mna_get_node_voltage(&solver, 2));

    // Calculate actual diode current
    double vd = solver.components[diode_idx].last_voltage;
    double diode_current, diode_cond;
    diode_func(&solver, diode_idx, vd, &diode_current, &diode_cond);
    printf("Diode Voltage: %.3f V\n", vd);
    printf("Diode Current: %.3f mA\n", diode_current * 1000);
    printf("Diode Conductance: %.3f S\n\n", diode_cond);

    // Perform AC analysis at different frequencies
    double frequencies[] = {10, 100, 1000, 10000, 100000};
    int num_freqs = sizeof(frequencies)/sizeof(frequencies[0]);

    printf("AC Analysis Results (Magnitude at Node 2):\n");
    printf("Frequency(Hz)\t|Vout|(V)\tPhase(°)\n");
    printf("----------------------------------------\n");

    for (int i = 0; i < num_freqs; i++) {
        double freq = frequencies[i];
        mna_solve_ac(&solver, freq);

        double complex v_out = mna_get_ac_node_voltage(&solver, 2);
        double magnitude = cabs(v_out);
        double phase = carg(v_out) * 180 / M_PI;

        printf("%.0f\t\t%.6f\t\t%.1f\n", freq, magnitude, phase);
    }

    return 0;
}
