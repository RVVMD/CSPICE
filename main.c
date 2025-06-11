#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mna_solver.h"

#define DIODE_IS 1e-14   // More realistic saturation current
#define DIODE_N   1.7    // Typical emission coefficient

void diode_model(MNASolver* solver, int comp_index,
                 double vd, double* Ic, double* G) {
    // Use the global VT constant from MNA solver
    double vt = MNA_VT * DIODE_N;

    double exponent = vd / vt;
    // Clamp exponent to avoid overflow
    exponent = (exponent > 700) ? 700 : (exponent < -700) ? -700 : exponent;

    *Ic = DIODE_IS * (exp(exponent) - 1.0);
    *G = DIODE_IS / vt * exp(exponent);

    // Clamp conductance to valid range
    *G = (*G < MNA_MIN_CONDUCTANCE) ? MNA_MIN_CONDUCTANCE :
         (*G > MNA_MAX_CONDUCTANCE) ? MNA_MAX_CONDUCTANCE : *G;
}

int main() {
    MNASolver solver;
    mna_init(&solver);
    double last_vd = 0.0;  // Declare and initialize here

    // Create circuit: V1 -- R1 -- D1 -- GND
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, 0, 1, 0.0);  // V1
    mna_add_component(&solver, MNA_RESISTOR, 1, 2, 1000.0);     // R1 = 1kΩ
    mna_add_custom_nonlinear(&solver, 0, 2, diode_model, NULL, 0.0);  // Diode

    FILE *csv = fopen("diode_sweep.csv", "w");
    fprintf(csv, "Vin,Vd,Id\n");

    const int points = 100;
    for (int i = 0; i <= points; i++) {
        double vin = 5.0 * i / points;  // 0-1V sweep

        // Reset solver system
        mna_reset_system(&solver);
        solver.components[0].value = vin;

        // Start from last solution (except first iteration)
        if (i > 0) {
            solver.components[2].last_voltage = last_vd;
        }

        if (mna_solve_dc(&solver)) {
            double vd = solver.components[2].last_voltage;
            last_vd = vd;  // Store for next iteration

            double id = 0;
            double g = 0;
            diode_model(&solver, 2, vd, &id, &g);
            fprintf(csv, "%.6f,%.6f,%.6e\n", vin, vd, id);
        } else {
            fprintf(csv, "%.6f,nan,nan\n", vin);
        }
    }

    fclose(csv);
    return 0;
}
