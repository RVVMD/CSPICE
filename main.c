#include <stdio.h>
#include <math.h>
#include "mna_solver.h"

// Nonlinear element wrappers with sign adjustments
void nonlinear_wrapper1(MNASolver* solver, int comp_index, double voltage, double current,
                        double* value1, double* value2, NonlinearType nl_type) {
    double V_physical = voltage;
    double In = 0.316 * (1 - exp(-0.0995 * V_physical));
    double g_physical = 0.316 * 0.0995 * exp(-0.0995 * V_physical);
    *value1 = In;
    *value2 = g_physical;
}

void nonlinear_wrapper2(MNASolver* solver, int comp_index, double voltage, double current,
                        double* value1, double* value2, NonlinearType nl_type) {
    double In = 0.025 * (exp(0.085 * voltage) - 0.025);
    double g = 0.025 * 0.085 * exp(0.085 * voltage);
    *value1 = In;
    *value2 = g;
}

void nonlinear_wrapper3(MNASolver* solver, int comp_index, double voltage, double current,
                        double* value1, double* value2, NonlinearType nl_type) {
    double In = 0.02467 * voltage - 0.001 * pow(voltage, 2) + 0.00001333 * pow(voltage, 3);
    double g = 0.02467 - 0.002 * voltage + 0.00003999 * pow(voltage, 2);
    *value1 = In;
    *value2 = g;
}

int main() {
    MNASolver solver;
    mna_init(&solver);

    // Create intermediate nodes (2,3,4) for voltage sources
    int node_a = 1;  // Main node a
    int node_b = 0;  // Ground (node b)
    int node2 = 2;   // Between E1 and NE1
    int node3 = 3;   // Between E2 and NE2
    int node4 = 4;   // Between E3 and NE3

    // Add components to the circuit
    // Branch 1: E1 (10V) and NE1
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, node_b, node2, 10.0);
    mna_add_custom_nonlinear(&solver, node_a, node2, NONLINEAR_RESISTOR,
                             nonlinear_wrapper1, NULL, 0.0, 0.0);

    // Branch 2: E2 (10V) and NE2
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, node3, node_b, 10.0);
    mna_add_custom_nonlinear(&solver, node_a, node3, NONLINEAR_RESISTOR,
                             nonlinear_wrapper2, NULL, 0.0, 0.0);

    // Branch 3: E3 (40V) and NE3
    mna_add_component(&solver, MNA_VOLTAGE_SOURCE, node_b, node4, 40.0);
    mna_add_custom_nonlinear(&solver, node_a, node4, NONLINEAR_RESISTOR,
                             nonlinear_wrapper3, NULL, 0.0, 0.0);

    // Solve DC operating point
    if (mna_solve_dc(&solver)) {
        printf("\nDC Solution Converged:\n");
        printf("---------------------------------\n");

        // Print node voltages
        printf("Node Voltages:\n");
        printf("Ua (node 1) = %.6f V\n", mna_get_node_voltage(&solver, node_a));
        printf("U2 (node 2) = %.6f V\n", mna_get_node_voltage(&solver, node2));
        printf("U3 (node 3) = %.6f V\n", mna_get_node_voltage(&solver, node3));
        printf("U4 (node 4) = %.6f V\n", mna_get_node_voltage(&solver, node4));

        // Print branch currents
        printf("\nBranch Currents:\n");
        int ne_count = 0;
        for (int i = 0; i < solver.num_components; i++) {
            Component* comp = &solver.components[i];
            if (comp->type == MNA_CUSTOM_NONLINEAR) {
                double current, conductance;
                comp->nonlinear_func(&solver, i, comp->last_voltage, 0,
                                    &current, &conductance, comp->nonlinear_type);
                printf("I_NE%d = %.6f A\n", ne_count + 1, current);
                ne_count++;
            }
        }

        // Calculate and print voltage source currents
        printf("\nVoltage Source Currents:\n");
        int vs_count = 0;
        for (int i = 0; i < solver.num_components; i++) {
            if (solver.components[i].type == MNA_VOLTAGE_SOURCE) {
                int idx = solver.max_node_index + vs_count;
                printf("I_E%d = %.6f A\n", vs_count + 1, solver.x[idx]);
                vs_count++;
            }
        }

        printf("---------------------------------\n");
    } else {
        printf("DC Solution Failed to Converge\n");
    }

    // After solving DC:
    double total_source_power = 0;
    double total_resistor_power = 0;

    // Calculate source powers
    int vs_count = 0;
    for (int i = 0; i < solver.num_components; i++) {
        if (solver.components[i].type == MNA_VOLTAGE_SOURCE) {
            int idx = solver.max_node_index + vs_count;
            double Is = solver.x[idx];
            double V = solver.components[i].value;
            double P = V * Is;
            total_source_power += P;
            vs_count++;
        }
    }

    // Calculate resistor powers
    for (int i = 0; i < solver.num_components; i++) {
        Component* comp = &solver.components[i];
        if (comp->type == MNA_CUSTOM_NONLINEAR) {
            double Ir, G;
            comp->nonlinear_func(&solver, i, comp->last_voltage, 0, &Ir, &G, comp->nonlinear_type);
            double P = comp->last_voltage * Ir;
            total_resistor_power += P;
        }
    }

    printf("Total source power: %.6f W\n", total_source_power);
    printf("Total resistor power: %.6f W\n", total_resistor_power);
    printf("Power balance: %.6f W\n", total_source_power + total_resistor_power);
    return 0;
}
