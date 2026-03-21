/*
 * Test for ideal transformer n-pole element
 * 
 * Circuit:
 *   V1 (10V DC) connected to primary winding (nodes 1-0)
 *   Transformer with turns ratio n=2 (V1/V2 = 2)
 *   Load resistor R1 (100 ohm) connected to secondary winding (nodes 2-0)
 * 
 * Expected results:
 *   V_primary = 10V
 *   V_secondary = V_primary / n = 10/2 = 5V
 *   I_secondary = V_secondary / R1 = 5/100 = 0.05A = 50mA
 *   I_primary = -I_secondary / n = -0.05/2 = -0.025A = -25mA
 */

#include "mna_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    MNASolver* solver = (MNASolver*)malloc(sizeof(MNASolver));
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver\n");
        return 1;
    }

    MNAStatus status = mna_init(solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to initialize solver: %d\n", status);
        free(solver);
        return 1;
    }

    /* Create nodes: 0 is ground, create nodes 1 and 2 */
    int node1 = mna_create_node(solver);  /* Primary positive */
    int node2 = mna_create_node(solver);  /* Secondary positive */
    printf("Created nodes: primary=%d, secondary=%d\n", node1, node2);

    /* Add 10V DC voltage source on primary */
    ComponentHandle vs_handle;
    status = mna_add_voltage_source(solver, node1, 0, 10.0, &vs_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add voltage source: %d\n", status);
        mna_destroy(solver);
        free(solver);
        return 1;
    }

    /* Add ideal transformer: primary (node1-0), secondary (node2-0), turns ratio n=2 */
    ComponentHandle xf_handle;
    status = mna_add_ideal_transformer(solver, node1, 0, node2, 0, 2.0, &xf_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add transformer: %d\n", status);
        mna_destroy(solver);
        free(solver);
        return 1;
    }

    /* Add 100 ohm load resistor on secondary */
    ComponentHandle r_handle;
    status = mna_add_resistor(solver, node2, 0, 100.0, &r_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add resistor: %d\n", status);
        mna_destroy(solver);
        free(solver);
        return 1;
    }

    /* Solve DC */
    status = mna_solve_dc(solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "DC solve failed: %d\n", status);
        mna_destroy(solver);
        free(solver);
        return 1;
    }

    /* Get results */
    double v_primary = mna_get_node_voltage(solver, node1);
    double v_secondary = mna_get_node_voltage(solver, node2);
    double i_vs = mna_get_component_current(solver, vs_handle);
    double i_r = mna_get_component_current(solver, r_handle);

    printf("\n=== Results ===\n");
    printf("V_primary   = %.6f V (expected: 10.0 V)\n", v_primary);
    printf("V_secondary = %.6f V (expected: 5.0 V)\n", v_secondary);
    printf("I_primary   = %.6f A (expected: -0.025 A)\n", i_vs);
    printf("I_secondary = %.6f A (expected: 0.05 A)\n", i_r);

    /* Verify turns ratio */
    double ratio;
    mna_get_transformer_turns_ratio(solver, xf_handle, &ratio);
    printf("\nTransformer turns ratio: %.2f\n", ratio);

    /* Check accuracy */
    double v_sec_expected = v_primary / ratio;
    double error_v = fabs(v_secondary - v_sec_expected) / v_sec_expected * 100.0;
    
    printf("\nVoltage ratio error: %.4f%%\n", error_v);

    mna_destroy(solver);
    free(solver);

    if (error_v < 1.0) {
        printf("\nTEST PASSED\n");
        return 0;
    } else {
        printf("\nTEST FAILED\n");
        return 1;
    }
}
