/*
 * Minimal Test Runner - Runs only basic tests
 */

#include "../mna_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char** argv) {
    printf("=== CSPICE Minimal Tests ===\n\n");
    
    /* Test 1: Resistor voltage divider */
    printf("Test 1: Resistor voltage divider...\n");
    {
        MNASolver solver;
        if (mna_init(&solver) != MNA_SUCCESS) {
            printf("  FAILED: mna_init\n");
            return 1;
        }
        
        int node0 = 0;
        int node1 = mna_create_node(&solver);
        int node2 = mna_create_node(&solver);
        
        ComponentHandle v1, r1, r2;
        mna_add_voltage_source(&solver, node1, node0, 10.0, &v1);
        mna_add_resistor(&solver, node1, node2, 2000.0, &r1);
        mna_add_resistor(&solver, node2, node0, 3000.0, &r2);
        
        if (mna_solve_dc(&solver) != MNA_SUCCESS) {
            printf("  FAILED: mna_solve_dc\n");
            return 1;
        }
        
        double v2 = mna_get_node_voltage(&solver, node2);
        printf("  V(node2) = %.6f (expected: 6.0)\n", v2);
        
        if (fabs(v2 - 6.0) > 0.001) {
            printf("  FAILED: Wrong voltage\n");
            return 1;
        }
        
        mna_destroy(&solver);
        printf("  PASSED\n");
    }
    
    /* Test 2: Transformer creation */
    printf("\nTest 2: Transformer creation...\n");
    {
        MNASolver solver;
        if (mna_init(&solver) != MNA_SUCCESS) {
            printf("  FAILED: mna_init\n");
            return 1;
        }
        
        int node0 = 0;
        int node1 = mna_create_node(&solver);
        int node2 = mna_create_node(&solver);
        
        ComponentHandle xfmr;
        MNAStatus status = mna_add_ideal_transformer(&solver, node1, node0, node2, node0, 2.0, &xfmr);
        if (status != MNA_SUCCESS) {
            printf("  FAILED: mna_add_ideal_transformer returned %d\n", status);
            mna_destroy(&solver);
            return 1;
        }
        printf("  Created transformer with handle %d\n", xfmr);
        
        double ratio;
        status = mna_get_transformer_turns_ratio(&solver, xfmr, &ratio);
        if (status != MNA_SUCCESS) {
            printf("  FAILED: mna_get_transformer_turns_ratio returned %d\n", status);
            mna_destroy(&solver);
            return 1;
        }
        printf("  Turns ratio = %.6f (expected: 2.0)\n", ratio);
        
        if (fabs(ratio - 2.0) > 0.001) {
            printf("  FAILED: Wrong ratio\n");
            mna_destroy(&solver);
            return 1;
        }
        
        mna_destroy(&solver);
        printf("  PASSED\n");
    }
    
    /* Test 3: AC analysis */
    printf("\nTest 3: AC analysis...\n");
    {
        MNASolver solver;
        if (mna_init(&solver) != MNA_SUCCESS) {
            printf("  FAILED: mna_init\n");
            return 1;
        }
        
        int node0 = 0;
        int node1 = mna_create_node(&solver);
        int node2 = mna_create_node(&solver);
        
        ComponentHandle v1, r1, c1;
        mna_add_voltage_source(&solver, node1, node0, 0.0, &v1);
        mna_set_ac_source(&solver, v1, 1.0, 0.0);
        mna_add_resistor(&solver, node1, node2, 1000.0, &r1);
        mna_add_capacitor(&solver, node2, node0, 1e-6, &c1);
        
        MNAStatus status = mna_solve_ac(&solver, 2.0 * M_PI * 159.15);
        if (status != MNA_SUCCESS) {
            printf("  FAILED: mna_solve_ac returned %d\n", status);
            mna_destroy(&solver);
            return 1;
        }
        
        double complex v_out = mna_get_ac_node_voltage(&solver, node2);
        double v_mag = cabs(v_out);
        printf("  |V(out)| = %.6f (expected: ~0.707)\n", v_mag);
        
        if (fabs(v_mag - 0.707) > 0.05) {
            printf("  FAILED: Wrong magnitude\n");
            mna_destroy(&solver);
            return 1;
        }
        
        mna_destroy(&solver);
        printf("  PASSED\n");
    }
    
    printf("\n=== All tests PASSED ===\n");
    return 0;
}
