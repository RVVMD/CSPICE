/**
 * Simple test - just voltage source, switch and resistor
 */

#include "mna_solver.h"
#include <stdio.h>

int main(void) {
    MNASolver* solver = (MNASolver*)calloc(1, sizeof(MNASolver));
    mna_init(solver);

    int n1 = mna_create_node(solver);
    int n2 = mna_create_node(solver);
    int gnd = 0;

    /* VS -- switch -- node -- resistor -- gnd */
    ComponentHandle vs_h, sw_h, r_h;
    mna_add_voltage_source(solver, n1, gnd, 10.0, &vs_h);
    mna_add_switch(solver, n1, n2, 0.001, &sw_h);
    mna_set_switch_state(solver, sw_h, 0);
    mna_add_resistor(solver, n2, gnd, 10.0, &r_h);

    printf("Testing simple circuit...\n");
    MNAStatus status = mna_solve_dc(solver);
    printf("DC result: %d\n", status);

    mna_init_transient(solver);
    
    FILE* f = fopen("simple_test.csv", "w");
    fprintf(f, "time,V_n1,V_n2,I_R\n");
    
    for (int i = 0; i < 100; i++) {
        if (solver->time >= 0.01) mna_set_switch_state(solver, sw_h, 1);
        
        double v1 = mna_get_node_voltage(solver, n1);
        double v2 = mna_get_node_voltage(solver, n2);
        double ir = mna_get_component_current(solver, r_h);
        fprintf(f, "%.6e,%.6e,%.6e,%.6e\n", solver->time, v1, v2, ir);
        
        status = mna_solve_transient_step(solver, 1e-4);
        if (status != MNA_SUCCESS) {
            printf("Transient failed at step %d: %d\n", i, status);
            break;
        }
    }
    fclose(f);
    printf("Simple test complete. Output: simple_test.csv\n");

    mna_destroy(solver);
    free(solver);
    return 0;
}
