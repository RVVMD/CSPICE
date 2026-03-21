/**
 * Test without coupled elements - simple circuit
 */

#include "mna_solver.h"
#include <stdio.h>

int main(void) {
    MNASolver* solver = (MNASolver*)calloc(1, sizeof(MNASolver));
    mna_init(solver);

    int n_pri = mna_create_node(solver);
    int n_sec = mna_create_node(solver);
    int gnd = 0;

    /* Voltage source on primary */
    ComponentHandle vs_h;
    mna_add_voltage_source(solver, n_pri, gnd, 10.0, &vs_h);

    /* Inductor on primary */
    ComponentHandle L_h;
    mna_add_inductor(solver, n_pri, gnd, 1.0, &L_h);

    /* Load on secondary (not connected to anything else - floating) */
    ComponentHandle load_h;
    mna_add_resistor(solver, n_sec, gnd, 10.0, &load_h);

    printf("Testing VS + L + floating resistor...\n");
    MNAStatus status = mna_solve_dc(solver);
    printf("DC result: %d\n", status);
    
    if (status == MNA_SUCCESS) {
        printf("V_pri = %.6f V\n", mna_get_node_voltage(solver, n_pri));
        printf("V_sec = %.6f V\n", mna_get_node_voltage(solver, n_sec));
    }

    mna_destroy(solver);
    free(solver);
    return 0;
}
