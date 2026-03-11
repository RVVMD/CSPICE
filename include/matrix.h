#ifndef MNA_MATRIX_H
#define MNA_MATRIX_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAT(solver, i, j) ((solver)->A[(i) * (solver)->matrix_cap_size + (j)])

bool mna_resize_components(MNASolver* solver, int required);
bool mna_resize_matrix(MNASolver* solver, int req_size);

void mna_reset_system(MNASolver* solver);
MNAStatus mna_solve_linear_system(MNASolver* solver, int size);

void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g);
void mna_stamp_current_source(MNASolver* solver, int node1, int node2, double current_val);
void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx);
void mna_ensure_ground_paths(MNASolver* solver);

static inline int mna_active_size(const MNASolver* solver) {
    return solver->max_node_index + solver->num_sources;
}

#ifdef __cplusplus
}
#endif

#endif /* MNA_MATRIX_H */
