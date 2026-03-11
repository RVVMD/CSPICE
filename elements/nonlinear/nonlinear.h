#ifndef MNA_ELEMENTS_NONLINEAR_H
#define MNA_ELEMENTS_NONLINEAR_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                                   NonlinearType nl_type,
                                   CustomNonlinearFunc func, void* user_data,
                                   double initial_value1, double initial_value2,
                                   ComponentHandle* handle);

MNAStatus mna_add_custom_n_pole(MNASolver* solver, const int* nodes, int num_nodes,
                                NPoleStampFunc stamp_func, void* user_data,
                                int num_branch_currents, ComponentHandle* handle);

double mna_compute_alpha(double rate, double prev_alpha);

void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc, int stage);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_NONLINEAR_H */
