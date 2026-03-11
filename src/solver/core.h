#ifndef MNA_SOLVER_CORE_H
#define MNA_SOLVER_CORE_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_init(MNASolver* solver);
void mna_destroy(MNASolver* solver);

MNAStatus mna_set_integration_method(MNASolver* solver, IntegrationMethod method);
MNAStatus mna_set_preserve_dc_state(MNASolver* solver, int enable);

int mna_create_node(MNASolver* solver);

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle);

MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2);

double mna_get_node_voltage(MNASolver* solver, int node);
double complex mna_get_ac_node_voltage(MNASolver* solver, int node);
double mna_get_component_current(MNASolver* solver, ComponentHandle handle);
double mna_get_component_voltage(MNASolver* solver, ComponentHandle handle);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_CORE_H */
