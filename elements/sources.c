#include "sources.h"
#include "../include/matrix.h"
#include <stdlib.h>

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle);

MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2,
                                 double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SOURCE, node1, node2, value,
                             SOURCE_VOLTAGE, handle);
}

MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2,
                                 double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SOURCE, node1, node2, value,
                             SOURCE_CURRENT, handle);
}

MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle,
                            double magnitude, double phase) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return MNA_INVALID_HANDLE;
    }
    solver->components[handle].ac_magnitude = magnitude;
    solver->components[handle].ac_phase = phase;
    return MNA_SUCCESS;
}
