#include "passive.h"
#include "../include/matrix.h"
#include <stdlib.h>

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle);
MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2);

MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2,
                           double value, ComponentHandle* handle) {
    if (value <= 0.0) return MNA_INVALID_PARAMETER;
    return mna_add_component(solver, MNA_RESISTOR, node1, node2, value,
                             SOURCE_CURRENT, handle);
}

MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2,
                            double value, ComponentHandle* handle) {
    if (value <= 0.0) return MNA_INVALID_PARAMETER;
    return mna_add_component(solver, MNA_CAPACITOR, node1, node2, value,
                             SOURCE_CURRENT, handle);
}

MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2,
                           double value, ComponentHandle* handle) {
    if (value <= 0.0) return MNA_INVALID_PARAMETER;
    return mna_add_component(solver, MNA_INDUCTOR, node1, node2, value,
                             SOURCE_CURRENT, handle);
}

MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2,
                         double value, ComponentHandle* handle) {
    if (value <= 0.0) return MNA_INVALID_PARAMETER;
    return mna_add_component(solver, MNA_SWITCH, node1, node2, value,
                             SOURCE_CURRENT, handle);
}

MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return MNA_INVALID_HANDLE;
    }

    Component* comp = &solver->components[handle];
    if (comp->type == MNA_SWITCH) {
        comp->state = state;

        for (int i = 0; i < solver->num_components; i++) {
            Component* c = &solver->components[i];
            if (c->type == MNA_CAPACITOR || c->type == MNA_INDUCTOR ||
                c->type == MNA_CUSTOM_NONLINEAR) {
                c->prev_voltage = c->last_voltage;
                c->prev_current = c->last_current;
                c->stage1_voltage = c->last_voltage;
                c->stage1_current = c->last_current;
            }
        }
        return MNA_SUCCESS;
    }
    return MNA_INVALID_PARAMETER;
}
