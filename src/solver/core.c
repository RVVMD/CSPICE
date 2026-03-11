#include "core.h"
#include "../include/matrix.h"
#include <stdlib.h>

MNAStatus mna_init(MNASolver* solver) {
    if (!solver) return MNA_INVALID_PARAMETER;

    memset(solver, 0, sizeof(MNASolver));
    solver->cap_components = MNA_INIT_CAPACITY;
    solver->matrix_cap_size = MNA_INIT_CAPACITY;
    solver->integration_method = MNA_INTEGRATION_TRBDF2;
    solver->preserve_dc_state = 0;

    solver->components = (Component*)calloc((size_t)solver->cap_components, sizeof(Component));
    solver->A = (double*)calloc((size_t)solver->matrix_cap_size * (size_t)solver->matrix_cap_size, sizeof(double));
    solver->b = (double*)calloc((size_t)solver->matrix_cap_size, sizeof(double));
    solver->x = (double*)calloc((size_t)solver->matrix_cap_size, sizeof(double));
    solver->ac_solution = (double complex*)calloc((size_t)solver->matrix_cap_size, sizeof(double complex));

    if (!solver->components || !solver->A || !solver->b || !solver->x || !solver->ac_solution) {
        mna_destroy(solver);
        return MNA_INSUFFICIENT_MEMORY;
    }

    return MNA_SUCCESS;
}

void mna_destroy(MNASolver* solver) {
    if (!solver) return;

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
            NPoleData* npole = comp->data.npole.npole_data;
            free(npole->nodes);
            free(npole->last_values);
            free(npole->branch_current_indices);
            free(npole);
        }
    }

    free(solver->components);
    free(solver->A);
    free(solver->b);
    free(solver->x);
    free(solver->ac_solution);

    memset(solver, 0, sizeof(MNASolver));
}

MNAStatus mna_set_integration_method(MNASolver* solver, IntegrationMethod method) {
    if (!solver) return MNA_INVALID_PARAMETER;
    solver->integration_method = method;
    return MNA_SUCCESS;
}

MNAStatus mna_set_preserve_dc_state(MNASolver* solver, int enable) {
    if (!solver) return MNA_INVALID_PARAMETER;
    solver->preserve_dc_state = (enable != 0);
    return MNA_SUCCESS;
}

int mna_create_node(MNASolver* solver) {
    if (!solver) return -1;

    int new_index = solver->max_node_index + 1;
    if (!mna_resize_matrix(solver, new_index + solver->num_sources)) {
        return -1;
    }

    solver->max_node_index = new_index;
    solver->num_nodes = (solver->num_nodes < new_index) ? new_index : solver->num_nodes;
    return new_index;
}

MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2) {
    if (!solver) return MNA_INVALID_HANDLE;
    if (node1 < 0 || node1 > solver->max_node_index ||
        node2 < 0 || node2 > solver->max_node_index) {
        return MNA_INVALID_NODE;
    }
    return MNA_SUCCESS;
}

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;

    MNAStatus status = mna_validate_nodes(solver, node1, node2);
    if (status != MNA_SUCCESS) return status;

    int new_count = solver->num_components + 1;
    if (!mna_resize_components(solver, new_count)) {
        return MNA_INSUFFICIENT_MEMORY;
    }

    Component comp;
    memset(&comp, 0, sizeof(comp));
    comp.type = type;
    comp.node1 = node1;
    comp.node2 = node2;
    comp.value = value;
    comp.state = 1;
    comp.source_type = src_type;
    comp.smoothed_alpha = 1.0;

    int index = solver->num_components++;
    solver->components[index] = comp;

    if (type == MNA_SOURCE && src_type == SOURCE_VOLTAGE) {
        if (!mna_resize_matrix(solver, solver->max_node_index + solver->num_sources + 1)) {
            return MNA_INSUFFICIENT_MEMORY;
        }
        solver->num_sources++;
    }

    if (handle) *handle = index;
    return MNA_SUCCESS;
}

double mna_get_node_voltage(MNASolver* solver, int node) {
    if (!solver || node <= 0 || node > solver->max_node_index) return 0.0;
    return solver->x[node-1];
}

double complex mna_get_ac_node_voltage(MNASolver* solver, int node) {
    if (!solver || node <= 0 || node > solver->max_node_index) return 0.0;
    return solver->ac_solution[node-1];
}

double mna_get_component_current(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    int n1 = comp->node1;
    int n2 = comp->node2;
    double v1 = (n1 > 0) ? mna_get_node_voltage(solver, n1) : 0.0;
    double v2 = (n2 > 0) ? mna_get_node_voltage(solver, n2) : 0.0;
    double v = v1 - v2;

    switch (comp->type) {
        case MNA_RESISTOR:
            return v / comp->value;

        case MNA_CAPACITOR:
            if (!solver->transient_initialized || solver->dt == 0.0) return 0.0;
            return comp->last_current;

        case MNA_INDUCTOR:
            return comp->last_current;

        case MNA_SOURCE:
            if (comp->source_type == SOURCE_VOLTAGE) {
                int vs_index = 0;
                for (int i = 0; i < handle; i++) {
                    if (solver->components[i].type == MNA_SOURCE &&
                        solver->components[i].source_type == SOURCE_VOLTAGE) {
                        vs_index++;
                    }
                }
                return solver->x[solver->max_node_index + vs_index];
            } else {
                return comp->value;
            }

        case MNA_SWITCH:
            return v / (comp->state ? comp->value : 1.0/MNA_MIN_CONDUCTANCE);

        case MNA_CUSTOM_NONLINEAR: {
            ComponentState state = {
                .voltage = v,
                .current = comp->last_current,
                .charge = comp->last_charge,
                .flux = comp->last_flux,
                .dt = solver->dt
            };
            double current, conductance;
            comp->nonlinear_func(&state, comp->user_data, &current, &conductance);
            return current;
        }

        case MNA_CUSTOM_NPOLE:
            return 0.0;

        default:
            return 0.0;
    }
}

double mna_get_component_voltage(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    int n1 = comp->node1;
    int n2 = comp->node2;
    double v1 = (n1 > 0) ? mna_get_node_voltage(solver, n1) : 0.0;
    double v2 = (n2 > 0) ? mna_get_node_voltage(solver, n2) : 0.0;
    return v1 - v2;
}
