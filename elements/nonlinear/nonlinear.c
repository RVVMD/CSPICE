#include "nonlinear.h"
#include "../../include/matrix.h"
#include <stdlib.h>

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle);
MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2);

MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                                   NonlinearType nl_type,
                                   CustomNonlinearFunc func, void* user_data,
                                   double initial_value1, double initial_value2,
                                   ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;

    MNAStatus status = mna_validate_nodes(solver, node1, node2);
    if (status != MNA_SUCCESS) return status;

    int new_count = solver->num_components + 1;
    if (!mna_resize_components(solver, new_count)) {
        return MNA_INSUFFICIENT_MEMORY;
    }

    Component comp;
    memset(&comp, 0, sizeof(comp));
    comp.type = MNA_CUSTOM_NONLINEAR;
    comp.node1 = node1;
    comp.node2 = node2;
    comp.is_nonlinear = true;
    comp.source_type = SOURCE_CURRENT;
    comp.nonlinear_type = nl_type;
    comp.nonlinear_func = func;
    comp.user_data = user_data;
    comp.smoothed_alpha = 1.0;

    if (nl_type == NONLINEAR_RESISTOR) {
        comp.last_voltage = initial_value1;
    } else if (nl_type == NONLINEAR_CAPACITOR) {
        comp.last_voltage = initial_value1;
        comp.last_charge = initial_value2;
    } else if (nl_type == NONLINEAR_INDUCTOR) {
        comp.last_current = initial_value1;
        comp.last_flux = initial_value2;
    }

    int index = solver->num_components++;
    solver->components[index] = comp;
    solver->num_nonlinear++;

    if (handle) *handle = index;
    return MNA_SUCCESS;
}

MNAStatus mna_add_custom_n_pole(MNASolver* solver, const int* nodes, int num_nodes,
                                NPoleStampFunc stamp_func, void* user_data,
                                int num_branch_currents, ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;

    for (int i = 0; i < num_nodes; i++) {
        if (nodes[i] < 0 || nodes[i] > solver->max_node_index) {
            return MNA_INVALID_NODE;
        }
    }

    int new_count = solver->num_components + 1;
    if (!mna_resize_components(solver, new_count)) {
        return MNA_INSUFFICIENT_MEMORY;
    }

    if (num_branch_currents > 0) {
        if (!mna_resize_matrix(solver, solver->max_node_index + solver->num_sources + num_branch_currents)) {
            return MNA_INSUFFICIENT_MEMORY;
        }
        solver->num_sources += num_branch_currents;
    }

    Component comp;
    memset(&comp, 0, sizeof(comp));
    comp.type = MNA_CUSTOM_NPOLE;
    comp.node1 = (num_nodes > 0) ? nodes[0] : 0;
    comp.node2 = (num_nodes > 1) ? nodes[1] : 0;
    comp.is_nonlinear = true;

    comp.data.npole.npole_data = (NPoleData*)malloc(sizeof(NPoleData));
    if (!comp.data.npole.npole_data) return MNA_INSUFFICIENT_MEMORY;

    comp.data.npole.npole_data->nodes = (int*)malloc((size_t)num_nodes * sizeof(int));
    if (!comp.data.npole.npole_data->nodes) {
        free(comp.data.npole.npole_data);
        return MNA_INSUFFICIENT_MEMORY;
    }

    comp.data.npole.npole_data->last_values = (double*)calloc((size_t)num_nodes, sizeof(double));
    if (!comp.data.npole.npole_data->last_values) {
        free(comp.data.npole.npole_data->nodes);
        free(comp.data.npole.npole_data);
        return MNA_INSUFFICIENT_MEMORY;
    }

    comp.data.npole.npole_data->num_branch_currents = num_branch_currents;
    if (num_branch_currents > 0) {
        comp.data.npole.npole_data->branch_current_indices =
            (int*)malloc((size_t)num_branch_currents * sizeof(int));
        if (!comp.data.npole.npole_data->branch_current_indices) {
            free(comp.data.npole.npole_data->last_values);
            free(comp.data.npole.npole_data->nodes);
            free(comp.data.npole.npole_data);
            return MNA_INSUFFICIENT_MEMORY;
        }
        for (int i = 0; i < num_branch_currents; i++) {
            comp.data.npole.npole_data->branch_current_indices[i] =
                solver->max_node_index + solver->num_sources - num_branch_currents + i;
        }
    } else {
        comp.data.npole.npole_data->branch_current_indices = NULL;
    }

    memcpy(comp.data.npole.npole_data->nodes, nodes, (size_t)num_nodes * sizeof(int));
    comp.data.npole.npole_data->num_nodes = num_nodes;
    comp.data.npole.npole_data->stamp_func = stamp_func;
    comp.data.npole.npole_data->user_data = user_data;

    int index = solver->num_components++;
    solver->components[index] = comp;
    solver->num_nonlinear++;

    if (handle) *handle = index;
    return MNA_SUCCESS;
}

double mna_compute_alpha(double rate, double prev_alpha) {
    double target_alpha = 1.0;

    if (rate >= MNA_DAMPING_THRESHOLD) {
        double excess = rate - MNA_DAMPING_THRESHOLD;
        target_alpha = exp(-MNA_DAMPING_STEEPNESS * excess);
    }

    double recovery_factor = 0.5;

    if (target_alpha < prev_alpha) {
        return target_alpha;
    } else {
        return prev_alpha + recovery_factor * (target_alpha - prev_alpha);
    }
}

void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc, int stage) {
    Component* comp = &solver->components[comp_index];
    int n1 = comp->node1;
    int n2 = comp->node2;

    ComponentState state = {
        .voltage = comp->last_voltage,
        .current = comp->last_current,
        .charge = comp->last_charge,
        .flux = comp->last_flux,
        .dt = solver->dt
    };

    double alpha = comp->smoothed_alpha;

    if (!is_dc && solver->dt > 0.0) {
        double rate = 0.0;
        if (comp->nonlinear_type == NONLINEAR_CAPACITOR ||
            comp->nonlinear_type == NONLINEAR_RESISTOR) {
            double dv = comp->last_voltage - comp->prev_voltage;
            rate = fabs(dv) / solver->dt;
        } else if (comp->nonlinear_type == NONLINEAR_INDUCTOR) {
            double di = comp->last_current - comp->prev_current;
            rate = fabs(di) / solver->dt;
        }

        alpha = mna_compute_alpha(rate, comp->smoothed_alpha);
        comp->smoothed_alpha = alpha;
    }

    switch (comp->nonlinear_type) {
        case NONLINEAR_RESISTOR: {
            double current, conductance;
            comp->nonlinear_func(&state, comp->user_data, &current, &conductance);

            if (conductance < MNA_MIN_CONDUCTANCE) conductance = MNA_MIN_CONDUCTANCE;
            if (conductance > MNA_MAX_CONDUCTANCE) conductance = MNA_MAX_CONDUCTANCE;

            mna_stamp_conductance(solver, n1, n2, conductance);
            mna_stamp_current_source(solver, n1, n2, (current - conductance * state.voltage));
            break;
        }
        case NONLINEAR_CAPACITOR: {
            if (is_dc) {
                mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
            } else {
                double q0, C0;
                comp->nonlinear_func(&state, comp->user_data, &q0, &C0);

                double G_eq, I_eq;
                double dt = solver->dt;
                double v_n = comp->last_voltage;
                double gamma = MNA_TRBDF2_GAMMA;

                if (stage == 1) {
                    double h1 = gamma * dt;
                    G_eq = (2.0 * C0) / h1;
                    I_eq = (G_eq * v_n);
                } else {
                    double v_ng = comp->stage1_voltage;
                    G_eq = (C0 * (2.0 - gamma)) / (dt * (1.0 - gamma));
                    I_eq = (C0 / (dt * gamma * (1.0 - gamma))) * v_ng -
                           (C0 * (1.0 - gamma) / (dt * gamma)) * v_n;
                }

                comp->trans_G_eq = G_eq;
                comp->trans_I_eq = I_eq;

                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, -I_eq);
            }
            break;
        }
        case NONLINEAR_INDUCTOR: {
            if (is_dc) {
                mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
            } else {
                double phi0, L0;
                comp->nonlinear_func(&state, comp->user_data, &phi0, &L0);

                double G_eq, I_eq;
                double dt = solver->dt;
                double i_n = comp->last_current;
                double gamma = MNA_TRBDF2_GAMMA;

                if (stage == 1) {
                    double h1 = gamma * dt;
                    G_eq = h1 / (2.0 * L0);
                    I_eq = i_n;
                } else {
                    double i_ng = comp->stage1_current;
                    G_eq = (dt * (1.0 - gamma)) / (L0 * (2.0 - gamma));
                    I_eq = (1.0 / (gamma * (2.0 - gamma))) * i_ng -
                           (((1.0 - gamma) * (1.0 - gamma)) / (gamma * (2.0 - gamma))) * i_n;
                }

                comp->trans_G_eq = G_eq;
                comp->trans_I_eq = I_eq;

                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
            }
            break;
        }
    }
}
