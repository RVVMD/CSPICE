#include "passive.h"
#include "../include/matrix.h"
#include "../src/solver/core.h"
#include <stdlib.h>
#include <math.h>

/*
 * Ideal Transformer N-pole Element
 * 
 * A 4-terminal ideal transformer with turns ratio n = N1/N2:
 *   - Primary winding: nodes p1, p2
 *   - Secondary winding: nodes s1, s2
 *   - Turns ratio: n (V1/V2 = n, I1*n = -I2)
 * 
 * MNA formulation adds 1 branch current variable i_x for the transformer coupling.
 * Equations:
 *   KCL at p1: +i_x = 0
 *   KCL at p2: -i_x = 0
 *   KCL at s1: -n*i_x = 0
 *   KCL at s2: +n*i_x = 0
 *   Constraint: V_p1 - V_p2 - n*(V_s1 - V_s2) = 0
 */

typedef struct {
    double turns_ratio;  /* n = N_primary / N_secondary */
    int branch_var_index; /* Index of the branch current variable in the MNA matrix */
} TransformerData;

static void mna_stamp_ideal_transformer(MNASolver* solver,
                                        const int* nodes,
                                        int num_nodes,
                                        void* user_data,
                                        double time,
                                        double dt,
                                        int stage) {
    if (num_nodes != 4 || !user_data) return;

    TransformerData* tdata = (TransformerData*)user_data;
    double n = tdata->turns_ratio;
    int branch_idx = tdata->branch_var_index;

    int p1 = nodes[0];  /* Primary positive */
    int p2 = nodes[1];  /* Primary negative */
    int s1 = nodes[2];  /* Secondary positive */
    int s2 = nodes[3];  /* Secondary negative */

    (void)time;
    (void)dt;
    (void)stage;

    /* Stamp B matrix contributions (branch current into KCL equations)
     * Primary: i_x flows from p1 to p2
     * Secondary: -n*i_x flows from s1 to s2 (current direction reversed due to transformer action)
     */
    if (p1 > 0) {
        MAT(solver, p1 - 1, branch_idx - 1) += 1.0;
    }
    if (p2 > 0) {
        MAT(solver, p2 - 1, branch_idx - 1) += -1.0;
    }
    if (s1 > 0) {
        MAT(solver, s1 - 1, branch_idx - 1) += -n;
    }
    if (s2 > 0) {
        MAT(solver, s2 - 1, branch_idx - 1) += n;
    }

    /* Stamp constraint equation row: V_p1 - V_p2 - n*V_s1 + n*V_s2 = 0
     * This goes in row branch_idx-1, columns for node voltages
     */
    if (p1 > 0) MAT(solver, branch_idx - 1, p1 - 1) += 1.0;
    if (p2 > 0) MAT(solver, branch_idx - 1, p2 - 1) += -1.0;
    if (s1 > 0) MAT(solver, branch_idx - 1, s1 - 1) += -n;
    if (s2 > 0) MAT(solver, branch_idx - 1, s2 - 1) += n;
    
    /* RHS for constraint equation is 0 (no independent source) */
    solver->b[branch_idx - 1] = 0.0;
}

MNAStatus mna_add_ideal_transformer(MNASolver* solver,
                                    int primary_p, int primary_n,
                                    int secondary_p, int secondary_n,
                                    double turns_ratio,
                                    ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;
    if (turns_ratio <= 0.0) return MNA_INVALID_PARAMETER;

    MNAStatus status = mna_validate_nodes(solver, primary_p, primary_n);
    if (status != MNA_SUCCESS) return status;
    status = mna_validate_nodes(solver, secondary_p, secondary_n);
    if (status != MNA_SUCCESS) return status;

    int new_count = solver->num_components + 1;
    if (!mna_resize_components(solver, new_count)) {
        return MNA_INSUFFICIENT_MEMORY;
    }

    /* Allocate transformer-specific data */
    TransformerData* tdata = (TransformerData*)malloc(sizeof(TransformerData));
    if (!tdata) return MNA_INSUFFICIENT_MEMORY;
    tdata->turns_ratio = turns_ratio;

    /* Allocate NPoleData for the 4-terminal element */
    NPoleData* npole = (NPoleData*)malloc(sizeof(NPoleData));
    if (!npole) {
        free(tdata);
        return MNA_INSUFFICIENT_MEMORY;
    }

    npole->nodes = (int*)malloc(4 * sizeof(int));
    npole->last_values = (double*)calloc(4, sizeof(double));
    npole->branch_current_indices = NULL;  /* We manage branch_var_index ourselves */

    if (!npole->nodes || !npole->last_values) {
        free(npole->nodes);
        free(npole->last_values);
        free(npole);
        free(tdata);
        return MNA_INSUFFICIENT_MEMORY;
    }

    npole->nodes[0] = primary_p;
    npole->nodes[1] = primary_n;
    npole->nodes[2] = secondary_p;
    npole->nodes[3] = secondary_n;
    npole->num_nodes = 4;
    npole->stamp_func = mna_stamp_ideal_transformer;
    npole->user_data = tdata;
    npole->num_branch_currents = 0;  /* We handle the branch variable directly */

    /* Need to expand matrix for the branch current variable (like a voltage source) */
    int branch_var_index = solver->max_node_index + solver->num_sources + 1;
    if (!mna_resize_matrix(solver, branch_var_index + 1)) {
        free(npole->nodes);
        free(npole->last_values);
        free(npole);
        free(tdata);
        return MNA_INSUFFICIENT_MEMORY;
    }
    tdata->branch_var_index = branch_var_index;
    solver->num_sources++;  /* Increment source count for the branch current variable */

    /* Create the component */
    Component comp;
    memset(&comp, 0, sizeof(comp));
    comp.type = MNA_CUSTOM_NPOLE;
    comp.node1 = primary_p;
    comp.node2 = primary_n;
    comp.value = turns_ratio;
    comp.state = 1;
    comp.source_type = SOURCE_CURRENT;
    comp.smoothed_alpha = 1.0;
    comp.data.npole.npole_data = npole;

    int index = solver->num_components++;
    solver->components[index] = comp;

    if (handle) *handle = index;
    return MNA_SUCCESS;
}

MNAStatus mna_set_transformer_turns_ratio(MNASolver* solver,
                                          ComponentHandle handle,
                                          double turns_ratio) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return MNA_INVALID_HANDLE;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) {
        return MNA_INVALID_PARAMETER;
    }

    if (turns_ratio <= 0.0) return MNA_INVALID_PARAMETER;

    TransformerData* tdata = (TransformerData*)comp->data.npole.npole_data->user_data;
    if (tdata) {
        tdata->turns_ratio = turns_ratio;
        comp->value = turns_ratio;
        return MNA_SUCCESS;
    }

    return MNA_INVALID_PARAMETER;
}

MNAStatus mna_get_transformer_turns_ratio(MNASolver* solver,
                                          ComponentHandle handle,
                                          double* turns_ratio) {
    if (!solver || handle < 0 || handle >= solver->num_components || !turns_ratio) {
        return MNA_INVALID_HANDLE;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) {
        return MNA_INVALID_PARAMETER;
    }

    TransformerData* tdata = (TransformerData*)comp->data.npole.npole_data->user_data;
    if (tdata) {
        *turns_ratio = tdata->turns_ratio;
        return MNA_SUCCESS;
    }

    return MNA_INVALID_PARAMETER;
}
