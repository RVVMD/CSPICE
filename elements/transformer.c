#include "transformer.h"
#include "../include/matrix.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

static void transformer_stamp_func(MNASolver* solver,
                                    const int* nodes,
                                    int num_nodes,
                                    void* user_data,
                                    double time,
                                    double dt);

static double compute_saturation_current(TransformerData* xf) {
    if (!xf->use_saturation || xf->Lm <= 0) {
        return xf->phi / xf->Lm;
    }

    double arg = xf->phi / (xf->Lm * xf->saturation_current);
    return xf->saturation_current * tanh(arg);
}

static double compute_incremental_inductance(TransformerData* xf) {
    if (!xf->use_saturation || xf->Lm <= 0) {
        return xf->Lm;
    }

    double arg = xf->phi / (xf->Lm * xf->saturation_current);
    double sech_sq = 1.0 / (cosh(arg) * cosh(arg));
    double di_dphi = (1.0 / xf->Lm) * sech_sq;

    if (fabs(di_dphi) < 1e-12) return xf->Lm * 100;
    return 1.0 / di_dphi;
}

static void transformer_stamp_func(MNASolver* solver,
                                    const int* nodes,
                                    int num_nodes,
                                    void* user_data,
                                    double time,
                                    double dt) {
    (void)time;
    (void)num_nodes;

    TransformerData* xf = (TransformerData*)user_data;
    double n = xf->turns_ratio;
    int branch_idx = xf->branch_current_index;

    int p1 = nodes[0];
    int p2 = nodes[1];
    int s1 = nodes[2];
    int s2 = nodes[3];

    double v_p1 = (p1 > 0) ? solver->x[p1 - 1] : 0.0;
    double v_p2 = (p2 > 0) ? solver->x[p2 - 1] : 0.0;
    double v_primary = v_p1 - v_p2;

    double i_mag_prev = xf->i_mag;
    if (dt > 0 && xf->Lm > 0) {
        xf->prev_phi = xf->phi;
        xf->phi = xf->phi + v_primary * dt;
        xf->i_mag = compute_saturation_current(xf);
    } else if (xf->Lm > 0) {
        xf->i_mag = xf->phi / xf->Lm;
    } else {
        xf->i_mag = 0.0;
    }

    bool is_transformer_mode = (dt > 0) || xf->is_ac_analysis || (xf->Lm > 0);

    if (!is_transformer_mode) {
        if (p1 > 0) MAT(solver, branch_idx, p1 - 1) = 1.0;
        if (p2 > 0) MAT(solver, branch_idx, p2 - 1) = -1.0;

        if (p1 > 0) MAT(solver, p1 - 1, branch_idx) = 1.0;
        if (p2 > 0) MAT(solver, p2 - 1, branch_idx) = -1.0;

    } else {
        if (p1 > 0) MAT(solver, branch_idx, p1 - 1) = -n;
        if (p2 > 0) MAT(solver, branch_idx, p2 - 1) = n;
        if (s1 > 0) MAT(solver, branch_idx, s1 - 1) = 1.0;
        if (s2 > 0) MAT(solver, branch_idx, s2 - 1) = -1.0;

        if (p1 > 0) MAT(solver, p1 - 1, branch_idx) = 1.0;
        if (p2 > 0) MAT(solver, p2 - 1, branch_idx) = -1.0;
        if (s1 > 0) MAT(solver, s1 - 1, branch_idx) = -1.0 / n;
        if (s2 > 0) MAT(solver, s2 - 1, branch_idx) = 1.0 / n;
    }

    if (xf->Lm > 0 && dt > 0) {
        double Gm = dt / xf->Lm;

        if (Gm > MNA_MAX_CONDUCTANCE) Gm = MNA_MAX_CONDUCTANCE;
        if (Gm < MNA_MIN_CONDUCTANCE) Gm = MNA_MIN_CONDUCTANCE;

        if (p1 > 0) {
            MAT(solver, p1 - 1, p1 - 1) += Gm;
            if (p2 > 0) MAT(solver, p1 - 1, p2 - 1) -= Gm;
        }
        if (p2 > 0) {
            if (p1 > 0) MAT(solver, p2 - 1, p1 - 1) -= Gm;
            MAT(solver, p2 - 1, p2 - 1) += Gm;
        }

        if (p1 > 0) solver->b[p1 - 1] -= i_mag_prev;
        if (p2 > 0) solver->b[p2 - 1] += i_mag_prev;
    }
}

MNAStatus mna_add_ideal_transformer(MNASolver* solver,
                                     int node_p1, int node_p2,
                                     int node_s1, int node_s2,
                                     double turns_ratio,
                                     ComponentHandle* handle) {
    return mna_add_voltage_transformer(solver, node_p1, node_p2,
                                        node_s1, node_s2,
                                        turns_ratio, 0.0, handle);
}

MNAStatus mna_add_voltage_transformer(MNASolver* solver,
                                       int node_p1, int node_p2,
                                       int node_s1, int node_s2,
                                       double turns_ratio,
                                       double Lm,
                                       ComponentHandle* handle) {
    return mna_add_transformer_with_saturation(solver, node_p1, node_p2,
                                                node_s1, node_s2,
                                                turns_ratio, Lm,
                                                0.0, 1.0, handle);
}

MNAStatus mna_add_transformer_with_saturation(MNASolver* solver,
                                               int node_p1, int node_p2,
                                               int node_s1, int node_s2,
                                               double turns_ratio,
                                               double Lm,
                                               double I_sat,
                                               double sat_factor,
                                               ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;

    int nodes[4] = {node_p1, node_p2, node_s1, node_s2};
    for (int i = 0; i < 4; i++) {
        if (nodes[i] < 0 || nodes[i] > solver->max_node_index) {
            return MNA_INVALID_NODE;
        }
    }

    if (turns_ratio <= 0.0) {
        return MNA_INVALID_PARAMETER;
    }

    TransformerData* xf = (TransformerData*)malloc(sizeof(TransformerData));
    if (!xf) return MNA_INSUFFICIENT_MEMORY;

    memset(xf, 0, sizeof(TransformerData));
    xf->turns_ratio = turns_ratio;
    xf->Lm = Lm;
    xf->i_mag = 0.0;
    xf->phi = 0.0;
    xf->prev_phi = 0.0;
    xf->Gm_eq = 0.0;
    xf->I_eq = 0.0;
    xf->branch_current_index = -1;
    xf->saturation_current = I_sat;
    xf->saturation_factor = (sat_factor > 0) ? sat_factor : 1.0;
    xf->use_saturation = (I_sat > 0 && Lm > 0);
    xf->is_ac_analysis = false;

    MNAStatus status = mna_add_custom_n_pole(solver, nodes, 4,
                                              transformer_stamp_func,
                                              xf, 1, handle);

    if (status != MNA_SUCCESS) {
        free(xf);
        return status;
    }

    Component* comp = &solver->components[*handle];
    if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
        xf->branch_current_index = comp->data.npole.npole_data->branch_current_indices[0];
    }

    return MNA_SUCCESS;
}

MNAStatus mna_add_current_transformer(MNASolver* solver,
                                       int node_p1, int node_p2,
                                       int node_s1, int node_s2,
                                       double turns_ratio,
                                       double burden,
                                       ComponentHandle* handle) {
    (void)burden;
    return mna_add_ideal_transformer(solver, node_p1, node_p2,
                                      node_s1, node_s2,
                                      turns_ratio, handle);
}

double mna_get_transformer_magnetizing_current(MNASolver* solver,
                                                ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return 0.0;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) {
        return 0.0;
    }

    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;

    if (xf->Lm > 0) {
        return xf->i_mag;
    } else {
        int branch_idx = xf->branch_current_index;
        if (branch_idx > 0 && branch_idx <= solver->max_node_index + solver->num_sources) {
            return solver->x[branch_idx - 1];
        }
        return 0.0;
    }
}

double mna_get_transformer_flux(MNASolver* solver,
                                 ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return 0.0;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) {
        return 0.0;
    }

    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    return xf->phi;
}

double mna_get_transformer_primary_voltage(MNASolver* solver,
                                            ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return 0.0;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) {
        return 0.0;
    }

    int* nodes = comp->data.npole.npole_data->nodes;
    int p1 = nodes[0];
    int p2 = nodes[1];

    double v_p1 = (p1 > 0) ? solver->x[p1 - 1] : 0.0;
    double v_p2 = (p2 > 0) ? solver->x[p2 - 1] : 0.0;

    return v_p1 - v_p2;
}

double mna_get_transformer_secondary_voltage(MNASolver* solver,
                                              ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return 0.0;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) {
        return 0.0;
    }

    int* nodes = comp->data.npole.npole_data->nodes;
    int s1 = nodes[2];
    int s2 = nodes[3];

    double v_s1 = (s1 > 0) ? solver->x[s1 - 1] : 0.0;
    double v_s2 = (s2 > 0) ? solver->x[s2 - 1] : 0.0;

    return v_s1 - v_s2;
}

int mna_get_transformer_branch_index(MNASolver* solver,
                                      ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) {
        return -1;
    }

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) {
        return -1;
    }

    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    return xf->branch_current_index;
}
