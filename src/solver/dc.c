#include "dc.h"
#include "../include/matrix.h"
#include "../elements/passive.h"
#include "../elements/sources.h"
#include "../elements/nonlinear/nonlinear.h"
#include "../elements/transformer.h"
#include <stdlib.h>

static MNAStatus mna_solve_dc_linear(MNASolver* solver) {
    mna_reset_system(solver);

    int source_count = 0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        const int n1 = comp->node1;
        const int n2 = comp->node2;

        switch (comp->type) {
            case MNA_RESISTOR:
                mna_stamp_conductance(solver, n1, n2, 1.0 / comp->value);
                break;

            case MNA_CAPACITOR:
                mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                break;

            case MNA_INDUCTOR:
                mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
                break;

            case MNA_SOURCE:
                if (comp->source_type == SOURCE_VOLTAGE) {
                    mna_stamp_voltage_source(solver, i, source_count++);
                } else {
                    mna_stamp_current_source(solver, n1, n2, comp->value);
                }
                break;

            case MNA_CUSTOM_NONLINEAR:
                mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                break;

            case MNA_CUSTOM_NPOLE: {
                NPoleData* npole = comp->data.npole.npole_data;
                if (npole) {
                    npole->stamp_func(solver, npole->nodes, npole->num_nodes,
                                      npole->user_data, 0.0, 0.0, 0);
                }
                break;
            }

            case MNA_SWITCH: {
                const double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }

            default:
                break;
        }
    }

    return mna_solve_linear_system(solver, mna_active_size(solver));
}

static MNAStatus mna_solve_dc_nonlinear(MNASolver* solver) {
    int matrix_size = mna_active_size(solver);

    double* orig_source_values = (double*)malloc((size_t)solver->num_components * sizeof(double));
    if (!orig_source_values) return MNA_INSUFFICIENT_MEMORY;

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];

        if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
            memset(comp->data.npole.npole_data->last_values, 0,
                   (size_t)comp->data.npole.npole_data->num_nodes * sizeof(double));
        }

        if (comp->type == MNA_SOURCE) {
            orig_source_values[i] = comp->value;
            comp->value = 0.0;
        }

        if (comp->type == MNA_CUSTOM_NONLINEAR) {
            comp->last_voltage = 0.0;
            comp->last_conductance = MNA_MIN_CONDUCTANCE;
        }
    }

    double current_factor = 0.0;
    const double min_step = 0.01;
    const double max_step = 0.2;
    double step_size = max_step;
    MNAStatus status = MNA_SUCCESS;
    int total_steps = 0;
    const int max_total_steps = 200;

    while (current_factor < 1.0 && total_steps < max_total_steps) {
        double next_factor = current_factor + step_size;
        if (next_factor > 1.0) next_factor = 1.0;

        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            if (comp->type == MNA_SOURCE) {
                comp->value = orig_source_values[i] * next_factor;
            }
        }

        int converged = 0;
        int iteration = 0;

        while (!converged && iteration < MNA_MAX_ITER) {
            mna_reset_system(solver);
            mna_ensure_ground_paths(solver);

            int source_count = 0;
            for (int i = 0; i < solver->num_components; i++) {
                Component* comp = &solver->components[i];
                const int n1 = comp->node1;
                const int n2 = comp->node2;

                switch (comp->type) {
                    case MNA_RESISTOR:
                        mna_stamp_conductance(solver, n1, n2, 1.0 / comp->value);
                        break;

                    case MNA_CAPACITOR:
                        mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                        break;

                    case MNA_INDUCTOR:
                        mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
                        break;

                    case MNA_SOURCE:
                        if (comp->source_type == SOURCE_VOLTAGE) {
                            mna_stamp_voltage_source(solver, i, source_count++);
                        } else {
                            mna_stamp_current_source(solver, n1, n2, comp->value);
                        }
                        break;

                    case MNA_CUSTOM_NONLINEAR:
                        mna_stamp_custom_nonlinear(solver, i, 1, 0);
                        break;

                    case MNA_CUSTOM_NPOLE: {
                        NPoleData* npole = comp->data.npole.npole_data;
                        if (npole) {
                            npole->stamp_func(solver, npole->nodes, npole->num_nodes,
                                              npole->user_data, 0.0, 0.0, 0);
                        }
                        break;
                    }

                    case MNA_SWITCH: {
                        const double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                        mna_stamp_conductance(solver, n1, n2, g);
                        break;
                    }

                    default:
                        break;
                }
            }

            status = mna_solve_linear_system(solver, matrix_size);
            if (status != MNA_SUCCESS) break;

            converged = 1;
            for (int i = 0; i < solver->num_components; i++) {
                Component* comp = &solver->components[i];

                if (comp->type == MNA_CUSTOM_NONLINEAR) {
                    const int n1 = comp->node1;
                    const int n2 = comp->node2;
                    const double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
                    const double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
                    const double new_voltage = v1 - v2;
                    const double voltage_diff = fabs(new_voltage - comp->last_voltage);
                    const double abs_tol = MNA_ABSTOL + MNA_RELTOL * (1.0 + fabs(new_voltage));

                    if (comp->nonlinear_type == NONLINEAR_RESISTOR) {
                        ComponentState state = {
                            .voltage = new_voltage,
                            .current = comp->last_current,
                            .charge = comp->last_charge,
                            .flux = comp->last_flux,
                            .dt = solver->dt
                        };
                        double val1, val2;
                        comp->nonlinear_func(&state, comp->user_data, &val1, &val2);
                        comp->last_conductance = (val2 < MNA_MIN_CONDUCTANCE) ? MNA_MIN_CONDUCTANCE :
                                                 (val2 > MNA_MAX_CONDUCTANCE) ? MNA_MAX_CONDUCTANCE : val2;
                    }

                    if (voltage_diff > 0.5 && iteration < 10) {
                        const double damping = 0.7;
                        comp->last_voltage = damping * comp->last_voltage + (1 - damping) * new_voltage;
                    } else {
                        comp->last_voltage = new_voltage;
                    }

                    if (voltage_diff > abs_tol) converged = 0;

                } else if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
                    NPoleData* npole = comp->data.npole.npole_data;
                    double max_diff = 0.0;

                    for (int j = 0; j < npole->num_nodes; j++) {
                        double current_val = (npole->nodes[j] > 0) ?
                                             solver->x[npole->nodes[j]-1] : 0.0;
                        double diff = fabs(current_val - npole->last_values[j]);
                        if (diff > max_diff) max_diff = diff;
                        npole->last_values[j] = current_val;
                    }

                    TransformerData* xf = (TransformerData*)npole->user_data;
                    if (xf && xf->Lm > 0) {
                        double flux_diff = fabs(xf->phi - xf->prev_phi);
                        if (flux_diff > max_diff) max_diff = flux_diff;
                    }

                    double abs_tol = MNA_ABSTOL + MNA_RELTOL * (1.0 + max_diff);
                    if (max_diff > abs_tol) converged = 0;
                }
            }

            if (status != MNA_SUCCESS || converged) break;
            iteration++;
        }

        if (status != MNA_SUCCESS || !converged) {
            step_size /= 2.0;
            if (step_size < min_step) {
                status = MNA_CONVERGENCE_FAILURE;
                break;
            }
            continue;
        }

        current_factor = next_factor;
        total_steps++;
        step_size *= 1.5;
        if (step_size > max_step) step_size = max_step;
    }

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->type == MNA_SOURCE) {
            comp->value = orig_source_values[i];
        }
    }

    free(orig_source_values);
    return status;
}

MNAStatus mna_solve_dc(MNASolver* solver) {
    if (!solver) return MNA_INVALID_HANDLE;

    const bool has_nonlinear = (solver->num_nonlinear > 0);

    if (!has_nonlinear) {
        return mna_solve_dc_linear(solver);
    } else {
        return mna_solve_dc_nonlinear(solver);
    }
}

void mna_solve_initial_conditions(MNASolver* solver) {
    if (!solver) return;

    int matrix_size = mna_active_size(solver);

    for (int i = 0; i < matrix_size; ++i) {
        memset(&MAT(solver, i, 0), 0, (size_t)matrix_size * sizeof(double));
    }
    memset(solver->b, 0, (size_t)matrix_size * sizeof(double));
    memset(solver->x, 0, (size_t)matrix_size * sizeof(double));

    int source_count = 0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        switch (comp->type) {
            case MNA_RESISTOR:
            case MNA_SWITCH: {
                double g = (comp->type == MNA_SWITCH && !comp->state) ?
                           MNA_MIN_CONDUCTANCE : 1.0 / comp->value;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }

            case MNA_CAPACITOR:
                mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
                break;

            case MNA_INDUCTOR:
                mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                break;

            case MNA_SOURCE:
                if (comp->source_type == SOURCE_VOLTAGE) {
                    mna_stamp_voltage_source(solver, i, source_count++);
                } else {
                    mna_stamp_current_source(solver, n1, n2, comp->value);
                }
                break;

            case MNA_CUSTOM_NONLINEAR:
                mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                break;

            case MNA_CUSTOM_NPOLE:
                break;

            default:
                break;
        }
    }

    mna_ensure_ground_paths(solver);
    MNAStatus status = mna_solve_linear_system(solver, matrix_size);
    if (status != MNA_SUCCESS) return;

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;
        double v1 = (n1 > 0) ? solver->x[n1 - 1] : 0.0;
        double v2 = (n2 > 0) ? solver->x[n2 - 1] : 0.0;
        double v_diff = v1 - v2;

        comp->prev_voltage = 0.0;
        comp->prev_current = 0.0;
        comp->stage1_voltage = 0.0;
        comp->stage1_current = 0.0;

        switch (comp->type) {
            case MNA_CAPACITOR:
                comp->last_voltage = 0.0;
                comp->last_current = MNA_MAX_CONDUCTANCE * v_diff;
                break;

            case MNA_INDUCTOR:
                comp->last_voltage = v_diff;
                comp->last_current = 0.0;
                break;

            case MNA_RESISTOR:
            case MNA_SOURCE:
            case MNA_SWITCH:
                comp->last_voltage = v_diff;
                comp->last_current = (comp->type == MNA_RESISTOR) ?
                                     v_diff / comp->value : 0.0;
                break;

            default:
                comp->last_voltage = 0.0;
                comp->last_current = 0.0;
                break;
        }

        comp->prev_voltage = comp->last_voltage;
        comp->prev_current = comp->last_current;
    }
}
