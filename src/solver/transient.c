#include "transient.h"
#include "../include/matrix.h"
#include "../elements/passive.h"
#include "../elements/sources.h"
#include "../elements/nonlinear/nonlinear.h"
#include <stdlib.h>

static void mna_solve_initial_conditions(MNASolver* solver) {
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
                {
                    NPoleData* npole = comp->data.npole.npole_data;
                    if (npole) {
                        npole->stamp_func(solver, npole->nodes, npole->num_nodes,
                                          npole->user_data, 0.0, 0.0, 0);
                    }
                }
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

void mna_init_transient(MNASolver* solver) {
    if (!solver) return;

    solver->time = 0.0;
    solver->transient_initialized = 1;

    if (solver->preserve_dc_state) {
        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            int n1 = comp->node1;
            int n2 = comp->node2;

            double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
            double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
            double v_diff = v1 - v2;

            comp->prev_voltage = v_diff;
            comp->stage1_voltage = v_diff;
            comp->last_voltage = v_diff;

            if (comp->type == MNA_CAPACITOR) {
                comp->last_current = 0.0;
                comp->prev_current = 0.0;
                comp->stage1_current = 0.0;
            }
            else if (comp->type == MNA_INDUCTOR) {
                comp->last_current = 0.0;
                comp->prev_current = 0.0;
                comp->stage1_current = 0.0;
            }
        }
    } else {
        mna_solve_initial_conditions(solver);
    }
}

MNAStatus mna_solve_transient_step(MNASolver* solver, double dt) {
    if (!solver) return MNA_INVALID_HANDLE;

    solver->dt = dt;
    int matrix_size = mna_active_size(solver);
    const int stages = 2;
    const double gamma = MNA_TRBDF2_GAMMA;

    double alpha_min = 1.0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->smoothed_alpha < alpha_min) {
            alpha_min = comp->smoothed_alpha;
        }
    }

    for (int stage = 1; stage <= stages; stage++) {
        mna_reset_system(solver);

        int source_count = 0;
        double h_eff = (stage == 1) ? (gamma * dt) : ((1.0 - gamma) * dt);

        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            int n1 = comp->node1;
            int n2 = comp->node2;

            switch (comp->type) {
                case MNA_CAPACITOR: {
                    double G_eq, I_eq;
                    double v_n = comp->last_voltage;

                    if (stage == 1) {
                        G_eq = (2.0 * comp->value) / h_eff;
                        I_eq = (G_eq * v_n) + (comp->last_current * alpha_min);
                    } else {
                        double v_ng = comp->stage1_voltage;
                        G_eq = (comp->value * (2.0 - gamma)) / (dt * (1.0 - gamma));
                        I_eq = (comp->value / (dt * gamma * (1.0 - gamma))) * v_ng -
                               (comp->value * (1.0 - gamma) / (dt * gamma)) * v_n;
                    }

                    comp->trans_G_eq = G_eq;
                    comp->trans_I_eq = I_eq;

                    mna_stamp_conductance(solver, n1, n2, G_eq);
                    mna_stamp_current_source(solver, n1, n2, -I_eq);
                    break;
                }

                case MNA_INDUCTOR: {
                    double G_eq, I_eq;
                    double i_n = comp->last_current;

                    if (stage == 1) {
                        G_eq = h_eff / (2.0 * comp->value);
                        I_eq = i_n + (G_eq * comp->last_voltage * alpha_min);
                    } else {
                        double i_ng = comp->stage1_current;
                        G_eq = (dt * (1.0 - gamma)) / (comp->value * (2.0 - gamma));
                        I_eq = (1.0 / (gamma * (2.0 - gamma))) * i_ng -
                               (((1.0 - gamma) * (1.0 - gamma)) / (gamma * (2.0 - gamma))) * i_n;
                    }

                    comp->trans_G_eq = G_eq;
                    comp->trans_I_eq = I_eq;

                    mna_stamp_conductance(solver, n1, n2, G_eq);
                    mna_stamp_current_source(solver, n1, n2, I_eq);
                    break;
                }

                case MNA_RESISTOR: {
                    double g = 1.0 / comp->value;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }

                case MNA_SOURCE: {
                    if (comp->source_type == SOURCE_VOLTAGE) {
                        mna_stamp_voltage_source(solver, i, source_count++);
                    } else {
                        mna_stamp_current_source(solver, n1, n2, comp->value);
                    }
                    break;
                }

                case MNA_SWITCH: {
                    double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }

                case MNA_CUSTOM_NONLINEAR: {
                    mna_stamp_custom_nonlinear(solver, i, 0, stage);
                    break;
                }

                case MNA_CUSTOM_NPOLE: {
                    NPoleData* npole = comp->data.npole.npole_data;
                    if (npole) {
                        npole->stamp_func(solver, npole->nodes, npole->num_nodes,
                                          npole->user_data, solver->time, dt, stage);
                    }
                    break;
                }

                default:
                    break;
            }
        }

        MNAStatus status = mna_solve_linear_system(solver, matrix_size);
        if (status != MNA_SUCCESS) return status;

        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            int n1 = comp->node1;
            int n2 = comp->node2;
            double v1 = (n1 > 0) ? solver->x[n1 - 1] : 0.0;
            double v2 = (n2 > 0) ? solver->x[n2 - 1] : 0.0;
            double v = v1 - v2;

            if (comp->type == MNA_CAPACITOR) {
                double i = comp->trans_G_eq * v - comp->trans_I_eq;
                if (stage == 1) {
                    comp->stage1_voltage = v;
                    comp->stage1_current = i;
                } else {
                    comp->last_voltage = v;
                    comp->last_current = i;
                }
            } else if (comp->type == MNA_INDUCTOR) {
                double i = comp->trans_G_eq * v + comp->trans_I_eq;
                if (stage == 1) {
                    comp->stage1_voltage = v;
                    comp->stage1_current = i;
                } else {
                    comp->last_voltage = v;
                    comp->last_current = i;
                }
            } else if (comp->type == MNA_CUSTOM_NPOLE) {
                /* Post-solve update */
            } else {
                if (stage == 1) {
                    comp->stage1_voltage = v;
                } else {
                    comp->last_voltage = v;
                }
            }
        }
    }

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        double rate = 0.0;

        if (comp->type == MNA_CAPACITOR) {
            double dv = comp->last_voltage - comp->prev_voltage;
            double v_ref = fmax(fabs(comp->last_voltage), fabs(comp->prev_voltage));
            rate = fabs(dv) / (v_ref + MNA_V_FLOOR);
        } else if (comp->type == MNA_INDUCTOR) {
            double di = comp->last_current - comp->prev_current;
            double i_ref = fmax(fabs(comp->last_current), fabs(comp->prev_current));
            rate = fabs(di) / (i_ref + MNA_I_FLOOR);
        } else if (comp->type == MNA_CUSTOM_NONLINEAR) {
            if (comp->nonlinear_type == NONLINEAR_INDUCTOR) {
                double di = comp->last_current - comp->prev_current;
                double i_ref = fmax(fabs(comp->last_current), fabs(comp->prev_current));
                rate = fabs(di) / (i_ref + MNA_I_FLOOR);
            } else {
                double dv = comp->last_voltage - comp->prev_voltage;
                double v_ref = fmax(fabs(comp->last_voltage), fabs(comp->prev_voltage));
                rate = fabs(dv) / (v_ref + MNA_V_FLOOR);
            }
        }

        double alpha = mna_compute_alpha(rate, comp->smoothed_alpha);
        comp->smoothed_alpha = alpha;
    }

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        comp->prev_voltage = comp->last_voltage;
        comp->prev_current = comp->last_current;
    }

    solver->time += dt;
    return MNA_SUCCESS;
}

double mna_get_time(MNASolver* solver) {
    if (!solver) return 0.0;
    return solver->time;
}
