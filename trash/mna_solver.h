#ifndef MNA_SOLVER_H
#define MNA_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdalign.h>
#include <string.h>
#include <stdbool.h>

// Optimized constants with better numerical properties
#define MNA_MAX_NODES 50
#define MNA_MAX_SOURCES 20
#define MNA_MAX_COMPONENTS 100
#define MNA_MAX_NONLINEAR 20
#define MNA_MAX_ITER 50
#define MNA_RELTOL 1e-6
#define MNA_ABSTOL 1e-9
#define MNA_VT 0.02585
#define MNA_MIN_CONDUCTANCE 1e-12
#define MNA_MAX_CONDUCTANCE 1e12
#define MNA_MAX_MATRIX_SIZE (MNA_MAX_NODES + MNA_MAX_SOURCES)
#define MNA_MATRIX_ELEMS (MNA_MAX_MATRIX_SIZE * MNA_MAX_MATRIX_SIZE)

// Precomputed constants
#define TWO_PI 6.28318530717958647692

typedef struct MNASolver MNASolver;

// Unified nonlinear element type
typedef enum {
    NONLINEAR_RESISTOR,
    NONLINEAR_CAPACITOR,
    NONLINEAR_INDUCTOR
} NonlinearType;

// Generalized nonlinear function interface
typedef void (*CustomNonlinearFunc)(MNASolver* solver, int comp_index,
                                  double voltage, double current,
                                  double* value1, double* value2,
                                  NonlinearType nl_type);

// Unified source type
typedef enum {
    SOURCE_VOLTAGE,
    SOURCE_CURRENT
} SourceType;

typedef enum {
    MNA_RESISTOR,
    MNA_CAPACITOR,
    MNA_INDUCTOR,
    MNA_SOURCE,        // Unified source component
    MNA_SWITCH,
    MNA_CUSTOM_NONLINEAR
} ComponentType;

typedef struct {
    ComponentType type;
    int node1;
    int node2;
    double value;        // DC value for sources/resistors/etc
    double ac_magnitude;
    double ac_phase;
    int state;
    int is_nonlinear;

    // Source-specific fields
    SourceType source_type;  // Only used when type == MNA_SOURCE

    // Nonlinear-specific fields
    NonlinearType nonlinear_type;
    double last_voltage;
    double last_current;
    double last_conductance;
    double last_charge;      // For capacitors
    double last_flux;        // For inductors
    CustomNonlinearFunc nonlinear_func;
    void* user_data;

    // Transient companion model
    double trans_G_eq;
    double trans_I_eq;
} Component;

struct MNASolver {
    int num_nodes;
    int num_components;
    int num_sources;      // Only counts voltage sources (for MNA variables)
    int num_nonlinear;
    int max_node_index;
    double time;
    double dt;            // Time step for transient analysis

    Component components[MNA_MAX_COMPONENTS];
    alignas(64) double A[MNA_MATRIX_ELEMS];  // 1D array for better cache locality
    double b[MNA_MAX_MATRIX_SIZE];
    double x[MNA_MAX_MATRIX_SIZE];
    double complex ac_solution[MNA_MAX_MATRIX_SIZE];
};

// Matrix access macros
#define MAT(solver, i, j) ((solver)->A[(i) * (solver->max_node_index + solver->num_sources) + (j)])

void mna_init(MNASolver* solver) {
    memset(solver, 0, sizeof(MNASolver));
}

void mna_add_resistor(MNASolver* solver, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_RESISTOR,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = 0,
        .source_type = SOURCE_CURRENT,
        .nonlinear_type = NONLINEAR_RESISTOR
    };

    solver->components[solver->num_components++] = comp;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_add_capacitor(MNASolver* solver, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_CAPACITOR,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = 0,
        .source_type = SOURCE_CURRENT,
        .nonlinear_type = NONLINEAR_CAPACITOR
    };

    solver->components[solver->num_components++] = comp;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_add_inductor(MNASolver* solver, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_INDUCTOR,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = 0,
        .source_type = SOURCE_CURRENT,
        .nonlinear_type = NONLINEAR_INDUCTOR
    };

    solver->components[solver->num_components++] = comp;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_add_voltage_source(MNASolver* solver, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_SOURCE,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = 0,
        .source_type = SOURCE_VOLTAGE,
        .nonlinear_type = NONLINEAR_RESISTOR
    };

    solver->components[solver->num_components++] = comp;
    solver->num_sources++;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_add_current_source(MNASolver* solver, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_SOURCE,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = 0,
        .source_type = SOURCE_CURRENT,
        .nonlinear_type = NONLINEAR_RESISTOR
    };

    solver->components[solver->num_components++] = comp;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_add_switch(MNASolver* solver, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_SWITCH,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = 0,
        .source_type = SOURCE_CURRENT,
        .nonlinear_type = NONLINEAR_RESISTOR
    };

    solver->components[solver->num_components++] = comp;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                            NonlinearType nl_type,
                            CustomNonlinearFunc func, void* user_data,
                            double initial_value1, double initial_value2) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_CUSTOM_NONLINEAR,
        .node1 = node1,
        .node2 = node2,
        .state = 1,
        .is_nonlinear = 1,
        .source_type = SOURCE_CURRENT,
        .nonlinear_type = nl_type,
        .nonlinear_func = func,
        .user_data = user_data
    };

    switch (nl_type) {
        case NONLINEAR_RESISTOR:
            comp.last_voltage = initial_value1;
            break;
        case NONLINEAR_CAPACITOR:
            comp.last_voltage = initial_value1;
            comp.last_charge = initial_value2;
            break;
        case NONLINEAR_INDUCTOR:
            comp.last_current = initial_value1;
            comp.last_flux = initial_value2;
            break;
    }

    solver->components[solver->num_components] = comp;
    solver->num_components++;
    solver->num_nonlinear++;
    solver->max_node_index = (node1 > solver->max_node_index) ? node1 : solver->max_node_index;
    solver->max_node_index = (node2 > solver->max_node_index) ? node2 : solver->max_node_index;
}

void mna_set_ac_source(MNASolver* solver, int comp_index, double magnitude, double phase) {
    if (comp_index < solver->num_components) {
        solver->components[comp_index].ac_magnitude = magnitude;
        solver->components[comp_index].ac_phase = phase;
    }
}

void mna_set_switch_state(MNASolver* solver, int comp_index, int state) {
    if (comp_index < solver->num_components &&
        solver->components[comp_index].type == MNA_SWITCH) {
        solver->components[comp_index].state = state;
    }
}

void mna_reset_system(MNASolver* solver) {
    int matrix_size = solver->max_node_index + solver->num_sources;
    int total_elems = matrix_size * matrix_size;

    for (int i = 0; i < total_elems; i++) {
        solver->A[i] = 0.0;
    }

    for (int i = 0; i < matrix_size; i++) {
        solver->b[i] = 0.0;
        solver->x[i] = 0.0;
    }
}

void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g) {
    if (node1 > 0) MAT(solver, node1-1, node1-1) += g;
    if (node2 > 0) MAT(solver, node2-1, node2-1) += g;
    if (node1 > 0 && node2 > 0) {
        MAT(solver, node1-1, node2-1) -= g;
        MAT(solver, node2-1, node1-1) -= g;
    }
}

void mna_stamp_current_source(MNASolver* solver, int node1, int node2, double current_val) {
    if (node1 > 0) solver->b[node1-1] -= current_val;
    if (node2 > 0) solver->b[node2-1] += current_val;
}

void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx) {
    Component* vs = &solver->components[comp_index];
    int v_index = solver->max_node_index + source_idx;
    int n1 = vs->node1;
    int n2 = vs->node2;

    if (n1 > 0) {
        MAT(solver, n1-1, v_index) = 1.0;
        MAT(solver, v_index, n1-1) = 1.0;
    }
    if (n2 > 0) {
        MAT(solver, n2-1, v_index) = -1.0;
        MAT(solver, v_index, n2-1) = -1.0;
    }

    solver->b[v_index] = vs->value;
}

void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc) {
    Component* comp = &solver->components[comp_index];
    int n1 = comp->node1;
    int n2 = comp->node2;

    switch (comp->nonlinear_type) {
        case NONLINEAR_RESISTOR: {
            double voltage = comp->last_voltage;
            double current, conductance;
            comp->nonlinear_func(solver, comp_index, voltage, 0,
                                &current, &conductance, NONLINEAR_RESISTOR);

            if (conductance < MNA_MIN_CONDUCTANCE) conductance = MNA_MIN_CONDUCTANCE;
            if (conductance > MNA_MAX_CONDUCTANCE) conductance = MNA_MAX_CONDUCTANCE;

            mna_stamp_conductance(solver, n1, n2, conductance);
            mna_stamp_current_source(solver, n1, n2, (current - conductance * voltage));
            break;
        }

        case NONLINEAR_CAPACITOR: {
            if (is_dc) {
                // DC behavior: open circuit
                mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
            } else {
                // Transient behavior: companion model
                double v0 = comp->last_voltage;
                double q0, C0;
                comp->nonlinear_func(solver, comp_index, v0, 0,
                                    &q0, &C0, NONLINEAR_CAPACITOR);

                double G_eq = C0 / solver->dt;
                double I_eq = (q0 - comp->last_charge) / solver->dt - G_eq * v0;

                comp->trans_G_eq = G_eq;
                comp->trans_I_eq = I_eq;

                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
            }
            break;
        }

        case NONLINEAR_INDUCTOR: {
            if (is_dc) {
                // DC behavior: short circuit
                mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
            } else {
                // Transient behavior: companion model
                double i0 = comp->last_current;
                double phi0, L0;
                comp->nonlinear_func(solver, comp_index, 0, i0,
                                    &phi0, &L0, NONLINEAR_INDUCTOR);

                double G_eq = solver->dt / L0;
                double I_eq = i0 + (phi0 - comp->last_flux) / L0;

                comp->trans_G_eq = G_eq;
                comp->trans_I_eq = I_eq;

                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
            }
            break;
        }
    }
}

int mna_solve_linear_system(MNASolver* solver, int size) {
    const double pivot_threshold = 1e-12;
    double max_val;
    int max_row;

    for (int pivot = 0; pivot < size; pivot++) {
        // Find pivot
        max_row = pivot;
        max_val = fabs(MAT(solver, pivot, pivot));
        for (int i = pivot + 1; i < size; i++) {
            double abs_val = fabs(MAT(solver, i, pivot));
            if (abs_val > max_val) {
                max_val = abs_val;
                max_row = i;
            }
        }

        // Singularity check with relative threshold
        if (max_val < pivot_threshold) {
            return 0;
        }

        // Swap rows if needed
        if (max_row != pivot) {
            for (int j = pivot; j < size; j++) {
                double temp = MAT(solver, pivot, j);
                MAT(solver, pivot, j) = MAT(solver, max_row, j);
                MAT(solver, max_row, j) = temp;
            }
            double temp = solver->b[pivot];
            solver->b[pivot] = solver->b[max_row];
            solver->b[max_row] = temp;
        }

        // Elimination
        const double pivot_inv = 1.0 / MAT(solver, pivot, pivot);
        for (int i = pivot + 1; i < size; i++) {
            double factor = MAT(solver, i, pivot) * pivot_inv;
            for (int j = pivot + 1; j < size; j++) {
                MAT(solver, i, j) -= factor * MAT(solver, pivot, j);
            }
            solver->b[i] -= factor * solver->b[pivot];
            MAT(solver, i, pivot) = 0.0;
        }
    }

    // Back substitution
    for (int i = size - 1; i >= 0; i--) {
        double sum = solver->b[i];
        for (int j = i + 1; j < size; j++) {
            sum -= MAT(solver, i, j) * solver->x[j];
        }
        solver->x[i] = sum / MAT(solver, i, i);
    }

    return 1;
}

int mna_solve_dc(MNASolver* solver) {
    const int matrix_size = solver->max_node_index + solver->num_sources;
    int converged = 0;
    int iteration = 0;
    double prev_voltages[MNA_MAX_NONLINEAR] = {0};
    const bool has_nonlinear = (solver->num_nonlinear > 0);
    int source_count;
    int nonlinear_count;

    // Initial solution
    mna_reset_system(solver);
    source_count = 0;
    nonlinear_count = 0;

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
                mna_stamp_custom_nonlinear(solver, i, 1);
                nonlinear_count++;
                break;
            case MNA_SWITCH: {
                const double g = comp->state ?
                    1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }
        }
    }

    if (!mna_solve_linear_system(solver, matrix_size)) {
        return 0;
    }

    // Early exit if no nonlinear components
    if (!has_nonlinear) return 1;

    // Store initial voltages for nonlinear components
    int nl_idx = 0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->type == MNA_CUSTOM_NONLINEAR) {
            int n1 = comp->node1;
            int n2 = comp->node2;
            double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
            double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
            comp->last_voltage = v1 - v2;
            prev_voltages[nl_idx++] = comp->last_voltage;
        }
    }

    // Nonlinear iteration
    while (!converged && iteration < MNA_MAX_ITER) {
        converged = 1;
        nl_idx = 0;

        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            if (comp->type != MNA_CUSTOM_NONLINEAR) continue;

            const int n1 = comp->node1;
            const int n2 = comp->node2;
            const double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
            const double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
            const double new_voltage = v1 - v2;
            const double voltage_diff = fabs(new_voltage - comp->last_voltage);
            const double abs_tol = MNA_ABSTOL + MNA_RELTOL * fabs(new_voltage);

            // Adaptive damping for convergence
            const double damping = (iteration < 5) ? 0.5 : 0.9;
            comp->last_voltage = damping * new_voltage +
                                (1 - damping) * comp->last_voltage;

            if (voltage_diff > abs_tol) {
                converged = 0;
            }

            // Update nonlinear parameters
            if (comp->nonlinear_type == NONLINEAR_RESISTOR) {
                double val1, val2;
                comp->nonlinear_func(solver, i, comp->last_voltage, 0,
                                   &val1, &val2, comp->nonlinear_type);
                comp->last_conductance =
                    (val2 < MNA_MIN_CONDUCTANCE) ? MNA_MIN_CONDUCTANCE :
                    (val2 > MNA_MAX_CONDUCTANCE) ? MNA_MAX_CONDUCTANCE : val2;
            }
            prev_voltages[nl_idx++] = comp->last_voltage;
        }

        if (converged) break;

        // Rebuild matrix with updated nonlinear components
        mna_reset_system(solver);
        source_count = 0;

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
                    mna_stamp_custom_nonlinear(solver, i, 1);
                    break;
                case MNA_SWITCH: {
                    const double g = comp->state ?
                        1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }
            }
        }

        if (!mna_solve_linear_system(solver, matrix_size)) {
            return 0;
        }

        iteration++;
    }

    return converged;
}

void mna_solve_ac(MNASolver* solver, double frequency) {
    double omega = TWO_PI * frequency;
    int matrix_size = solver->max_node_index + solver->num_sources;

    if (matrix_size > MNA_MAX_MATRIX_SIZE) {
        return;
    }

    // Use local arrays for complex arithmetic
    double complex A_complex[MNA_MAX_MATRIX_SIZE][MNA_MAX_MATRIX_SIZE] = {{0}};
    double complex b_complex[MNA_MAX_MATRIX_SIZE] = {0};
    double complex x_complex[MNA_MAX_MATRIX_SIZE] = {0};
    double complex admittance;
    int source_count = 0;

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        switch (comp->type) {
            case MNA_RESISTOR:
                admittance = 1.0 / comp->value;
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_CAPACITOR:
                admittance = I * omega * comp->value;
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_INDUCTOR:
                if (omega == 0) {
                    admittance = 1.0 / MNA_MIN_CONDUCTANCE;
                } else {
                    admittance = 1.0 / (I * omega * comp->value);
                }
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_SOURCE: {
                double complex source_val = comp->ac_magnitude *
                    (cos(comp->ac_phase) + I * sin(comp->ac_phase));

                if (comp->source_type == SOURCE_VOLTAGE) {
                    int v_index = solver->max_node_index + source_count;
                    if (n1 > 0) {
                        A_complex[n1-1][v_index] = 1.0;
                        A_complex[v_index][n1-1] = 1.0;
                    }
                    if (n2 > 0) {
                        A_complex[n2-1][v_index] = -1.0;
                        A_complex[v_index][n2-1] = -1.0;
                    }
                    b_complex[v_index] = source_val;
                    source_count++;
                } else {
                    if (n1 > 0) b_complex[n1-1] -= source_val;
                    if (n2 > 0) b_complex[n2-1] += source_val;
                }
                break;
            }

            case MNA_CUSTOM_NONLINEAR: {
                switch (comp->nonlinear_type) {
                    case NONLINEAR_RESISTOR:
                        admittance = comp->last_conductance;
                        break;
                    case NONLINEAR_CAPACITOR: {
                        double v0 = comp->last_voltage;
                        double q0, C0;
                        comp->nonlinear_func(solver, i, v0, 0,
                                           &q0, &C0, NONLINEAR_CAPACITOR);
                        admittance = I * omega * C0;
                        break;
                    }
                    case NONLINEAR_INDUCTOR: {
                        double i0 = comp->last_current;
                        double phi0, L0;
                        comp->nonlinear_func(solver, i, 0, i0,
                                           &phi0, &L0, NONLINEAR_INDUCTOR);
                        admittance = 1.0 / (I * omega * L0);
                        break;
                    }
                }

                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;
            }

            case MNA_SWITCH: {
                double g = comp->state ?
                    1.0 / comp->value :
                    MNA_MIN_CONDUCTANCE;

                if (n1 > 0) A_complex[n1-1][n1-1] += g;
                if (n2 > 0) A_complex[n2-1][n2-1] += g;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= g;
                    A_complex[n2-1][n1-1] -= g;
                }
                break;
            }
        }
    }

    // Solve complex linear system
    for (int pivot = 0; pivot < matrix_size; pivot++) {
        int max_row = pivot;
        double max_mag = cabs(A_complex[pivot][pivot]);
        for (int i = pivot + 1; i < matrix_size; i++) {
            double mag = cabs(A_complex[i][pivot]);
            if (mag > max_mag) {
                max_mag = mag;
                max_row = i;
            }
        }

        if (max_mag < MNA_MIN_CONDUCTANCE) {
            return;
        }

        if (max_row != pivot) {
            for (int j = pivot; j < matrix_size; j++) {
                double complex temp = A_complex[pivot][j];
                A_complex[pivot][j] = A_complex[max_row][j];
                A_complex[max_row][j] = temp;
            }
            double complex temp = b_complex[pivot];
            b_complex[pivot] = b_complex[max_row];
            b_complex[max_row] = temp;
        }

        for (int i = pivot + 1; i < matrix_size; i++) {
            double complex factor = A_complex[i][pivot] / A_complex[pivot][pivot];
            for (int j = pivot + 1; j < matrix_size; j++) {
                A_complex[i][j] -= factor * A_complex[pivot][j];
            }
            b_complex[i] -= factor * b_complex[pivot];
        }
    }

    // Back substitution
    for (int i = matrix_size - 1; i >= 0; i--) {
        x_complex[i] = b_complex[i];
        for (int j = i + 1; j < matrix_size; j++) {
            x_complex[i] -= A_complex[i][j] * x_complex[j];
        }
        x_complex[i] /= A_complex[i][i];
    }

    // Store solution
    for (int i = 0; i < matrix_size; i++) {
        solver->ac_solution[i] = x_complex[i];
    }
}

void mna_init_transient(MNASolver* solver) {
    solver->time = 0.0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        comp->last_voltage = 0.0;
        comp->last_current = 0.0;
        comp->last_charge = 0.0;
        comp->last_flux = 0.0;
    }
}

int mna_solve_transient_step(MNASolver* solver, double dt) {
    solver->dt = dt;
    int matrix_size = solver->max_node_index + solver->num_sources;
    mna_reset_system(solver);
    int source_count = 0;

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        switch (comp->type) {
            case MNA_RESISTOR: {
                double g = 1.0 / comp->value;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }

            case MNA_CAPACITOR: {
                double G_eq = comp->value / dt;
                double I_eq = -G_eq * comp->last_voltage;
                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
                break;
            }

            case MNA_INDUCTOR: {
                double G_eq = dt / comp->value;
                double I_eq = comp->last_current;
                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
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

            case MNA_CUSTOM_NONLINEAR: {
                mna_stamp_custom_nonlinear(solver, i, 0);
                break;
            }

            case MNA_SWITCH: {
                double g = comp->state ?
                    1.0 / comp->value :
                    MNA_MIN_CONDUCTANCE;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }
        }
    }

    if (!mna_solve_linear_system(solver, matrix_size)) {
        return 0;
    }

    // Update component states
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
        double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
        double v = v1 - v2;

        switch (comp->type) {
            case MNA_CAPACITOR:
                comp->last_current = (v - comp->last_voltage) * (comp->value / dt);
                comp->last_voltage = v;
                break;

            case MNA_INDUCTOR: {
                double di = v * (dt / comp->value);
                comp->last_current += di;
                break;
            }

            case MNA_CUSTOM_NONLINEAR:
                switch (comp->nonlinear_type) {
                    case NONLINEAR_RESISTOR:
                        comp->last_voltage = v;
                        break;

                    case NONLINEAR_CAPACITOR: {
                        double q, C;
                        comp->nonlinear_func(solver, i, v, 0,
                                           &q, &C, NONLINEAR_CAPACITOR);
                        comp->last_voltage = v;
                        comp->last_charge = q;
                        break;
                    }

                    case NONLINEAR_INDUCTOR: {
                        double i_new = comp->trans_I_eq + comp->trans_G_eq * v;
                        double phi, L;
                        comp->nonlinear_func(solver, i, 0, i_new,
                                           &phi, &L, NONLINEAR_INDUCTOR);
                        comp->last_current = i_new;
                        comp->last_flux = phi;
                        break;
                    }
                }
                break;

            default:
                break;
        }
    }

    solver->time += dt;
    return 1;
}

double mna_get_node_voltage(MNASolver* solver, int node) {
    if (node == 0) return 0.0;
    if (node > 0 && node <= solver->max_node_index) {
        return solver->x[node-1];
    }
    return 0.0;
}

double complex mna_get_ac_node_voltage(MNASolver* solver, int node) {
    if (node == 0) return 0.0;
    if (node > 0 && node <= solver->max_node_index) {
        return solver->ac_solution[node-1];
    }
    return 0.0;
}

#endif // MNA_SOLVER_H
