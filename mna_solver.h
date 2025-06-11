#ifndef MNA_SOLVER_H
#define MNA_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

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

// Forward declaration of MNASolver
typedef struct MNASolver MNASolver;

// Define function pointer type first
typedef void (*CustomNonlinearFunc)(MNASolver* solver, int comp_index,
                                  double voltage, double* current,
                                  double* conductance);

typedef enum {
    MNA_RESISTOR,
    MNA_CAPACITOR,
    MNA_INDUCTOR,
    MNA_VOLTAGE_SOURCE,
    MNA_CURRENT_SOURCE,
    MNA_SWITCH,
    MNA_CUSTOM_NONLINEAR
} ComponentType;

typedef struct {
    ComponentType type;
    int node1;
    int node2;
    double value;
    double ac_magnitude;
    double ac_phase;
    int state;
    int is_nonlinear;
    double last_voltage;
    double last_current;

    // Fields for custom nonlinear elements
    CustomNonlinearFunc nonlinear_func;
    void* user_data;
} Component;

struct MNASolver {
    int num_nodes;
    int num_components;
    int num_sources;
    int num_nonlinear;
    int max_node_index;
    double time;

    Component components[MNA_MAX_COMPONENTS];
    double A[MNA_MAX_NODES + MNA_MAX_SOURCES][MNA_MAX_NODES + MNA_MAX_SOURCES];
    double b[MNA_MAX_NODES + MNA_MAX_SOURCES];
    double x[MNA_MAX_NODES + MNA_MAX_SOURCES];
    double complex ac_solution[MNA_MAX_NODES + MNA_MAX_SOURCES];
};

// Rest of the function prototypes remain the same
void mna_init(MNASolver* solver);
void mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2, double value);
void mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                            CustomNonlinearFunc func, void* user_data,
                            double initial_voltage);
void mna_set_ac_source(MNASolver* solver, int comp_index, double magnitude, double phase);
void mna_set_switch_state(MNASolver* solver, int comp_index, int state);
void mna_reset_system(MNASolver* solver);
void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g);
void mna_stamp_current_source(MNASolver* solver, int node1, int node2, double current_val);
void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx);
void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index);
int mna_solve_linear_system(MNASolver* solver, int size);
int mna_solve_dc(MNASolver* solver);
void mna_solve_ac(MNASolver* solver, double frequency);
double mna_get_node_voltage(MNASolver* solver, int node);
double complex mna_get_ac_node_voltage(MNASolver* solver, int node);
void mna_init_transient(MNASolver* solver);
int mna_solve_transient_step(MNASolver* solver, double dt);

// Implementation
void mna_init(MNASolver* solver) {
    solver->num_nodes = 0;
    solver->num_components = 0;
    solver->num_sources = 0;
    solver->num_nonlinear = 0;
    solver->max_node_index = 0;
    solver->time = 0.0;
}

void mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2, double value) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = type,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .ac_magnitude = 0.0,
        .ac_phase = 0.0,
        .state = 1,
        .is_nonlinear = 0,
        .last_voltage = 0.0,
        .last_current = 0.0,
        .nonlinear_func = NULL,
        .user_data = NULL
    };

    solver->components[solver->num_components++] = comp;

    if (type == MNA_VOLTAGE_SOURCE) solver->num_sources++;
    if (node1 > solver->max_node_index) solver->max_node_index = node1;
    if (node2 > solver->max_node_index) solver->max_node_index = node2;
}

void mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                            CustomNonlinearFunc func, void* user_data,
                            double initial_voltage) {
    if (solver->num_components >= MNA_MAX_COMPONENTS) return;

    Component comp = {
        .type = MNA_CUSTOM_NONLINEAR,
        .node1 = node1,
        .node2 = node2,
        .value = 0.0,
        .ac_magnitude = 0.0,
        .ac_phase = 0.0,
        .state = 1,
        .is_nonlinear = 1,
        .last_voltage = initial_voltage,
        .last_current = 0.0,
        .nonlinear_func = func,
        .user_data = user_data
    };

    solver->components[solver->num_components] = comp;
    solver->num_components++;
    solver->num_nonlinear++;

    if (node1 > solver->max_node_index) solver->max_node_index = node1;
    if (node2 > solver->max_node_index) solver->max_node_index = node2;
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

    for (int i = 0; i < matrix_size; i++) {
        for (int j = 0; j < matrix_size; j++) {
            solver->A[i][j] = 0.0;
        }
        solver->b[i] = 0.0;
        solver->x[i] = 0.0;
    }
}

void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g) {
    if (node1 > 0) solver->A[node1-1][node1-1] += g;
    if (node2 > 0) solver->A[node2-1][node2-1] += g;
    if (node1 > 0 && node2 > 0) {
        solver->A[node1-1][node2-1] -= g;
        solver->A[node2-1][node1-1] -= g;
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
        solver->A[n1-1][v_index] = 1.0;
        solver->A[v_index][n1-1] = 1.0;
    }
    if (n2 > 0) {
        solver->A[n2-1][v_index] = -1.0;
        solver->A[v_index][n2-1] = -1.0;
    }

    solver->b[v_index] = vs->value;
}

void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index) {
    Component* comp = &solver->components[comp_index];
    int n1 = comp->node1;
    int n2 = comp->node2;

    double voltage = comp->last_voltage;
    double current, conductance;
    comp->nonlinear_func(solver, comp_index, voltage, &current, &conductance);

    // Ensure conductance stays within reasonable bounds
    if (conductance < MNA_MIN_CONDUCTANCE) conductance = MNA_MIN_CONDUCTANCE;
    if (conductance > MNA_MAX_CONDUCTANCE) conductance = MNA_MAX_CONDUCTANCE;

    // Stamp the Norton equivalent
    mna_stamp_conductance(solver, n1, n2, conductance);
    mna_stamp_current_source(solver, n1, n2, current - conductance * voltage);
}

int mna_solve_linear_system(MNASolver* solver, int size) {
    for (int pivot = 0; pivot < size; pivot++) {
        int max_row = pivot;
        for (int i = pivot + 1; i < size; i++) {
            if (fabs(solver->A[i][pivot]) > fabs(solver->A[max_row][pivot])) {
                max_row = i;
            }
        }

        if (fabs(solver->A[max_row][pivot]) < MNA_MIN_CONDUCTANCE) {
            return 0;  // Singular matrix
        }

        if (max_row != pivot) {
            for (int j = pivot; j < size; j++) {
                double temp = solver->A[pivot][j];
                solver->A[pivot][j] = solver->A[max_row][j];
                solver->A[max_row][j] = temp;
            }
            double temp = solver->b[pivot];
            solver->b[pivot] = solver->b[max_row];
            solver->b[max_row] = temp;
        }

        for (int i = pivot + 1; i < size; i++) {
            double factor = solver->A[i][pivot] / solver->A[pivot][pivot];
            for (int j = pivot + 1; j < size; j++) {
                solver->A[i][j] -= factor * solver->A[pivot][j];
            }
            solver->b[i] -= factor * solver->b[pivot];
            solver->A[i][pivot] = 0.0;
        }
    }

    for (int i = size - 1; i >= 0; i--) {
        solver->x[i] = solver->b[i];
        for (int j = i + 1; j < size; j++) {
            solver->x[i] -= solver->A[i][j] * solver->x[j];
        }
        solver->x[i] /= solver->A[i][i];
    }

    return 1;
}

int mna_solve_dc(MNASolver* solver) {
    int matrix_size = solver->max_node_index + solver->num_sources;
    int converged = 0;
    int iteration = 0;
    double prev_voltages[MNA_MAX_NONLINEAR] = {0};

    // Store initial nonlinear voltages for convergence tracking
    int nl_idx = 0;
    for (int i = 0; i < solver->num_components; i++) {
        if (solver->components[i].type == MNA_CUSTOM_NONLINEAR) {
            prev_voltages[nl_idx++] = solver->components[i].last_voltage;
        }
    }

    while (!converged && iteration < MNA_MAX_ITER) {
        mna_reset_system(solver);
        int source_count = 0;
        int nonlinear_count = 0;

        // Stamp all components
        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            int n1 = comp->node1;
            int n2 = comp->node2;

            switch (comp->type) {
                case MNA_RESISTOR:
                    mna_stamp_conductance(solver, n1, n2, 1.0 / comp->value);
                    break;

                case MNA_CAPACITOR:
                    // Open circuit for DC
                    break;

                case MNA_INDUCTOR:
                    // Small conductance approximation for DC short
                    mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                    break;

                case MNA_VOLTAGE_SOURCE:
                    mna_stamp_voltage_source(solver, i, source_count++);
                    break;

                case MNA_CURRENT_SOURCE:
                    mna_stamp_current_source(solver, n1, n2, comp->value);
                    break;

                case MNA_CUSTOM_NONLINEAR:
                    mna_stamp_custom_nonlinear(solver, i);
                    nonlinear_count++;
                    break;

                case MNA_SWITCH:
                    if (comp->state) {
                        mna_stamp_conductance(solver, n1, n2, 1.0 / comp->value);
                    } else {
                        mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                    }
                    break;
            }
        }

        // Solve linear system
        if (!mna_solve_linear_system(solver, matrix_size)) {
            return 0;  // Matrix solve failed
        }

        // Update nonlinear components and check convergence
        converged = 1;
        nl_idx = 0;
        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            if (comp->type == MNA_CUSTOM_NONLINEAR) {
                int n1 = comp->node1;
                int n2 = comp->node2;

                // Calculate new voltage
                double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
                double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
                double new_voltage = v1 - v2;

                // Calculate convergence criteria
                double voltage_diff = fabs(new_voltage - comp->last_voltage);
                double abs_tol = MNA_ABSTOL + MNA_RELTOL * fabs(new_voltage);

                // Update component state
                comp->last_voltage = new_voltage;

                // Check if any nonlinear component hasn't converged
                if (voltage_diff > abs_tol) {
                    converged = 0;
                }

                // Update previous voltages for damping
                prev_voltages[nl_idx++] = new_voltage;
            }
        }

        // Handle cases with no nonlinear components
        if (nonlinear_count == 0) {
            converged = 1;
        }

        iteration++;
    }

    return converged;
}

void mna_solve_ac(MNASolver* solver, double frequency) {
    double omega = 2 * M_PI * frequency;
    int matrix_size = solver->max_node_index + solver->num_sources;

    // Check if matrix size is within limits
    if (matrix_size > MNA_MAX_MATRIX_SIZE) {
        printf("Matrix size too large for AC analysis\n");
        return;
    }

    // Initialize complex matrix and vector
    double complex A_complex[MNA_MAX_MATRIX_SIZE][MNA_MAX_MATRIX_SIZE] = {{0}};
    double complex b_complex[MNA_MAX_MATRIX_SIZE] = {0};
    double complex x_complex[MNA_MAX_MATRIX_SIZE] = {0};

    int source_count = 0;  // For indexing voltage sources

    // Stamp components into complex matrix
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;
        double complex admittance;

        switch (comp->type) {
            case MNA_RESISTOR:
                // Stamp conductance (real part)
                admittance = 1.0 / comp->value;
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_CAPACITOR:
                // Stamp capacitive admittance (imaginary part)
                admittance = _Complex_I * omega * comp->value;
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_INDUCTOR:
                // Stamp inductive admittance (imaginary part)
                if (omega == 0) {
                    admittance = 1.0 / MNA_MIN_CONDUCTANCE;
                } else {
                    admittance = 1.0 / (_Complex_I * omega * comp->value);
                }
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_VOLTAGE_SOURCE: {
                // Convert to complex AC value
                double complex voltage = comp->ac_magnitude * (cos(comp->ac_phase) + _Complex_I * sin(comp->ac_phase));
                int v_index = solver->max_node_index + source_count;

                // Stamp voltage source constraints
                if (n1 > 0) {
                    A_complex[n1-1][v_index] = 1.0;
                    A_complex[v_index][n1-1] = 1.0;
                }
                if (n2 > 0) {
                    A_complex[n2-1][v_index] = -1.0;
                    A_complex[v_index][n2-1] = -1.0;
                }
                b_complex[v_index] = voltage;
                source_count++;
                break;
            }

            case MNA_CURRENT_SOURCE: {
                // Convert to complex AC value
                double complex current = comp->ac_magnitude * (cos(comp->ac_phase) + _Complex_I * sin(comp->ac_phase));
                if (n1 > 0) b_complex[n1-1] -= current;
                if (n2 > 0) b_complex[n2-1] += current;
                break;
            }

            case MNA_CUSTOM_NONLINEAR:
                // Skip nonlinear elements in AC analysis for now
                break;

            case MNA_SWITCH: {
                // Use conductance based on switch state
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

    // Solve complex linear system with Gaussian elimination
    int success = 1;
    for (int k = 0; k < matrix_size; k++) {
        // Find pivot
        int max_row = k;
        double max_mag = cabs(A_complex[k][k]);
        for (int i = k+1; i < matrix_size; i++) {
            double mag = cabs(A_complex[i][k]);
            if (mag > max_mag) {
                max_mag = mag;
                max_row = i;
            }
        }

        if (max_mag < MNA_MIN_CONDUCTANCE) {
            success = 0;
            break;
        }

        // Swap rows
        if (max_row != k) {
            for (int j = k; j < matrix_size; j++) {
                double complex temp = A_complex[k][j];
                A_complex[k][j] = A_complex[max_row][j];
                A_complex[max_row][j] = temp;
            }
            double complex temp = b_complex[k];
            b_complex[k] = b_complex[max_row];
            b_complex[max_row] = temp;
        }

        // Eliminate
        for (int i = k+1; i < matrix_size; i++) {
            double complex factor = A_complex[i][k] / A_complex[k][k];
            for (int j = k; j < matrix_size; j++) {
                A_complex[i][j] -= factor * A_complex[k][j];
            }
            b_complex[i] -= factor * b_complex[k];
        }
    }

    if (!success) {
        printf("AC solution failed at frequency %.0f Hz\n", frequency);
        return;
    }

    // Back substitution
    for (int i = matrix_size-1; i >= 0; i--) {
        x_complex[i] = b_complex[i];
        for (int j = i+1; j < matrix_size; j++) {
            x_complex[i] -= A_complex[i][j] * x_complex[j];
        }
        x_complex[i] /= A_complex[i][i];
    }

    // Store results
    for (int i = 0; i < matrix_size; i++) {
        solver->ac_solution[i] = x_complex[i];
    }
}

void mna_init_transient(MNASolver* solver) {
    solver->time = 0.0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        comp->last_voltage = 0.0;  // Reset capacitor voltages
        comp->last_current = 0.0;  // Reset inductor currents
    }
}

int mna_solve_transient_step(MNASolver* solver, double dt) {
    int matrix_size = solver->max_node_index + solver->num_sources;
    mna_reset_system(solver);
    int source_count = 0;

    // Stamp components
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
                double I_eq = G_eq * comp->last_voltage;
                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, -I_eq);
                break;
            }

            case MNA_INDUCTOR: {
                double G_eq = dt / comp->value;
                double I_eq = comp->last_current;
                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
                break;
            }

            case MNA_VOLTAGE_SOURCE: {
                mna_stamp_voltage_source(solver, i, source_count++);
                break;
            }

            case MNA_CURRENT_SOURCE: {
                mna_stamp_current_source(solver, n1, n2, comp->value);
                break;
            }

            case MNA_CUSTOM_NONLINEAR: {
                mna_stamp_custom_nonlinear(solver, i);
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

    // Solve the system
    if (!mna_solve_linear_system(solver, matrix_size)) {
        return 0;  // Matrix solve failed
    }

    // Update component states
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        if (comp->type == MNA_CAPACITOR) {
            double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
            double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
            comp->last_voltage = v1 - v2;
        }
        else if (comp->type == MNA_INDUCTOR) {
            double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
            double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
            double vl = v1 - v2;
            comp->last_current = comp->last_current + (dt / comp->value) * vl;
        }
        else if (comp->type == MNA_CUSTOM_NONLINEAR) {
            double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
            double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
            comp->last_voltage = v1 - v2;
        }
    }

    solver->time += dt;
    return 1;
}

double mna_get_node_voltage(MNASolver* solver, int node) {
    if (node == 0) return 0.0;  // Ground
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
