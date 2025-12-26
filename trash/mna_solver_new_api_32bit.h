#ifndef MNA_SOLVER_H
#define MNA_SOLVER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdalign.h>
#include <string.h>
#include <stdbool.h>

// Error codes
typedef enum {
    MNA_SUCCESS,
    MNA_COMPONENT_LIMIT_REACHED,
    MNA_MATRIX_SINGULAR,
    MNA_CONVERGENCE_FAILURE,
    MNA_INVALID_HANDLE,
    MNA_INVALID_NODE,
    MNA_INSUFFICIENT_MEMORY,
    MNA_INVALID_PARAMETER
} MNAStatus;

// Opaque component handle
typedef int ComponentHandle;

// Optimized constants
#define MNA_DEFAULT_MAX_NODES 50
#define MNA_DEFAULT_MAX_SOURCES 20
#define MNA_DEFAULT_MAX_COMPONENTS 100
#define MNA_MAX_NONLINEAR 20
#define MNA_MAX_ITER 50
#define MNA_RELTOL 1e-6f
#define MNA_ABSTOL 1e-9f
#define MNA_VT 0.02585f
#define MNA_MIN_CONDUCTANCE 1e-12f
#define MNA_MAX_CONDUCTANCE 1e12f
#define TWO_PI 6.28318530717958647692f
#define M_PI 3.14159265358979323846f

// Unified nonlinear element type
typedef enum {
    NONLINEAR_RESISTOR,
    NONLINEAR_CAPACITOR,
    NONLINEAR_INDUCTOR
} NonlinearType;

// Unified source type
typedef enum {
    SOURCE_VOLTAGE,
    SOURCE_CURRENT
} SourceType;

typedef enum {
    MNA_RESISTOR,
    MNA_CAPACITOR,
    MNA_INDUCTOR,
    MNA_SOURCE,
    MNA_SWITCH,
    MNA_CUSTOM_NONLINEAR
} ComponentType;

// Component state for nonlinear callbacks
typedef struct {
    float voltage;
    float current;
    float charge;
    float flux;
    float dt;
} ComponentState;

// Generalized nonlinear function interface
typedef void (*CustomNonlinearFunc)(const ComponentState* state, void* user_data,
                                  float* value1, float* value2);

typedef struct {
    ComponentType type;
    int node1;
    int node2;
    float value;        // DC value
    float ac_magnitude;
    float ac_phase;
    int state;
    bool is_nonlinear;

    // Source-specific fields
    SourceType source_type;

    // Nonlinear-specific fields
    NonlinearType nonlinear_type;
    float last_voltage;
    float last_current;
    float last_conductance;
    float last_charge;
    float last_flux;
    CustomNonlinearFunc nonlinear_func;
    void* user_data;

    // Transient companion model
    float trans_G_eq;
    float trans_I_eq;
} Component;

typedef struct MNASolver {
    int num_nodes;
    int num_components;
    int num_sources;
    int num_nonlinear;
    int max_node_index;
    float time;
    float dt;

    // Dynamic arrays
    Component* components;
    float* A;
    float* b;
    float* x;
    float complex* ac_solution;

    // Capacity tracking
    int cap_nodes;
    int cap_sources;
    int cap_components;
    int matrix_size;
} MNASolver;

// Matrix access macro
#define MAT(solver, i, j) ((solver)->A[(i) * (solver->matrix_size) + (j)])

MNAStatus mna_init_sized(MNASolver* solver, int max_nodes, int max_sources, int max_components) {
    memset(solver, 0, sizeof(MNASolver));

    // Allocate component array
    solver->components = calloc(max_components, sizeof(Component));
    if (!solver->components) return MNA_INSUFFICIENT_MEMORY;

    // Calculate matrix size and allocate
    solver->matrix_size = max_nodes + max_sources;
    size_t matrix_elems = solver->matrix_size * solver->matrix_size;

    solver->A = calloc(matrix_elems, sizeof(float));
    solver->b = calloc(solver->matrix_size, sizeof(float));
    solver->x = calloc(solver->matrix_size, sizeof(float));
    solver->ac_solution = calloc(solver->matrix_size, sizeof(float complex));

    if (!solver->A || !solver->b || !solver->x || !solver->ac_solution) {
        free(solver->components);
        free(solver->A);
        free(solver->b);
        free(solver->x);
        free(solver->ac_solution);
        return MNA_INSUFFICIENT_MEMORY;
    }

    solver->cap_nodes = max_nodes;
    solver->cap_sources = max_sources;
    solver->cap_components = max_components;
    return MNA_SUCCESS;
}

void mna_destroy(MNASolver* solver) {
    if (!solver) return;

    free(solver->components);
    free(solver->A);
    free(solver->b);
    free(solver->x);
    free(solver->ac_solution);

    memset(solver, 0, sizeof(MNASolver));
}

MNAStatus mna_init(MNASolver* solver) {
    return mna_init_sized(solver, MNA_DEFAULT_MAX_NODES,
                         MNA_DEFAULT_MAX_SOURCES,
                         MNA_DEFAULT_MAX_COMPONENTS);
}

int mna_create_node(MNASolver* solver) {
    if (solver->num_nodes >= solver->cap_nodes) return -1;
    return ++solver->max_node_index;
}

MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2) {
    if (node1 < 0 || node1 > solver->max_node_index ||
        node2 < 0 || node2 > solver->max_node_index) {
        return MNA_INVALID_NODE;
    }
    return MNA_SUCCESS;
}

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                          float value, SourceType src_type, ComponentHandle* handle) {
    MNAStatus status = mna_validate_nodes(solver, node1, node2);
    if (status != MNA_SUCCESS) return status;

    if (solver->num_components >= solver->cap_components) {
        return MNA_COMPONENT_LIMIT_REACHED;
    }

    Component comp = {
        .type = type,
        .node1 = node1,
        .node2 = node2,
        .value = value,
        .state = 1,
        .is_nonlinear = false,
        .source_type = src_type,
        .nonlinear_type = NONLINEAR_RESISTOR
    };

    int index = solver->num_components++;
    solver->components[index] = comp;

    // Update max node index
    if (node1 > solver->max_node_index) solver->max_node_index = node1;
    if (node2 > solver->max_node_index) solver->max_node_index = node2;

    if (type == MNA_SOURCE && src_type == SOURCE_VOLTAGE) {
        solver->num_sources++;
    }

    if (handle) *handle = index;
    return MNA_SUCCESS;
}

MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2, float value,
                         ComponentHandle* handle) {
    return mna_add_component(solver, MNA_RESISTOR, node1, node2, value,
                            SOURCE_CURRENT, handle);
}

MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2, float value,
                          ComponentHandle* handle) {
    return mna_add_component(solver, MNA_CAPACITOR, node1, node2, value,
                            SOURCE_CURRENT, handle);
}

MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2, float value,
                         ComponentHandle* handle) {
    return mna_add_component(solver, MNA_INDUCTOR, node1, node2, value,
                            SOURCE_CURRENT, handle);
}

MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2, float value,
                               ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SOURCE, node1, node2, value,
                            SOURCE_VOLTAGE, handle);
}

MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2, float value,
                               ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SOURCE, node1, node2, value,
                            SOURCE_CURRENT, handle);
}

MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2, float value,
                       ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SWITCH, node1, node2, value,
                            SOURCE_CURRENT, handle);
}

MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                                 NonlinearType nl_type,
                                 CustomNonlinearFunc func, void* user_data,
                                 float initial_value1, float initial_value2,
                                 ComponentHandle* handle) {
    MNAStatus status = mna_validate_nodes(solver, node1, node2);
    if (status != MNA_SUCCESS) return status;

    if (solver->num_components >= solver->cap_components) {
        return MNA_COMPONENT_LIMIT_REACHED;
    }

    Component comp = {
        .type = MNA_CUSTOM_NONLINEAR,
        .node1 = node1,
        .node2 = node2,
        .state = 1,
        .is_nonlinear = true,
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

    int index = solver->num_components;
    solver->components[index] = comp;
    solver->num_components++;
    solver->num_nonlinear++;

    // Update max node index
    if (node1 > solver->max_node_index) solver->max_node_index = node1;
    if (node2 > solver->max_node_index) solver->max_node_index = node2;

    if (handle) *handle = index;
    return MNA_SUCCESS;
}

MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle,
                          float magnitude, float phase) {
    if (handle < 0 || handle >= solver->num_components) {
        return MNA_INVALID_HANDLE;
    }

    solver->components[handle].ac_magnitude = magnitude;
    solver->components[handle].ac_phase = phase;
    return MNA_SUCCESS;
}

MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state) {
    if (handle < 0 || handle >= solver->num_components) {
        return MNA_INVALID_HANDLE;
    }

    Component* comp = &solver->components[handle];
    if (comp->type == MNA_SWITCH) {
        comp->state = state;
        return MNA_SUCCESS;
    }
    return MNA_INVALID_HANDLE;
}

void mna_reset_system(MNASolver* solver) {
    size_t matrix_elems = solver->matrix_size * solver->matrix_size;
    memset(solver->A, 0, matrix_elems * sizeof(float));
    memset(solver->b, 0, solver->matrix_size * sizeof(float));
    memset(solver->x, 0, solver->matrix_size * sizeof(float));
}

void mna_stamp_conductance(MNASolver* solver, int node1, int node2, float g) {
    if (node1 > 0) MAT(solver, node1-1, node1-1) += g;
    if (node2 > 0) MAT(solver, node2-1, node2-1) += g;
    if (node1 > 0 && node2 > 0) {
        MAT(solver, node1-1, node2-1) -= g;
        MAT(solver, node2-1, node1-1) -= g;
    }
}

void mna_stamp_current_source(MNASolver* solver, int node1, int node2, float current_val) {
    if (node1 > 0) solver->b[node1-1] -= current_val;
    if (node2 > 0) solver->b[node2-1] += current_val;
}

void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx) {
    Component* vs = &solver->components[comp_index];
    int v_index = solver->max_node_index + source_idx;
    int n1 = vs->node1;
    int n2 = vs->node2;

    if (n1 > 0) {
        MAT(solver, n1-1, v_index) = 1.0f;
        MAT(solver, v_index, n1-1) = 1.0f;
    }
    if (n2 > 0) {
        MAT(solver, n2-1, v_index) = -1.0f;
        MAT(solver, v_index, n2-1) = -1.0f;
    }

    solver->b[v_index] = vs->value;
}

void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc) {
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

    switch (comp->nonlinear_type) {
        case NONLINEAR_RESISTOR: {
            float current, conductance;
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
                float q0, C0;
                comp->nonlinear_func(&state, comp->user_data, &q0, &C0);

                float G_eq = C0 / solver->dt;
                float I_eq = (q0 - comp->last_charge) / solver->dt - G_eq * state.voltage;

                comp->trans_G_eq = G_eq;
                comp->trans_I_eq = I_eq;

                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
            }
            break;
        }

        case NONLINEAR_INDUCTOR: {
            if (is_dc) {
                mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
            } else {
                float phi0, L0;
                comp->nonlinear_func(&state, comp->user_data, &phi0, &L0);

                float G_eq = solver->dt / L0;
                float I_eq = state.current + (phi0 - comp->last_flux) / L0;

                comp->trans_G_eq = G_eq;
                comp->trans_I_eq = I_eq;

                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
            }
            break;
        }
    }
}

MNAStatus mna_solve_linear_system(MNASolver* solver, int size) {
    const float pivot_threshold = 1e-12f;
    float max_val;
    int max_row;

    for (int pivot = 0; pivot < size; pivot++) {
        // Find pivot
        max_row = pivot;
        max_val = fabsf(MAT(solver, pivot, pivot));
        for (int i = pivot + 1; i < size; i++) {
            float abs_val = fabsf(MAT(solver, i, pivot));
            if (abs_val > max_val) {
                max_val = abs_val;
                max_row = i;
            }
        }

        // Singularity check
        if (max_val < pivot_threshold) {
            return MNA_MATRIX_SINGULAR;
        }

        // Swap rows if needed
        if (max_row != pivot) {
            for (int j = pivot; j < size; j++) {
                float temp = MAT(solver, pivot, j);
                MAT(solver, pivot, j) = MAT(solver, max_row, j);
                MAT(solver, max_row, j) = temp;
            }
            float temp = solver->b[pivot];
            solver->b[pivot] = solver->b[max_row];
            solver->b[max_row] = temp;
        }

        // Elimination
        const float pivot_inv = 1.0f / MAT(solver, pivot, pivot);
        for (int i = pivot + 1; i < size; i++) {
            float factor = MAT(solver, i, pivot) * pivot_inv;
            for (int j = pivot + 1; j < size; j++) {
                MAT(solver, i, j) -= factor * MAT(solver, pivot, j);
            }
            solver->b[i] -= factor * solver->b[pivot];
            MAT(solver, i, pivot) = 0.0f;
        }
    }

    // Back substitution
    for (int i = size - 1; i >= 0; i--) {
        float sum = solver->b[i];
        for (int j = i + 1; j < size; j++) {
            sum -= MAT(solver, i, j) * solver->x[j];
        }
        solver->x[i] = sum / MAT(solver, i, i);
    }

    return MNA_SUCCESS;
}

MNAStatus mna_solve_dc(MNASolver* solver) {
    const int matrix_size = solver->max_node_index + solver->num_sources;
    const bool has_nonlinear = (solver->num_nonlinear > 0);

    if (!has_nonlinear) {
        // Standard linear solution for circuits without nonlinear components
        mna_reset_system(solver);
        int source_count = 0;

        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            const int n1 = comp->node1;
            const int n2 = comp->node2;

            switch (comp->type) {
                case MNA_RESISTOR:
                    mna_stamp_conductance(solver, n1, n2, 1.0f / comp->value);
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
                    // Shouldn't happen but handle safely
                    mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                    break;
                case MNA_SWITCH: {
                    const float g = comp->state ?
                        1.0f / comp->value : MNA_MIN_CONDUCTANCE;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }
            }
        }
        return mna_solve_linear_system(solver, matrix_size);
    }

    // Nonlinear circuit - implement source stepping with adaptive step size
    float* orig_source_values = (float*)malloc(solver->num_components * sizeof(float));
    if (!orig_source_values) return MNA_INSUFFICIENT_MEMORY;

    // Save original source values and set initial state to zero
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->type == MNA_SOURCE) {
            orig_source_values[i] = comp->value;
            comp->value = 0.0f;  // Start from zero
        }
        if (comp->type == MNA_CUSTOM_NONLINEAR) {
            comp->last_voltage = 0.0f;  // Reset to known state
            comp->last_conductance = MNA_MIN_CONDUCTANCE;
        }
    }

    // Source stepping parameters
    float current_factor = 0.0f;
    const float min_step = 0.01f;  // Minimum step size (1%)
    const float max_step = 0.2f;   // Maximum step size (20%)
    float step_size = max_step;
    MNAStatus status = MNA_SUCCESS;
    int total_steps = 0;
    const int max_total_steps = 200;

    while (current_factor < 1.0f && total_steps < max_total_steps) {
        float next_factor = current_factor + step_size;
        if (next_factor > 1.0f) next_factor = 1.0f;

        // Update sources to current factor
        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            if (comp->type == MNA_SOURCE) {
                comp->value = orig_source_values[i] * next_factor;
            }
        }

        // Solve at current source level
        mna_reset_system(solver);
        int source_count = 0;

        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            const int n1 = comp->node1;
            const int n2 = comp->node2;

            switch (comp->type) {
                case MNA_RESISTOR:
                    mna_stamp_conductance(solver, n1, n2, 1.0f / comp->value);
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
                    const float g = comp->state ?
                        1.0f / comp->value :
                        MNA_MIN_CONDUCTANCE;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }
            }
        }

        // Solve linear system to get initial guess
        status = mna_solve_linear_system(solver, matrix_size);
        if (status != MNA_SUCCESS) break;

        // Newton-Raphson iterations at current source level
        int converged = 0;
        int iteration = 0;

        while (!converged && iteration < MNA_MAX_ITER) {
            converged = 1;

            // Update nonlinear component states
            for (int i = 0; i < solver->num_components; i++) {
                Component* comp = &solver->components[i];
                if (comp->type != MNA_CUSTOM_NONLINEAR) continue;

                const int n1 = comp->node1;
                const int n2 = comp->node2;
                const float v1 = (n1 > 0) ? solver->x[n1-1] : 0.0f;
                const float v2 = (n2 > 0) ? solver->x[n2-1] : 0.0f;
                const float new_voltage = v1 - v2;
                const float voltage_diff = fabsf(new_voltage - comp->last_voltage);
                const float abs_tol = MNA_ABSTOL + MNA_RELTOL * fabsf(new_voltage);

                // Damping for convergence (0.5 for first 5 iterations, 0.9 after)
                const float damping = (iteration < 5) ? 0.5f : 0.9f;
                comp->last_voltage = damping * new_voltage +
                                    (1 - damping) * comp->last_voltage;

                // Update nonlinear parameters
                if (comp->nonlinear_type == NONLINEAR_RESISTOR) {
                    ComponentState state = {
                        .voltage = comp->last_voltage,
                        .current = comp->last_current,
                        .charge = comp->last_charge,
                        .flux = comp->last_flux,
                        .dt = solver->dt
                    };
                    float val1, val2;
                    comp->nonlinear_func(&state, comp->user_data, &val1, &val2);
                    comp->last_conductance =
                        (val2 < MNA_MIN_CONDUCTANCE) ? MNA_MIN_CONDUCTANCE :
                        (val2 > MNA_MAX_CONDUCTANCE) ? MNA_MAX_CONDUCTANCE : val2;
                }

                if (voltage_diff > abs_tol) {
                    converged = 0;
                }
            }

            if (converged) break;

            // Re-stamp nonlinear components with updated values
            mna_reset_system(solver);
            source_count = 0;

            for (int i = 0; i < solver->num_components; i++) {
                Component* comp = &solver->components[i];
                const int n1 = comp->node1;
                const int n2 = comp->node2;

                switch (comp->type) {
                    case MNA_RESISTOR:
                        mna_stamp_conductance(solver, n1, n2, 1.0f / comp->value);
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
                        const float g = comp->state ?
                            1.0f / comp->value : MNA_MIN_CONDUCTANCE;
                        mna_stamp_conductance(solver, n1, n2, g);
                        break;
                    }
                }
            }

            status = mna_solve_linear_system(solver, matrix_size);
            if (status != MNA_SUCCESS) break;

            iteration++;
        }

        if (status != MNA_SUCCESS || !converged) {
            // Reduce step size and retry
            step_size /= 2.0f;
            if (step_size < min_step) {
                status = MNA_CONVERGENCE_FAILURE;
                break;
            }
            continue;  // Retry with smaller step
        }

        // Successful step - move to next source level
        current_factor = next_factor;
        total_steps++;

        // Increase step size for next step (with upper bound)
        step_size *= 1.5f;
        if (step_size > max_step) step_size = max_step;
    }

    // Restore original source values
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->type == MNA_SOURCE) {
            comp->value = orig_source_values[i];
        }
    }
    free(orig_source_values);

    return status;
}

MNAStatus mna_solve_ac(MNASolver* solver, float frequency) {
    float omega = TWO_PI * frequency;
    int matrix_size = solver->max_node_index + solver->num_sources;

    if (matrix_size > solver->matrix_size) {
        return MNA_INVALID_PARAMETER;
    }

    // Use local arrays for complex arithmetic
    float complex (*A_complex)[matrix_size] = calloc(matrix_size * matrix_size, sizeof(float complex));
    float complex *b_complex = calloc(matrix_size, sizeof(float complex));
    float complex *x_complex = calloc(matrix_size, sizeof(float complex));

    if (!A_complex || !b_complex || !x_complex) {
        free(A_complex);
        free(b_complex);
        free(x_complex);
        return MNA_INSUFFICIENT_MEMORY;
    }

    float complex admittance;
    int source_count = 0;

    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        switch (comp->type) {
            case MNA_RESISTOR:
                admittance = 1.0f / comp->value;
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
                    admittance = 1.0f / MNA_MIN_CONDUCTANCE;
                } else {
                    admittance = 1.0f / (I * omega * comp->value);
                }
                if (n1 > 0) A_complex[n1-1][n1-1] += admittance;
                if (n2 > 0) A_complex[n2-1][n2-1] += admittance;
                if (n1 > 0 && n2 > 0) {
                    A_complex[n1-1][n2-1] -= admittance;
                    A_complex[n2-1][n1-1] -= admittance;
                }
                break;

            case MNA_SOURCE: {
                float complex source_val = comp->ac_magnitude *
                    (cosf(comp->ac_phase) + I * sinf(comp->ac_phase));

                if (comp->source_type == SOURCE_VOLTAGE) {
                    int v_index = solver->max_node_index + source_count;
                    if (n1 > 0) {
                        A_complex[n1-1][v_index] = 1.0f;
                        A_complex[v_index][n1-1] = 1.0f;
                    }
                    if (n2 > 0) {
                        A_complex[n2-1][v_index] = -1.0f;
                        A_complex[v_index][n2-1] = -1.0f;
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
                        ComponentState state = {
                            .voltage = comp->last_voltage,
                            .current = comp->last_current,
                            .charge = comp->last_charge,
                            .flux = comp->last_flux,
                            .dt = solver->dt
                        };
                        float q0, C0;
                        comp->nonlinear_func(&state, comp->user_data, &q0, &C0);
                        admittance = I * omega * C0;
                        break;
                    }
                    case NONLINEAR_INDUCTOR: {
                        ComponentState state = {
                            .voltage = comp->last_voltage,
                            .current = comp->last_current,
                            .charge = comp->last_charge,
                            .flux = comp->last_flux,
                            .dt = solver->dt
                        };
                        float phi0, L0;
                        comp->nonlinear_func(&state, comp->user_data, &phi0, &L0);
                        admittance = 1.0f / (I * omega * L0);
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
                float g = comp->state ?
                    1.0f / comp->value :
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
        float max_mag = cabsf(A_complex[pivot][pivot]);
        for (int i = pivot + 1; i < matrix_size; i++) {
            float mag = cabsf(A_complex[i][pivot]);
            if (mag > max_mag) {
                max_mag = mag;
                max_row = i;
            }
        }

        if (max_mag < MNA_MIN_CONDUCTANCE) {
            free(A_complex);
            free(b_complex);
            free(x_complex);
            return MNA_MATRIX_SINGULAR;
        }

        if (max_row != pivot) {
            for (int j = pivot; j < matrix_size; j++) {
                float complex temp = A_complex[pivot][j];
                A_complex[pivot][j] = A_complex[max_row][j];
                A_complex[max_row][j] = temp;
            }
            float complex temp = b_complex[pivot];
            b_complex[pivot] = b_complex[max_row];
            b_complex[max_row] = temp;
        }

        for (int i = pivot + 1; i < matrix_size; i++) {
            float complex factor = A_complex[i][pivot] / A_complex[pivot][pivot];
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

    free(A_complex);
    free(b_complex);
    free(x_complex);
    return MNA_SUCCESS;
}

void mna_init_transient(MNASolver* solver) {
    solver->time = 0.0f;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        comp->last_voltage = 0.0f;
        comp->last_current = 0.0f;
        comp->last_charge = 0.0f;
        comp->last_flux = 0.0f;
    }
}

MNAStatus mna_solve_transient_step(MNASolver* solver, float dt) {
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
                float g = 1.0f / comp->value;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }

            case MNA_CAPACITOR: {
                float G_eq = comp->value / dt;
                float I_eq = -G_eq * comp->last_voltage;
                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
                break;
            }

            case MNA_INDUCTOR: {
                float G_eq = dt / comp->value;
                float I_eq = comp->last_current;
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
                float g = comp->state ?
                    1.0f / comp->value :
                    MNA_MIN_CONDUCTANCE;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }
        }
    }

    MNAStatus status = mna_solve_linear_system(solver, matrix_size);
    if (status != MNA_SUCCESS) {
        return status;
    }

    // Update component states
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;

        float v1 = (n1 > 0) ? solver->x[n1-1] : 0.0f;
        float v2 = (n2 > 0) ? solver->x[n2-1] : 0.0f;
        float v = v1 - v2;

        switch (comp->type) {
            case MNA_CAPACITOR:
                comp->last_current = (v - comp->last_voltage) * (comp->value / dt);
                comp->last_voltage = v;
                break;

            case MNA_INDUCTOR: {
                float di = v * (dt / comp->value);
                comp->last_current += di;
                break;
            }

            case MNA_CUSTOM_NONLINEAR:
                switch (comp->nonlinear_type) {
                    case NONLINEAR_RESISTOR:
                        comp->last_voltage = v;
                        break;

                    case NONLINEAR_CAPACITOR: {
                        ComponentState state = {
                            .voltage = v,
                            .current = comp->last_current,
                            .charge = comp->last_charge,
                            .flux = comp->last_flux,
                            .dt = dt
                        };
                        float q, C;
                        comp->nonlinear_func(&state, comp->user_data, &q, &C);
                        comp->last_voltage = v;
                        comp->last_charge = q;
                        break;
                    }

                    case NONLINEAR_INDUCTOR: {
                        float i_new = comp->trans_I_eq + comp->trans_G_eq * v;
                        ComponentState state = {
                            .voltage = v,
                            .current = i_new,
                            .charge = comp->last_charge,
                            .flux = comp->last_flux,
                            .dt = dt
                        };
                        float phi, L;
                        comp->nonlinear_func(&state, comp->user_data, &phi, &L);
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
    return MNA_SUCCESS;
}

float mna_get_node_voltage(MNASolver* solver, int node) {
    if (node <= 0 || node > solver->max_node_index) return 0.0f;
    return solver->x[node-1];
}

float complex mna_get_ac_node_voltage(MNASolver* solver, int node) {
    if (node <= 0 || node > solver->max_node_index) return 0.0f;
    return solver->ac_solution[node-1];
}

float mna_get_component_current(MNASolver* solver, ComponentHandle handle) {
    if (handle < 0 || handle >= solver->num_components) return 0.0f;

    Component* comp = &solver->components[handle];
    int n1 = comp->node1;
    int n2 = comp->node2;

    float v1 = (n1 > 0) ? mna_get_node_voltage(solver, n1) : 0.0f;
    float v2 = (n2 > 0) ? mna_get_node_voltage(solver, n2) : 0.0f;
    float v = v1 - v2;

    switch (comp->type) {
        case MNA_RESISTOR:
            return v / comp->value;
        case MNA_CAPACITOR:
            return (v - comp->last_voltage) * (comp->value / solver->dt);
        case MNA_INDUCTOR:
            return comp->last_current;
        case MNA_SOURCE:
            if (comp->source_type == SOURCE_VOLTAGE) {
                // Find voltage source index
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
            return v / (comp->state ? comp->value : 1.0f/MNA_MIN_CONDUCTANCE);
        case MNA_CUSTOM_NONLINEAR: {
            ComponentState state = {
                .voltage = v,
                .current = comp->last_current,
                .charge = comp->last_charge,
                .flux = comp->last_flux,
                .dt = solver->dt
            };
            float current, conductance;
            comp->nonlinear_func(&state, comp->user_data, &current, &conductance);
            return current;
        }
        default:
            return 0.0f;
    }
}

float mna_get_component_voltage(MNASolver* solver, ComponentHandle handle) {
    if (handle < 0 || handle >= solver->num_components) return 0.0f;

    Component* comp = &solver->components[handle];
    int n1 = comp->node1;
    int n2 = comp->node2;

    float v1 = (n1 > 0) ? mna_get_node_voltage(solver, n1) : 0.0f;
    float v2 = (n2 > 0) ? mna_get_node_voltage(solver, n2) : 0.0f;
    return v1 - v2;
}

#endif // MNA_SOLVER_H
