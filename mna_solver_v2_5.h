#ifndef MNA_SOLVER_V2_5_H
#define MNA_SOLVER_V2_5_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdbool.h>

/* -------------------------------------------------------------------------- */
/*                                 Constants                                  */
/* -------------------------------------------------------------------------- */
#define MNA_MAX_ITER           50
#define MNA_RELTOL             1e-6
#define MNA_ABSTOL             1e-9
#define MNA_VT                 0.02585
#define MNA_MIN_CONDUCTANCE    1e-12
#define MNA_MAX_CONDUCTANCE    1e12
#define MNA_GROUND_CONDUCTANCE 1e-9
#define TWO_PI                 6.28318530717958647692
#define MNA_INIT_CAPACITY      4

/* Oscillation Detection Configuration */
#define MNA_EULER_COOLDOWN_STEPS 10      /* Steps to force Euler after oscillation */
#define MNA_OSCILLATION_TOL_REL  1e-6       /* Relative tolerance for oscillation */

/* -------------------------------------------------------------------------- */
/*                                 Types                                      */
/* -------------------------------------------------------------------------- */
typedef enum {
    MNA_SUCCESS,
    MNA_MATRIX_SINGULAR,
    MNA_CONVERGENCE_FAILURE,
    MNA_INVALID_HANDLE,
    MNA_INVALID_NODE,
    MNA_INSUFFICIENT_MEMORY,
    MNA_INVALID_PARAMETER,
    MNA_OSCILLATION_DETECTED
} MNAStatus;

typedef int ComponentHandle;

typedef enum {
    NONLINEAR_RESISTOR,
    NONLINEAR_CAPACITOR,
    NONLINEAR_INDUCTOR
} NonlinearType;

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
    MNA_CUSTOM_NONLINEAR,
    MNA_CUSTOM_NPOLE
} ComponentType;

typedef enum {
    MNA_INTEGRATION_TRAPEZOIDAL,
    MNA_INTEGRATION_EULER
} IntegrationMethod;

typedef struct {
    double voltage;
    double current;
    double charge;
    double flux;
    double dt;
} ComponentState;

typedef void (*CustomNonlinearFunc)(const ComponentState* state, void* user_data,
                                    double* value1, double* value2);

struct MNASolver;
typedef void (*NPoleStampFunc)(struct MNASolver* solver,
                               const int* nodes,
                               int num_nodes,
                               void* user_data,
                               double time,
                               double dt);

typedef struct {
    int* nodes;
    int num_nodes;
    NPoleStampFunc stamp_func;
    void* user_data;
    double* last_values;
} NPoleData;

typedef struct {
    ComponentType type;
    int node1;
    int node2;
    double value;
    double ac_magnitude;
    double ac_phase;
    int state;
    bool is_nonlinear;
    SourceType source_type;
    NonlinearType nonlinear_type;
    /* State variables for transient analysis */
    double last_voltage;
    double last_current;
    double prev_voltage;
    double prev_current;
    double last_conductance;
    double last_charge;
    double last_flux;
    CustomNonlinearFunc nonlinear_func;
    void* user_data;
    double trans_G_eq;
    double trans_I_eq;
    union {
        struct {
            NPoleData* npole_data;
        } npole;
    } data;
} Component;

typedef struct MNASolver {
    int num_nodes;
    int num_components;
    int num_sources;
    int num_nonlinear;
    int max_node_index;
    int transient_initialized;
    double time;
    double dt;
    IntegrationMethod integration_method;
    bool bypass_oscillation_check;
    int euler_cooldown_steps; /* Persistent cooldown counter */
    Component* components;
    double* A;
    double* b;
    double* x;
    double complex* ac_solution;
    int cap_components;
    int matrix_cap_size;
} MNASolver;

/* -------------------------------------------------------------------------- */
/*                             Function Prototypes                            */
/* -------------------------------------------------------------------------- */
void mna_destroy(MNASolver* solver);
MNAStatus mna_init(MNASolver* solver);
MNAStatus mna_set_integration_method(MNASolver* solver, IntegrationMethod method);
MNAStatus mna_set_oscillation_bypass(MNASolver* solver, bool enable);
int mna_create_node(MNASolver* solver);
MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle);
MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2, NonlinearType nl_type,
                                   CustomNonlinearFunc func, void* user_data,
                                   double initial_value1, double initial_value2, ComponentHandle* handle);
MNAStatus mna_add_custom_n_pole(MNASolver* solver, const int* nodes, int num_nodes,
                                NPoleStampFunc stamp_func, void* user_data, ComponentHandle* handle);
MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle, double magnitude, double phase);
MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state);
MNAStatus mna_solve_dc(MNASolver* solver);
MNAStatus mna_solve_ac(MNASolver* solver, double frequency);
void mna_init_transient(MNASolver* solver);
MNAStatus mna_solve_transient_step(MNASolver* solver, double dt);
double mna_get_node_voltage(MNASolver* solver, int node);
double complex mna_get_ac_node_voltage(MNASolver* solver, int node);
double mna_get_component_current(MNASolver* solver, ComponentHandle handle);
double mna_get_component_voltage(MNASolver* solver, ComponentHandle handle);

/* -------------------------------------------------------------------------- */
/*                             Internal Utilities                             */
/* -------------------------------------------------------------------------- */
static inline int mna_active_size(const MNASolver* s) {
    return s->max_node_index + s->num_sources;
}

#define MAT(solver, i, j) ((solver)->A[(i) * (solver)->matrix_cap_size + (j)])

static bool mna_resize_components(MNASolver* s, int required) {
    if (required <= s->cap_components) return true;
    int new_cap = s->cap_components ? s->cap_components : MNA_INIT_CAPACITY;
    while (new_cap < required) {
        new_cap = (int)(new_cap * 1.5) + 1;
    }
    Component* nxt = (Component*)realloc(s->components, (size_t)new_cap * sizeof(Component));
    if (!nxt) return false;
    if (new_cap > s->cap_components) {
        size_t added = (size_t)(new_cap - s->cap_components);
        memset(nxt + s->cap_components, 0, added * sizeof(Component));
    }
    s->components = nxt;
    s->cap_components = new_cap;
    return true;
}

static bool mna_resize_matrix(MNASolver* s, int req_size) {
    if (req_size <= s->matrix_cap_size) return true;
    int new_size = s->matrix_cap_size ? s->matrix_cap_size : MNA_INIT_CAPACITY;
    while (new_size < req_size) {
        new_size = (int)(new_size * 1.5) + 1;
    }
    size_t new_elems = (size_t)new_size * (size_t)new_size;
    size_t vec_elems = (size_t)new_size;
    double* A2 = (double*)calloc(new_elems, sizeof(double));
    double* b2 = (double*)calloc(vec_elems, sizeof(double));
    double* x2 = (double*)calloc(vec_elems, sizeof(double));
    double complex* ac2 = (double complex*)calloc(vec_elems, sizeof(double complex));
    if (!A2 || !b2 || !x2 || !ac2) {
        free(A2); free(b2); free(x2); free(ac2);
        return false;
    }
    int old_size = s->matrix_cap_size;
    if (s->A && old_size > 0) {
        for (int i = 0; i < old_size; ++i) {
            memcpy(A2 + (size_t)i * new_size, s->A + (size_t)i * old_size, (size_t)old_size * sizeof(double));
        }
        memcpy(b2, s->b, (size_t)old_size * sizeof(double));
        memcpy(x2, s->x, (size_t)old_size * sizeof(double));
        memcpy(ac2, s->ac_solution, (size_t)old_size * sizeof(double complex));
    }
    free(s->A); free(s->b); free(s->x); free(s->ac_solution);
    s->A = A2;
    s->b = b2;
    s->x = x2;
    s->ac_solution = ac2;
    s->matrix_cap_size = new_size;
    return true;
}

/* -------------------------------------------------------------------------- */
/*                             Implementation                                 */
/* -------------------------------------------------------------------------- */
MNAStatus mna_init(MNASolver* solver) {
    if (!solver) return MNA_INVALID_PARAMETER;
    memset(solver, 0, sizeof(MNASolver));
    solver->cap_components = MNA_INIT_CAPACITY;
    solver->matrix_cap_size = MNA_INIT_CAPACITY;
    solver->integration_method = MNA_INTEGRATION_TRAPEZOIDAL;
    solver->bypass_oscillation_check = false;
    solver->euler_cooldown_steps = 0;
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

MNAStatus mna_set_integration_method(MNASolver* solver, IntegrationMethod method) {
    if (!solver) return MNA_INVALID_PARAMETER;
    solver->integration_method = method;
    solver->euler_cooldown_steps = 0; /* Reset cooldown on manual change */
    return MNA_SUCCESS;
}

MNAStatus mna_set_oscillation_bypass(MNASolver* solver, bool enable) {
    if (!solver) return MNA_INVALID_PARAMETER;
    solver->bypass_oscillation_check = enable;
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

MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_RESISTOR, node1, node2, value, SOURCE_CURRENT, handle);
}

MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_CAPACITOR, node1, node2, value, SOURCE_CURRENT, handle);
}

MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_INDUCTOR, node1, node2, value, SOURCE_CURRENT, handle);
}

MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SOURCE, node1, node2, value, SOURCE_VOLTAGE, handle);
}

MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SOURCE, node1, node2, value, SOURCE_CURRENT, handle);
}

MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle) {
    return mna_add_component(solver, MNA_SWITCH, node1, node2, value, SOURCE_CURRENT, handle);
}

MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2, NonlinearType nl_type,
                                   CustomNonlinearFunc func, void* user_data,
                                   double initial_value1, double initial_value2, ComponentHandle* handle) {
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
    if (nl_type == NONLINEAR_RESISTOR) comp.last_voltage = initial_value1;
    else if (nl_type == NONLINEAR_CAPACITOR) {
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
                                NPoleStampFunc stamp_func, void* user_data, ComponentHandle* handle) {
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

MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle, double magnitude, double phase) {
    if (!solver || handle < 0 || handle >= solver->num_components) return MNA_INVALID_HANDLE;
    solver->components[handle].ac_magnitude = magnitude;
    solver->components[handle].ac_phase = phase;
    return MNA_SUCCESS;
}

MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state) {
    if (!solver || handle < 0 || handle >= solver->num_components) return MNA_INVALID_HANDLE;
    Component* comp = &solver->components[handle];
    if (comp->type == MNA_SWITCH) {
        comp->state = state;
        return MNA_SUCCESS;
    }
    return MNA_INVALID_PARAMETER;
}

static void mna_reset_system(MNASolver* solver) {
    int active = mna_active_size(solver);
    for (int i = 0; i < active; ++i) {
        memset(&MAT(solver, i, 0), 0, (size_t)active * sizeof(double));
    }
    memset(solver->b, 0, (size_t)active * sizeof(double));
    memset(solver->x, 0, (size_t)active * sizeof(double));
}

static void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g) {
    if (node1 > 0) MAT(solver, node1-1, node1-1) += g;
    if (node2 > 0) MAT(solver, node2-1, node2-1) += g;
    if (node1 > 0 && node2 > 0) {
        MAT(solver, node1-1, node2-1) -= g;
        MAT(solver, node2-1, node1-1) -= g;
    }
}

static void mna_stamp_current_source(MNASolver* solver, int node1, int node2, double current_val) {
    if (node1 > 0) solver->b[node1-1] -= current_val;
    if (node2 > 0) solver->b[node2-1] += current_val;
}

static void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx) {
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

static void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc, IntegrationMethod method) {
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
                if (method == MNA_INTEGRATION_TRAPEZOIDAL) {
                    G_eq = (2.0 * C0) / solver->dt;
                    I_eq = (G_eq * state.voltage) + state.current;
                } else {
                    G_eq = C0 / solver->dt;
                    I_eq = -G_eq * state.voltage;
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
                if (method == MNA_INTEGRATION_TRAPEZOIDAL) {
                    G_eq = solver->dt / (2.0 * L0);
                    I_eq = state.current + (G_eq * state.voltage);
                } else {
                    G_eq = solver->dt / L0;
                    I_eq = state.current;
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

static void mna_ensure_ground_paths(MNASolver* solver) {
    int* connected = (int*)calloc((size_t)(solver->max_node_index + 1), sizeof(int));
    if (!connected) return;
    connected[0] = 1;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->node1 == 0 && comp->node2 > 0) connected[comp->node2] = 1;
        else if (comp->node2 == 0 && comp->node1 > 0) connected[comp->node1] = 1;
    }
    for (int node = 1; node <= solver->max_node_index; node++) {
        if (!connected[node]) {
            mna_stamp_conductance(solver, node, 0, MNA_GROUND_CONDUCTANCE);
        }
    }
    free(connected);
}

static MNAStatus mna_solve_linear_system(MNASolver* solver, int size) {
    const double pivot_threshold = 1e-12;
    for (int pivot = 0; pivot < size; pivot++) {
        int max_row = pivot;
        double max_val = fabs(MAT(solver, pivot, pivot));
        for (int i = pivot + 1; i < size; i++) {
            double abs_val = fabs(MAT(solver, i, pivot));
            if (abs_val > max_val) {
                max_val = abs_val;
                max_row = i;
            }
        }
        if (max_val < pivot_threshold) return MNA_MATRIX_SINGULAR;
        if (max_row != pivot) {
            for (int j = pivot; j < size; j++) {
                double temp = MAT(solver, pivot, j);
                MAT(solver, pivot, j) = MAT(solver, max_row, j);
                MAT(solver, max_row, j) = temp;
            }
            double tempb = solver->b[pivot];
            solver->b[pivot] = solver->b[max_row];
            solver->b[max_row] = tempb;
        }
        const double pivot_inv = 1.0 / MAT(solver, pivot, pivot);
        for (int i = pivot + 1; i < size; i++) {
            double factor = MAT(solver, i, pivot) * pivot_inv;
            if (factor == 0.0) continue;
            for (int j = pivot + 1; j < size; j++) {
                MAT(solver, i, j) -= factor * MAT(solver, pivot, j);
            }
            solver->b[i] -= factor * solver->b[pivot];
            MAT(solver, i, pivot) = 0.0;
        }
    }
    for (int i = size - 1; i >= 0; i--) {
        double sum = solver->b[i];
        for (int j = i + 1; j < size; j++) {
            sum -= MAT(solver, i, j) * solver->x[j];
        }
        solver->x[i] = sum / MAT(solver, i, i);
    }
    return MNA_SUCCESS;
}

MNAStatus mna_solve_dc(MNASolver* solver) {
    if (!solver) return MNA_INVALID_HANDLE;
    const int matrix_size = mna_active_size(solver);
    const bool has_nonlinear = (solver->num_nonlinear > 0);
    if (!has_nonlinear) {
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
                    if (comp->source_type == SOURCE_VOLTAGE) mna_stamp_voltage_source(solver, i, source_count++);
                    else mna_stamp_current_source(solver, n1, n2, comp->value);
                    break;
                case MNA_CUSTOM_NONLINEAR:
                    mna_stamp_conductance(solver, n1, n2, MNA_MIN_CONDUCTANCE);
                    break;
                case MNA_CUSTOM_NPOLE: {
                    NPoleData* npole = comp->data.npole.npole_data;
                    if (npole) npole->stamp_func(solver, npole->nodes, npole->num_nodes, npole->user_data, 0.0, 0.0);
                    break;
                }
                case MNA_SWITCH: {
                    const double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }
                default: break;
            }
        }
        return mna_solve_linear_system(solver, matrix_size);
    }
    double* orig_source_values = (double*)malloc((size_t)solver->num_components * sizeof(double));
    if (!orig_source_values) return MNA_INSUFFICIENT_MEMORY;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
            memset(comp->data.npole.npole_data->last_values, 0, (size_t)comp->data.npole.npole_data->num_nodes * sizeof(double));
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
            if (comp->type == MNA_SOURCE) comp->value = orig_source_values[i] * next_factor;
        }
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
                    if (comp->source_type == SOURCE_VOLTAGE) mna_stamp_voltage_source(solver, i, source_count++);
                    else mna_stamp_current_source(solver, n1, n2, comp->value);
                    break;
                case MNA_CUSTOM_NONLINEAR:
                    mna_stamp_custom_nonlinear(solver, i, 1, solver->integration_method);
                    break;
                case MNA_CUSTOM_NPOLE: {
                    NPoleData* npole = comp->data.npole.npole_data;
                    if (npole) npole->stamp_func(solver, npole->nodes, npole->num_nodes, npole->user_data, 0.0, 0.0);
                    break;
                }
                case MNA_SWITCH: {
                    const double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                    mna_stamp_conductance(solver, n1, n2, g);
                    break;
                }
                default: break;
            }
        }
        status = mna_solve_linear_system(solver, matrix_size);
        if (status != MNA_SUCCESS) break;
        int converged = 0;
        int iteration = 0;
        while (!converged && iteration < MNA_MAX_ITER) {
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
                    const double abs_tol = MNA_ABSTOL + MNA_RELTOL * fabs(new_voltage);
                    const double damping = (iteration < 5) ? 0.5 : 0.9;
                    comp->last_voltage = damping * new_voltage + (1 - damping) * comp->last_voltage;
                    if (comp->nonlinear_type == NONLINEAR_RESISTOR) {
                        ComponentState state = { .voltage = comp->last_voltage, .current = comp->last_current,
                            .charge = comp->last_charge, .flux = comp->last_flux, .dt = solver->dt };
                        double val1, val2;
                        comp->nonlinear_func(&state, comp->user_data, &val1, &val2);
                        comp->last_conductance = (val2 < MNA_MIN_CONDUCTANCE) ? MNA_MIN_CONDUCTANCE :
                            (val2 > MNA_MAX_CONDUCTANCE) ? MNA_MAX_CONDUCTANCE : val2;
                    }
                    if (voltage_diff > abs_tol) converged = 0;
                } else if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
                    NPoleData* npole = comp->data.npole.npole_data;
                    double* current_values = (double*)malloc((size_t)npole->num_nodes * sizeof(double));
                    if (!current_values) { status = MNA_INSUFFICIENT_MEMORY; break; }
                    double max_diff = 0.0;
                    for (int j = 0; j < npole->num_nodes; j++) {
                        current_values[j] = (npole->nodes[j] > 0) ? solver->x[npole->nodes[j]-1] : 0.0;
                        double diff = fabs(current_values[j] - npole->last_values[j]);
                        if (diff > max_diff) max_diff = diff;
                        double damping = (iteration < 5) ? 0.5 : 0.9;
                        npole->last_values[j] = damping * current_values[j] + (1 - damping) * npole->last_values[j];
                    }
                    free(current_values);
                    if (status != MNA_SUCCESS) break;
                    double abs_tol = MNA_ABSTOL + MNA_RELTOL * max_diff;
                    if (max_diff > abs_tol) converged = 0;
                }
            }
            if (status != MNA_SUCCESS || converged) break;
            mna_reset_system(solver);
            mna_ensure_ground_paths(solver);
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
                        if (comp->source_type == SOURCE_VOLTAGE) mna_stamp_voltage_source(solver, i, source_count++);
                        else mna_stamp_current_source(solver, n1, n2, comp->value);
                        break;
                    case MNA_CUSTOM_NONLINEAR:
                        mna_stamp_custom_nonlinear(solver, i, 1, solver->integration_method);
                        break;
                    case MNA_CUSTOM_NPOLE: {
                        NPoleData* npole = comp->data.npole.npole_data;
                        if (npole) npole->stamp_func(solver, npole->nodes, npole->num_nodes, npole->user_data, 0.0, 0.0);
                        break;
                    }
                    case MNA_SWITCH: {
                        const double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                        mna_stamp_conductance(solver, n1, n2, g);
                        break;
                    }
                    default: break;
                }
            }
            status = mna_solve_linear_system(solver, matrix_size);
            if (status != MNA_SUCCESS) break;
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
        if (comp->type == MNA_SOURCE) comp->value = orig_source_values[i];
    }
    free(orig_source_values);
    return status;
}

MNAStatus mna_solve_ac(MNASolver* solver, double frequency) {
    if (!solver) return MNA_INVALID_HANDLE;
    double omega = TWO_PI * frequency;
    int matrix_size = mna_active_size(solver);
    if (matrix_size > solver->matrix_cap_size) return MNA_INVALID_PARAMETER;
    double complex* A_complex = (double complex*)calloc((size_t)matrix_size * (size_t)matrix_size, sizeof(double complex));
    double complex* b_complex = (double complex*)calloc((size_t)matrix_size, sizeof(double complex));
    double complex* x_complex = (double complex*)calloc((size_t)matrix_size, sizeof(double complex));
    if (!A_complex || !b_complex || !x_complex) {
        free(A_complex); free(b_complex); free(x_complex);
        return MNA_INSUFFICIENT_MEMORY;
    }
    int source_count = 0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;
        double complex admittance = 0.0;
        switch (comp->type) {
            case MNA_RESISTOR:
                admittance = 1.0 / comp->value;
                break;
            case MNA_CAPACITOR:
                admittance = I * omega * comp->value;
                break;
            case MNA_INDUCTOR:
                admittance = (omega == 0) ? 1.0 / MNA_MIN_CONDUCTANCE : 1.0 / (I * omega * comp->value);
                break;
            case MNA_SOURCE: {
                double complex source_val = comp->ac_magnitude * (cos(comp->ac_phase) + I * sin(comp->ac_phase));
                if (comp->source_type == SOURCE_VOLTAGE) {
                    int v_index = solver->max_node_index + source_count;
                    if (n1 > 0) {
                        A_complex[(n1-1) * matrix_size + v_index] = 1.0;
                        A_complex[v_index * matrix_size + (n1-1)] = 1.0;
                    }
                    if (n2 > 0) {
                        A_complex[(n2-1) * matrix_size + v_index] = -1.0;
                        A_complex[v_index * matrix_size + (n2-1)] = -1.0;
                    }
                    b_complex[v_index] = source_val;
                    source_count++;
                } else {
                    if (n1 > 0) b_complex[(n1-1)] -= source_val;
                    if (n2 > 0) b_complex[(n2-1)] += source_val;
                }
                continue;
            }
            case MNA_CUSTOM_NONLINEAR: {
                if (comp->nonlinear_type == NONLINEAR_RESISTOR) admittance = comp->last_conductance;
                else if (comp->nonlinear_type == NONLINEAR_CAPACITOR) {
                    ComponentState state = { .voltage = comp->last_voltage, .current = comp->last_current,
                        .charge = comp->last_charge, .flux = comp->last_flux, .dt = solver->dt };
                    double q0, C0;
                    comp->nonlinear_func(&state, comp->user_data, &q0, &C0);
                    admittance = I * omega * C0;
                } else if (comp->nonlinear_type == NONLINEAR_INDUCTOR) {
                    ComponentState state = { .voltage = comp->last_voltage, .current = comp->last_current,
                        .charge = comp->last_charge, .flux = comp->last_flux, .dt = solver->dt };
                    double phi0, L0;
                    comp->nonlinear_func(&state, comp->user_data, &phi0, &L0);
                    admittance = 1.0 / (I * omega * L0);
                }
                break;
            }
            case MNA_CUSTOM_NPOLE: {
                fprintf(stderr, "Warning: AC analysis for n-pole elements uses DC approximation\n");
                NPoleData* npole = comp->data.npole.npole_data;
                if (npole) npole->stamp_func(solver, npole->nodes, npole->num_nodes, npole->user_data, 0.0, 0.0);
                continue;
            }
            case MNA_SWITCH: {
                admittance = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                break;
            }
            default: continue;
        }
        if (n1 > 0) A_complex[(n1-1) * matrix_size + (n1-1)] += admittance;
        if (n2 > 0) A_complex[(n2-1) * matrix_size + (n2-1)] += admittance;
        if (n1 > 0 && n2 > 0) {
            A_complex[(n1-1) * matrix_size + (n2-1)] -= admittance;
            A_complex[(n2-1) * matrix_size + (n1-1)] -= admittance;
        }
    }
    for (int pivot = 0; pivot < matrix_size; pivot++) {
        int max_row = pivot;
        double max_mag = cabs(A_complex[pivot * matrix_size + pivot]);
        for (int i = pivot + 1; i < matrix_size; i++) {
            double mag = cabs(A_complex[i * matrix_size + pivot]);
            if (mag > max_mag) {
                max_mag = mag;
                max_row = i;
            }
        }
        if (max_mag < MNA_MIN_CONDUCTANCE) {
            free(A_complex); free(b_complex); free(x_complex);
            return MNA_MATRIX_SINGULAR;
        }
        if (max_row != pivot) {
            for (int j = pivot; j < matrix_size; j++) {
                double complex temp = A_complex[pivot * matrix_size + j];
                A_complex[pivot * matrix_size + j] = A_complex[max_row * matrix_size + j];
                A_complex[max_row * matrix_size + j] = temp;
            }
            double complex tempb = b_complex[pivot];
            b_complex[pivot] = b_complex[max_row];
            b_complex[max_row] = tempb;
        }
        for (int i = pivot + 1; i < matrix_size; i++) {
            double complex factor = A_complex[i * matrix_size + pivot] / A_complex[pivot * matrix_size + pivot];
            if (factor == 0.0) continue;
            for (int j = pivot + 1; j < matrix_size; j++) {
                A_complex[i * matrix_size + j] -= factor * A_complex[pivot * matrix_size + j];
            }
            b_complex[i] -= factor * b_complex[pivot];
            A_complex[i * matrix_size + pivot] = 0.0;
        }
    }
    for (int i = matrix_size - 1; i >= 0; i--) {
        x_complex[i] = b_complex[i];
        for (int j = i + 1; j < matrix_size; j++) {
            x_complex[i] -= A_complex[i * matrix_size + j] * x_complex[j];
        }
        x_complex[i] /= A_complex[i * matrix_size + i];
    }
    for (int i = 0; i < matrix_size; i++) solver->ac_solution[i] = x_complex[i];
    free(A_complex);
    free(b_complex);
    free(x_complex);
    return MNA_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/*                             Initial Condition Solve                        */
/* -------------------------------------------------------------------------- */
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
                double g = (comp->type == MNA_SWITCH && !comp->state) ? MNA_MIN_CONDUCTANCE : 1.0 / comp->value;
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
                comp->last_current = (comp->type == MNA_RESISTOR) ? v_diff / comp->value : 0.0;
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
    mna_solve_initial_conditions(solver);
}

MNAStatus mna_solve_transient_step(MNASolver* solver, double dt) {
    if (!solver) return MNA_INVALID_HANDLE;
    solver->dt = dt;
    int matrix_size = mna_active_size(solver);

    /* Determine method: Cooldown forces Euler */
    IntegrationMethod method = solver->integration_method;
    if (solver->euler_cooldown_steps > 0) {
        method = MNA_INTEGRATION_EULER;
    }

    bool fallback_triggered = false;

try_solve:
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
                double G_eq, I_eq;
                if (method == MNA_INTEGRATION_TRAPEZOIDAL) {
                    G_eq = (2.0 * comp->value) / dt;
                    I_eq = (G_eq * comp->last_voltage) + comp->last_current;
                    mna_stamp_conductance(solver, n1, n2, G_eq);
                    mna_stamp_current_source(solver, n1, n2, -I_eq);
                } else {
                    G_eq = comp->value / dt;
                    I_eq = -G_eq * comp->last_voltage;
                    mna_stamp_conductance(solver, n1, n2, G_eq);
                    mna_stamp_current_source(solver, n1, n2, I_eq);
                }
                break;
            }
            case MNA_INDUCTOR: {
                double G_eq, I_eq;
                if (method == MNA_INTEGRATION_TRAPEZOIDAL) {
                    G_eq = dt / (2.0 * comp->value);
                    I_eq = comp->last_current + (G_eq * comp->last_voltage);
                } else {
                    G_eq = dt / comp->value;
                    I_eq = comp->last_current;
                }
                mna_stamp_conductance(solver, n1, n2, G_eq);
                mna_stamp_current_source(solver, n1, n2, I_eq);
                break;
            }
            case MNA_SOURCE: {
                if (comp->source_type == SOURCE_VOLTAGE) mna_stamp_voltage_source(solver, i, source_count++);
                else mna_stamp_current_source(solver, n1, n2, comp->value);
                break;
            }
            case MNA_CUSTOM_NONLINEAR: {
                mna_stamp_custom_nonlinear(solver, i, 0, method);
                break;
            }
            case MNA_CUSTOM_NPOLE: {
                NPoleData* npole = comp->data.npole.npole_data;
                if (npole) npole->stamp_func(solver, npole->nodes, npole->num_nodes, npole->user_data, solver->time, dt);
                break;
            }
            case MNA_SWITCH: {
                double g = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                mna_stamp_conductance(solver, n1, n2, g);
                break;
            }
            default: break;
        }
    }
    MNAStatus status = mna_solve_linear_system(solver, matrix_size);
    if (status != MNA_SUCCESS) {
        if (method == MNA_INTEGRATION_TRAPEZOIDAL && !fallback_triggered && !solver->bypass_oscillation_check) {
            fprintf(stderr, "Warning: Linear solve failed. Retrying with Backward Euler.\n");
            method = MNA_INTEGRATION_EULER;
            fallback_triggered = true;
            goto try_solve;
        }
        return status;
    }

    /* Oscillation Detection (Only if Trapezoidal was attempted) */
    if (method == MNA_INTEGRATION_TRAPEZOIDAL && !fallback_triggered && !solver->bypass_oscillation_check && solver->time > dt) {
        bool oscillation = false;
        for (int i = 0; i < solver->num_components; i++) {
            Component* comp = &solver->components[i];
            if (comp->type == MNA_CAPACITOR || comp->type == MNA_INDUCTOR || comp->type == MNA_CUSTOM_NONLINEAR) {
                int n1 = comp->node1;
                int n2 = comp->node2;
                double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
                double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
                double v_new = v1 - v2;
                double delta_new = v_new - comp->last_voltage;
                double delta_old = comp->last_voltage - comp->prev_voltage;

                /* Relative tolerance check */
                double threshold = MNA_OSCILLATION_TOL_REL * (fabs(v_new) + 1.0);

                if (fabs(delta_new) > threshold && fabs(delta_old) > threshold) {
                    if (delta_new * delta_old < 0.0) {
                        oscillation = true;
                        break;
                    }
                }
            }
        }
        if (oscillation) {
            fprintf(stderr, "Warning: Trapezoidal oscillation detected. Retrying step with Backward Euler.\n");
            method = MNA_INTEGRATION_EULER;
            solver->euler_cooldown_steps = MNA_EULER_COOLDOWN_STEPS; /* Set persistent cooldown */
            fallback_triggered = true;
            goto try_solve;
        }
    }

    /* Update component states consistent with the method used */
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;
        double v1 = (n1 > 0) ? solver->x[n1-1] : 0.0;
        double v2 = (n2 > 0) ? solver->x[n2-1] : 0.0;
        double v = v1 - v2;

        comp->prev_voltage = comp->last_voltage;
        comp->prev_current = comp->last_current;

        switch (comp->type) {
            case MNA_CAPACITOR:
                comp->last_voltage = v;
                if (method == MNA_INTEGRATION_TRAPEZOIDAL) {
                    comp->last_current = ((2.0 * comp->value) / dt) * (v - comp->prev_voltage) - comp->last_current;
                } else {
                    comp->last_current = (comp->value / dt) * (v - comp->prev_voltage);
                }
                break;
            case MNA_INDUCTOR: {
                comp->last_voltage = v;
                if (method == MNA_INTEGRATION_TRAPEZOIDAL) {
                    comp->last_current += (dt / (2.0 * comp->value)) * (v + comp->prev_voltage);
                } else {
                    comp->last_current += (dt / comp->value) * v;
                }
                break;
            }
            case MNA_CUSTOM_NONLINEAR:
                switch (comp->nonlinear_type) {
                    case NONLINEAR_RESISTOR:
                        comp->last_voltage = v;
                        break;
                    case NONLINEAR_CAPACITOR: {
                        ComponentState state = { .voltage = v, .current = comp->last_current,
                            .charge = comp->last_charge, .flux = comp->last_flux, .dt = dt };
                        double q, C;
                        comp->nonlinear_func(&state, comp->user_data, &q, &C);
                        comp->last_current = (comp->trans_G_eq * (v - comp->prev_voltage)) - comp->last_current;
                        comp->last_voltage = v;
                        comp->last_charge = q;
                        break;
                    }
                    case NONLINEAR_INDUCTOR: {
                        double i_new = comp->trans_I_eq + comp->trans_G_eq * v;
                        ComponentState state = { .voltage = v, .current = i_new,
                            .charge = comp->last_charge, .flux = comp->last_flux, .dt = dt };
                        double phi, L;
                        comp->nonlinear_func(&state, comp->user_data, &phi, &L);
                        comp->last_current = i_new;
                        comp->last_flux = phi;
                        break;
                    }
                }
                break;
            case MNA_CUSTOM_NPOLE:
                if (comp->data.npole.npole_data) {
                    NPoleData* npole = comp->data.npole.npole_data;
                    for (int j = 0; j < npole->num_nodes; j++) {
                        npole->last_values[j] = (npole->nodes[j] > 0) ? solver->x[npole->nodes[j]-1] : 0.0;
                    }
                }
                break;
            default: break;
        }
    }

    /* Decrement cooldown for next step */
    if (solver->euler_cooldown_steps > 0) {
        solver->euler_cooldown_steps--;
    }

    solver->time += dt;
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
        case MNA_RESISTOR: return v / comp->value;
        case MNA_CAPACITOR:
            if (!solver->transient_initialized || solver->dt == 0.0) return 0.0;
            return (v - comp->prev_voltage) * (comp->value / solver->dt);
        case MNA_INDUCTOR: return comp->last_current;
        case MNA_SOURCE:
            if (comp->source_type == SOURCE_VOLTAGE) {
                int vs_index = 0;
                for (int i = 0; i < handle; i++) {
                    if (solver->components[i].type == MNA_SOURCE && solver->components[i].source_type == SOURCE_VOLTAGE) vs_index++;
                }
                return solver->x[solver->max_node_index + vs_index];
            } else return comp->value;
        case MNA_SWITCH: return v / (comp->state ? comp->value : 1.0/MNA_MIN_CONDUCTANCE);
        case MNA_CUSTOM_NONLINEAR: {
            ComponentState state = { .voltage = v, .current = comp->last_current,
                .charge = comp->last_charge, .flux = comp->last_flux, .dt = solver->dt };
            double current, conductance;
            comp->nonlinear_func(&state, comp->user_data, &current, &conductance);
            return current;
        }
        case MNA_CUSTOM_NPOLE: return 0.0;
        default: return 0.0;
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

#endif // MNA_SOLVER_V2_5_H
