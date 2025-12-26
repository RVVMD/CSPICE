#ifndef MNA_SOLVER_V2_1_H
#define MNA_SOLVER_V2_1_H

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

// Default initial capacities
#define MNA_INIT_NODES_CAP     16
#define MNA_INIT_SOURCES_CAP   8
#define MNA_INIT_COMPONENTS_CAP 64

#define MNA_MAX_NONLINEAR 20
#define MNA_MAX_ITER 50
#define MNA_RELTOL 1e-6
#define MNA_ABSTOL 1e-9
#define MNA_VT 0.02585
#define MNA_MIN_CONDUCTANCE 1e-12
#define MNA_MAX_CONDUCTANCE 1e12
#define TWO_PI 6.28318530717958647692

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
    double voltage;
    double current;
    double charge;
    double flux;
    double dt;
} ComponentState;

// Generalized nonlinear function interface
typedef void (*CustomNonlinearFunc)(const ComponentState* state, void* user_data,
                                  double* value1, double* value2);

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
    double last_voltage;
    double last_current;
    double last_conductance;
    double last_charge;
    double last_flux;
    CustomNonlinearFunc nonlinear_func;
    void* user_data;

    double trans_G_eq;
    double trans_I_eq;
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

    Component* components;
    double* A;
    double* b;
    double* x;
    double complex* ac_solution;

    int cap_nodes;
    int cap_sources;
    int cap_components;
    int matrix_cap_size;
} MNASolver;

// Matrix access macro
#define MAT(solver, i, j) ((solver)->A[(i) * (solver)->matrix_cap_size + (j)])

// Internal helpers that should remain inline
static inline int mna_active_size(const MNASolver* s) {
    return s->max_node_index + s->num_sources;
}

// ---- Public API ----
MNAStatus mna_init_sized(MNASolver* solver, int max_nodes, int max_sources, int max_components);
MNAStatus mna_init(MNASolver* solver);
void mna_destroy(MNASolver* solver);

int mna_create_node(MNASolver* solver);
MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2);

MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                          double value, SourceType src_type, ComponentHandle* handle);

MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2, double value, ComponentHandle* handle);
MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                                 NonlinearType nl_type,
                                 CustomNonlinearFunc func, void* user_data,
                                 double initial_value1, double initial_value2,
                                 ComponentHandle* handle);

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

#endif // MNA_SOLVER_V2_1_H
