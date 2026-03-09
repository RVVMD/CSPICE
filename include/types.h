#ifndef MNA_TYPES_H
#define MNA_TYPES_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdbool.h>

/* ============================================================================
 * Constants
 * ============================================================================ */
#define MNA_MAX_ITER           50
#define MNA_RELTOL             1e-6
#define MNA_ABSTOL             1e-9
#define MNA_VT                 0.02585
#define MNA_MIN_CONDUCTANCE    1e-12
#define MNA_MAX_CONDUCTANCE    1e12
#define MNA_GROUND_CONDUCTANCE 1e-9
#define TWO_PI                 6.283185307179586476925286766559
#define MNA_INIT_CAPACITY      4
#define MNA_TRBDF2_GAMMA       0.5857864376269049511983112757903
#define MNA_DAMPING_THRESHOLD  1
#define MNA_DAMPING_STEEPNESS  1e3
#define MNA_V_FLOOR            1e-9
#define MNA_I_FLOOR            1e-9

/* ============================================================================
 * Status Codes
 * ============================================================================ */
typedef enum {
    MNA_SUCCESS,
    MNA_MATRIX_SINGULAR,
    MNA_CONVERGENCE_FAILURE,
    MNA_INVALID_HANDLE,
    MNA_INVALID_NODE,
    MNA_INSUFFICIENT_MEMORY,
    MNA_INVALID_PARAMETER
} MNAStatus;

/* ============================================================================
 * Component Types
 * ============================================================================ */
typedef enum {
    MNA_RESISTOR,
    MNA_CAPACITOR,
    MNA_INDUCTOR,
    MNA_SOURCE,
    MNA_SWITCH,
    MNA_CUSTOM_NONLINEAR,
    MNA_CUSTOM_NPOLE
} ComponentType;

/* ============================================================================
 * Nonlinear Element Types
 * ============================================================================ */
typedef enum {
    NONLINEAR_RESISTOR,
    NONLINEAR_CAPACITOR,
    NONLINEAR_INDUCTOR
} NonlinearType;

/* ============================================================================
 * Source Types
 * ============================================================================ */
typedef enum {
    SOURCE_VOLTAGE,
    SOURCE_CURRENT
} SourceType;

/* ============================================================================
 * Integration Methods
 * ============================================================================ */
typedef enum {
    MNA_INTEGRATION_TRBDF2
} IntegrationMethod;

/* ============================================================================
 * Component State
 * ============================================================================ */
typedef struct {
    double voltage;
    double current;
    double charge;
    double flux;
    double dt;
} ComponentState;

/* ============================================================================
 * Forward Declarations
 * ============================================================================ */
struct MNASolver;
typedef struct MNASolver MNASolver;
typedef int ComponentHandle;

/* ============================================================================
 * Function Pointer Types
 * ============================================================================ */
typedef void (*CustomNonlinearFunc)(const ComponentState* state, void* user_data,
                                    double* value1, double* value2);

typedef void (*NPoleStampFunc)(struct MNASolver* solver,
                               const int* nodes,
                               int num_nodes,
                               void* user_data,
                               double time,
                               double dt);

/* ============================================================================
 * N-Pole Data Structure
 * ============================================================================ */
typedef struct {
    int* nodes;
    int num_nodes;
    NPoleStampFunc stamp_func;
    void* user_data;
    double* last_values;
} NPoleData;

/* ============================================================================
 * Component Structure
 * ============================================================================ */
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
    double prev_voltage;
    double prev_current;
    double stage2_voltage;
    double stage2_current;
    double stage1_voltage;
    double stage1_current;
    double last_conductance;
    double last_charge;
    double last_flux;
    double smoothed_alpha;
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

/* ============================================================================
 * MNA Solver Structure
 * ============================================================================ */
struct MNASolver {
    int num_nodes;
    int num_components;
    int num_sources;
    int num_nonlinear;
    int max_node_index;
    int transient_initialized;
    int preserve_dc_state;
    double time;
    double dt;
    IntegrationMethod integration_method;
    Component* components;
    double* A;
    double* b;
    double* x;
    double complex* ac_solution;
    int cap_components;
    int matrix_cap_size;
};

#endif /* MNA_TYPES_H */
