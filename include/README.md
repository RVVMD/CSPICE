# CSPICE Public API Headers

Core type definitions and public interfaces for the MNA solver library.

## Files

| File | Description |
|------|-------------|
| `types.h` | Core types, constants, and data structures |
| `matrix.h` | Matrix operations and access macros |

---

## types.h

### Status Codes

```c
typedef enum {
    MNA_SUCCESS,              // Operation completed successfully
    MNA_MATRIX_SINGULAR,      // System matrix is singular
    MNA_CONVERGENCE_FAILURE,  // Newton-Raphson did not converge
    MNA_INVALID_HANDLE,       // Invalid component handle
    MNA_INVALID_NODE,         // Invalid node index
    MNA_INSUFFICIENT_MEMORY,  // Memory allocation failed
    MNA_INVALID_PARAMETER     // Invalid parameter value
} MNAStatus;
```

### Component Types

```c
typedef enum {
    MNA_RESISTOR,
    MNA_CAPACITOR,
    MNA_INDUCTOR,
    MNA_SOURCE,
    MNA_SWITCH,
    MNA_CUSTOM_NONLINEAR,
    MNA_CUSTOM_NPOLE
} ComponentType;
```

### Constants

```c
#define MNA_MAX_ITER           50       // Max Newton-Raphson iterations
#define MNA_RELTOL             1e-6     // Relative tolerance
#define MNA_ABSTOL             1e-9     // Absolute tolerance
#define MNA_VT                 0.02585  // Thermal voltage (V) at 300K
#define MNA_TRBDF2_GAMMA       0.585786 // TR-BDF2 parameter (2-√2)
#define MNA_MIN_CONDUCTANCE    1e-12    // Minimum conductance (S)
#define MNA_MAX_CONDUCTANCE    1e12     // Maximum conductance (S)
```

### Core Structures

**Component State:**
```c
typedef struct {
    double voltage;    // Voltage across element
    double current;    // Current through element
    double charge;     // Stored charge (capacitors)
    double flux;       // Magnetic flux (inductors)
    double dt;         // Time step
} ComponentState;
```

**Solver Handle:**
```c
typedef struct MNASolver {
    int num_nodes;           // Number of nodes
    int num_components;      // Number of components
    int num_sources;         // Number of sources
    int num_nonlinear;       // Number of nonlinear elements
    double time;             // Current simulation time
    double dt;               // Current time step
    Component* components;   // Component array
    double* A;               // System matrix
    double* b;               // RHS vector
    double* x;               // Solution vector
} MNASolver;
```

**Component Handle:**
```c
typedef int ComponentHandle;  // Opaque handle to component
```

### Function Pointer Types

**Custom Nonlinear Element:**
```c
typedef void (*CustomNonlinearFunc)(
    const ComponentState* state,  // Current state
    void* user_data,              // User context
    double* value1,               // Output: current
    double* value2                // Output: conductance
);
```

**N-Pole Stamping:**
```c
typedef void (*NPoleStampFunc)(
    MNASolver* solver,
    const int* nodes,
    int num_nodes,
    void* user_data,
    double time,
    double dt
);
```

---

## matrix.h

### Matrix Access

```c
#define MAT(solver, i, j) ((solver)->A[(i) * (solver)->matrix_cap_size + (j)])
```

Access element A[i,j] in the system matrix.

### Memory Management

```c
bool mna_resize_components(MNASolver* solver, int required);
bool mna_resize_matrix(MNASolver* solver, int req_size);
```

Dynamically resize internal arrays. Called automatically when needed.

### Matrix Operations

```c
void mna_reset_system(MNASolver* solver);
```

Clear matrix A and vectors b, x for new stamping.

```c
MNAStatus mna_lu_solve(MNASolver* solver, int size);
```

Solve linear system using LU decomposition with partial pivoting.

---

## Usage Examples

### Basic Circuit Setup

```c
#include "types.h"
#include "matrix.h"

MNASolver solver;
mna_init(&solver);

// Create nodes
int node1 = mna_create_node(&solver);
int node2 = mna_create_node(&solver);

// Add components
ComponentHandle r1, c1;
mna_add_resistor(&solver, node1, node2, 1000.0, &r1);
mna_add_capacitor(&solver, node2, 0, 1e-6, &c1);

// Solve
mna_solve_dc(&solver);
double v = mna_get_node_voltage(&solver, node2);

mna_destroy(&solver);
```

### Custom Nonlinear Element

```c
void varistor(const ComponentState* state, void* user_data,
              double* current, double* conductance) {
    double v = state->voltage;
    double alpha = 30.0;  // Nonlinearity coefficient
    double vn = 100.0;    // Reference voltage
    
    *current = pow(fabs(v/vn), alpha) * (v > 0 ? 1 : -1);
    *conductance = (alpha/vn) * pow(fabs(v/vn), alpha-1);
}

// Usage
ComponentHandle mov;
mna_add_custom_nonlinear(&solver, n1, 0, varistor, NULL, &mov);
```

### N-Pole Element

```c
typedef struct {
    double L1, L2, M;  // Self and mutual inductances
} CoupledCoils;

void stamp_coupled_coils(MNASolver* solver, const int* nodes,
                         int num_nodes, void* user_data,
                         double time, double dt) {
    CoupledCoils* coils = (CoupledCoils*)user_data;
    
    // Stamp mutual inductance MNA equations
    // ...
}

// Usage
CoupledCoils params = {1.0, 1.0, 0.9};  // L1=L2=1H, k=0.9
int nodes[4] = {n1, n2, n3, n4};

mna_add_n_pole(&solver, nodes, 4, stamp_coupled_coils, &params, &handle);
```

---

## Thread Safety

The solver is **not thread-safe**. Each thread must use its own `MNASolver` instance:

```c
// OK: Separate solvers per thread
MNASolver solver1, solver2;
mna_init(&solver1);
mna_init(&solver2);

// NOT OK: Sharing solver between threads
// MNASolver shared_solver;
// mna_init(&shared_solver);  // Race conditions!
```

---

## Memory Management

### Initialization

Always initialize before use:

```c
MNASolver solver;
MNAStatus status = mna_init(&solver);
if (status != MNA_SUCCESS) {
    // Handle initialization failure
}
```

### Cleanup

Always destroy when done:

```c
mna_destroy(&solver);
```

This frees all allocated memory (components, matrix, vectors).

### State Preservation

For transient analysis following DC:

```c
// Compute DC operating point
mna_solve_dc(&solver);

// Preserve state for transient initialization
mna_set_preserve_dc_state(&solver, 1);

// Initialize and run transient
mna_init_transient(&solver);
```

---

## Error Handling Pattern

```c
MNAStatus status;

status = mna_add_resistor(&solver, n1, n2, value, &handle);
if (status == MNA_INVALID_NODE) {
    fprintf(stderr, "Invalid node index\n");
} else if (status == MNA_INSUFFICIENT_MEMORY) {
    fprintf(stderr, "Out of memory\n");
}

status = mna_solve_dc(&solver);
if (status == MNA_MATRIX_SINGULAR) {
    fprintf(stderr, "Circuit has floating nodes or short circuit\n");
} else if (status == MNA_CONVERGENCE_FAILURE) {
    fprintf(stderr, "Nonlinear iteration did not converge\n");
}
```

---

## Version Information

```c
#define MNA_SOLVER_VERSION_MAJOR 2
#define MNA_SOLVER_VERSION_MINOR 5
#define MNA_SOLVER_VERSION "2.5.0"
```

Check at runtime:

```c
printf("Using MNA Solver version %s\n", MNA_SOLVER_VERSION);
```

---

## Related Documentation

- [Solver Module](../src/solver/README.md) - Analysis implementation
- [Elements](../elements/README.md) - Component models
- [CLI](../src/cli/README.md) - Netlist format
