# CSPICE Solver Module

Implementation of Modified Nodal Analysis (MNA) with TR-BDF2 integration for circuit simulation.

## Architecture

```
solver/
├── core.c/h        # Solver lifecycle, matrix setup, node management
├── dc.c/h          # DC operating point analysis
├── ac.c/h          # AC small-signal analysis (frequency domain)
└── transient.c/h   # Transient analysis with TR-BDF2
```

## Modified Nodal Analysis

MNA extends classical nodal analysis to handle voltage sources and inductors without modification. The system is formulated as:

```
┌         ┐ ┌     ┐   ┌     ┐
│ G   B   │ │ v   │   │ i   │
│         │ │     │ = │     │
│ C   D   │ │ ix  │   │ e   │
└         ┘ └     ┘   └     ┘
```

Where:
- **G** - Conductance matrix (node equations)
- **B** - Voltage source incidence
- **C** - Current source incidence  
- **D** - Zero matrix (for independent sources)
- **v** - Node voltages (unknown)
- **ix** - Branch currents through voltage sources (unknown)
- **i** - Current source contributions
- **e** - Voltage source values

## Integration Method: TR-BDF2

The solver uses the **TR-BDF2** (Trapezoidal Rule with Backward Differentiation Formula) method, which combines:

1. **Trapezoidal Rule** (first half-step): A-stable, second-order accurate
2. **BDF2** (second half-step): L-stable, handles stiffness

This provides excellent stability for stiff power system transients while maintaining accuracy.

### Gamma Parameter

```
γ = 2 - √2 ≈ 0.5857864376269049511983112757903
```

### Two-Stage Integration

**Stage 1 (TR)** - from tₙ to tₙ₊ᵧ:
```
vₙ₊ᵧ = vₙ + (γ·dt/2) · [f(vₙ) + f(vₙ₊ᵧ)]
```

**Stage 2 (BDF2)** - from tₙ₊ᵧ to tₙ₊₁:
```
vₙ₊₁ = (1/γ(2-γ)) · vₙ₊ᵧ - ((1-γ)²/γ(2-γ)) · vₙ
```

## Analysis Types

### DC Analysis (`dc.c`)

Computes the operating point by solving the nonlinear DC equations:

```c
MNAStatus mna_solve_dc(MNASolver* solver);
```

**Algorithm:**
1. Stamp component contributions into matrix A and vector b
2. Solve linear system using Gaussian elimination with partial pivoting
3. For nonlinear circuits: Newton-Raphson iteration with adaptive damping

### AC Analysis (`ac.c`)

Small-signal frequency domain analysis:

```c
MNAStatus mna_solve_ac(MNASolver* solver, double frequency);
```

**Process:**
1. Linearize nonlinear elements around DC operating point
2. Convert to complex phasor domain
3. Solve complex linear system at specified frequency

### Transient Analysis (`transient.c`)

Time-domain simulation with energy storage elements:

```c
void mna_init_transient(MNASolver* solver);
MNAStatus mna_solve_transient_step(MNASolver* solver, double dt);
```

**Initialization:**
- Compute initial conditions for capacitors (voltage) and inductors (current)
- Optionally use DC operating point if `preserve_dc_state` is enabled

**Time Step:**
1. Stamp companion models for C and L using TR-BDF2 discretization
2. Handle nonlinear elements with Newton-Raphson iteration
3. Apply adaptive damping if convergence is slow

## Companion Models

Energy storage elements are replaced with equivalent conductance/current source models:

### Capacitor (TR-BDF2)
```
i_C = C · dv/dt  →  i_C = G_eq · v + I_eq

where:
G_eq = C / (γ · dt)
I_eq = history term from previous time step
```

### Inductor (TR-BDF2)
```
v_L = L · di/dt  →  v_L = R_eq · i + V_eq

where:
R_eq = L / (γ · dt)
V_eq = history term from previous time step
```

## Nonlinear Handling

Nonlinear elements (diodes, custom devices) use Newton-Raphson iteration:

```c
// Linearization at operating point (v_k, i_k)
i = i_k + g_k · (v - v_k)

where g_k = di/dv|_(v_k,i_k) is the small-signal conductance
```

**Adaptive Damping:**
When Newton-Raphson diverges, the update is scaled:
```
v_new = v_old + α · Δv

where α ∈ (0, 1] is reduced if residuals increase
```

## Matrix Operations

The solver uses dense matrix storage with Gaussian elimination:

```c
// Matrix access macro
#define MAT(solver, i, j) ((solver)->A[(i) * matrix_cap_size + (j)])

// LU decomposition with partial pivoting
// Forward/backward substitution
```

For large sparse systems (future optimization), consider:
- Sparse matrix formats (CSR/CSC)
- Iterative solvers (GMRES, BiCGSTAB)
- Preconditioners (ILU, Jacobi)

## API Reference

### Solver Lifecycle

```c
MNAStatus mna_init(MNASolver* solver);
void mna_destroy(MNASolver* solver);
```

### Configuration

```c
MNAStatus mna_set_integration_method(MNASolver* solver, IntegrationMethod method);
MNAStatus mna_set_preserve_dc_state(MNASolver* solver, int enable);
```

### Analysis Functions

```c
MNAStatus mna_solve_dc(MNASolver* solver);
MNAStatus mna_solve_ac(MNASolver* solver, double frequency);
void mna_init_transient(MNASolver* solver);
MNAStatus mna_solve_transient_step(MNASolver* solver, double dt);
```

### Solution Access

```c
double mna_get_node_voltage(MNASolver* solver, int node);
double complex mna_get_ac_node_voltage(MNASolver* solver, int node);
double mna_get_component_current(MNASolver* solver, ComponentHandle handle);
double mna_get_component_voltage(MNASolver* solver, ComponentHandle handle);
```

## Error Handling

Status codes defined in `types.h`:

| Code | Description |
|------|-------------|
| `MNA_SUCCESS` | Operation completed successfully |
| `MNA_MATRIX_SINGULAR` | System matrix is singular (check circuit topology) |
| `MNA_CONVERGENCE_FAILURE` | Newton-Raphson did not converge |
| `MNA_INVALID_HANDLE` | Invalid component handle |
| `MNA_INVALID_NODE` | Invalid node index |
| `MNA_INSUFFICIENT_MEMORY` | Memory allocation failed |

## Performance Considerations

1. **Matrix Size**: O(n²) storage, O(n³) solve time for dense matrices
2. **Time Steps**: Smaller dt increases accuracy but requires more steps
3. **Nonlinear Circuits**: Each Newton-Raphson iteration requires full matrix refactorization
4. **Memory**: Pre-allocate with estimated circuit size to avoid reallocations

## Future Enhancements

- [ ] Adaptive time-stepping based on truncation error
- [ ] Sparse matrix support for large circuits
- [ ] Parallel matrix operations
- [ ] Higher-order integration methods (BDF3-5)
- [ ] Frequency sweep for AC analysis

## References

1. C.W. Ho et al., "The Modified Nodal Approach to Network Analysis," IEEE TCS, 1975
2. J. Vlach and K. Singhal, "Computer Methods for Circuit Analysis and Design"
3. A. Protheroe et al., "TR-BDF2: A New Approach to Transient Integration"
