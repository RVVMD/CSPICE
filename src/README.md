# CSPICE Source Directory

Implementation source files for the MNA solver library.

## Structure

```
src/
├── matrix.c          # Matrix algebra and linear system solver
├── solver/           # Analysis engines
│   ├── core.c        # Solver lifecycle and utilities
│   ├── dc.c          # DC operating point analysis
│   ├── ac.c          # AC small-signal analysis
│   └── transient.c   # Time-domain transient analysis
└── cli/              # Command-line interface
    ├── netlist.c     # Netlist parser
    └── cspice_cli.c  # CLI main entry point
```

## Modules

### matrix.c

**Purpose:** Linear algebra operations for MNA system solution.

**Key Functions:**
- `mna_lu_solve()` - LU decomposition with partial pivoting
- `mna_resize_matrix()` - Dynamic matrix allocation
- `mna_reset_system()` - Clear matrix and vectors

**Algorithm:** Gaussian elimination with partial pivoting for numerical stability.

**See:** [include/README.md](../include/README.md) for matrix access macros.

---

### solver/

**Purpose:** Analysis engine implementations.

| Module | Function | Description |
|--------|----------|-------------|
| `core.c` | `mna_init()`, `mna_destroy()` | Solver lifecycle |
| `dc.c` | `mna_solve_dc()` | DC operating point |
| `ac.c` | `mna_solve_ac()` | Frequency domain analysis |
| `transient.c` | `mna_solve_transient_step()` | Time-domain simulation |

**See:** [solver/README.md](solver/README.md) for detailed documentation.

---

### cli/

**Purpose:** Command-line interface for netlist-based simulation.

**Components:**
- `netlist.c/h` - SPICE-like netlist parser
- `cspice_cli.c` - Main entry point with argument parsing

**Supported Netlist Syntax:**
```cir
* Comment
R1 n1 n2 1k        ; Resistor
V1 n1 0 DC 10      ; Voltage source
.analysis tran 1u 5m
.print v(n1)
.write output.csv
.end
```

**See:** [cli/README.md](cli/README.md) for user guide.

---

## Build Integration

Sources are compiled into the static and shared libraries:

```cmake
set(MNA_SOURCES
    src/matrix.c
    src/solver/core.c
    src/solver/dc.c
    src/solver/ac.c
    src/solver/transient.c
    elements/passive.c
    elements/sources.c
    elements/nonlinear/nonlinear.c
)
```

CLI is built as a separate executable:

```cmake
add_executable(cspice
    src/cli/netlist.c
    src/cli/cspice_cli.c
)
```

---

## Coding Conventions

### Documentation

All public functions use Doxygen-style comments:

```c
/**
 * @brief Solve DC operating point
 * @param solver The MNA solver
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_solve_dc(MNASolver* solver);
```

### Naming

- **Public API**: `mna_*` prefix (e.g., `mna_init`, `mna_solve_dc`)
- **Internal functions**: No prefix, file-local scope (e.g., `stamp_resistor`, `parse_value`)
- **Types**: `MNA*` prefix for enums/structs (e.g., `MNAStatus`, `MNASolver`)

### Error Handling

Return `MNAStatus` for all operations that can fail:

```c
MNAStatus status = mna_solve_dc(&solver);
if (status != MNA_SUCCESS) {
    // Handle error appropriately
}
```

### Memory Management

- Use `malloc`/`free` for dynamic allocation
- Always pair `mna_init()` with `mna_destroy()`
- Document ownership (caller vs. callee allocates)

---

## Testing

### Unit Tests (Planned)

Test files will be placed in `tests/`:

```
tests/
├── test_matrix.c      # Matrix operations
├── test_dc.c          # DC analysis
├── test_transient.c   # Transient analysis
└── test_netlist.c     # Netlist parser
```

### Regression Tests

Compare results against reference simulators:

```bash
# Run regression suite
ctest --test-dir build
```

---

## Debugging

### Build Types

```bash
# Debug build with symbols
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Common Issues

**Matrix Singular:**
- Check for floating nodes (no DC path to ground)
- Check for voltage source loops
- Add high-value resistors to floating nodes

**Convergence Failure:**
- Reduce time step for transient
- Check for discontinuous elements
- Enable damping in Newton-Raphson

**Memory Errors:**
- Check array bounds in stamping functions
- Verify `mna_resize_*()` calls before access
- Use valgrind for leak detection

---

## Performance Profiling

```bash
# Build with profiling
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DPROFILE=ON ..

# Run with gprof
./bin/cspice circuit.cir
gprof ./bin/cspice gmon.out > profile.txt
```

**Hotspots (typical):**
1. Matrix factorization (`mna_lu_solve`) - O(n³)
2. Newton-Raphson iteration for nonlinear circuits
3. Component stamping in large circuits

**Optimization Strategies:**
- Sparse matrix storage for large circuits
- Reuse LU factors when possible
- Parallel stamping for independent components

---

## Related Documentation

- [Main README](../README.md) - Project overview
- [Solver Module](solver/README.md) - MNA implementation
- [Elements](../elements/README.md) - Component models
- [Public API](../include/README.md) - Header documentation
- [CLI Guide](cli/README.md) - Netlist format
