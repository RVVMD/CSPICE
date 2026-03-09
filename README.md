# CSPICE - Circuit Simulation for Relay Protection and Automation

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]()
[![C Standard](https://img.shields.io/badge/C-C11-blue)]()

A lightweight, embeddable circuit simulation engine based on Modified Nodal Analysis (MNA) with TR-BDF2 integration, specifically designed for **Relay Protection and Automation (RZA)** applications.

## Features

- **MNA Solver** with TR-BDF2 integration method for stiff power systems
- **Nonlinear elements** with adaptive Newton-Raphson damping
- **Analysis types**: DC operating point, AC small-signal, Transient
- **Components**: R, L, C, V, I sources, switches, and power system elements
- **Export**: CSV and COMTRADE (planned) for external analysis
- **Pure C** (C11 standard) - no external dependencies except libm

## Quick Start

### Build

```bash
mkdir build && cd build
cmake ..
make
```

### Run Simulation

```bash
# Using the CLI
./bin/cspice ../examples/rc_circuit.cir

# Or use the library programmatically (see examples/)
```

### Example Netlist

```cir
* RC Circuit Transient Response
V1 1 0 DC 10
R1 1 2 1k
C1 2 0 1u

.analysis tran 10u 5m
.print v(2)
.write output.csv
.end
```

## Project Structure

```
CSPICE/
├── include/          # Public API headers
│   ├── types.h       # Core types and constants
│   └── matrix.h      # Matrix operations
├── src/              # Implementation
│   ├── solver/       # Analysis engines (DC, AC, Transient)
│   ├── cli/          # Command-line interface
│   └── matrix.c      # Matrix algebra
├── elements/         # Circuit components
│   ├── passive.c     # R, L, C, Switch
│   ├── sources.c     # Voltage/Current sources
│   └── nonlinear/    # Nonlinear elements
├── examples/         # Example netlists
├── backup/           # Reference implementations
└── tests/            # Unit tests (planned)
```

## Documentation

- **[CLI User Guide](src/cli/README.md)** - Netlist format and commands
- **[Library API](mna_solver.h)** - Programmatic usage
- **[Solver Architecture](src/solver/README.md)** - MNA implementation details
- **[Components](elements/README.md)** - Element models and equations

## Usage Modes

### 1. Command-Line Interface

Run simulations directly from netlist files:

```bash
./bin/cspice circuit.cir          # Run simulation
./bin/cspice -h                   # Show help
./bin/cspice -q circuit.cir       # Quiet mode
```

### 2. Library Integration

Link against `libmna_solver.a` or `libmna_solver.so`:

```c
#include "mna_solver.h"

MNASolver solver;
mna_init(&solver);

// Build circuit programmatically
int node1 = mna_create_node(&solver);
ComponentHandle r1;
mna_add_resistor(&solver, node1, 0, 1000.0, &r1);

// Run analysis
mna_solve_dc(&solver);
double v = mna_get_node_voltage(&solver, node1);

mna_destroy(&solver);
```

## Analysis Types

| Type | Command | Description |
|------|---------|-------------|
| DC | `.analysis dc` | Operating point calculation |
| AC | `.analysis ac <freq>` | Small-signal frequency response |
| Transient | `.analysis tran <dt> <tend>` | Time-domain simulation |

## Component Support

| Category | Components |
|----------|------------|
| Passive | Resistor, Capacitor, Inductor, Switch |
| Sources | DC Voltage, DC Current, AC Source |
| Nonlinear | Diode, Custom n-pole elements |
| **Planned (Phase 2)** | Transformers, CT/VT, Fault elements, Transmission lines |

## Roadmap

See [roadmap.txt](roadmap.txt) for the complete development plan.

### Phase 1: Foundation ✅
- [x] Modular code structure
- [x] CMake build system
- [x] Public C API
- [x] CLI netlist parser
- [x] Basic error logging

### Phase 2: Core RZA Components (Next)
- [ ] Mutual inductance and transformers
- [ ] Current/Voltage transformers with saturation
- [ ] Fault elements (1φ, 2φ, 3φ, ground)
- [ ] Time-controlled switches

### Phase 3: Interoperability
- [ ] COMTRADE export (IEEE C37.111)
- [ ] Standard netlist format enhancements

## Validation

Target accuracy: <1% deviation from reference simulators (LTspice, MATLAB/Simulink) for linear circuits.

## License

MIT License - See LICENSE file for details.

## Acknowledgments

- Modified Nodal Analysis: C.W. Ho et al., "The Modified Nodal Approach to Network Analysis"
- TR-BDF2 Method: A. Protheroe et al., "A new approach to transient integration"
- Target domain: Relay Protection and Automation (RZA) engineering
