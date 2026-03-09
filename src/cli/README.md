# CSPICE CLI User Guide

## Overview

The CSPICE command-line interface allows you to run circuit simulations using SPICE-like netlist files.

## Building

```bash
mkdir build && cd build
cmake ..
make
```

The `cspice` executable will be in `build/bin/`.

## Usage

```bash
cspice <netlist_file> [options]

Options:
  -h, --help     Show help message
  -v, --version  Show version information
  -q, --quiet    Suppress progress output
```

## Netlist Format

### Components

| Component | Prefix | Example | Description |
|-----------|--------|---------|-------------|
| Resistor | R | `R1 n1 n2 1k` | Resistor between nodes n1 and n2 |
| Capacitor | C | `C1 n1 0 1u` | Capacitor (node 0 is ground) |
| Inductor | L | `L1 n1 n2 1m` | Inductor |
| Voltage Source | V | `V1 n1 0 DC 10` | DC voltage source |
| Current Source | I | `I1 n1 n2 DC 1` | DC current source |
| Switch | S | `S1 n1 n2 0.001` | Switch (optional on-resistance) |

### Value Suffixes

| Suffix | Multiplier | Example |
|--------|------------|---------|
| k | 1e3 | `1k` = 1000 |
| m | 1e-3 | `1m` = 0.001 |
| u | 1e-6 | `1u` = 0.000001 |
| n | 1e-9 | `1n` = 1e-9 |
| p | 1e-12 | `1p` = 1e-12 |

### Analysis Commands

**DC Operating Point:**
```
.analysis dc
```

**AC Analysis:**
```
.analysis ac <frequency_hz>
```

**Transient Analysis:**
```
.analysis tran <time_step> <end_time>
```

### Output Commands

**Print to Console:**
```
.print v(<node>)      ; Node voltage
.print i(<component>) ; Component current
```

**Write to CSV:**
```
.write <filename.csv>
```

### Comments

Lines starting with `*` or `;` are comments. Inline comments after `;` are also supported.

## Examples

### DC Analysis - Resistor Divider

```cir
* Resistor Divider Circuit
V1 1 0 DC 12
R1 1 2 2k
R2 2 0 1k

.analysis dc
.print v(2)
.end
```

### Transient Analysis - RC Circuit

```cir
* RC Circuit Charging
V1 1 0 DC 10
R1 1 2 1k
C1 2 0 1u

.analysis tran 10u 5m
.print v(2)
.write rc_output.csv
.end
```

### Switch Operation

```cir
* Circuit with Switch
V1 1 0 DC 100
R1 1 2 10
S1 2 3 0.001
R2 3 0 50
L1 3 0 1m

.analysis tran 1u 10m
.print v(3)
.end
```

## Output Format

### Console Output

DC analysis prints voltages directly:
```
=== DC Analysis Results ===
V(2): 4.000000 V
V(3): 2.500000 V
```

Transient analysis shows progress and CSV file location:
```
=== Transient Analysis ===
Time step: 1.00e-05 s, End time: 5.00e-03 s, Steps: 500
Progress: 0% (t=0.000000 ms)
Progress: 10% (t=0.510000 ms)
...
Results written to: output.csv
```

### CSV Output

CSV files contain time and specified outputs:
```csv
time,V(2)
0.000000000e+00,0.000000000e+00
1.000000000e-05,9.950206325e-02
...
```

## Error Handling

Errors are reported with line numbers:
```
Error at line 5: Unknown component type
Error: Failed to parse netlist
```

Common errors:
- Unknown component type (unsupported prefix)
- Failed to create nodes (invalid node numbers)
- No analysis specified (missing `.analysis` command)
