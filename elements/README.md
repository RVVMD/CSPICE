# CSPICE Circuit Elements

Implementation of circuit components for Modified Nodal Analysis.

## Structure

```
elements/
├── passive.c/h         # R, L, C, Switch
├── sources.c/h         # Voltage and Current sources
└── nonlinear/
    ├── nonlinear.c/h   # Diodes and custom nonlinear elements
    └── ...             # Future: CT/VT saturation, arresters
```

## Component Interface

All elements implement a stamping interface that contributes to the MNA system:

```c
// Example: Resistor stamping
// Contributes to conductance matrix G:
//     node1    node2
// node1  +G      -G
// node2  -G      +G
// where G = 1/R
```

## Passive Elements (`passive.c`)

### Resistor

**API:**
```c
MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2, 
                           double value, ComponentHandle* handle);
```

**MNA Stamp:**
```
┌             ┐
│  1/R  -1/R  │
│ -1/R   1/R  │
└             ┘
```

**Equation:** `i = v/R = G·v`

---

### Capacitor

**API:**
```c
MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2,
                            double value, ComponentHandle* handle);
```

**Time-Domain Companion Model (TR-BDF2):**
```
i_C(t) = C · dv/dt ≈ G_eq · v(t) + I_eq

where:
G_eq = C / (γ · dt)
I_eq = -G_eq · v_history
```

**Stamp:**
```
Conductance: G_eq at (node1,node1) and (node2,node2)
Current:     I_eq injected at node1, -I_eq at node2
```

---

### Inductor

**API:**
```c
MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2,
                           double value, ComponentHandle* handle);
```

**Time-Domain Companion Model (TR-BDF2):**
```
v_L(t) = L · di/dt ≈ R_eq · i(t) + V_eq

where:
R_eq = L / (γ · dt)
V_eq = -R_eq · i_history
```

**Stamp:**
Requires additional current variable (MNA branch current):
```
┌             ┐ ┌     ┐   ┌         ┐
│  0     1    │ │ v1  │   │ 0       │
│ -1    R_eq  │ │ i_L │ = │ V_eq    │
└             ┘ └     ┘   └         ┘
```

---

### Switch

**API:**
```c
MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2,
                         double value, ComponentHandle* handle);
MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state);
```

**Model:**
- **Closed**: Small resistance R_on (default 1 mΩ)
- **Open**: Large resistance R_off (implemented as open circuit)

**Usage:**
```c
ComponentHandle sw;
mna_add_switch(&solver, n1, n2, 0.001, &sw);

// During transient simulation
mna_set_switch_state(&solver, sw, 1);  // Close at t = t_fault
```

---

## Sources (`sources.c`)

### DC Voltage Source

**API:**
```c
MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2,
                                 double value, ComponentHandle* handle);
```

**MNA Stamp:**
Requires branch current variable:
```
┌         ┐ ┌     ┐   ┌     ┐
│  0   1  │ │ v1  │   │ V_s │
│ -1   0  │ │ i_s │   │  0  │
└         ┘ └     ┘   └     ┘
```

**Equation:** `v1 - v2 = V_s`

---

### DC Current Source

**API:**
```c
MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2,
                                 double value, ComponentHandle* handle);
```

**MNA Stamp:**
Contributes to RHS vector:
```
b[node1] += I_s
b[node2] -= I_s
```

**Equation:** `i = I_s` (from node1 to node2)

---

### AC Source Parameters

**API:**
```c
MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle,
                            double magnitude, double phase);
```

**Usage:**
```c
ComponentHandle vs;
mna_add_voltage_source(&solver, n1, 0, 0.0, &vs);  // DC value = 0
mna_set_ac_source(&solver, vs, 230.0, 0.0);        // 230V RMS, 0° phase
```

---

## Nonlinear Elements (`nonlinear/`)

### Diode

**Model:** Shockley diode equation with series resistance

```
i = I_s · (exp(v/(n·V_T)) - 1)
```

**Newton-Raphson Linearization:**
```
i(v) ≈ i(v_k) + g_k · (v - v_k)

where g_k = di/dv|_vk = (I_s/(n·V_T)) · exp(v_k/(n·V_T))
```

**API:**
```c
MNAStatus mna_add_diode(MNASolver* solver, int node1, int node2,
                        double Is, double n, ComponentHandle* handle);
```

---

### Custom Nonlinear Element

**API:**
```c
typedef void (*CustomNonlinearFunc)(const ComponentState* state, void* user_data,
                                    double* value1, double* value2);

MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                                   CustomNonlinearFunc func, void* user_data,
                                   ComponentHandle* handle);
```

**Callback:**
The function receives the current state and must return:
- `value1`: Current through the element
- `value2`: Small-signal conductance (di/dv)

**Example:**
```c
void my_nonlinear(const ComponentState* state, void* user_data,
                  double* current, double* conductance) {
    double v = state->voltage;
    *current = v * v;           // i = v² (example)
    *conductance = 2.0 * v;     // di/dv = 2v
}
```

---

## Component State Tracking

Each component maintains state for transient analysis:

```c
typedef struct {
    ComponentType type;
    int node1, node2;
    double value;
    bool is_nonlinear;
    
    // Transient state
    double last_voltage;
    double last_current;
    double prev_voltage;
    double prev_current;
    double stage1_voltage;    // TR-BDF2 intermediate stage
    double stage1_current;
    double stage2_voltage;    // TR-BDF2 final stage
    double stage2_current;
    
    // Equivalent circuit for time-stepping
    double trans_G_eq;        // Equivalent conductance
    double trans_I_eq;        // Equivalent current source
} Component;
```

---

## Planned RZA Components (Phase 2)

### Current Transformer (CT)

**Model:** Ideal transformer with magnetizing branch and saturation

```
┌─────────────────┐
│  Ideal  │  Xm   │
│  N1:N2  │  Bsat │
└─────────────────┘
```

**Saturation Curve:**
```
Φ = f(H) with hysteresis (Jiles-Atherton model)
```

---

### Voltage Transformer (VT)

**Model:** Ideal transformer with accuracy class parameters

**Parameters:**
- Ratio (e.g., 110kV/110V)
- Burden (VA rating)
- Accuracy class (0.2, 0.5, 3P)

---

### Fault Element

**Types:**
- Single-line-to-ground (SLG)
- Line-to-line (LL)
- Double-line-to-ground (LLG)
- Three-phase (3P)

**Model:** Time-controlled resistance

```
R_fault(t) = {
    ∞,           t < t_fault
    R_transition, t ≥ t_fault
}
```

---

### Transmission Line

**Model:** Pi-section equivalent

```
     ┌───R/2───L/2───┐
-----│               │-----
     │               │
    C/2             C/2
     │               │
-----┴---------------┴-----
```

**API (planned):**
```c
MNAStatus mna_add_transmission_line(MNASolver* solver, 
                                    int node_from, int node_to,
                                    double R, double L, double C,
                                    int segments, ComponentHandle* handle);
```

---

## Value Suffixes

All component values support standard suffixes:

| Suffix | Multiplier | Example |
|--------|------------|---------|
| T | 1e12 | `1T` = 1 tera |
| G | 1e9 | `1G` = 1 giga |
| k | 1e3 | `1k` = 1 kilo |
| m | 1e-3 | `1m` = 1 milli |
| u | 1e-6 | `1u` = 1 micro |
| n | 1e-9 | `1n` = 1 nano |
| p | 1e-12 | `1p` = 1 pico |
| f | 1e-15 | `1f` = 1 femto |

---

## Error Handling

All element functions return `MNAStatus`:

```c
MNAStatus status = mna_add_resistor(&solver, n1, n2, 1000.0, &handle);
if (status != MNA_SUCCESS) {
    // Handle error
}
```

Common errors:
- `MNA_INVALID_NODE`: Node index out of range
- `MNA_INVALID_PARAMETER`: Negative resistance, zero capacitance, etc.
- `MNA_INSUFFICIENT_MEMORY`: Cannot allocate component slot

---

## References

1. J. Vlach and K. Singhal, "Computer Methods for Circuit Analysis and Design"
2. L.O. Chua and P.M. Lin, "Computer-Aided Analysis of Electronic Circuits"
3. IEEE C37 series for CT/VT modeling standards
