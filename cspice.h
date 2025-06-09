// cspice.h - Optimized SPICE-like Circuit Simulator Core
#ifndef CSPICE_H
#define CSPICE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define VT 0.02585635
#define DEFAULT_IS 1e-14
#define DEFAULT_N 1.0
#define MIN_CONDUCTANCE 1e-12
#define MIN_PIVOT_MAG 1e-18
#define DEFAULT_AC_FREQ 1000.0
#define MIN_INDUCTANCE 1e-12
#define MIN_DIVISOR 1e-30
#define DIODE_MAX_EXP_ARG 700.0

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double re;
    double im;
} Complex;

static Complex complex_add(Complex a, Complex b) {
    return (Complex){a.re + b.re, a.im + b.im};
}

static Complex complex_sub(Complex a, Complex b) {
    return (Complex){a.re - b.re, a.im - b.im};
}

static Complex complex_mul(Complex a, Complex b) {
    return (Complex){a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

static Complex complex_div(Complex a, Complex b) {
    double den = b.re * b.re + b.im * b.im;
    if (den < MIN_PIVOT_MAG * MIN_PIVOT_MAG) {  // Consistent magnitude check
        return (Complex){0.0, 0.0};
    }
    return (Complex){
        (a.re * b.re + a.im * b.im) / den,
        (a.im * b.re - a.re * b.im) / den
    };
}

static Complex complex_from_polar(double magnitude, double phase_rad) {
    return (Complex){magnitude * cos(phase_rad), magnitude * sin(phase_rad)};
}

static double complex_magnitude(Complex a) {
    return sqrt(a.re * a.re + a.im * a.im);
}

typedef enum {
    RESISTOR,
    CAPACITOR,
    INDUCTOR,
    VOLTAGE_SOURCE,
    CURRENT_SOURCE,
    DIODE
} ComponentType;

typedef enum {
    DC_SOURCE,
    AC_SOURCE,
    PULSE_SOURCE,
    SINE_SOURCE
} SourceType;

typedef struct {
    SourceType type;
    double dc_value;
    double ac_magnitude;
    double ac_phase;
    double ac_freq;
    double pulse[7];  // V1, V2, TD, TR, TF, PW, PER
    double sine[4];   // VO, VA, FREQ, PHASE
} SourceParams;

typedef struct {
    char name[20];
    ComponentType type;
    int node1;
    int node2;
    double value;
    SourceParams source;
    double state;     // Voltage for caps, current for inductors
    double Is;        // Diode saturation current
    double N;         // Diode emission coefficient
} Component;

typedef struct {
    Component *components;
    int component_count;
    int component_capacity;
    int node_count;
    double time_step;
    double rel_tol;
    double abs_tol;
    int max_iter;
    double current_time;
    double *node_voltages;
    int node_voltage_capacity;
} Circuit;

void init_circuit(Circuit *c, double rel_tol, double abs_tol, int max_iter, int initial_component_capacity, int initial_node_capacity);
void free_circuit(Circuit *c);
void set_fixed_time_step(Circuit *c, double step_size);
bool add_resistor(Circuit *c, const char* name, int n1, int n2, double value);
bool add_capacitor(Circuit *c, const char* name, int n1, int n2, double value, double initial_voltage);
bool add_inductor(Circuit *c, const char* name, int n1, int n2, double value, double initial_current);
bool add_voltage_source(Circuit *c, const char* name, int n1, int n2, SourceType stype, double params[]);
bool add_current_source(Circuit *c, const char* name, int n1, int n2, SourceType stype, double params[]);
bool add_diode(Circuit *c, const char* name, int anode_node, int cathode_node, double Is, double N);
int solve_dc(Circuit *c, double *node_voltages_out);
int solve_ac(Circuit *c, double frequency, double *node_voltages_re_out, double *node_voltages_im_out);
void transient_analysis(Circuit *c, double t_end, void (*output_callback)(double t, double *node_voltages));

// Implementation
void init_circuit(Circuit *c, double rel_tol, double abs_tol, int max_iter, int initial_component_capacity, int initial_node_capacity) {
    memset(c, 0, sizeof(Circuit));
    c->time_step = 1e-6;
    c->rel_tol = rel_tol > 0 ? rel_tol : 1e-6;
    c->abs_tol = abs_tol > 0 ? abs_tol : 1e-12;
    c->max_iter = max_iter > 0 ? max_iter : 100;

    c->component_capacity = (initial_component_capacity > 0) ? initial_component_capacity : 16;
    c->components = malloc(c->component_capacity * sizeof(Component));
    if (!c->components) exit(EXIT_FAILURE);

    c->node_voltage_capacity = (initial_node_capacity > 0) ? initial_node_capacity : 8;
    c->node_voltages = calloc(c->node_voltage_capacity, sizeof(double));
    if (!c->node_voltages) {
        free(c->components);
        exit(EXIT_FAILURE);
    }
}

void free_circuit(Circuit *c) {
    free(c->components);
    free(c->node_voltages);
    memset(c, 0, sizeof(Circuit));
}

void set_fixed_time_step(Circuit *c, double step_size) {
    if (step_size > 0) c->time_step = step_size;
}

static bool expand_components(Circuit *c) {
    if (c->component_count >= c->component_capacity) {
        int new_capacity = c->component_capacity * 2;
        Component *new_components = realloc(c->components, new_capacity * sizeof(Component));
        if (!new_components) return false;
        c->components = new_components;
        c->component_capacity = new_capacity;
    }
    return true;
}

static bool expand_node_voltages(Circuit *c, int required_nodes) {
    if (required_nodes >= c->node_voltage_capacity) {
        int new_capacity = c->node_voltage_capacity * 2;
        if (new_capacity <= required_nodes) new_capacity = required_nodes + 8;
        double *new_voltages = realloc(c->node_voltages, new_capacity * sizeof(double));
        if (!new_voltages) return false;
        memset(new_voltages + c->node_voltage_capacity, 0,
               (new_capacity - c->node_voltage_capacity) * sizeof(double));
        c->node_voltages = new_voltages;
        c->node_voltage_capacity = new_capacity;
    }
    return true;
}

bool add_resistor(Circuit *c, const char* name, int n1, int n2, double value) {
    if (value <= 0 || !expand_components(c)) return false;
    Component comp = {0};
    strncpy(comp.name, name, sizeof(comp.name) - 1);
    comp.type = RESISTOR;
    comp.node1 = n1;
    comp.node2 = n2;
    comp.value = value;
    c->components[c->component_count++] = comp;
    if (n1 > c->node_count) c->node_count = n1;
    if (n2 > c->node_count) c->node_count = n2;
    return expand_node_voltages(c, c->node_count);
}

bool add_capacitor(Circuit *c, const char* name, int n1, int n2, double value, double initial_voltage) {
    if (value < 0 || !expand_components(c)) return false;
    Component comp = {0};
    strncpy(comp.name, name, sizeof(comp.name) - 1);
    comp.type = CAPACITOR;
    comp.node1 = n1;
    comp.node2 = n2;
    comp.value = value;
    comp.state = initial_voltage;
    c->components[c->component_count++] = comp;
    if (n1 > c->node_count) c->node_count = n1;
    if (n2 > c->node_count) c->node_count = n2;
    return expand_node_voltages(c, c->node_count);
}

bool add_inductor(Circuit *c, const char* name, int n1, int n2, double value, double initial_current) {
    if (value <= 0 || !expand_components(c)) return false;
    Component comp = {0};
    strncpy(comp.name, name, sizeof(comp.name) - 1);
    comp.type = INDUCTOR;
    comp.node1 = n1;
    comp.node2 = n2;
    comp.value = value;
    comp.state = initial_current;
    c->components[c->component_count++] = comp;
    if (n1 > c->node_count) c->node_count = n1;
    if (n2 > c->node_count) c->node_count = n2;
    return expand_node_voltages(c, c->node_count);
}

bool add_voltage_source(Circuit *c, const char* name, int n1, int n2, SourceType stype, double params[]) {
    if (!expand_components(c)) return false;
    Component comp = {0};
    strncpy(comp.name, name, sizeof(comp.name) - 1);
    comp.type = VOLTAGE_SOURCE;
    comp.node1 = n1;
    comp.node2 = n2;
    comp.source.type = stype;

    switch (stype) {
        case DC_SOURCE:
            comp.source.dc_value = params[0];
            break;
        case AC_SOURCE:
            comp.source.ac_magnitude = params[0];
            comp.source.ac_phase = params[1];
            comp.source.ac_freq = (params[2] > 0) ? params[2] : DEFAULT_AC_FREQ;
            break;
        case PULSE_SOURCE:
            memcpy(comp.source.pulse, params, 7 * sizeof(double));
            if (comp.source.pulse[6] <= 0) comp.source.pulse[6] = 1e-9; // Default period
            break;
        case SINE_SOURCE:
            memcpy(comp.source.sine, params, 4 * sizeof(double));
            break;
    }
    c->components[c->component_count++] = comp;
    if (n1 > c->node_count) c->node_count = n1;
    if (n2 > c->node_count) c->node_count = n2;
    return expand_node_voltages(c, c->node_count);
}

bool add_current_source(Circuit *c, const char* name, int n1, int n2, SourceType stype, double params[]) {
    if (!expand_components(c)) return false;
    Component comp = {0};
    strncpy(comp.name, name, sizeof(comp.name) - 1);
    comp.type = CURRENT_SOURCE;
    comp.node1 = n1;
    comp.node2 = n2;
    comp.source.type = stype;

    switch (stype) {
        case DC_SOURCE:
            comp.source.dc_value = params[0];
            break;
        case AC_SOURCE:
            comp.source.ac_magnitude = params[0];
            comp.source.ac_phase = params[1];
            comp.source.ac_freq = (params[2] > 0) ? params[2] : DEFAULT_AC_FREQ;
            break;
        case PULSE_SOURCE:
            memcpy(comp.source.pulse, params, 7 * sizeof(double));
            if (comp.source.pulse[6] <= 0) comp.source.pulse[6] = 1e-9;
            break;
        case SINE_SOURCE:
            memcpy(comp.source.sine, params, 4 * sizeof(double));
            break;
    }
    c->components[c->component_count++] = comp;
    if (n1 > c->node_count) c->node_count = n1;
    if (n2 > c->node_count) c->node_count = n2;
    return expand_node_voltages(c, c->node_count);
}

bool add_diode(Circuit *c, const char* name, int anode_node, int cathode_node, double Is, double N) {
    if (!expand_components(c)) return false;
    Component comp = {0};
    strncpy(comp.name, name, sizeof(comp.name) - 1);
    comp.type = DIODE;
    comp.node1 = anode_node;
    comp.node2 = cathode_node;
    comp.Is = (Is > 0) ? Is : DEFAULT_IS;
    comp.N = (N > 0) ? N : DEFAULT_N;
    comp.state = 0.0;
    c->components[c->component_count++] = comp;
    if (anode_node > c->node_count) c->node_count = anode_node;
    if (cathode_node > c->node_count) c->node_count = cathode_node;
    return expand_node_voltages(c, c->node_count);
}

double get_source_value(Component *comp, double t) {
    SourceParams *s = &comp->source;
    switch (s->type) {
        case DC_SOURCE:
            return s->dc_value;
        case AC_SOURCE:
            return s->ac_magnitude * sin(2 * M_PI * s->ac_freq * t + s->ac_phase * M_PI / 180.0);
        case PULSE_SOURCE: {
            double V1 = s->pulse[0], V2 = s->pulse[1], TD = s->pulse[2];
            double TR = s->pulse[3], TF = s->pulse[4], PW = s->pulse[5], PER = s->pulse[6];
            if (PER <= 0) PER = 1e-9;
            double cycle_time = fmod(t - TD, PER);
            if (t < TD) return V1;
            if (cycle_time < TR) return V1 + (V2 - V1) * (cycle_time / TR);
            if (cycle_time < TR + PW) return V2;
            if (cycle_time < TR + PW + TF)
                return V2 + (V1 - V2) * ((cycle_time - TR - PW) / TF);
            return V1;
        }
        case SINE_SOURCE: {
            double VO = s->sine[0], VA = s->sine[1];
            double FREQ = s->sine[2], PHASE = s->sine[3];
            return VO + VA * sin(2 * M_PI * FREQ * t + PHASE * M_PI / 180.0);
        }
        default:
            return 0.0;
    }
}

static bool solve_real_linear_system(int size, double *A, double *b, double *x) {
    double *A_copy = malloc(size * size * sizeof(double));
    double *b_copy = malloc(size * sizeof(double));
    if (!A_copy || !b_copy) {
        free(A_copy); free(b_copy);
        return false;
    }

    memcpy(A_copy, A, size * size * sizeof(double));
    memcpy(b_copy, b, size * sizeof(double));

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < size; i++) {
        int pivot = i;
        for (int j = i + 1; j < size; j++) {
            if (fabs(A_copy[j*size + i]) > fabs(A_copy[pivot*size + i]))
                pivot = j;
        }

        if (pivot != i) {
            for (int k = 0; k < size; k++) {
                double tmp = A_copy[i*size + k];
                A_copy[i*size + k] = A_copy[pivot*size + k];
                A_copy[pivot*size + k] = tmp;
            }
            double tmp = b_copy[i];
            b_copy[i] = b_copy[pivot];
            b_copy[pivot] = tmp;
        }

        if (fabs(A_copy[i*size + i]) < MIN_PIVOT_MAG) {
            free(A_copy); free(b_copy);
            return false;
        }

        for (int j = i + 1; j < size; j++) {
            double factor = A_copy[j*size + i] / A_copy[i*size + i];
            for (int k = i; k < size; k++) {
                A_copy[j*size + k] -= factor * A_copy[i*size + k];
            }
            b_copy[j] -= factor * b_copy[i];
        }
    }

    // Back substitution
    for (int i = size - 1; i >= 0; i--) {
        x[i] = b_copy[i];
        for (int j = i + 1; j < size; j++) {
            x[i] -= A_copy[i*size + j] * x[j];
        }
        x[i] /= A_copy[i*size + i];
    }

    free(A_copy); free(b_copy);
    return true;
}

static bool solve_complex_linear_system(int size, Complex *A, Complex *b, Complex *x) {
    Complex *A_copy = malloc(size * size * sizeof(Complex));
    Complex *b_copy = malloc(size * sizeof(Complex));
    if (!A_copy || !b_copy) {
        free(A_copy); free(b_copy);
        return false;
    }

    memcpy(A_copy, A, size * size * sizeof(Complex));
    memcpy(b_copy, b, size * sizeof(Complex));

    for (int i = 0; i < size; i++) {
        int pivot = i;
        for (int j = i + 1; j < size; j++) {
            if (complex_magnitude(A_copy[j*size + i]) >
                complex_magnitude(A_copy[pivot*size + i])) {
                pivot = j;
            }
        }

        if (pivot != i) {
            for (int k = 0; k < size; k++) {
                Complex tmp = A_copy[i*size + k];
                A_copy[i*size + k] = A_copy[pivot*size + k];
                A_copy[pivot*size + k] = tmp;
            }
            Complex tmp = b_copy[i];
            b_copy[i] = b_copy[pivot];
            b_copy[pivot] = tmp;
        }

        if (complex_magnitude(A_copy[i*size + i]) < MIN_PIVOT_MAG) {
            free(A_copy); free(b_copy);
            return false;
        }

        for (int j = i + 1; j < size; j++) {
            Complex factor = complex_div(A_copy[j*size + i], A_copy[i*size + i]);
            for (int k = i; k < size; k++) {
                A_copy[j*size + k] = complex_sub(A_copy[j*size + k],
                                              complex_mul(factor, A_copy[i*size + k]));
            }
            b_copy[j] = complex_sub(b_copy[j], complex_mul(factor, b_copy[i]));
        }
    }

    // Back substitution
    for (int i = size - 1; i >= 0; i--) {
        x[i] = b_copy[i];
        for (int j = i + 1; j < size; j++) {
            x[i] = complex_sub(x[i], complex_mul(A_copy[i*size + j], x[j]));
        }
        x[i] = complex_div(x[i], A_copy[i*size + i]);
    }

    free(A_copy); free(b_copy);
    return true;
}

int solve_dc(Circuit *c, double *node_voltages_out) {
    int n_nodes = c->node_count;
    int num_voltage_sources = 0;

    // Count voltage sources and validate components
    for (int i = 0; i < c->component_count; i++) {
        if (c->components[i].type == VOLTAGE_SOURCE) num_voltage_sources++;
        if (c->components[i].type == INDUCTOR && c->components[i].value < MIN_INDUCTANCE) {
            c->components[i].value = MIN_INDUCTANCE;
        }
    }

    int matrix_size = n_nodes + num_voltage_sources;
    if (matrix_size == 0) return 1;  // Empty circuit

    double *G = calloc(matrix_size * matrix_size, sizeof(double));
    double *I = calloc(matrix_size, sizeof(double));
    double *x = calloc(matrix_size, sizeof(double));
    double *delta_x = calloc(matrix_size, sizeof(double));
    if (!G || !I || !x || !delta_x) {
        free(G); free(I); free(x); free(delta_x);
        return 0;
    }

    // Initialize with previous solution
    for (int i = 0; i < n_nodes; i++) {
        x[i] = c->node_voltages[i];
    }

    int converged = 0, iter = 0;
    while (iter < c->max_iter && !converged) {
        memset(G, 0, matrix_size * matrix_size * sizeof(double));
        memset(I, 0, matrix_size * sizeof(double));
        int vs_index = 0;

        for (int comp_index = 0; comp_index < c->component_count; comp_index++) {
            Component *comp = &c->components[comp_index];
            int n1 = comp->node1 - 1, n2 = comp->node2 - 1;
            double v1 = (n1 >= 0) ? x[n1] : 0.0;
            double v2 = (n2 >= 0) ? x[n2] : 0.0;
            double vd = v1 - v2;

            switch (comp->type) {
                case RESISTOR: {
                    double g = 1.0 / comp->value;
                    if (n1 >= 0) G[n1*matrix_size + n1] += g;
                    if (n2 >= 0) G[n2*matrix_size + n2] += g;
                    if (n1 >= 0 && n2 >= 0) {
                        G[n1*matrix_size + n2] -= g;
                        G[n2*matrix_size + n1] -= g;
                    }
                    break;
                }
                case INDUCTOR: {
                    double g = 1.0 / comp->value;
                    if (n1 >= 0) G[n1*matrix_size + n1] += g;
                    if (n2 >= 0) G[n2*matrix_size + n2] += g;
                    if (n1 >= 0 && n2 >= 0) {
                        G[n1*matrix_size + n2] -= g;
                        G[n2*matrix_size + n1] -= g;
                    }
                    break;
                }
                case VOLTAGE_SOURCE: {
                    int k = n_nodes + vs_index++;
                    if (n1 >= 0) {
                        G[n1*matrix_size + k] = 1.0;
                        G[k*matrix_size + n1] = 1.0;
                    }
                    if (n2 >= 0) {
                        G[n2*matrix_size + k] = -1.0;
                        G[k*matrix_size + n2] = -1.0;
                    }
                    I[k] = get_source_value(comp, 0.0);
                    break;
                }
                case CURRENT_SOURCE: {
                    double current = get_source_value(comp, 0.0);
                    if (n1 >= 0) I[n1] -= current;
                    if (n2 >= 0) I[n2] += current;
                    break;
                }
                case DIODE: {
                    double exp_term = (vd > DIODE_MAX_EXP_ARG * comp->N * VT) ? exp(DIODE_MAX_EXP_ARG) :
                                     (vd < -DIODE_MAX_EXP_ARG * comp->N * VT) ? 0.0 :
                                     exp(vd / (comp->N * VT));
                    double id = comp->Is * (exp_term - 1.0);
                    double gd = (comp->Is / (comp->N * VT)) * exp_term;
                    double ieq = id - gd * vd;
                    if (n1 >= 0) {
                        G[n1*matrix_size + n1] += gd;
                        I[n1] += ieq;
                    }
                    if (n2 >= 0) {
                        G[n2*matrix_size + n2] += gd;
                        I[n2] -= ieq;
                    }
                    if (n1 >= 0 && n2 >= 0) {
                        G[n1*matrix_size + n2] -= gd;
                        G[n2*matrix_size + n1] -= gd;
                    }
                    break;
                }
                default: break;
            }
        }

        // Form residual vector
        for (int row = 0; row < matrix_size; row++) {
            for (int col = 0; col < matrix_size; col++) {
                I[row] -= G[row*matrix_size + col] * x[col];
            }
        }

        // Solve linear system
        if (!solve_real_linear_system(matrix_size, G, I, delta_x)) {
            free(G); free(I); free(x); free(delta_x);
            return 0;
        }

        // Update solution and check convergence
        converged = 1;
        for (int i = 0; i < matrix_size; i++) {
            double prev = x[i];
            x[i] += delta_x[i];
            double abs_err = fabs(delta_x[i]);
            double rel_err = (fabs(prev) > MIN_DIVISOR) ? abs_err / fabs(prev) : abs_err;
            if (abs_err > c->abs_tol && rel_err > c->rel_tol) converged = 0;
        }
        iter++;
    }

    if (!converged) {
        free(G); free(I); free(x); free(delta_x);
        return 0;
    }

    // Update circuit state
    memcpy(node_voltages_out, x, n_nodes * sizeof(double));
    memcpy(c->node_voltages, node_voltages_out, n_nodes * sizeof(double));

    // Update diode voltages
    for (int i = 0; i < c->component_count; i++) {
        Component *comp = &c->components[i];
        if (comp->type == DIODE) {
            int n1 = comp->node1 - 1, n2 = comp->node2 - 1;
            comp->state = (n1 >= 0 ? x[n1] : 0.0) - (n2 >= 0 ? x[n2] : 0.0);
        }
    }

    free(G); free(I); free(x); free(delta_x);
    return 1;
}

int solve_ac(Circuit *c, double frequency, double *node_voltages_re_out, double *node_voltages_im_out) {
    int n_nodes = c->node_count;
    int num_voltage_sources = 0;

    // Count voltage sources
    for (int i = 0; i < c->component_count; i++) {
        if (c->components[i].type == VOLTAGE_SOURCE) num_voltage_sources++;
    }

    int matrix_size = n_nodes + num_voltage_sources;
    if (matrix_size == 0) return 1;  // Empty circuit

    Complex *Y = calloc(matrix_size * matrix_size, sizeof(Complex));
    Complex *J = calloc(matrix_size, sizeof(Complex));
    Complex *V = calloc(matrix_size, sizeof(Complex));
    double *dc_voltages = malloc(n_nodes * sizeof(double));
    if (!Y || !J || !V || !dc_voltages) {
        free(Y); free(J); free(V); free(dc_voltages);
        return 0;
    }

    // Solve DC operating point
    if (!solve_dc(c, dc_voltages)) {
        free(Y); free(J); free(V); free(dc_voltages);
        return 0;
    }

    double omega = 2 * M_PI * frequency;
    int vs_index = 0;

    // Build complex admittance matrix
    for (int comp_index = 0; comp_index < c->component_count; comp_index++) {
        Component *comp = &c->components[comp_index];
        int n1 = comp->node1 - 1, n2 = comp->node2 - 1;
        Complex admittance = {0};

        switch (comp->type) {
            case RESISTOR:
                admittance.re = 1.0 / comp->value;
                break;
            case CAPACITOR:
                admittance.im = omega * comp->value;
                break;
            case INDUCTOR:
                if (comp->value < MIN_INDUCTANCE) {
                    admittance.re = 1.0 / MIN_CONDUCTANCE;
                } else {
                    admittance.im = -1.0 / (omega * comp->value);
                }
                break;
            case VOLTAGE_SOURCE: {
                int k = n_nodes + vs_index++;
                if (n1 >= 0) {
                    Y[n1*matrix_size + k] = (Complex){1.0, 0.0};
                    Y[k*matrix_size + n1] = (Complex){1.0, 0.0};
                }
                if (n2 >= 0) {
                    Y[n2*matrix_size + k] = (Complex){-1.0, 0.0};
                    Y[k*matrix_size + n2] = (Complex){-1.0, 0.0};
                }
                if (comp->source.type == AC_SOURCE) {
                    J[k] = complex_from_polar(comp->source.ac_magnitude,
                                            comp->source.ac_phase * M_PI / 180.0);
                } else {
                    J[k].re = 0;  // DC sources are zero in AC analysis
                }
                break;
            }
            case CURRENT_SOURCE: {
                Complex current = {0};
                if (comp->source.type == AC_SOURCE) {
                    current = complex_from_polar(comp->source.ac_magnitude,
                                               comp->source.ac_phase * M_PI / 180.0);
                }
                if (n1 >= 0) J[n1] = complex_sub(J[n1], current);
                if (n2 >= 0) J[n2] = complex_add(J[n2], current);
                break;
            }
            case DIODE: {
                double vd = (n1 >= 0 ? dc_voltages[n1] : 0.0) - (n2 >= 0 ? dc_voltages[n2] : 0.0);
                double exp_term = (vd > DIODE_MAX_EXP_ARG * comp->N * VT) ? exp(DIODE_MAX_EXP_ARG) :
                                 (vd < -DIODE_MAX_EXP_ARG * comp->N * VT) ? 0.0 :
                                 exp(vd / (comp->N * VT));
                admittance.re = (comp->Is / (comp->N * VT)) * exp_term;
                break;
            }
            default: break;
        }

        // Stamp non-source components
        if (comp->type != VOLTAGE_SOURCE && comp->type != CURRENT_SOURCE) {
            if (n1 >= 0) {
                Y[n1*matrix_size + n1] = complex_add(Y[n1*matrix_size + n1], admittance);
            }
            if (n2 >= 0) {
                Y[n2*matrix_size + n2] = complex_add(Y[n2*matrix_size + n2], admittance);
            }
            if (n1 >= 0 && n2 >= 0) {
                Y[n1*matrix_size + n2] = complex_sub(Y[n1*matrix_size + n2], admittance);
                Y[n2*matrix_size + n1] = complex_sub(Y[n2*matrix_size + n1], admittance);
            }
        }
    }

    // Solve complex linear system
    if (!solve_complex_linear_system(matrix_size, Y, J, V)) {
        free(Y); free(J); free(V); free(dc_voltages);
        return 0;
    }

    // Copy results
    for (int i = 0; i < n_nodes; i++) {
        node_voltages_re_out[i] = V[i].re;
        node_voltages_im_out[i] = V[i].im;
    }

    free(Y); free(J); free(V); free(dc_voltages);
    return 1;
}

void transient_analysis(Circuit *c, double t_end, void (*output_callback)(double t, double *node_voltages)) {
    double t = 0.0, dt = c->time_step;
    int n_nodes = c->node_count;
    int num_voltage_sources = 0;

    // Count voltage sources
    for (int i = 0; i < c->component_count; i++) {
        if (c->components[i].type == VOLTAGE_SOURCE) num_voltage_sources++;
    }

    int matrix_size = n_nodes + num_voltage_sources;
    if (matrix_size == 0) return;  // Empty circuit

    double *G = calloc(matrix_size * matrix_size, sizeof(double));
    double *I = calloc(matrix_size, sizeof(double));
    double *x = calloc(matrix_size, sizeof(double));
    double *delta_x = calloc(matrix_size, sizeof(double));
    if (!G || !I || !x || !delta_x) return;

    // Use existing node voltages as initial condition
    memcpy(x, c->node_voltages, n_nodes * sizeof(double));
    if (output_callback) output_callback(t, c->node_voltages);

    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;
        c->current_time = t + dt;

        int converged = 0, iter = 0;
        while (iter < c->max_iter && !converged) {
            memset(G, 0, matrix_size * matrix_size * sizeof(double));
            memset(I, 0, matrix_size * sizeof(double));
            int vs_index = 0;

            for (int comp_index = 0; comp_index < c->component_count; comp_index++) {
                Component *comp = &c->components[comp_index];
                int n1 = comp->node1 - 1, n2 = comp->node2 - 1;
                double v1 = (n1 >= 0) ? x[n1] : 0.0;
                double v2 = (n2 >= 0) ? x[n2] : 0.0;
                double vd = v1 - v2;

                switch (comp->type) {
                    case RESISTOR: {
                        double g = 1.0 / comp->value;
                        if (n1 >= 0) G[n1*matrix_size + n1] += g;
                        if (n2 >= 0) G[n2*matrix_size + n2] += g;
                        if (n1 >= 0 && n2 >= 0) {
                            G[n1*matrix_size + n2] -= g;
                            G[n2*matrix_size + n1] -= g;
                        }
                        break;
                    }
                    case CAPACITOR: {
                        double geq = comp->value / dt;
                        double ieq = geq * comp->state;
                        if (n1 >= 0) {
                            G[n1*matrix_size + n1] += geq;
                            I[n1] += ieq;  // CHANGED FROM -= TO +=
                        }
                        if (n2 >= 0) {
                            G[n2*matrix_size + n2] += geq;
                            I[n2] -= ieq;  // CHANGED FROM += TO -=
                        }
                        if (n1 >= 0 && n2 >= 0) {
                            G[n1*matrix_size + n2] -= geq;
                            G[n2*matrix_size + n1] -= geq;
                        }
                        break;
                    }
                    case INDUCTOR: {
                        double geq = dt / comp->value;
                        double ieq = comp->state;
                        if (n1 >= 0) {
                            G[n1*matrix_size + n1] += geq;
                            I[n1] -= ieq;  // CHANGED FROM += TO -=
                        }
                        if (n2 >= 0) {
                            G[n2*matrix_size + n2] += geq;
                            I[n2] += ieq;  // CHANGED FROM -= TO +=
                        }
                        if (n1 >= 0 && n2 >= 0) {
                            G[n1*matrix_size + n2] -= geq;
                            G[n2*matrix_size + n1] -= geq;
                        }
                        break;
                    }
                    case VOLTAGE_SOURCE: {
                        int k = n_nodes + vs_index++;
                        if (n1 >= 0) {
                            G[n1*matrix_size + k] = 1.0;
                            G[k*matrix_size + n1] = 1.0;
                        }
                        if (n2 >= 0) {
                            G[n2*matrix_size + k] = -1.0;
                            G[k*matrix_size + n2] = -1.0;
                        }
                        I[k] = get_source_value(comp, c->current_time);
                        break;
                    }
                    case CURRENT_SOURCE: {
                        double current = get_source_value(comp, c->current_time);
                        if (n1 >= 0) I[n1] -= current;
                        if (n2 >= 0) I[n2] += current;
                        break;
                    }
                    case DIODE: {
                        double exp_term = (vd > DIODE_MAX_EXP_ARG * comp->N * VT) ? exp(DIODE_MAX_EXP_ARG) :
                                         (vd < -DIODE_MAX_EXP_ARG * comp->N * VT) ? 0.0 :
                                         exp(vd / (comp->N * VT));
                        double id = comp->Is * (exp_term - 1.0);
                        double gd = (comp->Is / (comp->N * VT)) * exp_term;
                        double ieq = id - gd * vd;
                        if (n1 >= 0) {
                            G[n1*matrix_size + n1] += gd;
                            I[n1] += ieq;
                        }
                        if (n2 >= 0) {
                            G[n2*matrix_size + n2] += gd;
                            I[n2] -= ieq;
                        }
                        if (n1 >= 0 && n2 >= 0) {
                            G[n1*matrix_size + n2] -= gd;
                            G[n2*matrix_size + n1] -= gd;
                        }
                        break;
                    }
                    default: break;
                }
            }

            // Form residual vector
            for (int row = 0; row < matrix_size; row++) {
                for (int col = 0; col < matrix_size; col++) {
                    I[row] -= G[row*matrix_size + col] * x[col];
                }
            }

            // Solve linear system
            if (!solve_real_linear_system(matrix_size, G, I, delta_x)) {
                free(G); free(I); free(x); free(delta_x);
                return;
            }

            // Update solution and check convergence
            converged = 1;
            for (int i = 0; i < matrix_size; i++) {
                double prev = x[i];
                x[i] += delta_x[i];
                double abs_err = fabs(delta_x[i]);
                double rel_err = (fabs(prev) > MIN_DIVISOR) ? abs_err / fabs(prev) : abs_err;
                if (abs_err > c->abs_tol && rel_err > c->rel_tol) converged = 0;
            }
            iter++;
        }

        // Update component states
        for (int i = 0; i < c->component_count; i++) {
            Component *comp = &c->components[i];
            int n1 = comp->node1 - 1, n2 = comp->node2 - 1;
            double v1 = (n1 >= 0) ? x[n1] : 0.0;
            double v2 = (n2 >= 0) ? x[n2] : 0.0;
            double vd = v1 - v2;

            if (comp->type == CAPACITOR) {
                comp->state = vd;  // Store capacitor voltage
            }
            else if (comp->type == INDUCTOR) {
                // Update inductor current: I = I_prev + (dt/L)*V
                comp->state += (dt / comp->value) * vd;
            }
            else if (comp->type == DIODE) {
                comp->state = vd;  // Store diode voltage
            }
        }

        // Update node voltages
        memcpy(c->node_voltages, x, n_nodes * sizeof(double));
        t += dt;
        if (output_callback) output_callback(t, c->node_voltages);
    }
    free(G); free(I); free(x); free(delta_x);
}

#endif
