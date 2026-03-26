#include "../include/types.h"
#include "transformer_sat.h"
#include "../include/matrix.h"
#include "../src/solver/core.h"
#include "nonlinear/magnetization.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Physical constants */
#define MU_0 (4.0 * M_PI * 1e-7)  /* Vacuum permeability [H/m] */

/**
 * TransformerSat - Internal structure for nonlinear transformer
 */
struct TransformerSat {
    /* Electrical parameters */
    double turns_ratio;          /* n = N_primary / N_secondary */
    double R1;                   /* Primary series resistance [ohm] */
    double R2;                   /* Secondary series resistance [ohm] */
    double L1_leak;              /* Primary leakage inductance [H] */
    double L2_leak;              /* Secondary leakage inductance [H] */
    double R_core;               /* Core loss resistance [ohm] (parallel) */
    
    /* Magnetic parameters */
    MagnetizationCurve* bh_curve; /* B-H curve for magnetization branch */
    int N_primary;                /* Primary turns (for geometry conversion) */
    
    /* State variables */
    double flux_linkage;          /* Core flux linkage [Wb] */
    double flux_prev;             /* Previous flux (for time integration) */
    double flux_stage1;           /* Stage 1 flux (TR-BDF2) */
    double magnetizing_current;   /* Current through magnetization branch [A] */
    double i_mag_prev;            /* Previous magnetizing current */
    double i_mag_stage1;          /* Stage 1 magnetizing current */
    double v_prev;                /* Previous primary voltage (for companion model) */
    double v_stage1;              /* Stage 1 primary voltage */
    
    /* N-pole data */
    int branch_var_ideal;         /* Branch current for ideal transformer constraint */
    int branch_var_L1;            /* Branch current for primary leakage inductance */
    int branch_var_L2;            /* Branch current for secondary leakage inductance */
    int component_index;          /* Component index for dynamic branch calculation */

    /* Terminal nodes */
    int nodes[4];                 /* p1, p2, s1, s2 (external terminals) */
    int node_internal_pri;        /* Internal primary magnetization node (after R1/L1) */
    int node_internal_sec;        /* Internal secondary node (before R2/L2) */
};

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * Calculate branch index for transformer at stamping time.
 * This matches how voltage sources calculate their indices.
 *
 * Each transformer has 1 branch: ideal constraint.
 * L1 uses a companion model without a branch variable.
 *
 * Branch indices are assigned in component order:
 * - Voltage sources: max_node_index + 0, max_node_index + 1, ...
 * - Transformers: interleaved with voltage sources based on component order
 */
static int get_transformer_branch_index(const MNASolver* solver, int comp_idx) {
    int branch_idx = solver->max_node_index;

    /* Count all branches (voltage sources + transformers) before this component */
    for (int i = 0; i < comp_idx; i++) {
        const Component* comp = &solver->components[i];
        if (comp->type == MNA_SOURCE && comp->source_type == SOURCE_VOLTAGE) {
            branch_idx++;  /* Voltage source has 1 branch */
        } else if (comp->type == MNA_CUSTOM_NPOLE) {
            branch_idx += 1;  /* Transformer has 1 branch */
        }
    }

    return branch_idx;
}

/**
 * Get magnetizing current from flux linkage using B-H curve
 */
static double get_magnetizing_current(TransformerSat* xf, double flux) {
    if (!xf->bh_curve || !mna_bh_curve_is_valid(xf->bh_curve)) {
        /* Linear model if no B-H curve: i = flux / L_m */
        double L_m = 1.0;  /* Default magnetizing inductance */
        return flux / L_m;
    }
    return mna_bh_curve_get_current(xf->bh_curve, flux);
}

/**
 * Get incremental magnetizing inductance d(flux)/di
 */
static double get_magnetizing_inductance(TransformerSat* xf, double flux) {
    if (!xf->bh_curve || !mna_bh_curve_is_valid(xf->bh_curve)) {
        return 1.0;  /* Default inductance */
    }
    return mna_bh_curve_get_incremental_inductance(xf->bh_curve, flux);
}

/**
 * Get flux from magnetizing current (inverse)
 */
static double get_flux_from_current(TransformerSat* xf, double i_mag) {
    if (!xf->bh_curve || !mna_bh_curve_is_valid(xf->bh_curve)) {
        return i_mag * 1.0;  /* Linear model */
    }
    return mna_bh_curve_get_flux(xf->bh_curve, i_mag);
}

/**
 * Stamp series RL branch for transient analysis
 * Uses companion model: G_eq in parallel with I_eq
 */
static void stamp_series_rl(MNASolver* solver, int n1, int n2,
                            double R, double L, double i_prev, double dt,
                            int branch_var, int is_dc) {
    if (is_dc) {
        /* DC: inductor is short, so total impedance is R */
        /* Don't use branch variable - just stamp conductance */
        if (R > 1e-12) {
            mna_stamp_conductance(solver, n1, n2, 1.0 / R);
        } else {
            mna_stamp_conductance(solver, n1, n2, MNA_MAX_CONDUCTANCE);
        }
    } else {
        /* Transient: use companion model */
        double G_eq, I_eq;
        double gamma = MNA_TRBDF2_GAMMA;

        if (L > 1e-12 && dt > 0) {
            /* RL branch: G_eq = dt/(2L) for trapezoidal, modified for TR-BDF2 */
            /* For series RL, we use the inductor companion model */
            double h1 = gamma * dt;
            G_eq = h1 / (2.0 * L);
            I_eq = i_prev;  /* History current source */

            /* Add series resistance to conductance */
            if (R > 1e-12) {
                /* Combined RL: G_eq = 1/(R + 2L/dt) */
                double Z_eq = R + 2.0 * L / h1;
                G_eq = 1.0 / Z_eq;
                I_eq = G_eq * (2.0 * L / h1) * i_prev;
            }
        } else if (R > 1e-12) {
            G_eq = 1.0 / R;
            I_eq = 0.0;
        } else {
            G_eq = MNA_MAX_CONDUCTANCE;
            I_eq = 0.0;
        }

        mna_stamp_conductance(solver, n1, n2, G_eq);
        mna_stamp_current_source(solver, n1, n2, -I_eq);
    }
}

/**
 * N-pole stamp function for transformer with saturation
 */
static void mna_stamp_transformer_sat(MNASolver* solver,
                                       const int* nodes,
                                       int num_nodes,
                                       void* user_data,
                                       double time,
                                       double dt,
                                       int stage) {
    if (num_nodes != 4 || !user_data) return;

    TransformerSat* xf = (TransformerSat*)user_data;
    int is_dc = (dt <= 0 || stage < 0);

    int p1 = nodes[0];  /* Primary positive (external) */
    int p2 = nodes[1];  /* Primary negative (external) */
    int s1 = nodes[2];  /* Secondary positive (external) */
    int s2 = nodes[3];  /* Secondary negative (external) */

    double n = xf->turns_ratio;

    /* Internal nodes:
     * - node_internal_pri (m1): magnetizing node after primary series impedance
     *
     * Topology:
     *   p1 (ext) --[R1/L1]-- m1 (int) --[L_mag]-- p2 (ext)
     *                           |
     *                           |-- Ideal Transformer Primary
     *                           |
     * Ideal Transformer Secondary: constrains V_m1/V_p2 = n * V_s1/V_s2
     *   (secondary leakage L2 is NOT implemented internally - add externally if needed)
     */
    int m1 = xf->node_internal_pri;  /* Primary internal node */
    int m2 = xf->node_internal_sec;  /* Secondary internal node - NOT USED for now */
    int p2_node = p2;  /* Primary return (usually ground) */
    int s2_node = s2;  /* Secondary return */

    /* If R1/L1 is not present, short p1 to m1 */
    int has_primary_rl = (xf->R1 > 1e-12 || xf->L1_leak > 1e-12);
    if (!has_primary_rl) {
        mna_stamp_conductance(solver, p1, m1, MNA_MAX_CONDUCTANCE);
    }

    /* Short m2 to s1 for now (L2_leak not implemented) */
    mna_stamp_conductance(solver, m2, s1, MNA_MAX_CONDUCTANCE);

    (void)time;

    /* ========================================================================
     * 1. Stamp ideal transformer constraint
     *    V_m1 - V_p2 = n * (V_s1 - V_s2)
     *    Note: Uses EXTERNAL secondary voltage (s1, s2), not internal.
     *    Secondary leakage inductance must be added externally if needed.
     * ======================================================================== */

    /* Use dynamic branch index calculation for consistency with voltage sources */
    int branch_ideal = get_transformer_branch_index(solver, xf->component_index);

    /* Stamp B matrix for ideal transformer constraint
     * The branch current i_ideal is defined as flowing from p2 to m1 (primary).
     * By dot convention: i_secondary = n * i_primary (current out of dotted secondary).
     */
    /* Primary KCL: +i_ideal flows into m1 from transformer */
    if (m1 > 0 && branch_ideal > 0) {
        MAT(solver, m1 - 1, branch_ideal) += 1.0;
    }
    if (p2_node > 0 && branch_ideal > 0) {
        MAT(solver, p2_node - 1, branch_ideal) += -1.0;
    }

    /* Secondary KCL: -n * i_ideal flows out of s1 into transformer
     * For ideal transformer: i_s = n * i_p (current out of dotted secondary)
     * This current leaves node s1 into the transformer secondary.
     */
    if (s1 > 0 && branch_ideal > 0) {
        MAT(solver, s1 - 1, branch_ideal) += -n;
    }
    if (s2_node > 0 && branch_ideal > 0) {
        MAT(solver, s2_node - 1, branch_ideal) += n;
    }

    /* Constraint equation: V_m1 - V_p2 - n*(V_s1 - V_s2) = 0 */
    if (branch_ideal > 0) {
        if (m1 > 0) {
            MAT(solver, branch_ideal, m1 - 1) += 1.0;
        }
        if (p2_node > 0) {
            MAT(solver, branch_ideal, p2_node - 1) += -1.0;
        }
        if (s1 > 0) {
            MAT(solver, branch_ideal, s1 - 1) += -n;
        }
        if (s2_node > 0) {
            MAT(solver, branch_ideal, s2_node - 1) += n;
        }
        /* Add small diagonal term for numerical stability during LU factorization */
        MAT(solver, branch_ideal, branch_ideal) += -1e-12;
    }
    
    /* ========================================================================
     * 2. Stamp magnetization branch (nonlinear inductor with B-H curve)
     *    Connected between internal node m1 and m2 (= p2)
     *    This is in parallel with the ideal transformer primary.
     * ======================================================================== */

    double flux = xf->flux_linkage;
    double i_mag = xf->magnetizing_current;
    double v_prev = xf->v_prev;

    if (is_dc) {
        /* DC: magnetizing branch modeling for transformer.
         * At DC, the transformer core flux is constant, so there's no induced voltage
         * across the magnetizing inductance. The magnetizing current is determined
         * by the core's B-H curve and the residual flux.
         *
         * For DC operating point calculation, we model the magnetizing branch as
         * a high resistance (representing the DC resistance of the magnetizing path)
         * in parallel with the core loss resistance.
         *
         * Use a very small conductance to avoid singularity while allowing proper
         * voltage transformation.
         */
        double G_mag = 1e-6;  /* Very high resistance for DC */
        if (xf->R_core > 1e-6) {
            G_mag += 1.0 / xf->R_core;
        }
        mna_stamp_conductance(solver, m1, p2_node, G_mag);
    } else {
        /* Transient: nonlinear inductor companion model
         * i = i(flux), flux_dot = v
         * Using TR-BDF2: i_new = G_eq * v_new + I_history
         * where I_history = i_old + G_eq * v_old (for trapezoidal)
         *
         * The magnetizing current flows from m1 to m2 through Lm || R_core
         */
        double L_inc = get_magnetizing_inductance(xf, flux);
        double gamma = MNA_TRBDF2_GAMMA;

        double G_eq, I_history;

        if (stage == 1) {
            /* Stage 1 of TR-BDF2 */
            double h1 = gamma * dt;
            G_eq = h1 / (2.0 * L_inc);
            /* I_history = i_prev + G_eq * v_prev */
            I_history = i_mag + G_eq * v_prev;
        } else {
            /* Stage 2 of TR-BDF2 */
            double h2 = (1.0 - gamma) * dt;
            G_eq = h2 / (2.0 * L_inc);
            /* I_history = i_stage1 + G_eq * v_stage1 */
            I_history = xf->i_mag_stage1 + G_eq * xf->v_stage1;
        }

        /* Add core loss resistance in parallel if present
         * R_core provides damping for inrush current decay
         */
        if (xf->R_core > 1e-6) {
            G_eq += 1.0 / xf->R_core;
        }

        /* Magnetizing branch is connected between m1 (internal primary) and p2 (primary return) */
        mna_stamp_conductance(solver, m1, p2_node, G_eq);
        /* I_history flows from m1 to p2, stamp as current leaving m1 */
        mna_stamp_current_source(solver, m1, p2_node, I_history);
    }

    /* ========================================================================
     * 3. Stamp series elements (leakage inductance and winding resistance)
     *    Primary series RL: between p1 (external) and m1 (internal)
     *    This ensures inrush current flows through R1, providing damping.
     *    Note: Uses companion model without branch variable.
     * ======================================================================== */

    /* Primary series RL: between p1 (external) and m1 (internal primary) */
    if (xf->R1 > 0 || xf->L1_leak > 0) {
        /* Use companion model without branch variable */
        stamp_series_rl(solver, p1, m1, xf->R1, xf->L1_leak,
                       xf->i_mag_prev, dt, 0, is_dc);
    }

    /* Secondary series RL: NOT IMPLEMENTED
     * The L2_leak should be in series with the ideal transformer secondary,
     * but proper implementation requires a more complex topology.
     * For now, L2_leak is ignored. Users can add external inductors if needed.
     */
    /* L2_leak stamping disabled */
}

/**
 * Update transformer state after solution
 */
static void update_transformer_state(TransformerSat* xf, MNASolver* solver,
                                     int stage, double dt) {
    if (!xf || !solver) return;
    /* This function is no longer used - state update is in mna_transformer_sat_update_state */
    (void)stage;
    (void)dt;
}

/* ============================================================================
 * Public API Implementation
 * ============================================================================ */

MNAStatus mna_add_transformer_sat(MNASolver* solver,
                                   int primary_p, int primary_n,
                                   int secondary_p, int secondary_n,
                                   double turns_ratio,
                                   ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;
    if (turns_ratio <= 0.0) return MNA_INVALID_PARAMETER;

    MNAStatus status = mna_validate_nodes(solver, primary_p, primary_n);
    if (status != MNA_SUCCESS) return status;
    status = mna_validate_nodes(solver, secondary_p, secondary_n);
    if (status != MNA_SUCCESS) return status;

    int new_count = solver->num_components + 1;
    if (!mna_resize_components(solver, new_count)) {
        return MNA_INSUFFICIENT_MEMORY;
    }

    /* Allocate transformer structure */
    TransformerSat* xf = (TransformerSat*)calloc(1, sizeof(TransformerSat));
    if (!xf) return MNA_INSUFFICIENT_MEMORY;

    xf->turns_ratio = turns_ratio;
    xf->R1 = 0.0;
    xf->R2 = 0.0;
    xf->L1_leak = 0.0;
    xf->L2_leak = 0.0;
    xf->R_core = 0.0;  /* No core loss by default */
    xf->bh_curve = NULL;
    xf->N_primary = 1;
    xf->flux_linkage = 0.0;
    xf->flux_prev = 0.0;
    xf->flux_stage1 = 0.0;
    xf->magnetizing_current = 0.0;
    xf->i_mag_prev = 0.0;
    xf->i_mag_stage1 = 0.0;
    xf->v_prev = 0.0;
    xf->v_stage1 = 0.0;
    
    xf->nodes[0] = primary_p;
    xf->nodes[1] = primary_n;
    xf->nodes[2] = secondary_p;
    xf->nodes[3] = secondary_n;

    /* Always allocate internal nodes for the transformer model:
     * - node_internal_pri: magnetizing node after primary series impedance
     * - node_internal_sec: internal secondary node before secondary series impedance
     * If R1/L1 are not present, node_internal_pri is shorted to p1.
     * If R2/L2 are not present, node_internal_sec is shorted to s1.
     */
    xf->node_internal_pri = mna_create_node(solver);
    if (xf->node_internal_pri <= 0) {
        free(xf);
        return MNA_INSUFFICIENT_MEMORY;
    }
    
    xf->node_internal_sec = mna_create_node(solver);
    if (xf->node_internal_sec <= 0) {
        free(xf);
        return MNA_INSUFFICIENT_MEMORY;
    }

    /* Allocate branch variables */
    /* branch_var_ideal: for ideal transformer constraint */
    /* branch_var_L1: for primary leakage inductance (if used) - ONLY for transient */
    /* branch_var_L2: for secondary leakage inductance (if used) */

    int base_index = solver->max_node_index + solver->num_sources + 1;
    /* Only allocate 1 branch per transformer: for the ideal constraint.
     * The L1 branch is NOT allocated because at DC it's not needed (inductor is short),
     * and for transient we can use a companion model without a branch variable.
     * This simplifies the MNA formulation and avoids singular matrices. */
    int branch_count = 1;

    /* Resize matrix to accommodate branch variables */
    if (!mna_resize_matrix(solver, base_index + branch_count)) {
        free(xf);
        return MNA_INSUFFICIENT_MEMORY;
    }

    xf->branch_var_ideal = base_index;
    xf->branch_var_L1 = 0;  /* Not used - companion model doesn't need branch variable */
    xf->branch_var_L2 = 0;  /* Not implemented */

    solver->num_sources += branch_count;
    
    /* Create n-pole component */
    int nodes[4] = {primary_p, primary_n, secondary_p, secondary_n};
    
    Component comp;
    memset(&comp, 0, sizeof(comp));
    comp.type = MNA_CUSTOM_NPOLE;
    comp.node1 = primary_p;
    comp.node2 = primary_n;
    comp.value = turns_ratio;
    comp.state = 1;
    comp.source_type = SOURCE_CURRENT;
    comp.is_nonlinear = true;
    comp.smoothed_alpha = 1.0;
    
    comp.data.npole.npole_data = (NPoleData*)malloc(sizeof(NPoleData));
    if (!comp.data.npole.npole_data) {
        free(xf);
        return MNA_INSUFFICIENT_MEMORY;
    }
    
    NPoleData* npole = comp.data.npole.npole_data;
    npole->nodes = (int*)malloc(4 * sizeof(int));
    npole->last_values = (double*)calloc(4, sizeof(double));
    npole->branch_current_indices = (int*)malloc(sizeof(int));
    
    if (!npole->nodes || !npole->last_values || !npole->branch_current_indices) {
        free(npole->nodes);
        free(npole->last_values);
        free(npole->branch_current_indices);
        free(npole);
        free(xf);
        return MNA_INSUFFICIENT_MEMORY;
    }
    
    memcpy(npole->nodes, nodes, 4 * sizeof(int));
    npole->num_nodes = 4;
    npole->stamp_func = mna_stamp_transformer_sat;
    npole->user_data = xf;
    npole->num_branch_currents = 1;
    npole->branch_current_indices[0] = xf->branch_var_ideal;
    
    int index = solver->num_components++;

    /* Store component index for dynamic branch calculation */
    xf->component_index = index;

    solver->components[index] = comp;
    solver->num_nonlinear++;
    
    if (handle) *handle = index;
    return MNA_SUCCESS;
}

TransformerSat* mna_get_transformer_sat(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return NULL;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return NULL;
    
    return (TransformerSat*)comp->data.npole.npole_data->user_data;
}

MNAStatus mna_transformer_sat_set_bh_curve(TransformerSat* xf, MagnetizationCurve* curve) {
    if (!xf || !curve) return MNA_INVALID_HANDLE;
    
    /* Free existing curve if present */
    if (xf->bh_curve) {
        mna_bh_curve_destroy(xf->bh_curve);
    }
    
    xf->bh_curve = curve;
    return MNA_SUCCESS;
}

MNAStatus mna_transformer_sat_set_primary_resistance(TransformerSat* xf, double resistance) {
    if (!xf) return MNA_INVALID_HANDLE;
    if (resistance < 0) return MNA_INVALID_PARAMETER;
    xf->R1 = resistance;
    return MNA_SUCCESS;
}

MNAStatus mna_transformer_sat_set_secondary_resistance(TransformerSat* xf, double resistance) {
    if (!xf) return MNA_INVALID_HANDLE;
    if (resistance < 0) return MNA_INVALID_PARAMETER;
    xf->R2 = resistance;
    return MNA_SUCCESS;
}

MNAStatus mna_transformer_sat_set_primary_leakage_inductance(TransformerSat* xf, double inductance) {
    if (!xf) return MNA_INVALID_HANDLE;
    if (inductance < 0) return MNA_INVALID_PARAMETER;
    xf->L1_leak = inductance;

    /* Allocate internal node if not already allocated and we now have series impedance */
    if ((xf->R1 > 1e-12 || xf->L1_leak > 1e-12) && xf->node_internal_pri <= 0) {
        /* Note: Can't allocate node here - we don't have access to solver */
        /* Node will be allocated lazily in stamp function or during init */
    }
    
    return MNA_SUCCESS;
}

MNAStatus mna_transformer_sat_set_secondary_leakage_inductance(TransformerSat* xf, double inductance) {
    if (!xf) return MNA_INVALID_HANDLE;
    if (inductance < 0) return MNA_INVALID_PARAMETER;
    xf->L2_leak = inductance;
    return MNA_SUCCESS;
}

MNAStatus mna_transformer_sat_set_core_loss_resistance(TransformerSat* xf, double resistance) {
    if (!xf) return MNA_INVALID_HANDLE;
    if (resistance < 0) return MNA_INVALID_PARAMETER;
    xf->R_core = resistance;
    return MNA_SUCCESS;
}

MNAStatus mna_transformer_sat_setup_langevin_core(TransformerSat* xf,
                                                   double B_sat_T,
                                                   double mu_r_initial,
                                                   double H_c_A_m,
                                                   double core_area_m2,
                                                   double path_length_m,
                                                   int N_primary) {
    if (!xf) return MNA_INVALID_HANDLE;
    if (B_sat_T <= 0 || mu_r_initial <= 0 || core_area_m2 <= 0 ||
        path_length_m <= 0 || N_primary <= 0) {
        return MNA_INVALID_PARAMETER;
    }

    /* Create B-H curve with Langevin model */
    MagnetizationCurve* curve = mna_bh_curve_create(BH_MODEL_ANALYTIC_LANGEVIN);
    if (!curve) return MNA_INSUFFICIENT_MEMORY;

    /* Set geometry */
    MNAStatus status = mna_bh_curve_set_geometry(curve, core_area_m2, path_length_m, N_primary);
    if (status != MNA_SUCCESS) {
        mna_bh_curve_destroy(curve);
        return status;
    }

    /* Set Langevin parameters */
    status = mna_bh_curve_set_langevin_params(curve, B_sat_T, mu_r_initial, H_c_A_m);
    if (status != MNA_SUCCESS) {
        mna_bh_curve_destroy(curve);
        return status;
    }

    xf->bh_curve = curve;
    xf->N_primary = N_primary;

    return MNA_SUCCESS;
}

double mna_transformer_sat_get_magnetizing_current(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;
    
    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf) return 0.0;
    
    return xf->magnetizing_current;
}

double mna_transformer_sat_get_flux_linkage(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;
    
    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf) return 0.0;
    
    return xf->flux_linkage;
}

double mna_transformer_sat_get_primary_current(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf) return 0.0;

    /* Use dynamic branch index calculation (0-indexed) */
    int idx = get_transformer_branch_index(solver, xf->component_index);
    if (idx <= 0) return 0.0;

    double i_ideal = solver->x[idx];
    return i_ideal + xf->magnetizing_current;
}

double mna_transformer_sat_get_secondary_current(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf) return 0.0;

    /* Use dynamic branch index calculation (0-indexed) */
    int idx = get_transformer_branch_index(solver, xf->component_index);
    if (idx <= 0) return 0.0;

    double i_ideal = solver->x[idx];

    /* Secondary current: i_s = n * i_p_ideal (current scales with turns ratio)
     * For step-down (n > 1), secondary current is LARGER than primary
     * The negative sign indicates current direction (dot convention)
     */
    return xf->turns_ratio * i_ideal;
}

double mna_transformer_sat_get_magnetizing_inductance(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf) return 0.0;

    return get_magnetizing_inductance(xf, xf->flux_linkage);
}

double mna_transformer_sat_get_B_field(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf || !xf->bh_curve || !xf->bh_curve->geometry_set) return 0.0;

    /* B = flux / (N * A_c) */
    return xf->flux_linkage / (xf->N_primary * xf->bh_curve->core_area);
}

double mna_transformer_sat_get_H_field(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf || !xf->bh_curve || !xf->bh_curve->geometry_set) return 0.0;

    /* H = N * i_mag / l_e */
    return xf->N_primary * xf->magnetizing_current / xf->bh_curve->magnetic_path_length;
}

/* Cleanup function to be called from mna_destroy */
void mna_transformer_sat_cleanup(Component* comp) {
    if (!comp || comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return;
    
    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (xf) {
        if (xf->bh_curve) {
            mna_bh_curve_destroy(xf->bh_curve);
        }
        free(xf);
    }
}

/* Initialize transformer state for transient analysis */
void mna_transformer_sat_init_state(Component* comp) {
    if (!comp || comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (xf) {
        xf->flux_linkage = 0.0;
        xf->flux_prev = 0.0;
        xf->flux_stage1 = 0.0;
        xf->magnetizing_current = 0.0;
        xf->i_mag_prev = 0.0;
        xf->i_mag_stage1 = 0.0;
        xf->v_prev = 0.0;
        xf->v_stage1 = 0.0;
        /* Don't reset node_internal - it's allocated once and persists */
    }
}

/* Update transformer state during transient step (called after each stage) */
void mna_transformer_sat_update_state(Component* comp, double* solver_x, double dt, int stage) {
    if (!comp || comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    NPoleData* npole = comp->data.npole.npole_data;
    if (!xf || !npole || npole->num_nodes < 2) return;
    if (dt <= 0 || stage < 0) return;

    /* Get voltage across magnetizing branch.
     * The magnetizing branch is connected between node_internal_pri and primary_n (p2).
     */
    double v_mag = 0.0;

    if (xf->node_internal_pri > 0) {
        /* Get voltage at internal primary node */
        if (xf->node_internal_pri > 0) v_mag += solver_x[xf->node_internal_pri - 1];
        if (xf->nodes[1] > 0) v_mag -= solver_x[xf->nodes[1] - 1];
    } else {
        /* Fallback: use terminal voltage directly */
        if (xf->nodes[0] > 0) v_mag += solver_x[xf->nodes[0] - 1];
        if (xf->nodes[1] > 0) v_mag -= solver_x[xf->nodes[1] - 1];
    }

    /* Update flux linkage using trapezoidal integration */
    double gamma = MNA_TRBDF2_GAMMA;
    if (stage == 1) {
        double h1 = gamma * dt;
        xf->flux_linkage = xf->flux_prev + (h1 / 2.0) * (xf->v_prev + v_mag);
        xf->v_stage1 = v_mag;
    } else {
        double h2 = (1.0 - gamma) * dt;
        xf->flux_linkage = xf->flux_stage1 + (h2 / 2.0) * (xf->v_stage1 + v_mag);
        xf->v_prev = v_mag;
    }

    /* Update magnetizing current from B-H curve */
    if (xf->bh_curve && mna_bh_curve_is_valid(xf->bh_curve)) {
        xf->magnetizing_current = mna_bh_curve_get_current(xf->bh_curve, xf->flux_linkage);
    } else {
        /* Linear model fallback */
        xf->magnetizing_current = xf->flux_linkage;
    }

    if (stage == 1) {
        xf->flux_stage1 = xf->flux_linkage;
        xf->i_mag_stage1 = xf->magnetizing_current;
    }
}

/* Finalize transformer state at end of time step */
void mna_transformer_sat_finalize_state(Component* comp) {
    if (!comp || comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (xf) {
        xf->flux_prev = xf->flux_linkage;
        xf->i_mag_prev = xf->magnetizing_current;
        xf->v_prev = xf->v_stage1;  /* Use stage 2 voltage as prev for next step */
    }
}
