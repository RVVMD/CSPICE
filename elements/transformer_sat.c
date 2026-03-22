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
    
    /* Terminal nodes */
    int nodes[4];                 /* p1, p2, s1, s2 */
};

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

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
        /* DC: inductor is short circuit */
        if (L > 0 && branch_var > 0) {
            /* Stamp inductor as short with branch current */
            if (n1 > 0) MAT(solver, n1 - 1, branch_var - 1) = 1.0;
            if (n2 > 0) MAT(solver, n2 - 1, branch_var - 1) = -1.0;
            if (branch_var > 0) {
                if (n1 > 0) MAT(solver, branch_var - 1, n1 - 1) = 1.0;
                if (n2 > 0) MAT(solver, branch_var - 1, n2 - 1) = -1.0;
            }
        } else if (R > 0) {
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
    int s1 = nodes[2];  /* Secondary positive */
    int s2 = nodes[3];  /* Secondary negative */
    
    double n = xf->turns_ratio;
    
    /* Internal magnetization node - same as primary terminal for now
     * Note: Series R1/L1 are stamped but have no effect since m1=p1
     * This is a known limitation - use external series R for now
     */
    int m1 = p1;  /* Magnetization branch connects to primary */
    int m2 = p2;
    
    (void)time;
    
    /* ========================================================================
     * 1. Stamp ideal transformer constraint
     *    V_pm1 - V_pm2 = n * (V_s1 - V_s2)
     *    where V_pm1, V_pm2 are voltages at internal primary nodes
     *    For simplicity, we assume magnetization branch is in parallel,
     *    so the ideal transformer sees the same voltage as the terminals
     *    minus the drop across series elements.
     * ======================================================================== */
    
    /* For the basic model, we'll use a simplified approach:
     * - Magnetization branch is directly across primary terminals
     * - Ideal transformer enforces V_primary = n * V_secondary
     * - Series R/L are handled separately
     */
    
    int branch_ideal = xf->branch_var_ideal;
    
    /* Stamp B matrix for ideal transformer constraint */
    /* Primary KCL: +i_ideal */
    if (p1 > 0 && branch_ideal > 0) {
        MAT(solver, p1 - 1, branch_ideal - 1) += 1.0;
    }
    if (p2 > 0 && branch_ideal > 0) {
        MAT(solver, p2 - 1, branch_ideal - 1) += -1.0;
    }
    
    /* Secondary KCL: -n * i_ideal (current scaled by turns ratio) */
    if (s1 > 0 && branch_ideal > 0) {
        MAT(solver, s1 - 1, branch_ideal - 1) += -n;
    }
    if (s2 > 0 && branch_ideal > 0) {
        MAT(solver, s2 - 1, branch_ideal - 1) += n;
    }
    
    /* Constraint equation: V_p1 - V_p2 - n*(V_s1 - V_s2) = 0 */
    if (branch_ideal > 0) {
        if (p1 > 0) MAT(solver, branch_ideal - 1, p1 - 1) += 1.0;
        if (p2 > 0) MAT(solver, branch_ideal - 1, p2 - 1) += -1.0;
        if (s1 > 0) MAT(solver, branch_ideal - 1, s1 - 1) += -n;
        if (s2 > 0) MAT(solver, branch_ideal - 1, s2 - 1) += n;
    }
    
    /* ========================================================================
     * 2. Stamp magnetization branch (nonlinear inductor with B-H curve)
     *    Connected in parallel with primary
     * ======================================================================== */
    
    double flux = xf->flux_linkage;
    double i_mag = xf->magnetizing_current;
    double v_prev = xf->v_prev;

    if (is_dc) {
        /* DC: magnetizing branch is a short (inductor at steady state)
         * But with saturation, we need to handle it carefully.
         * Use a large conductance for DC.
         */
        mna_stamp_conductance(solver, m1, m2, MNA_MAX_CONDUCTANCE);
    } else {
        /* Transient: nonlinear inductor companion model
         * i = i(flux), flux_dot = v
         * Using TR-BDF2: i_new = G_eq * v_new + I_history
         * where I_history = i_old + G_eq * v_old (for trapezoidal)
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

        /* Add core loss resistance in parallel if present */
        if (xf->R_core > 1e-6) {
            G_eq += 1.0 / xf->R_core;
        }

        mna_stamp_conductance(solver, m1, m2, G_eq);
        /* I_history flows from m1 to m2, so stamp +I_history at m1 (leaving) */
        mna_stamp_current_source(solver, m1, m2, I_history);
    }
    
    /* ========================================================================
     * 3. Stamp series elements (leakage inductance and winding resistance)
     *    These are handled as separate branches
     * ======================================================================== */
    
    /* Primary series RL */
    if (xf->R1 > 0 || xf->L1_leak > 0) {
        /* For now, stamp as simple series impedance */
        stamp_series_rl(solver, p1, m1, xf->R1, xf->L1_leak, 
                       xf->i_mag_prev, dt, xf->branch_var_L1, is_dc);
    }
    
    /* Secondary series RL */
    if (xf->R2 > 0 || xf->L2_leak > 0) {
        stamp_series_rl(solver, s1, s2, xf->R2, xf->L2_leak,
                       xf->i_mag_prev / n, dt, xf->branch_var_L2, is_dc);
    }
}

/**
 * Update transformer state after solution
 */
static void update_transformer_state(TransformerSat* xf, MNASolver* solver,
                                     int stage, double dt) {
    if (!xf || !solver) return;
    
    /* Get primary voltage */
    double v_p = 0.0;
    if (xf->nodes[0] > 0) v_p += solver->x[xf->nodes[0] - 1];
    if (xf->nodes[1] > 0) v_p -= solver->x[xf->nodes[1] - 1];
    
    /* Update flux linkage: d(flux)/dt = v */
    if (dt > 0 && stage >= 0) {
        double gamma = MNA_TRBDF2_GAMMA;
        
        if (stage == 1) {
            /* Stage 1: flux = flux_prev + (gamma*dt/2) * v */
            double h1 = gamma * dt;
            xf->flux_linkage = xf->flux_prev + (h1 / 2.0) * v_p;
        } else {
            /* Stage 2: flux = flux_stage1 + ((1-gamma)*dt/2) * v */
            double h2 = (1.0 - gamma) * dt;
            xf->flux_linkage = xf->flux_stage1 + (h2 / 2.0) * v_p;
        }
        
        /* Update magnetizing current from B-H curve */
        xf->magnetizing_current = get_magnetizing_current(xf, xf->flux_linkage);
        
        if (stage == 1) {
            xf->flux_stage1 = xf->flux_linkage;
            xf->i_mag_stage1 = xf->magnetizing_current;
        }
    }
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
    
    /* Allocate branch variables */
    /* branch_var_ideal: for ideal transformer constraint */
    /* branch_var_L1: for primary leakage inductance (if used) */
    /* branch_var_L2: for secondary leakage inductance (if used) */
    
    int base_index = solver->max_node_index + solver->num_sources + 1;
    
    /* We need 1 branch variable for the ideal transformer constraint */
    if (!mna_resize_matrix(solver, base_index + 1)) {
        free(xf);
        return MNA_INSUFFICIENT_MEMORY;
    }
    
    xf->branch_var_ideal = base_index;
    xf->branch_var_L1 = 0;  /* Not used by default */
    xf->branch_var_L2 = 0;  /* Not used by default */
    
    solver->num_sources++;  /* Account for branch variable */
    
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

MNAStatus mna_transformer_sat_setup_tanh_core(TransformerSat* xf,
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
    
    /* Create B-H curve with tanh model */
    MagnetizationCurve* curve = mna_bh_curve_create(BH_MODEL_ANALYTIC_TANH);
    if (!curve) return MNA_INSUFFICIENT_MEMORY;
    
    /* Set geometry */
    MNAStatus status = mna_bh_curve_set_geometry(curve, core_area_m2, path_length_m, N_primary);
    if (status != MNA_SUCCESS) {
        mna_bh_curve_destroy(curve);
        return status;
    }
    
    /* Set tanh parameters */
    status = mna_bh_curve_set_tanh_params(curve, B_sat_T, mu_r_initial, H_c_A_m);
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

    NPoleData* npole = comp->data.npole.npole_data;
    if (!npole || npole->num_branch_currents < 1) return 0.0;

    /* Primary current = ideal transformer branch current + magnetizing current */
    int idx = npole->branch_current_indices[0];
    if (idx <= 0) return 0.0;

    double i_ideal = solver->x[idx - 1];
    return i_ideal + xf->magnetizing_current;
}

double mna_transformer_sat_get_secondary_current(MNASolver* solver, ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return 0.0;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    if (!xf) return 0.0;

    NPoleData* npole = comp->data.npole.npole_data;
    if (!npole || npole->num_branch_currents < 1) return 0.0;

    /* Get ideal transformer branch current (before adding magnetizing current) */
    int idx = npole->branch_current_indices[0];
    if (idx <= 0) return 0.0;

    double i_ideal = solver->x[idx - 1];

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
    }
}

/* Update transformer state during transient step (called after each stage) */
void mna_transformer_sat_update_state(Component* comp, double* solver_x, double dt, int stage) {
    if (!comp || comp->type != MNA_CUSTOM_NPOLE || !comp->data.npole.npole_data) return;

    TransformerSat* xf = (TransformerSat*)comp->data.npole.npole_data->user_data;
    NPoleData* npole = comp->data.npole.npole_data;
    if (!xf || !npole || npole->num_nodes < 2) return;

    /* Get primary voltage from solver solution */
    double v_p = 0.0;
    if (npole->nodes[0] > 0) v_p += solver_x[npole->nodes[0] - 1];
    if (npole->nodes[1] > 0) v_p -= solver_x[npole->nodes[1] - 1];

    /* Update flux linkage: d(flux)/dt = v */
    double gamma = MNA_TRBDF2_GAMMA;
    if (stage == 1) {
        double h1 = gamma * dt;
        xf->flux_linkage = xf->flux_prev + (h1 / 2.0) * v_p;
        xf->v_stage1 = v_p;
    } else {
        double h2 = (1.0 - gamma) * dt;
        xf->flux_linkage = xf->flux_stage1 + (h2 / 2.0) * v_p;
        xf->v_prev = v_p;
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
