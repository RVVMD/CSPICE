#include "transformer.h"
#include "../include/types.h"
#include "../include/matrix.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* Physical constants */
#define MU_0 (4.0 * M_PI * 1e-7)  /* Vacuum permeability [H/m] */

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * Compute magnetizing current from flux using B-H curve or linear model
 */
static double compute_magnetizing_current(TransformerData* xf) {
    if (xf->bh_curve && mna_bh_curve_is_valid(xf->bh_curve)) {
        return mna_bh_curve_get_current(xf->bh_curve, xf->phi);
    }
    
    if (xf->Lm > 0) {
        return xf->phi / xf->Lm;
    }
    
    return 0.0;
}

/**
 * Compute incremental inductance L_inc = dφ/di
 */
static double compute_incremental_inductance(TransformerData* xf) {
    if (xf->bh_curve && mna_bh_curve_is_valid(xf->bh_curve)) {
        return mna_bh_curve_get_incremental_inductance(xf->bh_curve, xf->phi);
    }
    
    if (xf->Lm > 0) {
        return xf->Lm;
    }
    
    return 1e-6;  /* Small default */
}

/**
 * Compute core loss current
 */
static double compute_core_loss_current(TransformerData* xf, double v_primary, double freq) {
    double i_loss = 0.0;
    
    /* Eddy current: i = V/Rc */
    if (xf->Rc > 0 && fabs(v_primary) > 1e-12) {
        i_loss += v_primary / xf->Rc;
    }
    
    /* Hysteresis: simplified model */
    if (xf->hysteresis_coeff > 0 && freq > 0 && xf->core_area > 0) {
        double B = fabs(xf->phi) / xf->core_area;
        double P_hyst = xf->hysteresis_coeff * freq * pow(B, 1.6);
        if (fabs(v_primary) > 1e-9) {
            i_loss += P_hyst / fabs(v_primary) * (v_primary > 0 ? 1 : -1);
        }
    }
    
    return i_loss;
}

/**
 * Post-solve flux update for voltage transformer
 * 
 * This function is called AFTER the matrix solve to update flux state
 * using the newly computed voltages. This is critical for correct
 * transient integration.
 * 
 * TR-BDF2 flux integration:
 *   Stage 1: φ_{n+γ} = φ_n + (γh/2) * (v_n + v_{n+γ})
 *   Stage 2: φ_{n+1} = φ_{n+γ} + ((1-γ)h/2) * (v_{n+γ} + v_{n+1})
 */
void mna_transformer_update_flux(MNASolver* solver,
                                  Component* comp,
                                  int stage,
                                  double dt) {
    if (comp->type != MNA_CUSTOM_NPOLE) return;
    
    NPoleData* npole = comp->data.npole.npole_data;
    if (!npole) return;
    
    TransformerData* xf = (TransformerData*)npole->user_data;
    if (!xf || xf->mode != TRANSFORMER_MODE_VOLTAGE) return;
    
    /* Only update if we have magnetizing branch */
    bool has_magnetizing = (xf->Lm > 0 || xf->bh_curve != NULL);
    if (!has_magnetizing) return;
    
    int p1 = npole->nodes[0];
    int p2 = npole->nodes[1];
    
    /* Get primary voltage from solved values */
    double v_p1 = (p1 > 0) ? solver->x[p1 - 1] : 0.0;
    double v_p2 = (p2 > 0) ? solver->x[p2 - 1] : 0.0;
    double v_primary = v_p1 - v_p2;
    
    double gamma = MNA_TRBDF2_GAMMA;
    
    if (stage == 1) {
        /* Stage 1: t → t + γh */
        double h_eff = gamma * dt;

        /* Flux update: φ_{n+γ} = φ_n + (h_eff/2) * (v_n + v_{n+γ}) */
        /* v_n is approximated by prev_phi/dt, use stored value */
        double v_prev = (xf->prev_phi - xf->phi) / (solver->dt > 0 ? solver->dt : dt);
        xf->phi_stage1 = xf->phi + (h_eff / 2.0) * (v_prev + v_primary);

        /* Update magnetizing current from new flux */
        xf->i_mag_stage1 = compute_magnetizing_current(xf);

        /* Store stage 1 voltage for companion model */
        xf->v_primary_stage1 = v_primary;

    } else {
        /* Stage 2: t + γh → t + h */
        double h_eff = (1.0 - gamma) * dt;

        /* Flux update: φ_{n+1} = φ_{n+γ} + (h_eff/2) * (v_{n+γ} + v_{n+1}) */
        /* v_{n+γ} approximated from stage1 flux */
        double v_stage1 = (xf->phi_stage1 - xf->phi) / (gamma * solver->dt > 0 ? gamma * solver->dt : dt);
        xf->prev_phi = xf->phi;
        xf->phi = xf->phi_stage1 + (h_eff / 2.0) * (v_stage1 + v_primary);

        /* Update magnetizing current from new flux */
        xf->i_mag = compute_magnetizing_current(xf);

        /* Compute incremental inductance at new operating point */
        xf->L_incremental = compute_incremental_inductance(xf);

        /* Store v_{n+1} for next timestep's companion model */
        xf->v_primary_n = v_primary;
    }
}

/**
 * Voltage transformer stamp function
 *
 * Note: This function ONLY stamps the matrix using OLD state.
 * Flux update happens AFTER the solve via mna_transformer_update_flux().
 *
 * TR-BDF2 companion model for magnetizing branch (parallel inductor):
 *   Stage 1: G_eq = γh/(2L), I_eq = i_n + G_eq*v_n
 *   Stage 2: G_eq = (1-γ)h/(2L), I_eq = i_{n+γ} + G_eq*v_{n+γ}
 */
static void voltage_transformer_stamp(MNASolver* solver,
                                       const int* nodes,
                                       TransformerData* xf,
                                       double dt,
                                       int stage) {
    int p1 = nodes[0];
    int p2 = nodes[1];
    int s1 = nodes[2];
    int s2 = nodes[3];
    double n = xf->turns_ratio;
    int branch_idx = xf->branch_current_index;

    bool is_dc = (dt == 0);
    bool has_magnetizing = (xf->Lm > 0 || xf->bh_curve != NULL);

    /* ========================================================================
     * Ideal transformer constraints:
     * Vp - n*Vs = 0  (branch equation, where n = Np/Ns)
     * Ip + n*Is = 0  (current relationship from power conservation)
     *
     * For n = 9.58 (230V/24V step-down):
     *   Vs = Vp/n = Vp/9.58 (voltage step-down)
     *   Is = -n*Ip = -9.58*Ip (current step-up)
     *
     * MNA formulation:
     *   Branch current Ib represents the ideal transformer primary current.
     *   Magnetizing current (I_mag) flows in a parallel branch across primary.
     *   
     *   Primary KCL: Ib + I_mag = I_total_primary
     *   Secondary KCL: -n*Ib + I_load = 0
     * ======================================================================== */
    /* Branch equation: Vp - n*Vs = 0 */
    if (p1 > 0) MAT(solver, branch_idx, p1 - 1) = 1.0;
    if (p2 > 0) MAT(solver, branch_idx, p2 - 1) = -1.0;
    if (s1 > 0) MAT(solver, branch_idx, s1 - 1) = -n;
    if (s2 > 0) MAT(solver, branch_idx, s2 - 1) = n;

    /* Primary KCL: Ib = ideal transformer primary current */
    if (p1 > 0) MAT(solver, p1 - 1, branch_idx) = 1.0;
    if (p2 > 0) MAT(solver, p2 - 1, branch_idx) = -1.0;

    /* Secondary KCL: Is = -n*Ib (secondary current from ideal transformer) */
    if (s1 > 0) MAT(solver, s1 - 1, branch_idx) = -n;
    if (s2 > 0) MAT(solver, s2 - 1, branch_idx) = n;

    /* ========================================================================
     * Magnetizing branch (parallel across primary)
     * 
     * Represents the core magnetizing inductance.
     * Total primary current = Ib (ideal) + I_mag (magnetizing).
     * 
     * NOTE: During no-load (secondary open), the MNA formulation may show
     * non-zero secondary current due to coupling between magnetizing branch
     * and ideal transformer. This is a known limitation. Use I_s2 (load
     * current) for accurate secondary current measurement.
     * 
     * Uses TR-BDF2 companion model:
     *   Stage 1: G_eq = γh/(2L), I_eq = i_n + G_eq*v_n
     *   Stage 2: G_eq = (1-γ)h/(L*(2-γ)), I_eq = factor1*i_ng - factor2*i_n
     * ======================================================================== */
    if (has_magnetizing) {
        xf->L_incremental = compute_incremental_inductance(xf);

        /* Limit for numerical stability */
        if (xf->L_incremental < 1e-12) xf->L_incremental = 1e-12;
        if (xf->L_incremental > 1e6) xf->L_incremental = 1e6;

        if (is_dc) {
            /* DC: small conductance for convergence */
            double G_dc = 1.0 / (xf->L_incremental * 1e-6 + 1e-6);
            if (G_dc > MNA_MAX_CONDUCTANCE) G_dc = MNA_MAX_CONDUCTANCE;
            if (G_dc < MNA_MIN_CONDUCTANCE) G_dc = MNA_MIN_CONDUCTANCE;

            if (p1 > 0) {
                MAT(solver, p1 - 1, p1 - 1) += G_dc;
                if (p2 > 0) MAT(solver, p1 - 1, p2 - 1) -= G_dc;
            }
            if (p2 > 0) {
                if (p1 > 0) MAT(solver, p2 - 1, p1 - 1) -= G_dc;
                MAT(solver, p2 - 1, p2 - 1) += G_dc;
            }
        } else {
            /* Transient: TR-BDF2 integration
             * Use same companion model formulas as standard inductor in transient.c
             */
            double gamma = MNA_TRBDF2_GAMMA;
            double h_eff = (stage == 1) ? (gamma * dt) : ((1.0 - gamma) * dt);

            /* Get alpha_min for numerical damping (same as standard inductor) */
            double alpha_min = 1.0;

            /* Get history values based on stage */
            double i_n = xf->i_mag;           /* i_mag at t_n */
            double v_n = xf->v_primary_n;     /* v_primary at t_n */
            double i_ng, v_ng;                /* values at t_{n+γ} */

            if (stage == 1) {
                /* Stage 1: history is from t_n */
                i_ng = i_n;
                v_ng = v_n;
            } else {
                /* Stage 2: history includes t_{n+γ} */
                i_ng = xf->i_mag_stage1;
                v_ng = xf->v_primary_stage1;
            }

            double G_eq, I_eq;

            if (stage == 1) {
                /* Stage 1: G_eq = h_eff/(2L), I_eq = i_n + G_eq*v_n*alpha_min */
                G_eq = h_eff / (2.0 * xf->L_incremental);
                I_eq = i_n + G_eq * v_n * alpha_min;
            } else {
                /* Stage 2: Match standard inductor formula from transient.c
                 * G_eq = (dt*(1-γ)) / (L*(2-γ))
                 * I_eq = (1/(γ*(2-γ)))*i_ng - (((1-γ)²)/(γ*(2-γ)))*i_n
                 * Note: Standard inductor does NOT include voltage terms in I_eq
                 */
                G_eq = (dt * (1.0 - gamma)) / (xf->L_incremental * (2.0 - gamma));

                double factor1 = 1.0 / (gamma * (2.0 - gamma));
                double factor2 = ((1.0 - gamma) * (1.0 - gamma)) / (gamma * (2.0 - gamma));
                I_eq = factor1 * i_ng - factor2 * i_n;
            }

            if (G_eq > MNA_MAX_CONDUCTANCE) G_eq = MNA_MAX_CONDUCTANCE;
            if (G_eq < MNA_MIN_CONDUCTANCE) G_eq = MNA_MIN_CONDUCTANCE;

            xf->Gm_eq = G_eq;
            xf->I_eq = I_eq;

            /* Stamp G_eq across primary */
            if (p1 > 0) {
                MAT(solver, p1 - 1, p1 - 1) += G_eq;
                if (p2 > 0) MAT(solver, p1 - 1, p2 - 1) -= G_eq;
            }
            if (p2 > 0) {
                if (p1 > 0) MAT(solver, p2 - 1, p1 - 1) -= G_eq;
                MAT(solver, p2 - 1, p2 - 1) += G_eq;
            }

            /* Stamp current source (subtract because moving to RHS) */
            if (p1 > 0) solver->b[p1 - 1] -= I_eq;
            if (p2 > 0) solver->b[p2 - 1] += I_eq;
        }
    }
    
    /* ========================================================================
     * Core loss (parallel across primary)
     * ======================================================================== */
    if (xf->Rc > 0) {
        double Gc = 1.0 / xf->Rc;
        if (Gc > MNA_MAX_CONDUCTANCE) Gc = MNA_MAX_CONDUCTANCE;

        if (p1 > 0) {
            MAT(solver, p1 - 1, p1 - 1) += Gc;
            if (p2 > 0) MAT(solver, p1 - 1, p2 - 1) -= Gc;
        }
        if (p2 > 0) {
            if (p1 > 0) MAT(solver, p2 - 1, p1 - 1) -= Gc;
            MAT(solver, p2 - 1, p2 - 1) += Gc;
        }
    }

    /* Hysteresis loss current */
    if (xf->hysteresis_coeff > 0 && !is_dc) {
        double freq = 50.0;  /* Default power frequency */
        double v_p1 = (p1 > 0) ? solver->x[p1 - 1] : 0.0;
        double v_p2 = (p2 > 0) ? solver->x[p2 - 1] : 0.0;
        double v_primary = v_p1 - v_p2;

        xf->i_core_loss = compute_core_loss_current(xf, v_primary, freq);

        if (p1 > 0) solver->b[p1 - 1] -= xf->i_core_loss;
        if (p2 > 0) solver->b[p2 - 1] += xf->i_core_loss;
    }
}

/**
 * Current transformer stamp function
 */
static void current_transformer_stamp(MNASolver* solver,
                                       const int* nodes,
                                       TransformerData* xf,
                                       double dt,
                                       int stage) {
    (void)stage;  /* CT mode uses simpler integration */
    int s1 = nodes[2];
    int s2 = nodes[3];
    int branch_idx = xf->branch_current_index;
    int burden_idx = xf->burden_branch_index;
    
    double Np = (double)xf->N_primary;
    double Ns = (double)xf->N_secondary;
    
    bool is_dc = (dt == 0);
    bool has_magnetizing = (xf->Lm > 0 || xf->bh_curve != NULL);
    bool has_burden = (xf->burden_resistance > 0);
    
    /* Get secondary voltage */
    double v_s1 = (s1 > 0) ? solver->x[s1 - 1] : 0.0;
    double v_s2 = (s2 > 0) ? solver->x[s2 - 1] : 0.0;
    double v_secondary = v_s1 - v_s2;
    
    /* Flux update from secondary voltage (for CT mode) */
    if (!is_dc && dt > 0 && has_magnetizing) {
        xf->prev_phi = xf->phi;
        xf->phi = xf->phi + v_secondary * dt;
        xf->i_mag = compute_magnetizing_current(xf);
    } else if (has_magnetizing) {
        xf->i_mag = xf->Lm > 0 ? xf->phi / xf->Lm : 0.0;
    }
    
    /* Burden resistor stamp */
    if (has_burden) {
        if (s1 > 0) MAT(solver, burden_idx, s1 - 1) = 1.0;
        if (s2 > 0) MAT(solver, burden_idx, s2 - 1) = -1.0;
        MAT(solver, burden_idx, burden_idx) = -xf->burden_resistance;
        
        if (s1 > 0) MAT(solver, s1 - 1, burden_idx) = 1.0;
        if (s2 > 0) MAT(solver, s2 - 1, burden_idx) = -1.0;
    }
    
    /* Magnetizing branch */
    if (has_magnetizing) {
        xf->L_incremental = compute_incremental_inductance(xf);
        
        if (is_dc) {
            MAT(solver, branch_idx, branch_idx) = 1.0;
        } else {
            double gamma = MNA_TRBDF2_GAMMA;
            double G_eq = (dt * (1.0 - gamma) / 2.0) / xf->L_incremental;
            
            if (G_eq > MNA_MAX_CONDUCTANCE) G_eq = MNA_MAX_CONDUCTANCE;
            if (G_eq < MNA_MIN_CONDUCTANCE) G_eq = MNA_MIN_CONDUCTANCE;
            
            xf->Gm_eq = G_eq;
            xf->I_eq = G_eq * xf->prev_phi + xf->i_mag;
            
            if (s1 > 0) {
                MAT(solver, s1 - 1, s1 - 1) += G_eq;
                if (s2 > 0) MAT(solver, s1 - 1, s2 - 1) -= G_eq;
            }
            if (s2 > 0) {
                if (s1 > 0) MAT(solver, s2 - 1, s1 - 1) -= G_eq;
                MAT(solver, s2 - 1, s2 - 1) += G_eq;
            }
            
            if (s1 > 0) solver->b[s1 - 1] -= xf->I_eq;
            if (s2 > 0) solver->b[s2 - 1] += xf->I_eq;
        }
    }
    
    /* Core loss */
    if (xf->Rc > 0) {
        double Gc = 1.0 / xf->Rc;
        if (Gc > MNA_MAX_CONDUCTANCE) Gc = MNA_MAX_CONDUCTANCE;
        
        if (s1 > 0) {
            MAT(solver, s1 - 1, s1 - 1) += Gc;
            if (s2 > 0) MAT(solver, s1 - 1, s2 - 1) -= Gc;
        }
        if (s2 > 0) {
            if (s1 > 0) MAT(solver, s2 - 1, s1 - 1) -= Gc;
            MAT(solver, s2 - 1, s2 - 1) += Gc;
        }
    }
}

/**
 * Unified transformer stamp function
 */
static void transformer_stamp_func(MNASolver* solver,
                                    const int* nodes,
                                    int num_nodes,
                                    void* user_data,
                                    double time,
                                    double dt,
                                    int stage) {
    (void)time;
    (void)num_nodes;

    TransformerData* xf = (TransformerData*)user_data;

    if (xf->mode == TRANSFORMER_MODE_VOLTAGE) {
        voltage_transformer_stamp(solver, nodes, xf, dt, stage);
    } else {
        current_transformer_stamp(solver, nodes, xf, dt, stage);
    }
}

/* ============================================================================
 * Transformer Creation
 * ============================================================================ */

MNAStatus mna_add_transformer(MNASolver* solver,
                               int node_p1, int node_p2,
                               int node_s1, int node_s2,
                               const TransformerConfig* config,
                               ComponentHandle* handle) {
    if (!solver) return MNA_INVALID_HANDLE;
    if (!config) return MNA_INVALID_PARAMETER;
    if (config->turns_ratio <= 0.0) return MNA_INVALID_PARAMETER;
    
    int nodes[4] = {node_p1, node_p2, node_s1, node_s2};
    for (int i = 0; i < 4; i++) {
        if (nodes[i] < 0 || nodes[i] > solver->max_node_index) {
            return MNA_INVALID_NODE;
        }
    }
    
    /* Allocate transformer data */
    TransformerData* xf = (TransformerData*)calloc(1, sizeof(TransformerData));
    if (!xf) return MNA_INSUFFICIENT_MEMORY;
    
    /* Copy configuration */
    xf->mode = config->mode;
    xf->turns_ratio = config->turns_ratio;
    xf->Lm = config->Lm;
    xf->bh_curve = config->bh_curve;
    xf->Rc = config->Rc;
    xf->hysteresis_coeff = config->hysteresis_coeff;
    xf->eddy_current_coeff = config->eddy_current_coeff;
    xf->core_area = config->core_area;
    xf->core_path_length = config->magnetic_path_length;
    xf->N_primary = config->N_primary;
    xf->N_secondary = config->N_secondary;
    xf->geometry_set = (config->core_area > 0 && config->magnetic_path_length > 0 &&
                        config->N_primary > 0 && config->N_secondary > 0);
    xf->burden_resistance = config->burden_resistance;
    xf->is_ac_analysis = false;
    
    /* Initialize state */
    xf->phi = config->initial_flux;
    xf->prev_phi = config->initial_flux;
    xf->phi_stage1 = config->initial_flux;
    xf->i_mag = config->initial_current;
    xf->i_mag_stage1 = config->initial_current;
    xf->Gm_eq = 0.0;
    xf->I_eq = 0.0;
    xf->v_primary_n = 0.0;
    xf->v_primary_stage1 = 0.0;
    xf->L_incremental = config->Lm > 0 ? config->Lm : 1e-6;
    xf->branch_current_index = -1;
    xf->burden_branch_index = -1;
    xf->sensed_primary_current = 0.0;
    
    /* Determine number of branch currents needed */
    int num_branch_currents = 1;
    if (config->mode == TRANSFORMER_MODE_CURRENT && config->burden_resistance > 0) {
        num_branch_currents = 2;
    }
    
    /* Create n-pole component */
    MNAStatus status = mna_add_custom_n_pole(solver, nodes, 4,
                                              transformer_stamp_func,
                                              xf, num_branch_currents, handle);
    
    if (status != MNA_SUCCESS) {
        free(xf);
        return status;
    }
    
    /* Store branch indices */
    Component* comp = &solver->components[*handle];
    if (comp->type == MNA_CUSTOM_NPOLE && comp->data.npole.npole_data) {
        xf->branch_current_index = comp->data.npole.npole_data->branch_current_indices[0];
        if (num_branch_currents >= 2) {
            xf->burden_branch_index = comp->data.npole.npole_data->branch_current_indices[1];
        }
    }
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Query Functions - Universal
 * ============================================================================ */

double mna_get_transformer_magnetizing_current(MNASolver* solver,
                                                ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    return xf->i_mag;
}

double mna_get_transformer_flux(MNASolver* solver,
                                 ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    return xf->phi;
}

double mna_get_transformer_primary_voltage(MNASolver* solver,
                                            ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    int* nodes = comp->data.npole.npole_data->nodes;
    int p1 = nodes[0];
    int p2 = nodes[1];
    
    double v_p1 = (p1 > 0) ? solver->x[p1 - 1] : 0.0;
    double v_p2 = (p2 > 0) ? solver->x[p2 - 1] : 0.0;
    
    return v_p1 - v_p2;
}

double mna_get_transformer_secondary_voltage(MNASolver* solver,
                                              ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    int* nodes = comp->data.npole.npole_data->nodes;
    int s1 = nodes[2];
    int s2 = nodes[3];
    
    double v_s1 = (s1 > 0) ? solver->x[s1 - 1] : 0.0;
    double v_s2 = (s2 > 0) ? solver->x[s2 - 1] : 0.0;
    
    return v_s1 - v_s2;
}

/* ============================================================================
 * Query Functions - Currents
 * ============================================================================ */

double mna_get_transformer_primary_current(MNASolver* solver,
                                            ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;

    NPoleData* npole = comp->data.npole.npole_data;
    if (!npole || npole->num_branch_currents < 1) return 0.0;

    TransformerData* xf = (TransformerData*)npole->user_data;
    if (!xf) return 0.0;

    if (xf->mode == TRANSFORMER_MODE_VOLTAGE) {
        /* Voltage transformer:
         * Branch current variable = primary current (Ip)
         * Total primary current = Ip + i_mag (load + magnetizing)
         * Note: i_mag flows in parallel with the ideal transformer primary
         */
        int branch_idx = npole->branch_current_indices[0];
        double i_branch = solver->x[branch_idx - 1];

        /* Branch current is the ideal transformer primary current
         * Total primary current includes magnetizing current
         */
        return i_branch + xf->i_mag;
    } else {
        /* Current transformer: return sensed primary current */
        return xf->sensed_primary_current;
    }
}

double mna_get_transformer_secondary_current(MNASolver* solver,
                                              ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;

    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;

    NPoleData* npole = comp->data.npole.npole_data;
    if (!npole || npole->num_branch_currents < 1) return 0.0;

    TransformerData* xf = (TransformerData*)npole->user_data;
    if (!xf) return 0.0;

    if (xf->mode == TRANSFORMER_MODE_VOLTAGE) {
        /* Voltage transformer:
         * Branch current = primary current (Ip)
         * Secondary current: Is = -n * Ip (from ideal transformer relationship)
         */
        int branch_idx = npole->branch_current_indices[0];
        double i_primary = solver->x[branch_idx - 1];
        double n = xf->turns_ratio;

        return -n * i_primary;
    } else {
        /* Current transformer: use dedicated CT function */
        return mna_get_ct_secondary_current(solver, handle);
    }
}

/* ============================================================================
 * Query Functions - Core Magnetics
 * ============================================================================ */

double mna_get_transformer_B_field(MNASolver* solver,
                                    ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    if (!xf->geometry_set || xf->core_area <= 0 || xf->N_secondary <= 0) return 0.0;
    
    int N = (xf->mode == TRANSFORMER_MODE_CURRENT) ? xf->N_secondary : xf->N_primary;
    return xf->phi / ((double)N * xf->core_area);
}

double mna_get_transformer_H_field(MNASolver* solver,
                                    ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    if (!xf->geometry_set || xf->core_path_length <= 0) return 0.0;
    
    int N = (xf->mode == TRANSFORMER_MODE_CURRENT) ? xf->N_secondary : xf->N_primary;
    return (double)N * xf->i_mag / xf->core_path_length;
}

double mna_get_transformer_relative_permeability(MNASolver* solver,
                                                  ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 1.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 1.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    double B = mna_get_transformer_B_field(solver, handle);
    double H = mna_get_transformer_H_field(solver, handle);
    
    if (fabs(H) < 1e-10) return 1.0;
    
    return B / (MU_0 * H);
}

/* ============================================================================
 * Query Functions - Loss
 * ============================================================================ */

double mna_get_transformer_core_loss(MNASolver* solver,
                                      ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    double v = (xf->mode == TRANSFORMER_MODE_VOLTAGE) ?
               mna_get_transformer_primary_voltage(solver, handle) :
               mna_get_transformer_secondary_voltage(solver, handle);
    
    double P_eddy = (xf->Rc > 0) ? (v * v) / xf->Rc : 0.0;
    double P_hyst = v * xf->i_core_loss;
    
    return P_eddy + P_hyst;
}

double mna_get_transformer_incremental_inductance(MNASolver* solver,
                                                   ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    return xf->L_incremental;
}

/* ============================================================================
 * Query Functions - Current Transformer Specific
 * ============================================================================ */

double mna_get_ct_secondary_current(MNASolver* solver,
                                     ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    if (xf->mode != TRANSFORMER_MODE_CURRENT) return 0.0;
    
    if (xf->burden_branch_index > 0) {
        return solver->x[xf->burden_branch_index - 1];
    }
    
    return -xf->i_mag / (double)xf->N_secondary;
}

double mna_get_ct_burden_voltage(MNASolver* solver,
                                  ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    if (xf->mode != TRANSFORMER_MODE_CURRENT) return 0.0;
    
    double Is = mna_get_ct_secondary_current(solver, handle);
    return Is * xf->burden_resistance;
}

double mna_get_ct_primary_current(MNASolver* solver,
                                   ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    if (xf->mode != TRANSFORMER_MODE_CURRENT) return 0.0;
    
    return xf->sensed_primary_current;
}

double mna_get_ct_ratio_error(MNASolver* solver,
                               ComponentHandle handle) {
    if (!solver || handle < 0 || handle >= solver->num_components) return 0.0;
    
    Component* comp = &solver->components[handle];
    if (comp->type != MNA_CUSTOM_NPOLE) return 0.0;
    
    TransformerData* xf = (TransformerData*)comp->data.npole.npole_data->user_data;
    
    if (xf->mode != TRANSFORMER_MODE_CURRENT) return 0.0;
    if (xf->N_primary <= 0 || xf->N_secondary <= 0) return 0.0;
    if (fabs(xf->sensed_primary_current) < 1e-12) return 0.0;
    
    double Is = mna_get_ct_secondary_current(solver, handle);
    double ideal_ratio = (double)xf->N_secondary / (double)xf->N_primary;
    double actual_ratio = -Is / xf->sensed_primary_current;
    
    return (actual_ratio - ideal_ratio) / ideal_ratio;
}
