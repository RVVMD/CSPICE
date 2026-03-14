#ifndef MNA_ELEMENTS_TRANSFORMER_H
#define MNA_ELEMENTS_TRANSFORMER_H

#include "../include/types.h"
#include "nonlinear/nonlinear.h"
#include "nonlinear/magnetization.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Transformer operating mode
 */
typedef enum {
    TRANSFORMER_MODE_VOLTAGE,   /* Voltage transformer (primary connected to voltage source) */
    TRANSFORMER_MODE_CURRENT    /* Current transformer (primary is conductor through toroid) */
} TransformerMode;

/**
 * TransformerConfig - Configuration for unified transformer model
 * 
 * This struct supports both voltage and current transformer modes with
 * optional saturation (B-H curve) and core loss modeling.
 */
typedef struct {
    /* === Required === */
    TransformerMode mode;        /* Operating mode */
    double turns_ratio;          /* N = Np/Ns (for CT: Ns/Np where Np=1) */
    
    /* === Magnetizing Branch === */
    double Lm;                   /* Magnetizing inductance [H] (0 = no magnetizing branch) */
    MagnetizationCurve* bh_curve; /* B-H curve for saturation (NULL = linear) */
    
    /* === Core Loss === */
    double Rc;                   /* Core loss resistance [Ω] parallel across primary (0 = none) */
    double hysteresis_coeff;     /* Hysteresis loss coefficient */
    double eddy_current_coeff;   /* Eddy current loss coefficient */
    
    /* === Core Geometry (required for B-H curve) === */
    double core_area;            /* A_c: Cross-sectional area [m²] */
    double magnetic_path_length; /* l_e: Effective magnetic path length [m] */
    int N_primary;               /* Primary turns (1 for CT) */
    int N_secondary;             /* Secondary turns */
    
    /* === CT-Specific === */
    double burden_resistance;    /* Burden resistor [Ω] across secondary (CT mode) */
    
    /* === Initial Conditions === */
    double initial_flux;         /* Initial flux linkage [Wb] */
    double initial_current;      /* Initial magnetizing current [A] */
    
} TransformerConfig;

/**
 * TransformerData - Internal transformer state
 */
typedef struct {
    /* === Configuration === */
    TransformerMode mode;        /* Operating mode */
    double turns_ratio;          /* N = Np/Ns */
    double Lm;                   /* Magnetizing inductance [H] */
    MagnetizationCurve* bh_curve; /* B-H curve */
    
    /* === MNA Indices === */
    int branch_current_index;    /* Branch current index for primary/ideal constraint */
    int burden_branch_index;     /* Branch current index for burden (CT mode) */
    
    /* === State Variables === */
    double phi;                  /* Flux linkage [Wb] */
    double prev_phi;             /* Previous flux [Wb] */
    double phi_stage1;           /* Flux at stage 1 (TR-BDF2) */
    double i_mag;                /* Magnetizing current [A] */
    double i_mag_stage1;         /* Current at stage 1 (TR-BDF2) */
    
    /* === Companion Model === */
    double Gm_eq;                /* Equivalent conductance [S] */
    double I_eq;                 /* Equivalent current source [A] */
    double v_primary_n;          /* Primary voltage at t_n [V] */
    double v_primary_stage1;     /* Primary voltage at t_{n+γ} [V] */
    
    /* === Core Loss === */
    double Rc;                   /* Core loss resistance [Ω] */
    double i_core_loss;          /* Core loss current [A] */
    double hysteresis_coeff;     /* Hysteresis coefficient */
    double eddy_current_coeff;   /* Eddy current coefficient */
    
    /* === Core Geometry === */
    double core_area;            /* [m²] */
    double core_path_length;     /* [m] */
    int N_primary;               /* Primary turns */
    int N_secondary;             /* Secondary turns */
    bool geometry_set;           /* Geometry configured flag */
    
    /* === CT Specific === */
    double burden_resistance;    /* Burden resistor [Ω] */
    double sensed_primary_current; /* Primary current (CT mode) [A] */
    
    /* === Numerical === */
    double L_incremental;        /* dφ/di at operating point [H] */
    bool is_ac_analysis;         /* AC analysis flag */
    
} TransformerData;

/* ============================================================================
 * Unified Transformer Creation
 * ============================================================================ */

/**
 * Create a transformer (voltage or current mode)
 * 
 * @param solver   MNA solver instance
 * @param node_p1  Primary + node (for CT: unused, primary is external conductor)
 * @param node_p2  Primary - node (for CT: unused)
 * @param node_s1  Secondary + node
 * @param node_s2  Secondary - node
 * @param config   Transformer configuration
 * @param handle   Output: component handle
 * @return MNAStatus
 */
MNAStatus mna_add_transformer(MNASolver* solver,
                               int node_p1, int node_p2,
                               int node_s1, int node_s2,
                               const TransformerConfig* config,
                               ComponentHandle* handle);

/* ============================================================================
 * Query Functions - Universal
 * ============================================================================ */

/**
 * Get magnetizing current
 */
double mna_get_transformer_magnetizing_current(MNASolver* solver,
                                                ComponentHandle handle);

/**
 * Get flux linkage [Wb]
 */
double mna_get_transformer_flux(MNASolver* solver,
                                 ComponentHandle handle);

/**
 * Get primary voltage [V]
 */
double mna_get_transformer_primary_voltage(MNASolver* solver,
                                            ComponentHandle handle);

/**
 * Get secondary voltage [V]
 */
double mna_get_transformer_secondary_voltage(MNASolver* solver,
                                              ComponentHandle handle);

/* ============================================================================
 * Query Functions - Currents
 * ============================================================================ */

/**
 * Get primary current [A]
 * 
 * For voltage transformer: Returns total primary current including
 * both reflected load current and magnetizing current.
 * For current transformer: Returns sensed primary current.
 */
double mna_get_transformer_primary_current(MNASolver* solver,
                                            ComponentHandle handle);

/**
 * Get secondary current [A]
 * 
 * For voltage transformer: Returns secondary winding current.
 * For current transformer: Returns secondary current (same as mna_get_ct_secondary_current).
 */
double mna_get_transformer_secondary_current(MNASolver* solver,
                                              ComponentHandle handle);

/* ============================================================================
 * Query Functions - Core Magnetics (requires geometry_set)
 * ============================================================================ */

/**
 * Get B field in core [T]
 */
double mna_get_transformer_B_field(MNASolver* solver,
                                    ComponentHandle handle);

/**
 * Get H field in core [A/m]
 */
double mna_get_transformer_H_field(MNASolver* solver,
                                    ComponentHandle handle);

/**
 * Get relative permeability μ_r = B/(μ₀H)
 */
double mna_get_transformer_relative_permeability(MNASolver* solver,
                                                  ComponentHandle handle);

/* ============================================================================
 * Query Functions - Loss
 * ============================================================================ */

/**
 * Get core loss power [W]
 */
double mna_get_transformer_core_loss(MNASolver* solver,
                                      ComponentHandle handle);

/**
 * Get incremental inductance L_inc = dφ/di [H]
 */
double mna_get_transformer_incremental_inductance(MNASolver* solver,
                                                   ComponentHandle handle);

/* ============================================================================
 * Query Functions - Current Transformer Specific
 * ============================================================================ */

/**
 * Get CT secondary current [A]
 */
double mna_get_ct_secondary_current(MNASolver* solver,
                                     ComponentHandle handle);

/**
 * Get CT burden voltage [V]
 */
double mna_get_ct_burden_voltage(MNASolver* solver,
                                  ComponentHandle handle);

/**
 * Get CT sensed primary current [A]
 */
double mna_get_ct_primary_current(MNASolver* solver,
                                   ComponentHandle handle);

/**
 * Get CT ratio error (actual_ratio - ideal_ratio) / ideal_ratio
 */
double mna_get_ct_ratio_error(MNASolver* solver,
                               ComponentHandle handle);

/* ============================================================================
 * Internal Functions (used by solver)
 * ============================================================================ */

/**
 * Post-solve flux update for voltage transformers
 * Called by transient solver after each matrix solve
 * @internal - requires types.h to be included
 */
void mna_transformer_update_flux(MNASolver* solver,
                                  Component* comp,
                                  int stage,
                                  double dt);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_TRANSFORMER_H */
