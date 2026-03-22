#ifndef MNA_ELEMENTS_TRANSFORMER_SAT_H
#define MNA_ELEMENTS_TRANSFORMER_SAT_H

#include "../include/types.h"
#include "nonlinear/magnetization.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * TransformerSat - Nonlinear transformer model with saturation and magnetization branch
 * 
 * This model represents a real transformer with:
 * - Ideal transformer coupling (turns ratio n = N1/N2)
 * - Magnetization branch (nonlinear inductor) in parallel with primary
 * - Optional series resistance and leakage inductance on both windings
 * 
 * The magnetization branch uses a B-H curve (tanh model by default) to capture
 * core saturation effects, enabling simulation of inrush currents.
 * 
 * Circuit topology:
 * 
 *   Primary side:                    Secondary side:
 *   p1 o----[R1]----[L1]----+----o s1
 *                           |
 *   p2 o--------------------+----o s2
 *                           |
 *                    [Lmagnetizing]
 *                    (nonlinear)
 *                           |
 *   (internal node m) ------+
 * 
 * Terminals:
 *   - Terminal 0: p1 (primary positive)
 *   - Terminal 1: p2 (primary negative)
 *   - Terminal 2: s1 (secondary positive)
 *   - Terminal 3: s2 (secondary negative)
 */

typedef struct TransformerSat TransformerSat;

/* ============================================================================
 * Creation and Configuration
 * ============================================================================ */

/**
 * Create a new nonlinear transformer with saturation
 * 
 * @param solver         MNA solver instance
 * @param primary_p      Primary winding positive node
 * @param primary_n      Primary winding negative node
 * @param secondary_p    Secondary winding positive node
 * @param secondary_n    Secondary winding negative node
 * @param turns_ratio    Turns ratio n = N_primary/N_secondary
 * @param handle         Output component handle
 * @return MNAStatus
 */
MNAStatus mna_add_transformer_sat(MNASolver* solver,
                                   int primary_p, int primary_n,
                                   int secondary_p, int secondary_n,
                                   double turns_ratio,
                                   ComponentHandle* handle);

/**
 * Get the internal transformer structure for configuration
 * 
 * @param solver  MNA solver instance
 * @param handle  Component handle from mna_add_transformer_sat
 * @return Pointer to TransformerSat structure, or NULL if invalid
 */
TransformerSat* mna_get_transformer_sat(MNASolver* solver, ComponentHandle handle);

/**
 * Set the magnetization B-H curve for the transformer
 * 
 * The B-H curve defines the nonlinear magnetization characteristic.
 * Use mna_bh_curve_create() and mna_bh_curve_set_tanh_params() to create.
 * 
 * @param xf     Transformer instance
 * @param curve  B-H curve (transformer takes ownership, will free on destroy)
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_set_bh_curve(TransformerSat* xf, MagnetizationCurve* curve);

/**
 * Set primary winding series resistance
 * 
 * @param xf         Transformer instance
 * @param resistance Series resistance [ohm]
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_set_primary_resistance(TransformerSat* xf, double resistance);

/**
 * Set secondary winding series resistance
 * 
 * @param xf         Transformer instance
 * @param resistance Series resistance [ohm]
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_set_secondary_resistance(TransformerSat* xf, double resistance);

/**
 * Set primary winding leakage inductance
 * 
 * @param xf         Transformer instance
 * @param inductance Leakage inductance [H]
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_set_primary_leakage_inductance(TransformerSat* xf, double inductance);

/**
 * Set secondary winding leakage inductance
 * 
 * @param xf         Transformer instance
 * @param inductance Leakage inductance [H]
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_set_secondary_leakage_inductance(TransformerSat* xf, double inductance);

/**
 * Set core loss resistance (parallel with magnetization branch)
 * 
 * Models eddy current and hysteresis losses in the core.
 * 
 * @param xf         Transformer instance
 * @param resistance Core loss resistance [ohm] (referenced to primary)
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_set_core_loss_resistance(TransformerSat* xf, double resistance);

/* ============================================================================
 * Langevin B-H Curve Convenience Setup
 * ============================================================================ */

/**
 * Configure transformer with physics-based Langevin B-H model
 *
 * The Langevin function: L(x) = coth(x) - 1/x
 * B(H) = B_sat * L(H / H_c)
 *
 * This model naturally saturates without artificial linear terms,
 * producing sharper saturation knees and narrower inrush current spikes.
 *
 * @param xf              Transformer instance
 * @param B_sat_T         Saturation flux density [T]
 * @param mu_r_initial    Initial relative permeability (at H->0)
 * @param H_c_A_m         Critical field scale [A/m] (0 = auto-compute)
 * @param core_area_m2    Core cross-sectional area [m^2]
 * @param path_length_m   Magnetic path length [m]
 * @param N_primary       Number of primary turns
 * @return MNAStatus
 */
MNAStatus mna_transformer_sat_setup_langevin_core(TransformerSat* xf,
                                                   double B_sat_T,
                                                   double mu_r_initial,
                                                   double H_c_A_m,
                                                   double core_area_m2,
                                                   double path_length_m,
                                                   int N_primary);

/* ============================================================================
 * State Queries
 * ============================================================================ */

/**
 * Get magnetizing current (current through magnetization branch)
 * 
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Magnetizing current [A]
 */
double mna_transformer_sat_get_magnetizing_current(MNASolver* solver, ComponentHandle handle);

/**
 * Get flux linkage in the core
 * 
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Flux linkage [Wb]
 */
double mna_transformer_sat_get_flux_linkage(MNASolver* solver, ComponentHandle handle);

/**
 * Get primary winding current
 * 
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Primary current [A]
 */
double mna_transformer_sat_get_primary_current(MNASolver* solver, ComponentHandle handle);

/**
 * Get secondary winding current
 * 
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Secondary current [A]
 */
double mna_transformer_sat_get_secondary_current(MNASolver* solver, ComponentHandle handle);

/**
 * Get instantaneous inductance of magnetization branch
 * 
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Incremental inductance [H]
 */
double mna_transformer_sat_get_magnetizing_inductance(MNASolver* solver, ComponentHandle handle);

/**
 * Get magnetic flux density B in the core
 *
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Flux density B [T]
 */
double mna_transformer_sat_get_B_field(MNASolver* solver, ComponentHandle handle);

/**
 * Get magnetic field intensity H in the core
 *
 * @param solver  MNA solver instance
 * @param handle  Component handle
 * @return Field intensity H [A/m]
 */
double mna_transformer_sat_get_H_field(MNASolver* solver, ComponentHandle handle);

/* Internal cleanup function - called from mna_destroy */
void mna_transformer_sat_cleanup(Component* comp);

/* Internal state update functions - called from transient solver */
void mna_transformer_sat_init_state(Component* comp);
void mna_transformer_sat_update_state(Component* comp, double* solver_x, double dt, int stage);
void mna_transformer_sat_finalize_state(Component* comp);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_TRANSFORMER_SAT_H */
