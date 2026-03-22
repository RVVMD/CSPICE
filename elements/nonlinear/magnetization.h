#ifndef MNA_ELEMENTS_MAGNETIZATION_H
#define MNA_ELEMENTS_MAGNETIZATION_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * B-H Curve Model Types
 *
 * ANALYTIC_LANGEVIN: Physics-based Langevin function (coth(x) - 1/x)
 *   - Most accurate for transformer core saturation
 *   - Natural saturation behavior (no artificial linear term)
 *   - Sharper knee than tanh model
 *
 * PIECEWISE_LINEAR: Tabulated data with linear interpolation
 *   - Fast and robust
 *   - Can represent any B-H curve from measured data
 */
typedef enum {
    BH_MODEL_ANALYTIC_LANGEVIN,
    BH_MODEL_PIECEWISE_LINEAR
} BHModelType;

/**
 * MagnetizationCurve - Represents the B-H characteristic of a magnetic core
 * 
 * This structure supports multiple models for representing the nonlinear
 * relationship between magnetic field intensity H [A/m] and magnetic flux
 * density B [T]. The curve can be created from:
 * - Analytic equations (tanh model)
 * - Measured data points (piecewise interpolation)
 * 
 * Core geometry parameters convert between:
 * - B (flux density) <-> φ (flux linkage): φ = B * A_c * N
 * - H (field intensity) <-> i (current): i = H * l_e / N
 */
typedef struct {
    /* Core geometry - required for unit conversions */
    double core_area;            /* A_c: Cross-sectional area [m²] */
    double magnetic_path_length; /* l_e: Effective magnetic path length [m] */
    int N_primary;               /* Primary winding turns */
    bool geometry_set;           /* Flag: has geometry been configured? */
    
    /* B-H model type */
    BHModelType model_type;

    /* === Analytic Langevin Model Parameters ===
     * B(H) = B_sat * L(H/H_c) where L(x) = coth(x) - 1/x
     * This model naturally saturates without artificial linear terms
     */
    struct {
        double B_sat;            /* Saturation flux density [T] */
        double H_c;              /* Critical field scale [A/m] */
        double B_rem;            /* Remanence flux density [T] (optional) */
    } langevin_params;
    
    /* === Piecewise Model Data ===
     * Arrays of (H, B) data points sorted by H
     */
    struct {
        double* H_points;        /* H field values [A/m] */
        double* B_points;        /* B field values [T] */
        int num_points;          /* Number of data points */
        int allocated_size;      /* Allocated array capacity */
    } pw_params;
    
    /* === Hysteresis State (for future extension) ===
     * Currently unused, reserved for Jiles-Atherton model
     */
    struct {
        double B_last;           /* Last B value for hysteresis tracking */
        double H_last;           /* Last H value */
        double hysteresis_loss;  /* Accumulated hysteresis loss [J] */
    } hysteresis_state;
    
} MagnetizationCurve;

/* ============================================================================
 * Creation and Destruction
 * ============================================================================ */

/**
 * Create a new B-H curve with specified model type
 * 
 * @param model_type  The type of B-H model to use
 * @return Pointer to new MagnetizationCurve, or NULL on failure
 */
MagnetizationCurve* mna_bh_curve_create(BHModelType model_type);

/**
 * Destroy a B-H curve and free all associated memory
 * 
 * @param curve  The curve to destroy (can be NULL)
 */
void mna_bh_curve_destroy(MagnetizationCurve* curve);

/**
 * Set core geometry parameters for B-H curve
 * 
 * These parameters are required to convert between:
 * - Magnetic quantities (B [T], H [A/m])
 * - Circuit quantities (φ [Wb], i [A])
 * 
 * @param curve               The B-H curve to configure
 * @param core_area_m2        Core cross-sectional area [m²]
 * @param magnetic_path_length_m  Effective magnetic path length [m]
 * @param N_primary           Number of primary winding turns
 * @return MNAStatus
 */
MNAStatus mna_bh_curve_set_geometry(MagnetizationCurve* curve,
                                     double core_area_m2,
                                     double magnetic_path_length_m,
                                     int N_primary);

/* ============================================================================
 * Analytic Langevin Model Configuration
 * ============================================================================ */

/**
 * Configure analytic Langevin model from physical parameters
 *
 * The Langevin function: L(x) = coth(x) - 1/x
 * B(H) = B_sat * L(H / H_c)
 *
 * Initial permeability: μ_initial = B_sat / (3 * H_c)
 * So H_c = B_sat / (3 * μ_0 * μ_r_initial)
 *
 * @param curve           The B-H curve (must be BH_MODEL_ANALYTIC_LANGEVIN)
 * @param B_sat_T         Saturation flux density [T]
 * @param mu_r_initial    Initial relative permeability (at H→0)
 * @param H_c_A_m         Critical field scale [A/m] (optional, 0 = auto)
 * @return MNAStatus
 */
MNAStatus mna_bh_curve_set_langevin_params(MagnetizationCurve* curve,
                                            double B_sat_T,
                                            double mu_r_initial,
                                            double H_c_A_m);

/* ============================================================================
 * Piecewise Model: Data Point Management
 * ============================================================================ */

/**
 * Add a data point to piecewise B-H curve
 * 
 * Points should be added in order of increasing H. The function will
 * automatically sort if needed. For cubic spline, coefficients are
 * computed on-demand when first queried.
 * 
 * @param curve     The B-H curve (must be PIECEWISE_LINEAR or PIECEWISE_CUBIC)
 * @param H_A_m     Magnetic field intensity [A/m]
 * @param B_Tesla   Magnetic flux density [T]
 * @return MNAStatus
 */
MNAStatus mna_bh_curve_add_point(MagnetizationCurve* curve,
                                  double H_A_m, double B_Tesla);

/**
 * Add multiple data points at once
 * 
 * @param curve       The B-H curve
 * @param H_values    Array of H values [A/m]
 * @param B_values    Array of B values [T]
 * @param num_points  Number of points
 * @return MNAStatus
 */
MNAStatus mna_bh_curve_add_points(MagnetizationCurve* curve,
                                   const double* H_values,
                                   const double* B_values,
                                   int num_points);

/**
 * Clear all data points from piecewise curve
 * 
 * @param curve  The B-H curve to clear
 * @return MNAStatus
 */
MNAStatus mna_bh_curve_clear_points(MagnetizationCurve* curve);

/* ============================================================================
 * Core Queries: B-H Relationship
 * ============================================================================ */

/**
 * Get B field from H field using the configured model
 * 
 * @param curve   The B-H curve
 * @param H_A_m   Magnetic field intensity [A/m]
 * @return B field [T]
 */
double mna_bh_curve_get_B(const MagnetizationCurve* curve, double H_A_m);

/**
 * Get H field from B field (inverse relationship)
 * 
 * Uses interpolation or Newton-Raphson as needed.
 * 
 * @param curve    The B-H curve
 * @param B_Tesla  Magnetic flux density [T]
 * @return H field [A/m]
 */
double mna_bh_curve_get_H(const MagnetizationCurve* curve, double B_Tesla);

/**
 * Get incremental permeability dB/dH at given operating point
 * 
 * This is essential for Newton-Raphson convergence in MNA.
 * 
 * @param curve   The B-H curve
 * @param H_A_m   Magnetic field intensity [A/m]
 * @return dB/dH [T·m/A = H/m]
 */
double mna_bh_curve_get_dB_dH(const MagnetizationCurve* curve, double H_A_m);

/* ============================================================================
 * Circuit-Level Queries (with geometry conversion)
 * ============================================================================ */

/**
 * Get magnetizing current from flux linkage
 * 
 * Converts φ → B → H → i using core geometry.
 * 
 * @param curve      The B-H curve
 * @param flux_Wb    Flux linkage [Wb]
 * @return Magnetizing current [A]
 */
double mna_bh_curve_get_current(const MagnetizationCurve* curve,
                                 double flux_Wb);

/**
 * Get flux linkage from magnetizing current
 * 
 * Converts i → H → B → φ using core geometry.
 * 
 * @param curve       The B-H curve
 * @param current_A   Magnetizing current [A]
 * @return Flux linkage [Wb]
 */
double mna_bh_curve_get_flux(const MagnetizationCurve* curve,
                              double current_A);

/**
 * Get incremental inductance dφ/di at operating point
 * 
 * L_inc = dφ/di = N² * A_c / l_e * (dB/dH)
 * 
 * This is critical for TR-BDF2 integration and Newton-Raphson.
 * 
 * @param curve      The B-H curve
 * @param flux_Wb    Current flux linkage [Wb]
 * @return Incremental inductance [H]
 */
double mna_bh_curve_get_incremental_inductance(const MagnetizationCurve* curve,
                                                double flux_Wb);

/**
 * Get B field from flux linkage
 * 
 * B = φ / (N * A_c)
 * 
 * @param curve      The B-H curve
 * @param flux_Wb    Flux linkage [Wb]
 * @return B field [T]
 */
double mna_bh_curve_B_from_flux(const MagnetizationCurve* curve,
                                 double flux_Wb);

/**
 * Get flux linkage from B field
 * 
 * φ = B * N * A_c
 * 
 * @param curve      The B-H curve
 * @param B_Tesla    B field [T]
 * @return Flux linkage [Wb]
 */
double mna_bh_curve_flux_from_B(const MagnetizationCurve* curve,
                                 double B_Tesla);

/**
 * Get H field from magnetizing current
 * 
 * H = N * i / l_e
 * 
 * @param curve      The B-H curve
 * @param current_A  Current [A]
 * @return H field [A/m]
 */
double mna_bh_curve_H_from_current(const MagnetizationCurve* curve,
                                    double current_A);

/**
 * Get magnetizing current from H field
 * 
 * i = H * l_e / N
 * 
 * @param curve      The B-H curve
 * @param H_A_m      H field [A/m]
 * @return Current [A]
 */
double mna_bh_curve_current_from_H(const MagnetizationCurve* curve,
                                    double H_A_m);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Check if curve is properly configured
 * 
 * @param curve  The B-H curve to check
 * @return true if ready for use, false otherwise
 */
bool mna_bh_curve_is_valid(const MagnetizationCurve* curve);

/**
 * Get the number of data points (for piecewise models)
 * 
 * @param curve  The B-H curve
 * @return Number of points, or 0 for analytic models
 */
int mna_bh_curve_get_point_count(const MagnetizationCurve* curve);

/**
 * Get a specific data point (for piecewise models)
 * 
 * @param curve     The B-H curve
 * @param index     Point index (0 to count-1)
 * @param H_out     Output: H value [A/m]
 * @param B_out     Output: B value [T]
 * @return MNAStatus
 */
MNAStatus mna_bh_curve_get_point(const MagnetizationCurve* curve,
                                  int index, double* H_out, double* B_out);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_MAGNETIZATION_H */
