#ifndef MNA_ELEMENTS_NONLINEAR_H
#define MNA_ELEMENTS_NONLINEAR_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Custom Nonlinear Element
 * ============================================================================ */

/**
 * @brief Add a custom nonlinear element to the circuit
 * 
 * The nonlinear function should compute either:
 * - NONLINEAR_RESISTOR: current and conductance (I-V characteristic)
 * - NONLINEAR_CAPACITOR: charge and capacitance (Q-V characteristic)
 * - NONLINEAR_INDUCTOR: flux and inductance (Phi-I characteristic)
 * 
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param nl_type Type of nonlinear element
 * @param func User-provided function to compute element values
 * @param user_data User data passed to the function
 * @param initial_value1 Initial voltage (for resistor/capacitor) or current (for inductor)
 * @param initial_value2 Initial charge (for capacitor) or flux (for inductor)
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2, 
                                   NonlinearType nl_type,
                                   CustomNonlinearFunc func, void* user_data,
                                   double initial_value1, double initial_value2, 
                                   ComponentHandle* handle);

/* ============================================================================
 * Custom N-Pole Element
 * ============================================================================ */

/**
 * @brief Add a custom n-pole element to the circuit
 * 
 * The stamp function is called during matrix assembly to stamp
 * the element's contribution to the MNA matrix.
 * 
 * @param solver The MNA solver
 * @param nodes Array of node indices (1-based, 0 is ground)
 * @param num_nodes Number of nodes
 * @param stamp_func User-provided stamping function
 * @param user_data User data passed to the stamping function
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_custom_n_pole(MNASolver* solver, const int* nodes, int num_nodes,
                                NPoleStampFunc stamp_func, void* user_data, 
                                ComponentHandle* handle);

/* ============================================================================
 * Nonlinear Element Utilities
 * ============================================================================ */

/**
 * @brief Compute adaptive damping factor based on rate of change
 * @param rate Rate of change (dv/dt or di/dt normalized)
 * @param prev_alpha Previous alpha value
 * @return Computed alpha damping factor
 */
double mna_compute_alpha(double rate, double prev_alpha);

/* ============================================================================
 * Internal Stamping Functions (for solver use)
 * ============================================================================ */

/**
 * @brief Stamp a custom nonlinear element (internal use)
 * @param solver The MNA solver
 * @param comp_index Component index
 * @param is_dc Non-zero for DC analysis, zero for transient
 * @param stage Integration stage (1 or 2 for TRBDF2)
 */
void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc, int stage);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_NONLINEAR_H */
