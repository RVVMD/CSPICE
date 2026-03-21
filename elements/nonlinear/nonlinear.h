#ifndef MNA_ELEMENTS_NONLINEAR_H
#define MNA_ELEMENTS_NONLINEAR_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_add_custom_nonlinear(MNASolver* solver, int node1, int node2,
                                   NonlinearType nl_type,
                                   CustomNonlinearFunc func, void* user_data,
                                   double initial_value1, double initial_value2,
                                   ComponentHandle* handle);

MNAStatus mna_add_custom_n_pole(MNASolver* solver, const int* nodes, int num_nodes,
                                NPoleStampFunc stamp_func, void* user_data,
                                int num_branch_currents, ComponentHandle* handle);

double mna_compute_alpha(double rate, double prev_alpha);

void mna_stamp_custom_nonlinear(MNASolver* solver, int comp_index, int is_dc, int stage);

/* ============================================================================
 * N-Pole Element Query Utilities
 * 
 * Safe extraction of voltages and currents from n-pole elements.
 * These utilities provide a consistent API for accessing n-pole state
 * without directly accessing solver internals.
 * ============================================================================ */

/**
 * Get voltage at a specific terminal of an n-pole element
 * 
 * @param solver         MNA solver instance
 * @param handle         Component handle for the n-pole element
 * @param terminal_index Index of the terminal (0-based, must be < num_nodes)
 * @return Node voltage at the terminal, or 0.0 if terminal is ground or invalid
 */
double mna_get_npole_node_voltage(MNASolver* solver, ComponentHandle handle,
                                   int terminal_index);

/**
 * Get voltage difference between two terminals of an n-pole element
 * 
 * @param solver            MNA solver instance
 * @param handle            Component handle for the n-pole element
 * @param positive_terminal Index of the positive terminal (0-based)
 * @param negative_terminal Index of the negative terminal (0-based)
 * @return Voltage at positive_terminal minus voltage at negative_terminal
 */
double mna_get_npole_terminal_voltage_diff(MNASolver* solver, ComponentHandle handle,
                                            int positive_terminal, int negative_terminal);

/**
 * Get branch current for an n-pole element
 * 
 * N-pole elements can have multiple branch currents
 * one for the ideal constraint). This function retrieves the current at
 * the specified branch index.
 * 
 * @param solver        MNA solver instance
 * @param handle        Component handle for the n-pole element
 * @param branch_index  Index of the branch current (0-based, must be < num_branch_currents)
 * @return Branch current value, or 0.0 if invalid
 */
double mna_get_npole_branch_current(MNASolver* solver, ComponentHandle handle,
                                     int branch_index);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_NONLINEAR_H */
