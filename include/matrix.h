#ifndef MNA_MATRIX_H
#define MNA_MATRIX_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Matrix Access Macro
 * ============================================================================ */
#define MAT(solver, i, j) ((solver)->A[(i) * (solver)->matrix_cap_size + (j)])

/* ============================================================================
 * Matrix Memory Management
 * ============================================================================ */

/**
 * @brief Resize the component array to accommodate more components
 * @param solver The MNA solver
 * @param required Minimum required capacity
 * @return true if successful, false otherwise
 */
bool mna_resize_components(MNASolver* solver, int required);

/**
 * @brief Resize the matrix arrays to accommodate larger systems
 * @param solver The MNA solver
 * @param req_size Required matrix size
 * @return true if successful, false otherwise
 */
bool mna_resize_matrix(MNASolver* solver, int req_size);

/* ============================================================================
 * Matrix Operations
 * ============================================================================ */

/**
 * @brief Reset the matrix system (clear A, b, and x arrays)
 * @param solver The MNA solver
 */
void mna_reset_system(MNASolver* solver);

/**
 * @brief Solve a linear system using Gaussian elimination with partial pivoting
 * @param solver The MNA solver
 * @param size The size of the matrix
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_solve_linear_system(MNASolver* solver, int size);

/* ============================================================================
 * Matrix Stamping Functions
 * ============================================================================ */

/**
 * @brief Stamp a conductance into the matrix
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param g Conductance value
 */
void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g);

/**
 * @brief Stamp a current source into the RHS vector
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param current_val Current value (positive from node1 to node2)
 */
void mna_stamp_current_source(MNASolver* solver, int node1, int node2, double current_val);

/**
 * @brief Stamp a voltage source into the matrix and RHS
 * @param solver The MNA solver
 * @param comp_index Index of the voltage source component
 * @param source_idx Index of this voltage source among all voltage sources
 */
void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx);

/**
 * @brief Ensure all nodes have a DC path to ground
 * @param solver The MNA solver
 */
void mna_ensure_ground_paths(MNASolver* solver);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * @brief Get the active size of the matrix
 * @param solver The MNA solver
 * @return Active matrix size
 */
static inline int mna_active_size(const MNASolver* solver) {
    return solver->max_node_index + solver->num_sources;
}

#ifdef __cplusplus
}
#endif

#endif /* MNA_MATRIX_H */
