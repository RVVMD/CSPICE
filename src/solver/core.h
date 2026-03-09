#ifndef MNA_SOLVER_CORE_H
#define MNA_SOLVER_CORE_H

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Solver Lifecycle
 * ============================================================================ */

/**
 * @brief Initialize the MNA solver
 * @param solver The MNA solver to initialize
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_init(MNASolver* solver);

/**
 * @brief Destroy the MNA solver and free all resources
 * @param solver The MNA solver to destroy
 */
void mna_destroy(MNASolver* solver);

/* ============================================================================
 * Solver Configuration
 * ============================================================================ */

/**
 * @brief Set the integration method for transient analysis
 * @param solver The MNA solver
 * @param method Integration method (currently only TRBDF2 supported)
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_set_integration_method(MNASolver* solver, IntegrationMethod method);

/**
 * @brief Enable or disable DC state preservation for transient analysis
 * @param solver The MNA solver
 * @param enable Non-zero to preserve DC state, zero to compute initial conditions
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_set_preserve_dc_state(MNASolver* solver, int enable);

/* ============================================================================
 * Node Management
 * ============================================================================ */

/**
 * @brief Create a new node in the circuit
 * @param solver The MNA solver
 * @return Node index (1-based), or -1 on failure
 */
int mna_create_node(MNASolver* solver);

/* ============================================================================
 * Component Management (Internal)
 * ============================================================================ */

/**
 * @brief Add a component to the circuit (internal use)
 * @param solver The MNA solver
 * @param type Component type
 * @param node1 First node index
 * @param node2 Second node index
 * @param value Component value
 * @param src_type Source type (for sources)
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_component(MNASolver* solver, ComponentType type, int node1, int node2,
                            double value, SourceType src_type, ComponentHandle* handle);

/**
 * @brief Validate node indices
 * @param solver The MNA solver
 * @param node1 First node index
 * @param node2 Second node index
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_validate_nodes(MNASolver* solver, int node1, int node2);

/* ============================================================================
 * Solution Access
 * ============================================================================ */

/**
 * @brief Get the DC voltage at a node
 * @param solver The MNA solver
 * @param node Node index (1-based)
 * @return Node voltage
 */
double mna_get_node_voltage(MNASolver* solver, int node);

/**
 * @brief Get the AC voltage at a node (complex)
 * @param solver The MNA solver
 * @param node Node index (1-based)
 * @return Complex node voltage
 */
double complex mna_get_ac_node_voltage(MNASolver* solver, int node);

/**
 * @brief Get the current through a component
 * @param solver The MNA solver
 * @param handle Component handle
 * @return Component current
 */
double mna_get_component_current(MNASolver* solver, ComponentHandle handle);

/**
 * @brief Get the voltage across a component
 * @param solver The MNA solver
 * @param handle Component handle
 * @return Component voltage
 */
double mna_get_component_voltage(MNASolver* solver, ComponentHandle handle);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_CORE_H */
