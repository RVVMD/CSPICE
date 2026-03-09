#ifndef MNA_ELEMENTS_SOURCES_H
#define MNA_ELEMENTS_SOURCES_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Voltage Source
 * ============================================================================ */

/**
 * @brief Add a DC voltage source to the circuit
 * @param solver The MNA solver
 * @param node1 Positive node index (1-based, 0 is ground)
 * @param node2 Negative node index (1-based, 0 is ground)
 * @param value Voltage in volts
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2, 
                                 double value, ComponentHandle* handle);

/* ============================================================================
 * Current Source
 * ============================================================================ */

/**
 * @brief Add a DC current source to the circuit
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param value Current in amperes (positive from node1 to node2)
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2, 
                                 double value, ComponentHandle* handle);

/* ============================================================================
 * AC Source Configuration
 * ============================================================================ */

/**
 * @brief Set AC parameters for a source
 * @param solver The MNA solver
 * @param handle Component handle of the source
 * @param magnitude AC magnitude
 * @param phase AC phase in radians
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle, 
                            double magnitude, double phase);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_SOURCES_H */
