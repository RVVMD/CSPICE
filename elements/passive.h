#ifndef MNA_ELEMENTS_PASSIVE_H
#define MNA_ELEMENTS_PASSIVE_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Resistor
 * ============================================================================ */

/**
 * @brief Add a resistor to the circuit
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param value Resistance in ohms
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2, 
                           double value, ComponentHandle* handle);

/* ============================================================================
 * Capacitor
 * ============================================================================ */

/**
 * @brief Add a capacitor to the circuit
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param value Capacitance in farads
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2, 
                            double value, ComponentHandle* handle);

/* ============================================================================
 * Inductor
 * ============================================================================ */

/**
 * @brief Add an inductor to the circuit
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param value Inductance in henries
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2, 
                           double value, ComponentHandle* handle);

/* ============================================================================
 * Switch
 * ============================================================================ */

/**
 * @brief Add a switch to the circuit
 * @param solver The MNA solver
 * @param node1 First node index (1-based, 0 is ground)
 * @param node2 Second node index (1-based, 0 is ground)
 * @param value On-state resistance in ohms
 * @param handle Output parameter for component handle
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2, 
                         double value, ComponentHandle* handle);

/**
 * @brief Set the state of a switch
 * @param solver The MNA solver
 * @param handle Component handle of the switch
 * @param state New state (1 = closed, 0 = open)
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_PASSIVE_H */
