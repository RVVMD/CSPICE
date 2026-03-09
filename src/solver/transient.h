#ifndef MNA_SOLVER_TRANSIENT_H
#define MNA_SOLVER_TRANSIENT_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Transient Analysis
 * ============================================================================ */

/**
 * @brief Initialize transient analysis
 * 
 * Sets up initial conditions for energy storage elements.
 * If preserve_dc_state is enabled, uses the DC operating point.
 * Otherwise, computes initial conditions with capacitors as
 * open circuits and inductors as short circuits.
 * 
 * @param solver The MNA solver
 */
void mna_init_transient(MNASolver* solver);

/**
 * @brief Solve one time step of transient analysis
 * 
 * Uses TRBDF2 (Trapezoidal Rule with Backward Differentiation Formula)
 * integration method with adaptive damping for nonlinear elements.
 * 
 * @param solver The MNA solver
 * @param dt Time step size in seconds
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_solve_transient_step(MNASolver* solver, double dt);

/**
 * @brief Get the current simulation time
 * @param solver The MNA solver
 * @return Current time in seconds
 */
double mna_get_time(MNASolver* solver);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_TRANSIENT_H */
