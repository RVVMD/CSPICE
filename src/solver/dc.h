#ifndef MNA_SOLVER_DC_H
#define MNA_SOLVER_DC_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * DC Analysis
 * ============================================================================ */

/**
 * @brief Perform DC operating point analysis
 * 
 * For linear circuits, solves the MNA system directly.
 * For nonlinear circuits, uses Newton-Raphson iteration with
 * adaptive source stepping for convergence.
 * 
 * @param solver The MNA solver
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_solve_dc(MNASolver* solver);

/* ============================================================================
 * Initial Conditions for Transient Analysis
 * ============================================================================ */

/**
 * @brief Solve initial conditions for transient analysis
 * 
 * Computes the initial state of energy storage elements
 * (capacitors and inductors) before transient analysis begins.
 * 
 * @param solver The MNA solver
 */
void mna_solve_initial_conditions(MNASolver* solver);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_DC_H */
