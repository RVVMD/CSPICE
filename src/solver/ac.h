#ifndef MNA_SOLVER_AC_H
#define MNA_SOLVER_AC_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * AC (Small-Signal) Analysis
 * ============================================================================ */

/**
 * @brief Perform AC small-signal analysis at a given frequency
 * 
 * Solves the complex MNA system for linearized circuit around
 * the DC operating point. Must call mna_solve_dc first for
 * circuits with nonlinear elements.
 * 
 * @param solver The MNA solver
 * @param frequency Frequency in Hz
 * @return MNAStatus indicating success or failure
 */
MNAStatus mna_solve_ac(MNASolver* solver, double frequency);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_AC_H */
