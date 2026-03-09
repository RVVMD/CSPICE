#ifndef MNA_SOLVER_H
#define MNA_SOLVER_H

/**
 * @file mna_solver.h
 * @brief MNA Circuit Solver - Main Public API
 * 
 * A modified nodal analysis circuit solver supporting:
 * - DC operating point analysis
 * - AC small-signal analysis
 * - Transient analysis with TRBDF2 integration
 * - Nonlinear elements with adaptive damping
 * - Custom n-pole elements
 */

#include "types.h"
#include "matrix.h"
#include "solver/core.h"
#include "elements/passive.h"
#include "elements/sources.h"
#include "elements/nonlinear/nonlinear.h"
#include "src/solver/dc.h"
#include "src/solver/ac.h"
#include "src/solver/transient.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Version Information
 * ============================================================================ */
#define MNA_SOLVER_VERSION_MAJOR 2
#define MNA_SOLVER_VERSION_MINOR 5
#define MNA_SOLVER_VERSION "2.5.0"

/* ============================================================================
 * Quick Start Example
 * ============================================================================
 * 
 * MNASolver solver;
 * ComponentHandle r1, c1, v1;
 * 
 * // Initialize
 * mna_init(&solver);
 * 
 * // Create nodes (node 0 is always ground)
 * int n1 = mna_create_node(&solver);
 * int n2 = mna_create_node(&solver);
 * 
 * // Add components
 * mna_add_voltage_source(&solver, 0, n1, 5.0, &v1);
 * mna_add_resistor(&solver, n1, n2, 1000.0, &r1);
 * mna_add_capacitor(&solver, n2, 0, 1e-6, &c1);
 * 
 * // DC analysis
 * mna_solve_dc(&solver);
 * printf("V(n2) = %f V\n", mna_get_node_voltage(&solver, n2));
 * 
 * // AC analysis
 * mna_set_ac_source(&solver, v1, 1.0, 0.0);
 * mna_solve_ac(&solver, 1000.0);
 * 
 * // Transient analysis
 * mna_init_transient(&solver);
 * for (int i = 0; i < 1000; i++) {
 *     mna_solve_transient_step(&solver, 1e-6);
 *     printf("t=%f V(n2)=%f\n", mna_get_time(&solver), 
 *            mna_get_node_voltage(&solver, n2));
 * }
 * 
 * // Cleanup
 * mna_destroy(&solver);
 * ============================================================================ */

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_H */
