#ifndef MNA_SOLVER_DC_H
#define MNA_SOLVER_DC_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_solve_dc(MNASolver* solver);
void mna_solve_initial_conditions(MNASolver* solver);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_DC_H */
