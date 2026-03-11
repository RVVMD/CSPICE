#ifndef MNA_SOLVER_TRANSIENT_H
#define MNA_SOLVER_TRANSIENT_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

void mna_init_transient(MNASolver* solver);
MNAStatus mna_solve_transient_step(MNASolver* solver, double dt);
double mna_get_time(MNASolver* solver);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_TRANSIENT_H */
