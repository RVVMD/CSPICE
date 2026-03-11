#ifndef MNA_SOLVER_AC_H
#define MNA_SOLVER_AC_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_solve_ac(MNASolver* solver, double frequency);

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_AC_H */
