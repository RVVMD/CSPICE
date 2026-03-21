#ifndef MNA_SOLVER_H
#define MNA_SOLVER_H

#include "types.h"
#include "matrix.h"
#include "solver/core.h"
#include "elements/passive.h"
#include "elements/sources.h"
#include "elements/transformer.h"
#include "elements/nonlinear/nonlinear.h"
#include "elements/nonlinear/magnetization.h"
#include "src/solver/dc.h"
#include "src/solver/ac.h"
#include "src/solver/transient.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MNA_SOLVER_VERSION_MAJOR 2
#define MNA_SOLVER_VERSION_MINOR 5
#define MNA_SOLVER_VERSION "2.5.0"

#ifdef __cplusplus
}
#endif

#endif /* MNA_SOLVER_H */
