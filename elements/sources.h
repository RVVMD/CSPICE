#ifndef MNA_ELEMENTS_SOURCES_H
#define MNA_ELEMENTS_SOURCES_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_add_voltage_source(MNASolver* solver, int node1, int node2,
                                 double value, ComponentHandle* handle);

MNAStatus mna_add_current_source(MNASolver* solver, int node1, int node2,
                                 double value, ComponentHandle* handle);

MNAStatus mna_set_ac_source(MNASolver* solver, ComponentHandle handle,
                            double magnitude, double phase);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_SOURCES_H */
