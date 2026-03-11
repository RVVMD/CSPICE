#ifndef MNA_ELEMENTS_PASSIVE_H
#define MNA_ELEMENTS_PASSIVE_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

MNAStatus mna_add_resistor(MNASolver* solver, int node1, int node2,
                           double value, ComponentHandle* handle);

MNAStatus mna_add_capacitor(MNASolver* solver, int node1, int node2,
                            double value, ComponentHandle* handle);

MNAStatus mna_add_inductor(MNASolver* solver, int node1, int node2,
                           double value, ComponentHandle* handle);

MNAStatus mna_add_switch(MNASolver* solver, int node1, int node2,
                         double value, ComponentHandle* handle);

MNAStatus mna_set_switch_state(MNASolver* solver, ComponentHandle handle, int state);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_PASSIVE_H */
