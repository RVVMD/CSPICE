#ifndef MNA_ELEMENTS_TRANSFORMER_H
#define MNA_ELEMENTS_TRANSFORMER_H

#include "../include/types.h"
#include "nonlinear/nonlinear.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double turns_ratio;
    double Lm;
    double i_mag;
    double phi;
    double prev_phi;
    double Gm_eq;
    double I_eq;
    int branch_current_index;
    double saturation_current;
    double saturation_factor;
    bool use_saturation;
    bool is_ac_analysis;
} TransformerData;

MNAStatus mna_add_ideal_transformer(MNASolver* solver,
                                     int node_p1, int node_p2,
                                     int node_s1, int node_s2,
                                     double turns_ratio,
                                     ComponentHandle* handle);

MNAStatus mna_add_voltage_transformer(MNASolver* solver,
                                       int node_p1, int node_p2,
                                       int node_s1, int node_s2,
                                       double turns_ratio,
                                       double Lm,
                                       ComponentHandle* handle);

MNAStatus mna_add_transformer_with_saturation(MNASolver* solver,
                                               int node_p1, int node_p2,
                                               int node_s1, int node_s2,
                                               double turns_ratio,
                                               double Lm,
                                               double I_sat,
                                               double sat_factor,
                                               ComponentHandle* handle);

MNAStatus mna_add_current_transformer(MNASolver* solver,
                                       int node_p1, int node_p2,
                                       int node_s1, int node_s2,
                                       double turns_ratio,
                                       double burden,
                                       ComponentHandle* handle);

double mna_get_transformer_magnetizing_current(MNASolver* solver,
                                                ComponentHandle handle);

double mna_get_transformer_flux(MNASolver* solver,
                                 ComponentHandle handle);

double mna_get_transformer_primary_voltage(MNASolver* solver,
                                            ComponentHandle handle);

double mna_get_transformer_secondary_voltage(MNASolver* solver,
                                              ComponentHandle handle);

int mna_get_transformer_branch_index(MNASolver* solver,
                                      ComponentHandle handle);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_TRANSFORMER_H */
