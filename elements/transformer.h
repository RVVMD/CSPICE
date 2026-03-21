#ifndef MNA_ELEMENTS_TRANSFORMER_H
#define MNA_ELEMENTS_TRANSFORMER_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Add an ideal transformer to the circuit
 * 
 * @param solver MNA solver instance
 * @param primary_p Primary winding positive node
 * @param primary_n Primary winding negative node
 * @param secondary_p Secondary winding positive node
 * @param secondary_n Secondary winding negative node
 * @param turns_ratio Turns ratio n = N_primary/N_secondary (V_primary/V_secondary = n)
 * @param handle Output component handle
 * @return MNAStatus MNA_SUCCESS on success, error code otherwise
 * 
 * The ideal transformer enforces:
 *   V_primary = n * V_secondary
 *   I_primary = -I_secondary / n
 */
MNAStatus mna_add_ideal_transformer(MNASolver* solver,
                                    int primary_p, int primary_n,
                                    int secondary_p, int secondary_n,
                                    double turns_ratio,
                                    ComponentHandle* handle);

/**
 * @brief Set the turns ratio of an existing transformer
 * 
 * @param solver MNA solver instance
 * @param handle Component handle
 * @param turns_ratio New turns ratio (must be > 0)
 * @return MNAStatus MNA_SUCCESS on success, error code otherwise
 */
MNAStatus mna_set_transformer_turns_ratio(MNASolver* solver,
                                          ComponentHandle handle,
                                          double turns_ratio);

/**
 * @brief Get the turns ratio of a transformer
 * 
 * @param solver MNA solver instance
 * @param handle Component handle
 * @param turns_ratio Output turns ratio
 * @return MNAStatus MNA_SUCCESS on success, error code otherwise
 */
MNAStatus mna_get_transformer_turns_ratio(MNASolver* solver,
                                          ComponentHandle handle,
                                          double* turns_ratio);

#ifdef __cplusplus
}
#endif

#endif /* MNA_ELEMENTS_TRANSFORMER_H */
