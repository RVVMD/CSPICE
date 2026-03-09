#include "ac.h"
#include "../include/matrix.h"
#include "../elements/nonlinear/nonlinear.h"
#include <stdlib.h>

/* ============================================================================
 * AC Analysis Implementation
 * ============================================================================ */

MNAStatus mna_solve_ac(MNASolver* solver, double frequency) {
    if (!solver) return MNA_INVALID_HANDLE;
    
    double omega = TWO_PI * frequency;
    int matrix_size = mna_active_size(solver);
    
    if (matrix_size > solver->matrix_cap_size) return MNA_INVALID_PARAMETER;
    
    /* Allocate complex matrices */
    double complex* A_complex = (double complex*)calloc((size_t)matrix_size * 
                                         (size_t)matrix_size, sizeof(double complex));
    double complex* b_complex = (double complex*)calloc((size_t)matrix_size, 
                                         sizeof(double complex));
    double complex* x_complex = (double complex*)calloc((size_t)matrix_size, 
                                         sizeof(double complex));
    
    if (!A_complex || !b_complex || !x_complex) {
        free(A_complex); free(b_complex); free(x_complex);
        return MNA_INSUFFICIENT_MEMORY;
    }
    
    /* Stamp components */
    int source_count = 0;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;
        double complex admittance = 0.0;
        
        switch (comp->type) {
            case MNA_RESISTOR:
                admittance = 1.0 / comp->value;
                break;
                
            case MNA_CAPACITOR:
                admittance = I * omega * comp->value;
                break;
                
            case MNA_INDUCTOR:
                admittance = (omega == 0) ? 1.0 / MNA_MIN_CONDUCTANCE : 
                             1.0 / (I * omega * comp->value);
                break;
                
            case MNA_SOURCE: {
                double complex source_val = comp->ac_magnitude * 
                    (cos(comp->ac_phase) + I * sin(comp->ac_phase));
                
                if (comp->source_type == SOURCE_VOLTAGE) {
                    int v_index = solver->max_node_index + source_count;
                    
                    if (n1 > 0) {
                        A_complex[(n1-1) * matrix_size + v_index] = 1.0;
                        A_complex[v_index * matrix_size + (n1-1)] = 1.0;
                    }
                    if (n2 > 0) {
                        A_complex[(n2-1) * matrix_size + v_index] = -1.0;
                        A_complex[v_index * matrix_size + (n2-1)] = -1.0;
                    }
                    b_complex[v_index] = source_val;
                    source_count++;
                } else {
                    if (n1 > 0) b_complex[n1-1] -= source_val;
                    if (n2 > 0) b_complex[n2-1] += source_val;
                }
                continue;
            }
            
            case MNA_CUSTOM_NONLINEAR: {
                if (comp->nonlinear_type == NONLINEAR_RESISTOR) {
                    admittance = comp->last_conductance;
                } else if (comp->nonlinear_type == NONLINEAR_CAPACITOR) {
                    ComponentState state = { 
                        .voltage = comp->last_voltage, 
                        .current = comp->last_current,
                        .charge = comp->last_charge, 
                        .flux = comp->last_flux, 
                        .dt = solver->dt 
                    };
                    double q0, C0;
                    comp->nonlinear_func(&state, comp->user_data, &q0, &C0);
                    admittance = I * omega * C0;
                } else if (comp->nonlinear_type == NONLINEAR_INDUCTOR) {
                    ComponentState state = { 
                        .voltage = comp->last_voltage, 
                        .current = comp->last_current,
                        .charge = comp->last_charge, 
                        .flux = comp->last_flux, 
                        .dt = solver->dt 
                    };
                    double phi0, L0;
                    comp->nonlinear_func(&state, comp->user_data, &phi0, &L0);
                    admittance = 1.0 / (I * omega * L0);
                }
                break;
            }
            
            case MNA_CUSTOM_NPOLE: {
                NPoleData* npole = comp->data.npole.npole_data;
                if (npole) {
                    npole->stamp_func(solver, npole->nodes, npole->num_nodes, 
                                      npole->user_data, 0.0, 0.0);
                }
                continue;
            }
            
            case MNA_SWITCH:
                admittance = comp->state ? 1.0 / comp->value : MNA_MIN_CONDUCTANCE;
                break;
                
            default:
                continue;
        }
        
        /* Stamp admittance */
        if (n1 > 0) A_complex[(n1-1) * matrix_size + (n1-1)] += admittance;
        if (n2 > 0) A_complex[(n2-1) * matrix_size + (n2-1)] += admittance;
        if (n1 > 0 && n2 > 0) {
            A_complex[(n1-1) * matrix_size + (n2-1)] -= admittance;
            A_complex[(n2-1) * matrix_size + (n1-1)] -= admittance;
        }
    }
    
    /* Solve complex linear system using Gaussian elimination */
    for (int pivot = 0; pivot < matrix_size; pivot++) {
        /* Partial pivoting */
        int max_row = pivot;
        double max_mag = cabs(A_complex[pivot * matrix_size + pivot]);
        
        for (int i = pivot + 1; i < matrix_size; i++) {
            double mag = cabs(A_complex[i * matrix_size + pivot]);
            if (mag > max_mag) {
                max_mag = mag;
                max_row = i;
            }
        }
        
        if (max_mag < MNA_MIN_CONDUCTANCE) {
            free(A_complex); free(b_complex); free(x_complex);
            return MNA_MATRIX_SINGULAR;
        }
        
        /* Swap rows */
        if (max_row != pivot) {
            for (int j = pivot; j < matrix_size; j++) {
                double complex temp = A_complex[pivot * matrix_size + j];
                A_complex[pivot * matrix_size + j] = A_complex[max_row * matrix_size + j];
                A_complex[max_row * matrix_size + j] = temp;
            }
            double complex tempb = b_complex[pivot];
            b_complex[pivot] = b_complex[max_row];
            b_complex[max_row] = tempb;
        }
        
        /* Elimination */
        for (int i = pivot + 1; i < matrix_size; i++) {
            double complex factor = A_complex[i * matrix_size + pivot] / 
                                    A_complex[pivot * matrix_size + pivot];
            if (factor == 0.0) continue;
            
            for (int j = pivot + 1; j < matrix_size; j++) {
                A_complex[i * matrix_size + j] -= factor * A_complex[pivot * matrix_size + j];
            }
            b_complex[i] -= factor * b_complex[pivot];
            A_complex[i * matrix_size + pivot] = 0.0;
        }
    }
    
    /* Back substitution */
    for (int i = matrix_size - 1; i >= 0; i--) {
        x_complex[i] = b_complex[i];
        for (int j = i + 1; j < matrix_size; j++) {
            x_complex[i] -= A_complex[i * matrix_size + j] * x_complex[j];
        }
        x_complex[i] /= A_complex[i * matrix_size + i];
    }
    
    /* Copy solution */
    for (int i = 0; i < matrix_size; i++) {
        solver->ac_solution[i] = x_complex[i];
    }
    
    free(A_complex);
    free(b_complex);
    free(x_complex);
    
    return MNA_SUCCESS;
}
