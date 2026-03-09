#include "matrix.h"
#include <stdlib.h>

/* ============================================================================
 * Matrix Memory Management
 * ============================================================================ */

bool mna_resize_components(MNASolver* solver, int required) {
    if (required <= solver->cap_components) return true;
    
    int new_cap = solver->cap_components ? solver->cap_components : MNA_INIT_CAPACITY;
    while (new_cap < required) {
        new_cap = (int)(new_cap * 1.5) + 1;
    }
    
    Component* nxt = (Component*)realloc(solver->components, (size_t)new_cap * sizeof(Component));
    if (!nxt) return false;
    
    if (new_cap > solver->cap_components) {
        size_t added = (size_t)(new_cap - solver->cap_components);
        memset(nxt + solver->cap_components, 0, added * sizeof(Component));
    }
    
    solver->components = nxt;
    solver->cap_components = new_cap;
    return true;
}

bool mna_resize_matrix(MNASolver* solver, int req_size) {
    if (req_size <= solver->matrix_cap_size) return true;
    
    int new_size = solver->matrix_cap_size ? solver->matrix_cap_size : MNA_INIT_CAPACITY;
    while (new_size < req_size) {
        new_size = (int)(new_size * 1.5) + 1;
    }
    
    size_t new_elems = (size_t)new_size * (size_t)new_size;
    size_t vec_elems = (size_t)new_size;
    
    double* A2 = (double*)calloc(new_elems, sizeof(double));
    double* b2 = (double*)calloc(vec_elems, sizeof(double));
    double* x2 = (double*)calloc(vec_elems, sizeof(double));
    double complex* ac2 = (double complex*)calloc(vec_elems, sizeof(double complex));
    
    if (!A2 || !b2 || !x2 || !ac2) {
        free(A2); free(b2); free(x2); free(ac2);
        return false;
    }
    
    int old_size = solver->matrix_cap_size;
    if (solver->A && old_size > 0) {
        for (int i = 0; i < old_size; ++i) {
            memcpy(A2 + (size_t)i * new_size, 
                   solver->A + (size_t)i * old_size, 
                   (size_t)old_size * sizeof(double));
        }
        memcpy(b2, solver->b, (size_t)old_size * sizeof(double));
        memcpy(x2, solver->x, (size_t)old_size * sizeof(double));
        memcpy(ac2, solver->ac_solution, (size_t)old_size * sizeof(double complex));
    }
    
    free(solver->A); free(solver->b); free(solver->x); free(solver->ac_solution);
    
    solver->A = A2;
    solver->b = b2;
    solver->x = x2;
    solver->ac_solution = ac2;
    solver->matrix_cap_size = new_size;
    return true;
}

/* ============================================================================
 * Matrix Operations
 * ============================================================================ */

void mna_reset_system(MNASolver* solver) {
    int active = mna_active_size(solver);
    for (int i = 0; i < active; ++i) {
        memset(&MAT(solver, i, 0), 0, (size_t)active * sizeof(double));
    }
    memset(solver->b, 0, (size_t)active * sizeof(double));
    memset(solver->x, 0, (size_t)active * sizeof(double));
}

MNAStatus mna_solve_linear_system(MNASolver* solver, int size) {
    const double pivot_threshold = 1e-12;
    
    for (int pivot = 0; pivot < size; pivot++) {
        /* Partial pivoting */
        int max_row = pivot;
        double max_val = fabs(MAT(solver, pivot, pivot));
        
        for (int i = pivot + 1; i < size; i++) {
            double abs_val = fabs(MAT(solver, i, pivot));
            if (abs_val > max_val) {
                max_val = abs_val;
                max_row = i;
            }
        }
        
        if (max_val < pivot_threshold) {
            return MNA_MATRIX_SINGULAR;
        }
        
        /* Swap rows if needed */
        if (max_row != pivot) {
            for (int j = pivot; j < size; j++) {
                double temp = MAT(solver, pivot, j);
                MAT(solver, pivot, j) = MAT(solver, max_row, j);
                MAT(solver, max_row, j) = temp;
            }
            double tempb = solver->b[pivot];
            solver->b[pivot] = solver->b[max_row];
            solver->b[max_row] = tempb;
        }
        
        /* Elimination */
        const double pivot_inv = 1.0 / MAT(solver, pivot, pivot);
        for (int i = pivot + 1; i < size; i++) {
            double factor = MAT(solver, i, pivot) * pivot_inv;
            if (factor == 0.0) continue;
            
            for (int j = pivot + 1; j < size; j++) {
                MAT(solver, i, j) -= factor * MAT(solver, pivot, j);
            }
            solver->b[i] -= factor * solver->b[pivot];
            MAT(solver, i, pivot) = 0.0;
        }
    }
    
    /* Back substitution */
    for (int i = size - 1; i >= 0; i--) {
        double sum = solver->b[i];
        for (int j = i + 1; j < size; j++) {
            sum -= MAT(solver, i, j) * solver->x[j];
        }
        solver->x[i] = sum / MAT(solver, i, i);
    }
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Matrix Stamping Functions
 * ============================================================================ */

void mna_stamp_conductance(MNASolver* solver, int node1, int node2, double g) {
    if (node1 > 0) MAT(solver, node1-1, node1-1) += g;
    if (node2 > 0) MAT(solver, node2-1, node2-1) += g;
    if (node1 > 0 && node2 > 0) {
        MAT(solver, node1-1, node2-1) -= g;
        MAT(solver, node2-1, node1-1) -= g;
    }
}

void mna_stamp_current_source(MNASolver* solver, int node1, int node2, double current_val) {
    if (node1 > 0) solver->b[node1-1] -= current_val;
    if (node2 > 0) solver->b[node2-1] += current_val;
}

void mna_stamp_voltage_source(MNASolver* solver, int comp_index, int source_idx) {
    Component* vs = &solver->components[comp_index];
    int v_index = solver->max_node_index + source_idx;
    int n1 = vs->node1;
    int n2 = vs->node2;
    
    if (n1 > 0) {
        MAT(solver, n1-1, v_index) = 1.0;
        MAT(solver, v_index, n1-1) = 1.0;
    }
    if (n2 > 0) {
        MAT(solver, n2-1, v_index) = -1.0;
        MAT(solver, v_index, n2-1) = -1.0;
    }
    solver->b[v_index] = vs->value;
}

void mna_ensure_ground_paths(MNASolver* solver) {
    int* connected = (int*)calloc((size_t)(solver->max_node_index + 1), sizeof(int));
    if (!connected) return;
    
    connected[0] = 1;
    for (int i = 0; i < solver->num_components; i++) {
        Component* comp = &solver->components[i];
        if (comp->node1 == 0 && comp->node2 > 0) connected[comp->node2] = 1;
        else if (comp->node2 == 0 && comp->node1 > 0) connected[comp->node1] = 1;
    }
    
    for (int node = 1; node <= solver->max_node_index; node++) {
        if (!connected[node]) {
            mna_stamp_conductance(solver, node, 0, MNA_GROUND_CONDUCTANCE);
        }
    }
    
    free(connected);
}
