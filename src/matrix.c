#include "matrix.h"
#include "../include/types.h"
#include <stdlib.h>
#include <string.h>

bool mna_resize_components(MNASolver* s, int required) {
    if (required <= s->cap_components) return true;
    int new_cap = s->cap_components ? s->cap_components : MNA_INIT_CAPACITY;
    while (new_cap < required) {
        new_cap = (int)(new_cap * 1.5) + 1;
    }
    Component* nxt = (Component*)realloc(s->components, (size_t)new_cap * sizeof(Component));
    if (!nxt) return false;
    if (new_cap > s->cap_components) {
        size_t added = (size_t)(new_cap - s->cap_components);
        memset(nxt + s->cap_components, 0, added * sizeof(Component));
    }
    s->components = nxt;
    s->cap_components = new_cap;
    return true;
}

bool mna_resize_matrix(MNASolver* s, int req_size) {
    if (req_size <= s->matrix_cap_size) return true;
    int new_size = s->matrix_cap_size ? s->matrix_cap_size : MNA_INIT_CAPACITY;
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
    int old_size = s->matrix_cap_size;
    if (s->A && old_size > 0) {
        for (int i = 0; i < old_size; ++i) {
            memcpy(A2 + (size_t)i * new_size, s->A + (size_t)i * old_size, (size_t)old_size * sizeof(double));
        }
        memcpy(b2, s->b, (size_t)old_size * sizeof(double));
        memcpy(x2, s->x, (size_t)old_size * sizeof(double));
        memcpy(ac2, s->ac_solution, (size_t)old_size * sizeof(double complex));
    }
    free(s->A); free(s->b); free(s->x); free(s->ac_solution);
    s->A = A2;
    s->b = b2;
    s->x = x2;
    s->ac_solution = ac2;
    s->matrix_cap_size = new_size;
    return true;
}

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
    int n1 = vs->node1;
    int n2 = vs->node2;
    int v_index = solver->max_node_index + source_idx;

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

void mna_reset_system(MNASolver* solver) {
    if (!solver) return;

    int size = mna_active_size(solver);
    if (size > solver->matrix_cap_size) size = solver->matrix_cap_size;

    for (int i = 0; i < size; i++) {
        memset(&MAT(solver, i, 0), 0, (size_t)size * sizeof(double));
    }
    memset(solver->b, 0, (size_t)size * sizeof(double));
}

MNAStatus mna_solve_linear_system(MNASolver* solver, int size) {
    if (!solver || size <= 0) return MNA_INVALID_PARAMETER;

    for (int pivot = 0; pivot < size; pivot++) {
        int max_row = pivot;
        double max_val = fabs(MAT(solver, pivot, pivot));

        for (int i = pivot + 1; i < size; i++) {
            double val = fabs(MAT(solver, i, pivot));
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        if (max_val < MNA_MIN_CONDUCTANCE) {
            return MNA_MATRIX_SINGULAR;
        }

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

        for (int i = pivot + 1; i < size; i++) {
            double factor = MAT(solver, i, pivot) / MAT(solver, pivot, pivot);
            if (factor == 0.0) continue;

            for (int j = pivot + 1; j < size; j++) {
                MAT(solver, i, j) -= factor * MAT(solver, pivot, j);
            }
            solver->b[i] -= factor * solver->b[pivot];
            MAT(solver, i, pivot) = 0.0;
        }
    }

    for (int i = size - 1; i >= 0; i--) {
        solver->x[i] = solver->b[i];
        for (int j = i + 1; j < size; j++) {
            solver->x[i] -= MAT(solver, i, j) * solver->x[j];
        }
        solver->x[i] /= MAT(solver, i, i);
    }

    return MNA_SUCCESS;
}
