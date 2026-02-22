#include "mna_solver_v2_5.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int* nodes;
    ComponentHandle* series_handles;
    int segments;
} TransmissionLine;

TransmissionLine add_transmission_line(MNASolver* solver, int input_node, int output_node,
                                      double R_total, double L_total, double C_total, int segments) {
    TransmissionLine line = {0};
    if (segments <= 0) return line;
    line.nodes = (int*)malloc((segments + 1) * sizeof(int));
    line.series_handles = (ComponentHandle*)malloc(segments * sizeof(ComponentHandle));
    if (!line.nodes || !line.series_handles) {
        free(line.nodes);
        free(line.series_handles);
        line.nodes = NULL;
        return line;
    }
    line.segments = segments;
    line.nodes[0] = input_node;
    double R_seg = R_total / segments;
    double L_seg = L_total / segments;
    double C_seg = C_total / segments;
    int current_node = input_node;
    for (int i = 0; i < segments; i++) {
        int next_node = (i == segments - 1) ? output_node : mna_create_node(solver);
        int mid_node = mna_create_node(solver);
        ComponentHandle r_handle, l_handle;
        mna_add_resistor(solver, current_node, mid_node, R_seg, &r_handle);
        mna_add_inductor(solver, mid_node, next_node, L_seg, &l_handle);
        mna_add_capacitor(solver, next_node, 0, C_seg, NULL);
        line.nodes[i + 1] = next_node;
        line.series_handles[i] = l_handle;
        current_node = next_node;
    }
    return line;
}

int main() {
    MNASolver solver;
    double Z0_line1 = 350.0;
    double v_line1 = 3.0e8;
    double L_per_m_line1 = Z0_line1 / v_line1;
    double C_per_m_line1 = 1.0 / (Z0_line1 * v_line1);
    double length_line1 = 100e3;
    double R_total_line1 = 0.01 * length_line1/1000;
    double L_total_line1 = L_per_m_line1 * length_line1;
    double C_total_line1 = C_per_m_line1 * length_line1;
    double Z0_line23 = 75.0;
    double v_line23 = 1.5e8;
    double L_per_m_line23 = Z0_line23 / v_line23;
    double C_per_m_line23 = 1.0 / (Z0_line23 * v_line23);
    double length_line2 = 60e3;
    double R_total_line2 = 0.01 * length_line2/1000;
    double L_total_line2 = L_per_m_line23 * length_line2;
    double C_total_line2 = C_per_m_line23 * length_line2;
    double length_line3 = 30e3;
    double R_total_line3 = 0.01 * length_line3/1000;
    double L_total_line3 = L_per_m_line23 * length_line3;
    double C_total_line3 = C_per_m_line23 * length_line3;
    int source_node = mna_create_node(&solver);
    int line1_end = mna_create_node(&solver);
    int line2_start = mna_create_node(&solver);
    int line2_end = mna_create_node(&solver);
    int inductor_node = mna_create_node(&solver);
    int line3_end = mna_create_node(&solver);
    mna_add_voltage_source(&solver, source_node, 0, 35000.0, NULL);
    ComponentHandle switch_handle;
    mna_add_switch(&solver, line1_end, line2_start, 0.001, &switch_handle);
    mna_set_switch_state(&solver, switch_handle, 0);
    TransmissionLine line1 = add_transmission_line(&solver, source_node, line1_end,
                                                  R_total_line1, L_total_line1, C_total_line1, 400);
    TransmissionLine line2 = add_transmission_line(&solver, line2_start, line2_end,
                                                  R_total_line2, L_total_line2, C_total_line2, 240);
    ComponentHandle load_ind_handle;
    mna_add_inductor(&solver, line2_end, inductor_node, 0.1, &load_ind_handle);
    mna_add_resistor(&solver, inductor_node, 0, 100.0, NULL);
    TransmissionLine line3 = add_transmission_line(&solver, inductor_node, line3_end,
                                                  R_total_line3, L_total_line3, C_total_line3, 120);
    mna_add_capacitor(&solver, line3_end, 0, 1e-6, NULL);
    MNAStatus status = mna_solve_dc(&solver);
    if (status != MNA_SUCCESS) {
        free(line1.nodes); free(line1.series_handles);
        free(line2.nodes); free(line2.series_handles);
        free(line3.nodes); free(line3.series_handles);
        mna_destroy(&solver);
        return 1;
    }
    solver.transient_initialized = 1;
    solver.time = 0.0;
    for (int i = 0; i < solver.num_components; i++) {
        Component* comp = &solver.components[i];
        int n1 = comp->node1;
        int n2 = comp->node2;
        double v1 = (n1 > 0) ? solver.x[n1-1] : 0.0;
        double v2 = (n2 > 0) ? solver.x[n2-1] : 0.0;
        double voltage = v1 - v2;
        switch (comp->type) {
            case MNA_CAPACITOR:
                comp->last_voltage = voltage;
                comp->last_current = 0.0;
                break;
            case MNA_INDUCTOR:
                comp->last_current = 0.0;
                break;
        }
    }
    mna_set_switch_state(&solver, switch_handle, 1);
    FILE* csv_v = fopen("voltage_results.csv", "w");
    FILE* csv_i = fopen("current_results.csv", "w");
    if (!csv_v || !csv_i) {
        fprintf(stderr, "Failed to create CSV files\n");
        free(line1.nodes); free(line1.series_handles);
        free(line2.nodes); free(line2.series_handles);
        free(line3.nodes); free(line3.series_handles);
        mna_destroy(&solver);
        return 1;
    }
    fprintf(csv_v, "time");
    fprintf(csv_i, "time");
    for (int i = 0; i <= 400; i++) fprintf(csv_v, ",line1_%.3f", (double)i/400);
    for (int i = 0; i <= 240; i++) fprintf(csv_v, ",line2_%.3f", (double)i/240);
    for (int i = 0; i <= 120; i++) fprintf(csv_v, ",line3_%.3f", (double)i/120);
    for (int i = 0; i < 400; i++) fprintf(csv_i, ",line1_I_%.3f", (double)i/400);
    for (int i = 0; i < 240; i++) fprintf(csv_i, ",line2_I_%.3f", (double)i/240);
    for (int i = 0; i < 120; i++) fprintf(csv_i, ",line3_I_%.3f", (double)i/120);
    fprintf(csv_v, "\n");
    fprintf(csv_i, "\n");
    fprintf(csv_v, "0.0");
    fprintf(csv_i, "0.0");
    for (int i = 0; i <= 400; i++) fprintf(csv_v, ",%.6f", mna_get_node_voltage(&solver, line1.nodes[i]));
    for (int i = 0; i <= 240; i++) fprintf(csv_v, ",%.6f", mna_get_node_voltage(&solver, line2.nodes[i]));
    for (int i = 0; i <= 120; i++) fprintf(csv_v, ",%.6f", mna_get_node_voltage(&solver, line3.nodes[i]));
    for (int i = 0; i < 400; i++) fprintf(csv_i, ",%.6f", mna_get_component_current(&solver, line1.series_handles[i]));
    for (int i = 0; i < 240; i++) fprintf(csv_i, ",%.6f", mna_get_component_current(&solver, line2.series_handles[i]));
    for (int i = 0; i < 120; i++) fprintf(csv_i, ",%.6f", mna_get_component_current(&solver, line3.series_handles[i]));
    fprintf(csv_v, "\n");
    fprintf(csv_i, "\n");
    double t_end = 1e-3;
    double dt = 1e-6;
    int steps = (int)(t_end / dt);
    for (int step = 1; step <= steps; step++) {
        double t = step * dt;
        status = mna_solve_transient_step(&solver, dt);
        if (status != MNA_SUCCESS) break;
        fprintf(csv_v, "%.9f", t);
        fprintf(csv_i, "%.9f", t);
        for (int i = 0; i <= 400; i++) fprintf(csv_v, ",%.6f", mna_get_node_voltage(&solver, line1.nodes[i]));
        for (int i = 0; i <= 240; i++) fprintf(csv_v, ",%.6f", mna_get_node_voltage(&solver, line2.nodes[i]));
        for (int i = 0; i <= 120; i++) fprintf(csv_v, ",%.6f", mna_get_node_voltage(&solver, line3.nodes[i]));
        for (int i = 0; i < 400; i++) fprintf(csv_i, ",%.6f", mna_get_component_current(&solver, line1.series_handles[i]));
        for (int i = 0; i < 240; i++) fprintf(csv_i, ",%.6f", mna_get_component_current(&solver, line2.series_handles[i]));
        for (int i = 0; i < 120; i++) fprintf(csv_i, ",%.6f", mna_get_component_current(&solver, line3.series_handles[i]));
        fprintf(csv_v, "\n");
        fprintf(csv_i, "\n");
        if (step % (steps/10) == 0) {
            printf("Progress: %.0f%% (t=%.6fs)\n", (double)step/steps*100, t);
        }
    }
    fclose(csv_v);
    fclose(csv_i);
    free(line1.nodes); free(line1.series_handles);
    free(line2.nodes); free(line2.series_handles);
    free(line3.nodes); free(line3.series_handles);
    mna_destroy(&solver);
    return 0;
}
