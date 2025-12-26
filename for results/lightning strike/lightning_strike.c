#include "mna_solver_v2.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// ПАРАМЕТРЫ СИМУЛЯЦИИ
const int num_sections = 1000;
const double line_length = 100000.0;  // 100 км
const double L0 = 0.98e-6;    // Гн/м
const double C0 = 8.6e-12;    // Ф/м
const double R0 = 0.08;       // Ом/м
const double I_peak = 30000.0; // 30 кА
const double total_time = 300e-6; // 200 мкс
const int strike_pos = num_sections / 2.5;
const int save_interval = 10;
const double termination_resistance = 338.0; // Волновое сопротивление

double lightning_current(double t) {
    double tau1 = 8.0e-6, tau2 = 20.0e-6;
    return I_peak * (exp(-t/tau1) - exp(-t/tau2));
}

int main() {
    double section_length = line_length / num_sections;
    double L_section = L0 * section_length;
    double C_section = C0 * section_length;
    double R_section = R0 * section_length;
    double wave_velocity = 1.0 / sqrt(L0 * C0);
    double dt = section_length / (2.0 * wave_velocity);
    int num_time_steps = (int)(total_time / dt) + 1;

    MNASolver solver;
    mna_init_sized(&solver, num_sections + 1, 1, num_sections * 3 + 3);

    int* nodes = malloc((num_sections + 1) * sizeof(int));
    for (int i = 0; i <= num_sections; i++)
        nodes[i] = mna_create_node(&solver);

    // Добавление компонентов
    for (int i = 0; i < num_sections; i++) {
        ComponentHandle ind, cap1, cap2;
        mna_add_inductor(&solver, nodes[i], nodes[i+1], L_section, &ind);

        if (i == 0) mna_add_capacitor(&solver, nodes[i], 0, C_section, &cap1);
        else mna_add_capacitor(&solver, nodes[i], 0, C_section/2, &cap1);

        if (i == num_sections - 1) mna_add_capacitor(&solver, nodes[i+1], 0, C_section, &cap2);
        else mna_add_capacitor(&solver, nodes[i+1], 0, C_section/2, &cap2);

        if (R_section > 0) {
            ComponentHandle r;
            mna_add_resistor(&solver, nodes[i], nodes[i+1], R_section, &r);
        }
    }

    // Добавление согласующих сопротивлений на концах линии
    ComponentHandle r_start, r_end;
    mna_add_resistor(&solver, nodes[0], 0, termination_resistance, &r_start);
    mna_add_resistor(&solver, nodes[num_sections], 0, termination_resistance, &r_end);

    ComponentHandle lightning_source;
    mna_add_current_source(&solver, nodes[strike_pos], 0, 0.0, &lightning_source);

    mna_init_transient(&solver);
    double* voltages = malloc((num_sections + 1) * sizeof(double));
    FILE* f = fopen("results.csv", "w");

    // Заголовок
    fprintf(f, "Time(us)");
    for (int i = 0; i <= num_sections; i++)
        fprintf(f, ",V_node_%d(kV)", i);
    fprintf(f, "\n");

    // Основной цикл
    for (int step = 0; step < num_time_steps; step++) {
        double t = step * dt;
        solver.components[lightning_source].value = lightning_current(t);

        mna_solve_transient_step(&solver, dt);

        for (int i = 0; i <= num_sections; i++)
            voltages[i] = mna_get_node_voltage(&solver, nodes[i]);

        // Сохранение результатов
        if (step % save_interval == 0) {
            fprintf(f, "%.3f", t * 1e6);
            for (int i = 0; i <= num_sections; i++)
                fprintf(f, ",%.3f", voltages[i] / 1000.0);
            fprintf(f, "\n");
            fflush(f);
        }
    }

    printf("Результаты сохранены в results.csv\n");
    fclose(f);
    free(voltages);
    free(nodes);
    mna_destroy(&solver);
    return 0;
}
