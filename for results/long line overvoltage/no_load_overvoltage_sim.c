#include "mna_solver_v2.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define LINE_LENGTH_KM 3000.0
#define NUM_SEGMENTS 1000
#define VOLTAGE_AMP 220000.0
#define FREQUENCY 50.0
#define SOURCE_FREQ (2 * M_PI * FREQUENCY)
#define R_PER_KM 0.05
#define L_PER_KM 0.001
#define C_PER_KM 10e-9
#define LOAD_RESISTANCE 220.0

int main() {
    MNASolver solver;
    if (mna_init(&solver) != MNA_SUCCESS) return 1;

    // Создание узлов
    int total_nodes = NUM_SEGMENTS * 2 + 1;
    int* node_ids = (int*)malloc(total_nodes * sizeof(int));
    if (!node_ids) return 1;

    node_ids[0] = 0;
    for (int i = 1; i < total_nodes; i++) {
        node_ids[i] = mna_create_node(&solver);
    }

    // Источник напряжения
    int source_node = node_ids[1];
    ComponentHandle source_handle;
    mna_add_voltage_source(&solver, source_node, 0, 0, &source_handle);
    mna_set_ac_source(&solver, source_handle, VOLTAGE_AMP, 0.0);

    // Параметры сегментов
    double segment_length = LINE_LENGTH_KM / NUM_SEGMENTS;
    double R_segment = R_PER_KM * segment_length;
    double L_segment = L_PER_KM * segment_length;
    double C_segment = C_PER_KM * segment_length;

    // Создание Т-образных сегментов
    for (int i = 0; i < NUM_SEGMENTS; i++) {
        int node_start = (i == 0) ? 0 : node_ids[i * 2];
        int node_mid = node_ids[i * 2 + 1];
        int node_end = node_ids[i * 2 + 2];

        ComponentHandle r1_handle, r2_handle, l1_handle, l2_handle;
        mna_add_resistor(&solver, node_start, node_mid, R_segment / 2, &r1_handle);
        mna_add_resistor(&solver, node_mid, node_end, R_segment / 2, &r2_handle);
        mna_add_inductor(&solver, node_start, node_mid, L_segment / 2, &l1_handle);
        mna_add_inductor(&solver, node_mid, node_end, L_segment / 2, &l2_handle);

        ComponentHandle c_handle;
        mna_add_capacitor(&solver, node_mid, 0, C_segment, &c_handle);
    }

    // Нагрузка с ключом
    int end_node = node_ids[NUM_SEGMENTS * 2];
    ComponentHandle switch_handle;
    mna_add_switch(&solver, end_node, 0, LOAD_RESISTANCE, &switch_handle);
    mna_set_switch_state(&solver, switch_handle, 1);

    // Инициализация и настройка_transient анализ
    mna_init_transient(&solver);
    mna_solve_ac(&solver, FREQUENCY);

    int matrix_size = mna_active_size(&solver);
    for (int i = 0; i < matrix_size; i++) {
        solver.x[i] = creal(solver.ac_solution[i]);
    }

    // Параметры симуляции
    double dt = 1e-4;
    double total_time = 0.15;
    double switch_time = 0.1;
    int steps = (int)(total_time / dt);
    int switch_step = (int)(switch_time / dt);

    // Период сохранения кадров для анимации (каждые 100 шагов = 0.5 мс)
    int frame_interval = 5;
    int num_frames = steps / frame_interval + 1;

    // Выделение памяти для кадров анимации
    double* frame_times = (double*)malloc(num_frames * sizeof(double));
    double** voltage_frames = (double**)malloc(num_frames * sizeof(double*));

    for (int f = 0; f < num_frames; f++) {
        voltage_frames[f] = (double*)malloc(NUM_SEGMENTS * sizeof(double));
    }

    // Основной цикл симуляции
    int frame_count = 0;
    for (int step = 0; step < steps; step++) {
        double t = step * dt;

        // Обновление напряжения источника
        double voltage_value = VOLTAGE_RMS * sqrt(2) * sin(SOURCE_FREQ * t);
        solver.components[source_handle].value = voltage_value;

        // Отключение нагрузки
        if (step == switch_step) {
            mna_set_switch_state(&solver, switch_handle, 0);
        }

        // Шаг симуляции
        mna_solve_transient_step(&solver, dt);

        // Сохранение кадра для анимации
        if (step % frame_interval == 0) {
            frame_times[frame_count] = t;
            for (int i = 0; i < NUM_SEGMENTS; i++) {
                int node_idx = i * 2 + 1; // Берем промежуточные узлы сегментов
                voltage_frames[frame_count][i] = mna_get_node_voltage(&solver, node_ids[node_idx]);
            }
            frame_count++;
        }

        // Прогресс-бар (только проценты)
        if (step % (steps / 100) == 0) {
            int percent = (step * 100) / steps;
            printf("\r%d%%", percent);
            fflush(stdout);
        }
    }
    printf("\r100%%\n");

    // Сохранение данных для анимации
    FILE* anim_file = fopen("animation_data.csv", "w");
    if (anim_file) {
        fprintf(anim_file, "frame,time");
        for (int i = 0; i < NUM_SEGMENTS; i++) {
            fprintf(anim_file, ",node_%d", i);
        }
        fprintf(anim_file, "\n");

        for (int f = 0; f < frame_count; f++) {
            fprintf(anim_file, "%d,%.9f", f, frame_times[f]);
            for (int i = 0; i < NUM_SEGMENTS; i++) {
                fprintf(anim_file, ",%.2f", voltage_frames[f][i]);
            }
            fprintf(anim_file, "\n");
        }
        fclose(anim_file);
    }

    // Очистка
    free(node_ids);
    for (int f = 0; f < num_frames; f++) {
        free(voltage_frames[f]);
    }
    free(voltage_frames);
    free(frame_times);
    mna_destroy(&solver);

    return 0;
}
