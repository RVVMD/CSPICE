/*
 * Transformer energization test - CSV output
 *
 * Circuit:
 *   V1 (AC source) --- Switch --- R_winding --- Transformer --- R_load
 *                                      |
 *                                   R_core (core loss)
 *
 * Realistic 120V, 100VA transformer parameters:
 * - Primary: 120V, 0.83A rated
 * - Secondary: 12V, 8.3A rated
 * - Winding resistance: ~0.5Ω primary (referred)
 * - Leakage inductance: ~5mH primary (referred)
 * - Core loss: ~2W (R_core ≈ 7kΩ)
 *
 * Outputs CSV data for plotting. Demonstrates inrush current
 * when transformer is energized at zero crossing.
 */

#include "mna_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Circuit parameters - 120V, 100VA transformer */
#define FREQ_HZ         60.0        /* 60 Hz mains */
#define V_RMS           120.0       /* Primary RMS voltage */
#define V_PEAK          (V_RMS * sqrt(2.0))
#define TURNS_RATIO     10.0        /* 10:1 step-down (120V -> 12V) */
#define R_LOAD          1.44        /* Secondary: 12V/8.3A = 1.44Ω (100VA load) */
#define R_LOAD_NOLOAD   10000.0     /* No-load: very light load for comparison */

/* Transformer winding parameters (referred to primary) */
#define R_WINDING       0.35        /* Primary winding resistance [ohm] */
#define L_LEAKAGE       0.005       /* Primary leakage inductance [H] */

/* Transformer core parameters - typical silicon steel */
#define B_SAT           1.6         /* Saturation flux density [T] */
#define B_REMANENT      1.2         /* Remanent flux density [T] */
#define MU_R_INITIAL    3000        /* Initial relative permeability */
#define CORE_AREA       4e-4        /* Core area [m^2] = 4 cm^2 */
#define PATH_LENGTH     0.25        /* Magnetic path length [m] = 25 cm */
#define N_PRIMARY       300         /* Primary turns (for ~1.6T at 120V) */

/* Core loss (eddy current + hysteresis) */
#define P_CORE_LOSS     2.0         /* Core loss [W] at rated voltage */
#define R_CORE          1e99  /* ~7200 ohm */

/* Simulation parameters */
#define T_FINAL         2.0         /* 2.0 seconds (120 cycles at 60Hz) */
#define DT              5e-5        /* 50 us time step (20kHz sample rate) */
#define PRINT_INTERVAL  20          /* Print every 1ms (1kHz output rate) */
#define SWITCH_CLOSE_TIME 0.00833   /* Switch closes at 8.33ms (1/2 cycle) */

int main(int argc, char** argv) {
    /* Optional: set initial phase for inrush demonstration */
    double initial_phase = 0.0;  /* 0 = zero crossing (max inrush), PI/2 = peak (min inrush) */
    if (argc > 1) {
        initial_phase = atof(argv[1]);
    }

    MNASolver* solver = (MNASolver*)malloc(sizeof(MNASolver));
    if (!solver) {
        fprintf(stderr, "Failed to allocate solver\n");
        return 1;
    }

    MNAStatus status = mna_init(solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to initialize solver: %d\n", status);
        free(solver);
        return 1;
    }

    /* Create nodes */
    int node_src = mna_create_node(solver);    /* Voltage source */
    int node_sw = mna_create_node(solver);     /* Switch output */
    int node_w = mna_create_node(solver);      /* After winding R */
    int node_l = mna_create_node(solver);      /* After leakage L (transformer primary) */
    int node_sec = mna_create_node(solver);    /* Transformer secondary */

    /* Add AC voltage source */
    ComponentHandle vs_handle;
    status = mna_add_voltage_source(solver, node_src, 0, V_PEAK, &vs_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add voltage source: %d\n", status);
        goto cleanup;
    }

    /* Add switch (initially open) */
    ComponentHandle sw_handle;
    status = mna_add_switch(solver, node_src, node_sw, 1e-12, &sw_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add switch: %d\n", status);
        goto cleanup;
    }
    mna_set_switch_state(solver, sw_handle, 0);  /* Start open */

    /* Add winding resistance */
    ComponentHandle r_w_handle;
    status = mna_add_resistor(solver, node_sw, node_w, R_WINDING, &r_w_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add winding resistance: %d\n", status);
        goto cleanup;
    }

    /* Add leakage inductance */
    ComponentHandle l_lk_handle;
    status = mna_add_inductor(solver, node_w, node_l, L_LEAKAGE, &l_lk_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add leakage inductance: %d\n", status);
        goto cleanup;
    }

    /* Add transformer with saturation */
    ComponentHandle xf_handle;
    status = mna_add_transformer_sat(solver, node_l, 0, node_sec, 0,
                                      TURNS_RATIO, &xf_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add transformer: %d\n", status);
        goto cleanup;
    }

    /* Configure transformer core */
    TransformerSat* xf = mna_get_transformer_sat(solver, xf_handle);
    if (!xf) {
        fprintf(stderr, "Failed to get transformer structure\n");
        goto cleanup;
    }

    status = mna_transformer_sat_setup_tanh_core(xf, B_SAT, MU_R_INITIAL, 0.0,
                                                  CORE_AREA, PATH_LENGTH, N_PRIMARY);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to setup transformer core: %d\n", status);
        goto cleanup;
    }

    /* Add core loss resistance in parallel with primary (models eddy current + hysteresis) */
    ComponentHandle r_core_handle;
    status = mna_add_resistor(solver, node_l, 0, R_CORE, &r_core_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add core loss resistor: %d\n", status);
        goto cleanup;
    }

    /* Add load resistor on secondary */
    ComponentHandle r_handle;
    status = mna_add_resistor(solver, node_sec, 0, R_LOAD, &r_handle);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Failed to add resistor: %d\n", status);
        goto cleanup;
    }

    /* Print CSV header to stdout */
    printf("time_s,v_source,v_primary,v_secondary,i_primary,i_secondary,i_mag,flux,L_mag,switch_state\n");

    /* Run transient simulation */
    double t = 0.0;
    double omega = 2.0 * M_PI * FREQ_HZ;
    int step = 0;
    int switch_closed = 0;

    while (t < T_FINAL) {
        /* Close switch at 50ms */
        if (!switch_closed && t >= SWITCH_CLOSE_TIME) {
            mna_set_switch_state(solver, sw_handle, 1);  /* Close switch */
            switch_closed = 1;
        }

        /* Update AC source voltage with initial phase - MUST be before solve */
        double v_src = V_PEAK * sin(omega * t + initial_phase);
        solver->components[vs_handle].value = v_src;

        /* Solve transient step */
        status = mna_solve_transient_step(solver, DT);
        if (status != MNA_SUCCESS) {
            fprintf(stderr, "Transient solve failed at t=%g: %d\n", t, status);
            break;
        }

        /* Get values - v_primary is at transformer primary (after R_winding and L_leakage) */
        double v_source = mna_get_node_voltage(solver, node_l);
        double v_primary = v_source;
        double v_secondary = mna_get_node_voltage(solver, node_sec);
        double i_primary = mna_transformer_sat_get_primary_current(solver, xf_handle);
        double i_secondary = mna_transformer_sat_get_secondary_current(solver, xf_handle);
        double i_mag = mna_transformer_sat_get_magnetizing_current(solver, xf_handle);
        double flux = mna_transformer_sat_get_flux_linkage(solver, xf_handle);
        double L_mag = mna_transformer_sat_get_magnetizing_inductance(solver, xf_handle);

        /* Output CSV data (every PRINT_INTERVAL steps = 1ms) */
        if (step % PRINT_INTERVAL == 0) {
            printf("%.8f,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.8f,%.6f,%d\n",
                   t, v_source, v_primary, v_secondary,
                   i_primary, i_secondary, i_mag, flux, L_mag, switch_closed);
        }

        t += DT;
        step++;
    }

    /* Print summary to stderr */
    fprintf(stderr, "# Simulation completed: %.1f ms, %d steps\n", t * 1000.0, step);
    fprintf(stderr, "# Initial phase: %.3f rad (%.1f deg)\n", initial_phase, initial_phase * 180.0 / M_PI);

cleanup:
    mna_destroy(solver);
    free(solver);
    return (status == MNA_SUCCESS) ? 0 : 1;
}
