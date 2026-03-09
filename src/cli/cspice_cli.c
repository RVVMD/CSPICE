/**
 * @file cspice_cli.c
 * @brief CSPICE Command Line Interface for netlist-based circuit simulation
 *
 * Usage:
 *   cspice <netlist_file> [options]
 *   cspice -h | --help
 *
 * Netlist Format (SPICE-like):
 *   * Comment
 *   R1 n1 n2 1k      ; Resistor
 *   C1 n1 0 1u       ; Capacitor
 *   L1 n1 n2 1m      ; Inductor
 *   V1 n1 0 DC 10    ; Voltage source
 *   I1 n1 n2 DC 1    ; Current source
 *   S1 n1 n2         ; Switch
 *   .analysis dc     ; DC analysis
 *   .analysis ac 50  ; AC analysis at 50 Hz
 *   .analysis tran 1u 1m  ; Transient: step=1us, end=1ms
 *   .print v(n1)     ; Print voltage at node 1
 *   .write out.csv   ; Write results to CSV
 *   .end
 */

#include "mna_solver.h"
#include "netlist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_help(const char* prog_name) {
    printf("CSPICE - Circuit Simulation for Relay Protection\n");
    printf("Version %s\n\n", MNA_SOLVER_VERSION);
    printf("Usage: %s <netlist_file> [options]\n\n", prog_name);
    printf("Options:\n");
    printf("  -h, --help     Show this help message\n");
    printf("  -v, --version  Show version information\n");
    printf("  -q, --quiet    Suppress progress output\n");
    printf("\nNetlist Format Examples:\n");
    printf("  R1 n1 n2 1k        Resistor between nodes n1 and n2\n");
    printf("  C1 n1 0 1u         Capacitor to ground\n");
    printf("  V1 n1 0 DC 10      DC voltage source\n");
    printf("  .analysis tran 1u 1m  Transient analysis\n");
    printf("  .write output.csv  Export results\n");
    printf("\nComponent Prefixes:\n");
    printf("  R - Resistor       V - Voltage source\n");
    printf("  C - Capacitor      I - Current source\n");
    printf("  L - Inductor       S - Switch\n");
    printf("\nAnalysis Commands:\n");
    printf("  .analysis dc           DC operating point\n");
    printf("  .analysis ac <freq>    AC analysis at frequency (Hz)\n");
    printf("  .analysis tran <dt> <tend>  Transient analysis\n");
    printf("\nOutput Commands:\n");
    printf("  .print v(<node>)       Print node voltage\n");
    printf("  .print i(<component>)  Print component current\n");
    printf("  .write <file.csv>      Write results to file\n");
}

static void print_version(void) {
    printf("CSPICE Version %s\n", MNA_SOLVER_VERSION);
    printf("Modified Nodal Analysis Solver with TR-BDF2 Integration\n");
    printf("Target: Relay Protection and Automation (RZA)\n");
}

int main(int argc, char* argv[]) {
    const char* netlist_file = NULL;
    int quiet = 0;
    
    /* Parse command line arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            print_version();
            return 0;
        } else if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--quiet") == 0) {
            quiet = 1;
        } else if (argv[i][0] != '-') {
            netlist_file = argv[i];
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            fprintf(stderr, "Use -h for help.\n");
            return 1;
        }
    }
    
    if (!netlist_file) {
        fprintf(stderr, "Error: No netlist file specified\n");
        fprintf(stderr, "Use -h for help.\n");
        return 1;
    }
    
    /* Initialize solver */
    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize solver\n");
        return 1;
    }
    
    /* Initialize netlist context */
    NetlistContext ctx;
    status = netlist_init(&ctx, &solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize netlist parser\n");
        mna_destroy(&solver);
        return 1;
    }
    
    /* Parse netlist file */
    if (!quiet) {
        printf("Parsing netlist: %s\n", netlist_file);
    }
    
    status = netlist_parse_file(&ctx, netlist_file);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Error: Failed to parse netlist\n");
        netlist_destroy(&ctx);
        mna_destroy(&solver);
        return 1;
    }
    
    /* Run analysis */
    if (!quiet) {
        printf("\n");
    }
    
    status = netlist_run_analysis(&ctx);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Error: Analysis failed\n");
        netlist_destroy(&ctx);
        mna_destroy(&solver);
        return 1;
    }
    
    /* Cleanup */
    netlist_destroy(&ctx);
    mna_destroy(&solver);
    
    if (!quiet) {
        printf("\nSimulation completed successfully.\n");
    }
    
    return 0;
}
