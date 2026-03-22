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
    printf("  TR1 n1 n2 n3 n4 10 Ideal transformer (n1-n2: primary, n3-n4: secondary)\n");
    printf("  TRSAT1 n1 n2 n3 n4 XFMR1  Saturated transformer with model\n");
    printf("  .analysis tran 1u 1m  Transient analysis\n");
    printf("  .write output.csv  Export results\n");
    printf("\nComponent Prefixes:\n");
    printf("  R - Resistor       V - Voltage source     TR - Transformer\n");
    printf("  C - Capacitor      I - Current source     TRSAT - Saturated transformer\n");
    printf("  L - Inductor       S - Switch             D - Diode\n");
    printf("\nAnalysis Commands:\n");
    printf("  .analysis dc           DC operating point\n");
    printf("  .analysis ac <freq>    AC analysis at frequency (Hz)\n");
    printf("  .analysis tran <dt> <tend>  Transient analysis\n");
    printf("\nModel Definitions:\n");
    printf("  .model DNAME D (Is=1e-14 n=1.0 BV=100)  Diode model\n");
    printf("  .model XFMR1 TRSAT (V1_nom=220 V2_nom=12 f_nom=50\n");
    printf("                       Acore=4e-4 lpath=0.25 Np=300\n");
    printf("                       Bsat=1.6 Mur=3000\n");
    printf("                       S_nom=100 P_core=2.5 P_cu=4.0 u_k=0.10)  Transformer\n");
    printf("\nTransformer Model Parameters (ГОСТ hybrid style):\n");
    printf("  Required - Electrical: V1_nom, V2_nom, f_nom\n");
    printf("  Required - Core: Acore, lpath, Np, Bsat, Mur\n");
    printf("  Optional - ГОСТ: S_nom, P_core, P_cu, u_k (auto-computes R/L)\n");
    printf("  Optional - Direct: R1, R2, L1, L2, Rcore\n");
    printf("\nOutput Commands:\n");
    printf("  .print v(<node>)       Print node voltage\n");
    printf("  .print i(<component>)  Print component current\n");
    printf("  .print i_pri(<xfmr>)   Print transformer primary current\n");
    printf("  .print i_sec(<xfmr>)   Print transformer secondary current\n");
    printf("  .print i_mag(<xfmr>)   Print transformer magnetizing current\n");
    printf("  .print flux(<xfmr>)    Print transformer flux linkage [Wb]\n");
    printf("  .print L_mag(<xfmr>)   Print transformer magnetizing inductance [H]\n");
    printf("  .print B(<xfmr>)       Print magnetic flux density [T]\n");
    printf("  .print H(<xfmr>)       Print magnetic field intensity [A/m]\n");
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

    MNASolver solver;
    MNAStatus status = mna_init(&solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize solver\n");
        return 1;
    }

    NetlistContext ctx;
    status = netlist_init(&ctx, &solver);
    if (status != MNA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize netlist parser\n");
        mna_destroy(&solver);
        return 1;
    }

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

    netlist_destroy(&ctx);
    mna_destroy(&solver);

    if (!quiet) {
        printf("\nSimulation completed successfully.\n");
    }

    return 0;
}
