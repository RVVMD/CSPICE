#ifndef MNA_CLI_NETLIST_H
#define MNA_CLI_NETLIST_H

#include "../../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MNA_NETLIST_MAX_LINE 512
#define MNA_NETLIST_MAX_NAME 32
#define MNA_NETLIST_MAX_NODES 256
#define MNA_NETLIST_MAX_OUTPUTS 64

/* ============================================================================
 * Analysis Types
 * ============================================================================ */
typedef enum {
    MNA_ANALYSIS_NONE,
    MNA_ANALYSIS_DC,
    MNA_ANALYSIS_AC,
    MNA_ANALYSIS_TRAN
} AnalysisType;

/* ============================================================================
 * Output Variable Types
 * ============================================================================ */
typedef enum {
    MNA_OUTPUT_VOLTAGE,
    MNA_OUTPUT_CURRENT
} OutputType;

/* ============================================================================
 * Output Variable Specification
 * ============================================================================ */
typedef struct {
    OutputType type;
    char name[MNA_NETLIST_MAX_NAME];  /* Component name for currents */
    int node;                          /* Node for voltages */
} OutputVar;

/* ============================================================================
 * Analysis Configuration
 * ============================================================================ */
typedef struct {
    AnalysisType type;
    /* DC */
    int dc_enabled;
    /* AC */
    double ac_freq;
    int ac_enabled;
    /* Transient */
    double tran_dt;
    double tran_end;
    int tran_enabled;
} AnalysisConfig;

/* ============================================================================
 * Netlist Parser Context
 * ============================================================================ */
typedef struct {
    MNASolver* solver;
    AnalysisConfig analysis;
    OutputVar outputs[MNA_NETLIST_MAX_OUTPUTS];
    int num_outputs;
    char output_file[256];
    int write_enabled;
    /* Node mapping: netlist node index -> solver node index */
    int node_map[MNA_NETLIST_MAX_NODES];
    int max_node_index;
    /* Error tracking */
    int error_line;
    char error_msg[256];
} NetlistContext;

/* ============================================================================
 * Parser Functions
 * ============================================================================ */

/**
 * @brief Initialize a netlist context
 * @param ctx Netlist context to initialize
 * @param solver MNA solver instance
 * @return MNAStatus indicating success or failure
 */
MNAStatus netlist_init(NetlistContext* ctx, MNASolver* solver);

/**
 * @brief Parse and execute a netlist file
 * @param ctx Netlist context
 * @param filename Path to the netlist file
 * @return MNAStatus indicating success or failure
 */
MNAStatus netlist_parse_file(NetlistContext* ctx, const char* filename);

/**
 * @brief Parse a netlist from a string
 * @param ctx Netlist context
 * @param content Netlist content as string
 * @return MNAStatus indicating success or failure
 */
MNAStatus netlist_parse_string(NetlistContext* ctx, const char* content);

/**
 * @brief Run the configured analysis
 * @param ctx Netlist context
 * @return MNAStatus indicating success or failure
 */
MNAStatus netlist_run_analysis(NetlistContext* ctx);

/**
 * @brief Free resources in netlist context
 * @param ctx Netlist context
 */
void netlist_destroy(NetlistContext* ctx);

#ifdef __cplusplus
}
#endif

#endif /* MNA_CLI_NETLIST_H */
