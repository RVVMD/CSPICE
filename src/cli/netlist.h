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
#define MNA_NETLIST_MAX_SWITCH_EVENTS 32
#define MNA_NETLIST_MAX_MODELS 16

typedef struct {
    double time;
    char name[MNA_NETLIST_MAX_NAME];
    int state;
    int component_idx;
    int valid;
} SwitchEvent;

typedef struct {
    char name[MNA_NETLIST_MAX_NAME];
    double Is;
    double n;
    double BV;
    double IBV;
    int valid;
} DiodeModel;

typedef enum {
    MNA_ANALYSIS_NONE,
    MNA_ANALYSIS_DC,
    MNA_ANALYSIS_AC,
    MNA_ANALYSIS_TRAN
} AnalysisType;

typedef enum {
    MNA_OUTPUT_VOLTAGE,
    MNA_OUTPUT_CURRENT
} OutputType;

typedef struct {
    OutputType type;
    char name[MNA_NETLIST_MAX_NAME];
    int node;
} OutputVar;

typedef struct {
    AnalysisType type;
    int dc_enabled;
    double ac_freq;
    int ac_enabled;
    double tran_dt;
    double tran_end;
    int tran_enabled;
    double sin_freq;
    int sin_source_idx;
    int sin_source_valid;
    SwitchEvent switch_events[MNA_NETLIST_MAX_SWITCH_EVENTS];
    int num_switch_events;
} AnalysisConfig;

typedef struct {
    MNASolver* solver;
    AnalysisConfig analysis;
    OutputVar outputs[MNA_NETLIST_MAX_OUTPUTS];
    int num_outputs;
    char output_file[256];
    int write_enabled;
    int node_map[MNA_NETLIST_MAX_NODES];
    int max_node_index;
    char switch_names[MNA_NETLIST_MAX_SWITCH_EVENTS][MNA_NETLIST_MAX_NAME];
    int switch_indices[MNA_NETLIST_MAX_SWITCH_EVENTS];
    int num_switches;
    DiodeModel diode_models[MNA_NETLIST_MAX_MODELS];
    int num_diode_models;
    int error_line;
    char error_msg[256];
} NetlistContext;

MNAStatus netlist_init(NetlistContext* ctx, MNASolver* solver);
MNAStatus netlist_parse_file(NetlistContext* ctx, const char* filename);
MNAStatus netlist_parse_string(NetlistContext* ctx, const char* content);
MNAStatus netlist_run_analysis(NetlistContext* ctx);
void netlist_destroy(NetlistContext* ctx);

#ifdef __cplusplus
}
#endif

#endif /* MNA_CLI_NETLIST_H */
