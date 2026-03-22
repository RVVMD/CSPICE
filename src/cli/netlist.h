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
#define MNA_NETLIST_MAX_TRANSFORMER_MODELS 16
#define MNA_NETLIST_MAX_COMPONENT_NAMES 128

typedef struct {
    char name[MNA_NETLIST_MAX_NAME];
    int handle;
    int valid;
} ComponentNameMap;

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

typedef struct {
    char name[MNA_NETLIST_MAX_NAME];
    int valid;
    int component_handle;  /* Handle of the instantiated component */
    
    /* Required: Electrical nominals (ГОСТ style) */
    double V1_nom;        /* primary voltage [V] */
    double V2_nom;        /* secondary voltage [V] */
    double f_nom;         /* frequency [Hz] */
    
    /* Required: Core geometry (for saturation) */
    double Acore;         /* core area [m^2] */
    double lpath;         /* magnetic path length [m] */
    int Np;               /* primary turns */
    
    /* Optional: Losses (for auto-compute) */
    double S_nom;         /* apparent power [VA], 0 = not provided */
    double P_core;        /* core loss [W] */
    double P_cu;          /* copper loss [W] */
    double u_k;           /* impedance voltage [pu] */
    
    /* Optional: Direct parameters (override auto-compute) */
    double R1;            /* primary resistance [ohm] */
    double R2;            /* secondary resistance [ohm] */
    double L1;            /* primary leakage inductance [H] */
    double L2;            /* secondary leakage inductance [H] */
    double Rcore;         /* core loss resistance [ohm] */
    
    /* Required for saturation: B-H curve */
    double Bsat;          /* saturation flux density [T] */
    double Mur;           /* initial relative permeability */
    double Brem;          /* remanent flux density [T] */
    
    /* Computed turns ratio */
    double turns_ratio;   /* n = V1_nom / V2_nom */
    
    /* Flags for which parameters were explicitly set */
    int has_S_nom;
    int has_P_core;
    int has_P_cu;
    int has_u_k;
    int has_R1;
    int has_R2;
    int has_L1;
    int has_L2;
    int has_Rcore;
    int has_Brem;
} TransformerModel;

typedef enum {
    MNA_ANALYSIS_NONE,
    MNA_ANALYSIS_DC,
    MNA_ANALYSIS_AC,
    MNA_ANALYSIS_TRAN
} AnalysisType;

typedef enum {
    MNA_OUTPUT_VOLTAGE,
    MNA_OUTPUT_CURRENT,
    MNA_OUTPUT_I_PRI,
    MNA_OUTPUT_I_SEC,
    MNA_OUTPUT_I_MAG,
    MNA_OUTPUT_FLUX,
    MNA_OUTPUT_L_MAG,
    MNA_OUTPUT_B_FIELD,
    MNA_OUTPUT_H_FIELD
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
    TransformerModel transformer_models[MNA_NETLIST_MAX_TRANSFORMER_MODELS];
    int num_transformer_models;
    ComponentNameMap component_names[MNA_NETLIST_MAX_COMPONENT_NAMES];
    int num_component_names;
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
