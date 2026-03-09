#include "netlist.h"
#include "../solver/core.h"
#include "../solver/dc.h"
#include "../solver/ac.h"
#include "../solver/transient.h"
#include "../../elements/passive.h"
#include "../../elements/sources.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

static void skip_whitespace(const char** p) {
    while (**p && isspace((unsigned char)**p)) (*p)++;
}

static void parse_token(const char** p, char* token, int max_len) {
    skip_whitespace(p);
    int i = 0;
    while (**p && !isspace((unsigned char)**p) && **p != ';' && **p != '\n' && i < max_len - 1) {
        token[i++] = **p;
        (*p)++;
    }
    token[i] = '\0';
}

static double parse_value(const char** p) {
    skip_whitespace(p);
    char num_str[64];
    int i = 0;
    
    while (**p && (isdigit((unsigned char)**p) || **p == '.' || **p == '-' || **p == '+' || 
                   **p == 'e' || **p == 'E') && i < 63) {
        num_str[i++] = **p;
        (*p)++;
    }
    num_str[i] = '\0';
    
    double value = atof(num_str);
    
    /* Parse scale suffix */
    char suffix = **p;
    if (suffix) {
        switch (tolower(suffix)) {
            case 'k': value *= 1e3; (*p)++; break;
            case 'm': value *= 1e-3; (*p)++; break;
            case 'u': value *= 1e-6; (*p)++; break;
            case 'n': value *= 1e-9; (*p)++; break;
            case 'p': value *= 1e-12; (*p)++; break;
            case 'f': value *= 1e-15; (*p)++; break;
            case 'g': value *= 1e9; (*p)++; break;
            case 't': value *= 1e12; (*p)++; break;
        }
    }
    
    return value;
}

static void to_upper(char* str) {
    for (int i = 0; str[i]; i++) {
        str[i] = toupper((unsigned char)str[i]);
    }
}

static MNAStatus report_error(NetlistContext* ctx, int line, const char* msg) {
    ctx->error_line = line;
    strncpy(ctx->error_msg, msg, sizeof(ctx->error_msg) - 1);
    ctx->error_msg[sizeof(ctx->error_msg) - 1] = '\0';
    fprintf(stderr, "Error at line %d: %s\n", line, msg);
    return MNA_INVALID_PARAMETER;
}

/* ============================================================================
 * Node Management
 * ============================================================================ */

static int get_or_create_node(NetlistContext* ctx, int netlist_node) {
    if (netlist_node < 0 || netlist_node >= MNA_NETLIST_MAX_NODES) {
        return -1;
    }
    
    if (netlist_node == 0) {
        return 0; /* Ground */
    }
    
    if (ctx->node_map[netlist_node] == 0 && netlist_node != 0) {
        /* Create new node */
        int solver_node = mna_create_node(ctx->solver);
        if (solver_node < 0) {
            return -1;
        }
        ctx->node_map[netlist_node] = solver_node;
        if (solver_node > ctx->max_node_index) {
            ctx->max_node_index = solver_node;
        }
    }
    
    return ctx->node_map[netlist_node];
}

/* ============================================================================
 * Component Parsing
 * ============================================================================ */

static MNAStatus parse_resistor(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64];
    int node1, node2;
    double value;
    ComponentHandle handle;
    const char* p = line + 1; /* Skip 'R' */
    
    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    value = parse_value(&p);
    
    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);
    
    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for resistor");
    }
    
    if (mna_add_resistor(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add resistor");
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_capacitor(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64];
    int node1, node2;
    double value;
    ComponentHandle handle;
    const char* p = line + 1; /* Skip 'C' */
    
    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    value = parse_value(&p);
    
    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);
    
    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for capacitor");
    }
    
    if (mna_add_capacitor(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add capacitor");
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_inductor(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64];
    int node1, node2;
    double value;
    ComponentHandle handle;
    const char* p = line + 1; /* Skip 'L' */
    
    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    value = parse_value(&p);
    
    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);
    
    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for inductor");
    }
    
    if (mna_add_inductor(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add inductor");
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_voltage_source(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64], type[32];
    int node1, node2;
    double value = 0;
    ComponentHandle handle;
    const char* p = line + 1; /* Skip 'V' */
    
    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    parse_token(&p, type, sizeof(type));
    to_upper(type);
    
    if (strcmp(type, "DC") == 0) {
        value = parse_value(&p);
    } else if (strcmp(type, "AC") == 0) {
        /* AC source - parse magnitude and optionally phase */
        value = parse_value(&p);
        /* Phase would require extended API support */
    } else {
        /* Assume DC by default */
        p -= strlen(type);
        value = parse_value(&p);
    }
    
    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);
    
    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for voltage source");
    }
    
    if (mna_add_voltage_source(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add voltage source");
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_current_source(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64], type[32];
    int node1, node2;
    double value = 0;
    ComponentHandle handle;
    const char* p = line + 1; /* Skip 'I' */
    
    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    parse_token(&p, type, sizeof(type));
    to_upper(type);
    
    if (strcmp(type, "DC") == 0) {
        value = parse_value(&p);
    } else {
        p -= strlen(type);
        value = parse_value(&p);
    }
    
    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);
    
    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for current source");
    }
    
    if (mna_add_current_source(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add current source");
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_switch(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64];
    int node1, node2;
    double value = 0.001; /* Default on-resistance */
    ComponentHandle handle;
    const char* p = line + 1; /* Skip 'S' */
    
    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    
    /* Check for optional resistance value */
    skip_whitespace(&p);
    if (*p && *p != ';' && *p != '\n') {
        value = parse_value(&p);
    }
    
    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);
    
    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for switch");
    }
    
    if (mna_add_switch(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add switch");
    }
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Command Parsing
 * ============================================================================ */

static MNAStatus parse_analysis_cmd(NetlistContext* ctx, const char* line, int line_num) {
    char cmd[32];
    const char* p = line;
    
    parse_token(&p, cmd, sizeof(cmd));
    to_upper(cmd);
    
    if (strcmp(cmd, "DC") == 0) {
        ctx->analysis.dc_enabled = 1;
        ctx->analysis.type = MNA_ANALYSIS_DC;
    } else if (strcmp(cmd, "AC") == 0) {
        ctx->analysis.ac_enabled = 1;
        ctx->analysis.type = MNA_ANALYSIS_AC;
        ctx->analysis.ac_freq = parse_value(&p);
    } else if (strcmp(cmd, "TRAN") == 0) {
        ctx->analysis.tran_enabled = 1;
        ctx->analysis.type = MNA_ANALYSIS_TRAN;
        ctx->analysis.tran_dt = parse_value(&p);
        ctx->analysis.tran_end = parse_value(&p);
    } else {
        return report_error(ctx, line_num, "Unknown analysis type");
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_print_cmd(NetlistContext* ctx, const char* line, int line_num) {
    const char* p = line;
    char token[64];
    
    while (*p && *p != ';' && *p != '\n') {
        parse_token(&p, token, sizeof(token));
        if (strlen(token) == 0) break;
        
        if (ctx->num_outputs >= MNA_NETLIST_MAX_OUTPUTS) {
            return report_error(ctx, line_num, "Too many output variables");
        }
        
        OutputVar* out = &ctx->outputs[ctx->num_outputs];
        
        if (tolower(token[0]) == 'v') {
            out->type = MNA_OUTPUT_VOLTAGE;
            /* Parse node from v(node) */
            char* paren = strchr(token, '(');
            if (paren) {
                out->node = atoi(paren + 1);
            } else {
                return report_error(ctx, line_num, "Invalid voltage specification");
            }
        } else if (tolower(token[0]) == 'i') {
            out->type = MNA_OUTPUT_CURRENT;
            /* Parse component name from i(name) */
            char* paren = strchr(token, '(');
            if (paren) {
                char* end = strchr(paren, ')');
                if (end) {
                    int len = end - paren - 1;
                    if (len >= MNA_NETLIST_MAX_NAME) len = MNA_NETLIST_MAX_NAME - 1;
                    strncpy(out->name, paren + 1, len);
                    out->name[len] = '\0';
                } else {
                    return report_error(ctx, line_num, "Invalid current specification");
                }
            } else {
                return report_error(ctx, line_num, "Invalid current specification");
            }
        }
        
        ctx->num_outputs++;
    }
    
    return MNA_SUCCESS;
}

static MNAStatus parse_write_cmd(NetlistContext* ctx, const char* line, int line_num) {
    const char* p = line;
    parse_token(&p, ctx->output_file, sizeof(ctx->output_file));
    
    if (strlen(ctx->output_file) == 0) {
        return report_error(ctx, line_num, "No output file specified");
    }
    
    ctx->write_enabled = 1;
    return MNA_SUCCESS;
}

/* ============================================================================
 * Line Parser
 * ============================================================================ */

static MNAStatus parse_line(NetlistContext* ctx, const char* line, int line_num) {
    const char* p = line;
    
    /* Skip whitespace */
    skip_whitespace(&p);
    
    /* Skip empty lines and comments */
    if (*p == '\0' || *p == '*' || *p == ';' || *p == '#') {
        return MNA_SUCCESS;
    }
    
    /* Remove inline comments */
    char clean_line[MNA_NETLIST_MAX_LINE];
    strncpy(clean_line, p, sizeof(clean_line) - 1);
    clean_line[sizeof(clean_line) - 1] = '\0';
    
    char* comment = strchr(clean_line, ';');
    if (comment) *comment = '\0';
    comment = strchr(clean_line, '\n');
    if (comment) *comment = '\0';
    
    /* Trim trailing whitespace */
    int len = strlen(clean_line);
    while (len > 0 && isspace((unsigned char)clean_line[len - 1])) {
        clean_line[--len] = '\0';
    }
    
    if (strlen(clean_line) == 0) {
        return MNA_SUCCESS;
    }
    
    /* Check for dot commands */
    if (clean_line[0] == '.') {
        char cmd[32];
        const char* cmd_p = clean_line + 1;
        parse_token(&cmd_p, cmd, sizeof(cmd));
        to_upper(cmd);
        
        if (strcmp(cmd, "ANALYSIS") == 0) {
            return parse_analysis_cmd(ctx, cmd_p, line_num);
        } else if (strcmp(cmd, "PRINT") == 0) {
            return parse_print_cmd(ctx, cmd_p, line_num);
        } else if (strcmp(cmd, "WRITE") == 0) {
            return parse_write_cmd(ctx, cmd_p, line_num);
        } else if (strcmp(cmd, "END") == 0) {
            return MNA_SUCCESS; /* Will stop parsing */
        } else if (strcmp(cmd, "INCLUDE") == 0 || strcmp(cmd, "LIB") == 0) {
            return report_error(ctx, line_num, ".INCLUDE and .LIB not yet implemented");
        } else {
            return report_error(ctx, line_num, "Unknown command");
        }
    }
    
    /* Parse components */
    char first = toupper((unsigned char)clean_line[0]);
    switch (first) {
        case 'R':
            return parse_resistor(ctx, clean_line, line_num);
        case 'C':
            return parse_capacitor(ctx, clean_line, line_num);
        case 'L':
            return parse_inductor(ctx, clean_line, line_num);
        case 'V':
            return parse_voltage_source(ctx, clean_line, line_num);
        case 'I':
            return parse_current_source(ctx, clean_line, line_num);
        case 'S':
            return parse_switch(ctx, clean_line, line_num);
        default:
            return report_error(ctx, line_num, "Unknown component type");
    }
}

/* ============================================================================
 * Public API
 * ============================================================================ */

MNAStatus netlist_init(NetlistContext* ctx, MNASolver* solver) {
    memset(ctx, 0, sizeof(NetlistContext));
    ctx->solver = solver;
    ctx->analysis.type = MNA_ANALYSIS_NONE;
    ctx->num_outputs = 0;
    ctx->write_enabled = 0;
    ctx->error_line = 0;
    ctx->error_msg[0] = '\0';
    
    /* Initialize node 0 (ground) */
    ctx->node_map[0] = 0;
    
    return MNA_SUCCESS;
}

MNAStatus netlist_parse_file(NetlistContext* ctx, const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open file '%s'\n", filename);
        return MNA_INVALID_PARAMETER;
    }
    
    char line[MNA_NETLIST_MAX_LINE];
    int line_num = 0;
    MNAStatus status = MNA_SUCCESS;
    
    while (fgets(line, sizeof(line), f)) {
        line_num++;
        status = parse_line(ctx, line, line_num);
        if (status != MNA_SUCCESS) {
            fclose(f);
            return status;
        }
        
        /* Check for .end command */
        if (strstr(line, ".END") || strstr(line, ".end")) {
            break;
        }
    }
    
    fclose(f);
    return status;
}

MNAStatus netlist_parse_string(NetlistContext* ctx, const char* content) {
    char* copy = strdup(content);
    if (!copy) {
        return MNA_INSUFFICIENT_MEMORY;
    }
    
    char* line = strtok(copy, "\n");
    int line_num = 0;
    MNAStatus status = MNA_SUCCESS;
    
    while (line) {
        line_num++;
        status = parse_line(ctx, line, line_num);
        if (status != MNA_SUCCESS) {
            free(copy);
            return status;
        }
        
        if (strstr(line, ".END") || strstr(line, ".end")) {
            break;
        }
        
        line = strtok(NULL, "\n");
    }
    
    free(copy);
    return status;
}

MNAStatus netlist_run_analysis(NetlistContext* ctx) {
    if (ctx->analysis.type == MNA_ANALYSIS_NONE) {
        fprintf(stderr, "Warning: No analysis specified\n");
        return MNA_SUCCESS;
    }
    
    MNAStatus status;
    FILE* out_file = NULL;
    
    if (ctx->write_enabled) {
        out_file = fopen(ctx->output_file, "w");
        if (!out_file) {
            fprintf(stderr, "Error: Cannot create output file '%s'\n", ctx->output_file);
            return MNA_INVALID_PARAMETER;
        }
    }
    
    /* DC Analysis */
    if (ctx->analysis.dc_enabled || ctx->analysis.type == MNA_ANALYSIS_DC) {
        status = mna_solve_dc(ctx->solver);
        if (status != MNA_SUCCESS) {
            fprintf(stderr, "Error: DC analysis failed\n");
            if (out_file) fclose(out_file);
            return status;
        }
        
        printf("=== DC Analysis Results ===\n");
        if (out_file) {
            fprintf(out_file, "=== DC Analysis Results ===\n");
        }
        
        for (int i = 0; i < ctx->num_outputs; i++) {
            OutputVar* out = &ctx->outputs[i];
            if (out->type == MNA_OUTPUT_VOLTAGE) {
                int node = ctx->node_map[out->node];
                double v = mna_get_node_voltage(ctx->solver, node);
                printf("V(%d): %.6f V\n", out->node, v);
                if (out_file) {
                    fprintf(out_file, "V(%d): %.6f V\n", out->node, v);
                }
            }
        }
        printf("\n");
        if (out_file) fprintf(out_file, "\n");
    }
    
    /* AC Analysis */
    if (ctx->analysis.ac_enabled || ctx->analysis.type == MNA_ANALYSIS_AC) {
        status = mna_solve_ac(ctx->solver, ctx->analysis.ac_freq);
        if (status != MNA_SUCCESS) {
            fprintf(stderr, "Error: AC analysis failed\n");
            if (out_file) fclose(out_file);
            return status;
        }
        
        printf("=== AC Analysis Results (f = %.2f Hz) ===\n", ctx->analysis.ac_freq);
        if (out_file) {
            fprintf(out_file, "=== AC Analysis Results (f = %.2f Hz) ===\n", ctx->analysis.ac_freq);
        }
        
        for (int i = 0; i < ctx->num_outputs; i++) {
            OutputVar* out = &ctx->outputs[i];
            if (out->type == MNA_OUTPUT_VOLTAGE) {
                int node = ctx->node_map[out->node];
                double complex v = mna_get_ac_node_voltage(ctx->solver, node);
                double mag = cabs(v);
                double phase = carg(v) * 180.0 / M_PI;
                printf("V(%d): %.6f < %.2f deg\n", out->node, mag, phase);
                if (out_file) {
                    fprintf(out_file, "V(%d): %.6f < %.2f deg\n", out->node, mag, phase);
                }
            }
        }
        printf("\n");
        if (out_file) fprintf(out_file, "\n");
    }
    
    /* Transient Analysis */
    if (ctx->analysis.tran_enabled || ctx->analysis.type == MNA_ANALYSIS_TRAN) {
        double dt = ctx->analysis.tran_dt;
        double t_end = ctx->analysis.tran_end;
        int steps = (int)(t_end / dt);
        
        /* Initialize transient */
        mna_init_transient(ctx->solver);
        
        printf("=== Transient Analysis ===\n");
        printf("Time step: %.2e s, End time: %.2e s, Steps: %d\n", dt, t_end, steps);
        
        if (out_file) {
            fprintf(out_file, "time");
            for (int i = 0; i < ctx->num_outputs; i++) {
                OutputVar* out = &ctx->outputs[i];
                if (out->type == MNA_OUTPUT_VOLTAGE) {
                    fprintf(out_file, ",V(%d)", out->node);
                } else {
                    fprintf(out_file, ",I(%s)", out->name);
                }
            }
            fprintf(out_file, "\n");
        }
        
        for (int step = 0; step <= steps; step++) {
            double t = step * dt;
            
            if (step > 0) {
                status = mna_solve_transient_step(ctx->solver, dt);
                if (status != MNA_SUCCESS) {
                    fprintf(stderr, "Error: Transient step %d failed\n", step);
                    if (out_file) fclose(out_file);
                    return status;
                }
            }
            
            /* Output results */
            if (out_file) {
                fprintf(out_file, "%.9e", t);
                for (int i = 0; i < ctx->num_outputs; i++) {
                    OutputVar* out = &ctx->outputs[i];
                    if (out->type == MNA_OUTPUT_VOLTAGE) {
                        int node = ctx->node_map[out->node];
                        double v = mna_get_node_voltage(ctx->solver, node);
                        fprintf(out_file, ",%.9e", v);
                    }
                }
                fprintf(out_file, "\n");
            }
            
            /* Progress indicator */
            if (step % (steps / 10 + 1) == 0) {
                printf("Progress: %.0f%% (t=%.6f ms)\n", (double)step / steps * 100, t * 1000);
            }
        }
        
        if (out_file) {
            fclose(out_file);
            printf("Results written to: %s\n", ctx->output_file);
        }
    }
    
    return MNA_SUCCESS;
}

void netlist_destroy(NetlistContext* ctx) {
    /* Nothing to free - solver is owned by caller */
    (void)ctx;
}
