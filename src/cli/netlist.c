#include "netlist.h"
#include "../solver/core.h"
#include "../solver/dc.h"
#include "../solver/ac.h"
#include "../solver/transient.h"
#include "../../elements/passive.h"
#include "../../elements/sources.h"
#include "../../elements/nonlinear/nonlinear.h"
#include "../../elements/transformer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static void skip_whitespace(const char** p) {
    while (**p && isspace((unsigned char)**p)) (*p)++;
}

static void parse_token(const char** p, char* token, int max_len) {
    skip_whitespace(p);
    int i = 0;
    while (**p && !isspace((unsigned char)**p) && **p != ';' && **p != '\n' && **p != '(' && i < max_len - 1) {
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

static int get_or_create_node(NetlistContext* ctx, int netlist_node) {
    if (netlist_node < 0 || netlist_node >= MNA_NETLIST_MAX_NODES) {
        return -1;
    }

    if (netlist_node == 0) {
        return 0;
    }

    if (ctx->node_map[netlist_node] == 0 && netlist_node != 0) {
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

static MNAStatus parse_resistor(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64];
    int node1, node2;
    double value;
    ComponentHandle handle;
    const char* p = line + 1;

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
    const char* p = line + 1;

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
    const char* p = line + 1;

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
    const char* p = line + 1;

    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);
    parse_token(&p, type, sizeof(type));
    to_upper(type);

    char* paren = strchr(type, '(');
    if (paren) *paren = '\0';

    if (strcmp(type, "DC") == 0) {
        value = parse_value(&p);
    } else if (strcmp(type, "AC") == 0) {
        double ac_mag = parse_value(&p);
        double ac_phase = 0.0;
        skip_whitespace(&p);
        if (*p && *p != ';' && *p != '\n') {
            ac_phase = parse_value(&p);
        }
        value = 0.0;

        int n1 = get_or_create_node(ctx, node1);
        int n2 = get_or_create_node(ctx, node2);

        if (n1 < 0 || n2 < 0) {
            return report_error(ctx, line_num, "Failed to create nodes for voltage source");
        }

        if (mna_add_voltage_source(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
            return report_error(ctx, line_num, "Failed to add voltage source");
        }

        mna_set_ac_source(ctx->solver, handle, ac_mag, ac_phase);
        return MNA_SUCCESS;
    } else if (strcmp(type, "SIN") == 0) {
        skip_whitespace(&p);
        if (*p == '(') p++;

        double offset = parse_value(&p);
        double peak = parse_value(&p);
        double freq = parse_value(&p);
        double phase = 0.0;
        skip_whitespace(&p);
        if (*p && *p != ';' && *p != '\n' && *p != ')') {
            phase = parse_value(&p);
        }

        value = offset;

        int n1 = get_or_create_node(ctx, node1);
        int n2 = get_or_create_node(ctx, node2);

        if (n1 < 0 || n2 < 0) {
            return report_error(ctx, line_num, "Failed to create nodes for voltage source");
        }

        if (mna_add_voltage_source(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
            return report_error(ctx, line_num, "Failed to add voltage source");
        }

        ctx->analysis.sin_freq = freq;
        ctx->analysis.sin_source_idx = handle;
        ctx->analysis.sin_source_valid = 1;

        ctx->solver->components[handle].ac_magnitude = peak;
        ctx->solver->components[handle].ac_phase = phase * M_PI / 180.0;

        return MNA_SUCCESS;
    } else {
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
    const char* p = line + 1;

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
    double value = 0.001;
    int initial_state = 0;
    ComponentHandle handle;
    const char* p = line + 1;

    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);

    skip_whitespace(&p);
    if (*p && *p != ';' && *p != '\n') {
        value = parse_value(&p);
    }

    skip_whitespace(&p);
    if (*p && isdigit((unsigned char)*p)) {
        initial_state = atoi(p);
    }

    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);

    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for switch");
    }

    if (mna_add_switch(ctx->solver, n1, n2, value, &handle) != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add switch");
    }

    mna_set_switch_state(ctx->solver, handle, initial_state);

    if (ctx->num_switches < MNA_NETLIST_MAX_SWITCH_EVENTS) {
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-truncation"
#endif
        snprintf(ctx->switch_names[ctx->num_switches], MNA_NETLIST_MAX_NAME, "S%s", name);
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
        ctx->switch_indices[ctx->num_switches] = handle;
        ctx->num_switches++;
    }

    return MNA_SUCCESS;
}

static MNAStatus parse_switch_cmd(NetlistContext* ctx, const char* line, int line_num) {
    const char* p = line;
    char token[64];
    double time;
    char sw_name[MNA_NETLIST_MAX_NAME];
    int state;

    skip_whitespace(&p);

    time = parse_value(&p);
    parse_token(&p, sw_name, sizeof(sw_name));
    parse_token(&p, token, sizeof(token));
    state = atoi(token);

    if (ctx->analysis.num_switch_events >= MNA_NETLIST_MAX_SWITCH_EVENTS) {
        return report_error(ctx, line_num, "Too many switch events");
    }

    SwitchEvent* evt = &ctx->analysis.switch_events[ctx->analysis.num_switch_events];
    evt->time = time;
    strncpy(evt->name, sw_name, MNA_NETLIST_MAX_NAME - 1);
    evt->name[MNA_NETLIST_MAX_NAME - 1] = '\0';
    evt->state = state ? 1 : 0;
    evt->valid = 1;
    evt->component_idx = -1;

    ctx->analysis.num_switch_events++;

    return MNA_SUCCESS;
}

static MNAStatus parse_transformer(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64];
    int node_p1, node_p2, node_s1, node_s2;
    double ratio, Lm = 0.0, Isat = 0.0, sat_factor = 1.0;
    ComponentHandle handle;
    const char* p = line + 1;

    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node_p1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node_p2 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node_s1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node_s2 = atoi(token);

    skip_whitespace(&p);
    if (!*p || *p == ';' || *p == '\n') {
        return report_error(ctx, line_num, "Transformer requires turns ratio");
    }
    ratio = parse_value(&p);

    skip_whitespace(&p);
    if (*p && *p != ';' && *p != '\n' && isdigit((unsigned char)*p)) {
        Lm = parse_value(&p);

        skip_whitespace(&p);
        if (*p && *p != ';' && *p != '\n' && isdigit((unsigned char)*p)) {
            Isat = parse_value(&p);

            skip_whitespace(&p);
            if (*p && *p != ';' && *p != '\n' && isdigit((unsigned char)*p)) {
                sat_factor = parse_value(&p);
            }
        }
    }

    int n_p1 = get_or_create_node(ctx, node_p1);
    int n_p2 = get_or_create_node(ctx, node_p2);
    int n_s1 = get_or_create_node(ctx, node_s1);
    int n_s2 = get_or_create_node(ctx, node_s2);

    if (n_p1 < 0 || n_p2 < 0 || n_s1 < 0 || n_s2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for transformer");
    }

    MNAStatus status;
    if (Lm > 0 && Isat > 0) {
        status = mna_add_transformer_with_saturation(ctx->solver,
                                                      n_p1, n_p2, n_s1, n_s2,
                                                      ratio, Lm, Isat, sat_factor,
                                                      &handle);
    } else if (Lm > 0) {
        status = mna_add_voltage_transformer(ctx->solver,
                                              n_p1, n_p2, n_s1, n_s2,
                                              ratio, Lm, &handle);
    } else {
        status = mna_add_ideal_transformer(ctx->solver,
                                            n_p1, n_p2, n_s1, n_s2,
                                            ratio, &handle);
    }

    if (status != MNA_SUCCESS) {
        return report_error(ctx, line_num, "Failed to add transformer");
    }

    return MNA_SUCCESS;
}

static void diode_iv_function(const ComponentState* state, void* user_data,
                              double* current, double* conductance) {
    double* params = (double*)user_data;
    double Is = params[0];
    double n = params[1];
    double BV = params[2];
    double IBV = params[3];

    double v = state->voltage;
    double vt = MNA_VT;

    if (BV > 0 && v < -BV) {
        double ibv_actual = (IBV > 0) ? IBV : 0.01;
        double Rz = 10.0;

        double v_over = -v - BV;
        *current = -(ibv_actual + v_over / Rz);
        *conductance = 1.0 / Rz;
    } else {
        double exp_arg = v / (n * vt);
        if (exp_arg > 700) exp_arg = 700;
        if (exp_arg < -700) exp_arg = -700;

        *current = Is * (exp(exp_arg) - 1.0);
        *conductance = Is / (n * vt) * exp(exp_arg);
    }

    if (*conductance < MNA_MIN_CONDUCTANCE) {
        *conductance = MNA_MIN_CONDUCTANCE;
    }
    if (*conductance > MNA_MAX_CONDUCTANCE) {
        *conductance = MNA_MAX_CONDUCTANCE;
    }
}

static MNAStatus parse_diode(NetlistContext* ctx, const char* line, int line_num) {
    char name[MNA_NETLIST_MAX_NAME], token[64], model_name[MNA_NETLIST_MAX_NAME];
    int node1, node2;
    ComponentHandle handle;
    const char* p = line + 1;

    parse_token(&p, name, sizeof(name));
    parse_token(&p, token, sizeof(token));
    node1 = atoi(token);
    parse_token(&p, token, sizeof(token));
    node2 = atoi(token);

    parse_token(&p, model_name, sizeof(model_name));
    to_upper(model_name);

    int n1 = get_or_create_node(ctx, node1);
    int n2 = get_or_create_node(ctx, node2);

    if (n1 < 0 || n2 < 0) {
        return report_error(ctx, line_num, "Failed to create nodes for diode");
    }

    DiodeModel* model = NULL;
    for (int i = 0; i < ctx->num_diode_models; i++) {
        char upper_name[MNA_NETLIST_MAX_NAME];
        strncpy(upper_name, ctx->diode_models[i].name, MNA_NETLIST_MAX_NAME - 1);
        upper_name[MNA_NETLIST_MAX_NAME - 1] = '\0';
        to_upper(upper_name);
        if (strcmp(upper_name, model_name) == 0) {
            model = &ctx->diode_models[i];
            break;
        }
    }

    if (!model) {
        static DiodeModel default_model = {
            .name = "DEFAULT",
            .Is = 1e-14,
            .n = 1.0,
            .BV = 0,
            .IBV = 0,
            .valid = 1
        };
        model = &default_model;
    }

    double* diode_params = (double*)malloc(4 * sizeof(double));
    if (!diode_params) {
        return report_error(ctx, line_num, "Failed to allocate diode parameters");
    }
    diode_params[0] = model->Is;
    diode_params[1] = model->n;
    diode_params[2] = model->BV;
    diode_params[3] = model->IBV;

    if (mna_add_custom_nonlinear(ctx->solver, n1, n2, NONLINEAR_RESISTOR,
                                  diode_iv_function, diode_params, 0.0, 0.0, &handle) != MNA_SUCCESS) {
        free(diode_params);
        return report_error(ctx, line_num, "Failed to add diode");
    }

    return MNA_SUCCESS;
}

static MNAStatus parse_model_cmd(NetlistContext* ctx, const char* line, int line_num) {
    char model_name[MNA_NETLIST_MAX_NAME], model_type[32];
    const char* p = line;

    parse_token(&p, model_name, sizeof(model_name));
    parse_token(&p, model_type, sizeof(model_type));
    to_upper(model_type);

    if (strcmp(model_type, "D") != 0) {
        return MNA_SUCCESS;
    }

    if (ctx->num_diode_models >= MNA_NETLIST_MAX_MODELS) {
        return report_error(ctx, line_num, "Too many diode models");
    }

    DiodeModel* model = &ctx->diode_models[ctx->num_diode_models];
    strncpy(model->name, model_name, MNA_NETLIST_MAX_NAME - 1);
    model->name[MNA_NETLIST_MAX_NAME - 1] = '\0';

    model->Is = 1e-14;
    model->n = 1.0;
    model->BV = 0;
    model->IBV = 0;
    model->valid = 1;

    skip_whitespace(&p);
    if (*p == '(') p++;

    char param_name[32];

    while (*p && *p != ';' && *p != '\n') {
        skip_whitespace(&p);
        if (*p == '\0' || *p == ';' || *p == '\n') break;

        int i = 0;
        while (*p && isalpha((unsigned char)*p) && i < 31) {
            param_name[i++] = *p++;
        }
        param_name[i] = '\0';

        if (i == 0) break;

        if (*p == '=') p++;

        double param_value = parse_value(&p);

        if (strcmp(param_name, "Is") == 0 || strcmp(param_name, "IS") == 0) {
            model->Is = param_value;
        } else if (strcmp(param_name, "n") == 0 || strcmp(param_name, "N") == 0) {
            model->n = param_value;
        } else if (strcmp(param_name, "BV") == 0 || strcmp(param_name, "bv") == 0) {
            model->BV = param_value;
        } else if (strcmp(param_name, "IBV") == 0 || strcmp(param_name, "ibv") == 0) {
            model->IBV = param_value;
        }
    }

    ctx->num_diode_models++;
    return MNA_SUCCESS;
}

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
        skip_whitespace(&p);

        int i = 0;
        while (*p && !isspace((unsigned char)*p) && *p != ';' && *p != '\n' && i < 63) {
            token[i++] = *p;
            p++;
        }
        token[i] = '\0';

        if (strlen(token) == 0) break;

        if (ctx->num_outputs >= MNA_NETLIST_MAX_OUTPUTS) {
            return report_error(ctx, line_num, "Too many output variables");
        }

        OutputVar* out = &ctx->outputs[ctx->num_outputs];

        if (tolower(token[0]) == 'v') {
            out->type = MNA_OUTPUT_VOLTAGE;

            char* paren = strchr(token, '(');
            if (paren) {
                out->node = atoi(paren + 1);
            } else {
                return report_error(ctx, line_num, "Invalid voltage specification");
            }
        } else if (tolower(token[0]) == 'i') {
            out->type = MNA_OUTPUT_CURRENT;

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

static MNAStatus parse_line(NetlistContext* ctx, const char* line, int line_num) {
    const char* p = line;

    skip_whitespace(&p);

    if (*p == '\0' || *p == '*' || *p == ';' || *p == '#') {
        return MNA_SUCCESS;
    }

    char clean_line[MNA_NETLIST_MAX_LINE];
    strncpy(clean_line, p, sizeof(clean_line) - 1);
    clean_line[sizeof(clean_line) - 1] = '\0';

    char* comment = strchr(clean_line, ';');
    if (comment) *comment = '\0';
    comment = strchr(clean_line, '\n');
    if (comment) *comment = '\0';

    int len = strlen(clean_line);
    while (len > 0 && isspace((unsigned char)clean_line[len - 1])) {
        clean_line[--len] = '\0';
    }

    if (strlen(clean_line) == 0) {
        return MNA_SUCCESS;
    }

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
        } else if (strcmp(cmd, "SWITCH_AT") == 0 || strcmp(cmd, "SW") == 0) {
            return parse_switch_cmd(ctx, cmd_p, line_num);
        } else if (strcmp(cmd, "MODEL") == 0) {
            return parse_model_cmd(ctx, cmd_p, line_num);
        } else if (strcmp(cmd, "END") == 0) {
            return MNA_SUCCESS;
        } else if (strcmp(cmd, "INCLUDE") == 0 || strcmp(cmd, "LIB") == 0) {
            return report_error(ctx, line_num, ".INCLUDE and .LIB not yet implemented");
        } else {
            return report_error(ctx, line_num, "Unknown command");
        }
    }

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
        case 'D':
            return parse_diode(ctx, clean_line, line_num);
        case 'T':
            return parse_transformer(ctx, clean_line, line_num);
        default:
            return report_error(ctx, line_num, "Unknown component type");
    }
}

MNAStatus netlist_init(NetlistContext* ctx, MNASolver* solver) {
    memset(ctx, 0, sizeof(NetlistContext));
    ctx->solver = solver;
    ctx->analysis.type = MNA_ANALYSIS_NONE;
    ctx->num_outputs = 0;
    ctx->write_enabled = 0;
    ctx->error_line = 0;
    ctx->error_msg[0] = '\0';
    ctx->num_diode_models = 0;

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

    if (ctx->analysis.tran_enabled || ctx->analysis.type == MNA_ANALYSIS_TRAN) {
        double dt = ctx->analysis.tran_dt;
        double t_end = ctx->analysis.tran_end;
        int steps = (int)(t_end / dt);

        for (int i = 0; i < ctx->analysis.num_switch_events; i++) {
            SwitchEvent* evt = &ctx->analysis.switch_events[i];
            evt->component_idx = -1;
            for (int j = 0; j < ctx->num_switches; j++) {
                if (strcmp(evt->name, ctx->switch_names[j]) == 0) {
                    evt->component_idx = ctx->switch_indices[j];
                    break;
                }
            }
            if (evt->component_idx < 0) {
                fprintf(stderr, "Warning: Switch '%s' not found for event at t=%.3f ms\n",
                        evt->name, evt->time * 1000.0);
            }
        }

        mna_init_transient(ctx->solver);

        printf("=== Transient Analysis ===\n");
        printf("Time step: %.2e s, End time: %.2e s, Steps: %d\n", dt, t_end, steps);
        if (ctx->analysis.num_switch_events > 0) {
            printf("Switch events: %d\n", ctx->analysis.num_switch_events);
        }

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

        int next_switch_event = 0;

        for (int step = 0; step <= steps; step++) {
            double t = step * dt;

            while (next_switch_event < ctx->analysis.num_switch_events &&
                   ctx->analysis.switch_events[next_switch_event].time <= t) {
                SwitchEvent* evt = &ctx->analysis.switch_events[next_switch_event];
                if (evt->valid && evt->component_idx >= 0) {
                    mna_set_switch_state(ctx->solver, evt->component_idx, evt->state);
                    printf("  [t=%.3f ms] Switch '%s' -> %s\n",
                           t * 1000.0, evt->name, evt->state ? "CLOSED" : "OPEN");
                }
                next_switch_event++;
            }

            if (ctx->analysis.sin_source_valid && step > 0) {
                double omega = 2.0 * M_PI * ctx->analysis.sin_freq;
                double peak = ctx->solver->components[ctx->analysis.sin_source_idx].ac_magnitude;
                double phase = ctx->solver->components[ctx->analysis.sin_source_idx].ac_phase;
                double v_sin = peak * sin(omega * t + phase);
                ctx->solver->components[ctx->analysis.sin_source_idx].value = v_sin;
            }

            if (step > 0) {
                status = mna_solve_transient_step(ctx->solver, dt);
                if (status != MNA_SUCCESS) {
                    fprintf(stderr, "Error: Transient step %d failed\n", step);
                    if (out_file) fclose(out_file);
                    return status;
                }
            }

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
    (void)ctx;
}
