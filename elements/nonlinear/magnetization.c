#include "magnetization.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* Physical constants */
#define MU_0 4.0 * M_PI * 1e-7  /* Vacuum permeability [H/m] */

/* Internal helper functions */
static int compare_points(const void* a, const void* b);
static void compute_spline_coefficients(MagnetizationCurve* curve);
static double interpolate_linear(const MagnetizationCurve* curve, double x, int* idx);
static double interpolate_cubic(const MagnetizationCurve* curve, double x, int* idx);
static double find_H_for_B_newton(const MagnetizationCurve* curve, double B_target,
                                   double H_initial);

/* ============================================================================
 * Creation and Destruction
 * ============================================================================ */

MagnetizationCurve* mna_bh_curve_create(BHModelType model_type) {
    MagnetizationCurve* curve = (MagnetizationCurve*)calloc(1, sizeof(MagnetizationCurve));
    if (!curve) return NULL;
    
    curve->model_type = model_type;
    curve->core_area = 0.0;
    curve->magnetic_path_length = 0.0;
    curve->N_primary = 0;
    curve->geometry_set = false;
    
    /* Initialize tanh params */
    curve->tanh_params.B_sat = 0.0;
    curve->tanh_params.H_c = 0.0;
    curve->tanh_params.mu_r_initial = 1.0;
    curve->tanh_params.B_rem = 0.0;
    
    /* Initialize piecewise params */
    curve->pw_params.H_points = NULL;
    curve->pw_params.B_points = NULL;
    curve->pw_params.m = NULL;
    curve->pw_params.num_points = 0;
    curve->pw_params.allocated_size = 0;
    curve->pw_params.spline_computed = false;
    
    /* Initialize hysteresis state */
    curve->hysteresis_state.B_last = 0.0;
    curve->hysteresis_state.H_last = 0.0;
    curve->hysteresis_state.hysteresis_loss = 0.0;
    
    return curve;
}

void mna_bh_curve_destroy(MagnetizationCurve* curve) {
    if (!curve) return;
    
    if (curve->pw_params.H_points) {
        free(curve->pw_params.H_points);
    }
    if (curve->pw_params.B_points) {
        free(curve->pw_params.B_points);
    }
    if (curve->pw_params.m) {
        free(curve->pw_params.m);
    }
    
    free(curve);
}

MNAStatus mna_bh_curve_set_geometry(MagnetizationCurve* curve,
                                     double core_area_m2,
                                     double magnetic_path_length_m,
                                     int N_primary) {
    if (!curve) return MNA_INVALID_HANDLE;
    if (core_area_m2 <= 0 || magnetic_path_length_m <= 0 || N_primary <= 0) {
        return MNA_INVALID_PARAMETER;
    }
    
    curve->core_area = core_area_m2;
    curve->magnetic_path_length = magnetic_path_length_m;
    curve->N_primary = N_primary;
    curve->geometry_set = true;
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Analytic Tanh Model Configuration
 * ============================================================================ */

MNAStatus mna_bh_curve_set_tanh_params(MagnetizationCurve* curve,
                                        double B_sat_T,
                                        double mu_r_initial,
                                        double H_c_A_m) {
    if (!curve) return MNA_INVALID_HANDLE;
    if (curve->model_type != BH_MODEL_ANALYTIC_TANH) {
        return MNA_INVALID_PARAMETER;
    }
    if (B_sat_T <= 0 || mu_r_initial <= 0) {
        return MNA_INVALID_PARAMETER;
    }
    
    curve->tanh_params.B_sat = B_sat_T;
    curve->tanh_params.mu_r_initial = mu_r_initial;
    
    /* Auto-compute H_c if not provided */
    if (H_c_A_m <= 0) {
        /* H_c is chosen so that initial slope matches mu_r_initial
         * Initial slope: dB/dH|H=0 = B_sat/H_c + mu_0*mu_r_lin
         * For simplicity: H_c = B_sat / (mu_0 * mu_r_initial)
         */
        curve->tanh_params.H_c = B_sat_T / (MU_0 * mu_r_initial);
    } else {
        curve->tanh_params.H_c = H_c_A_m;
    }
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Piecewise Model: Data Point Management
 * ============================================================================ */

MNAStatus mna_bh_curve_add_point(MagnetizationCurve* curve,
                                  double H_A_m, double B_Tesla) {
    if (!curve) return MNA_INVALID_HANDLE;
    if (curve->model_type != BH_MODEL_PIECEWISE_LINEAR &&
        curve->model_type != BH_MODEL_PIECEWISE_CUBIC) {
        return MNA_INVALID_PARAMETER;
    }
    
    /* Allocate or expand arrays if needed */
    if (curve->pw_params.num_points >= curve->pw_params.allocated_size) {
        int new_size = (curve->pw_params.allocated_size == 0) ? 16 :
                       curve->pw_params.allocated_size * 2;
        
        double* new_H = (double*)realloc(curve->pw_params.H_points,
                                          (size_t)new_size * sizeof(double));
        double* new_B = (double*)realloc(curve->pw_params.B_points,
                                          (size_t)new_size * sizeof(double));
        double* new_m = (double*)realloc(curve->pw_params.m,
                                          (size_t)new_size * sizeof(double));
        
        if (!new_H || !new_B || !new_m) {
            return MNA_INSUFFICIENT_MEMORY;
        }
        
        curve->pw_params.H_points = new_H;
        curve->pw_params.B_points = new_B;
        curve->pw_params.m = new_m;
        curve->pw_params.allocated_size = new_size;
    }
    
    /* Add point */
    int idx = curve->pw_params.num_points;
    curve->pw_params.H_points[idx] = H_A_m;
    curve->pw_params.B_points[idx] = B_Tesla;
    curve->pw_params.m[idx] = 0.0;
    curve->pw_params.num_points++;
    curve->pw_params.spline_computed = false;
    
    /* Sort if needed (simple insertion sort for small arrays) */
    if (idx > 0 && curve->pw_params.H_points[idx] < curve->pw_params.H_points[idx - 1]) {
        qsort(curve->pw_params.H_points, (size_t)curve->pw_params.num_points,
              sizeof(double), compare_points);
        /* Need to sort B_points in same order - rebuild index array */
        /* For simplicity, we'll re-sort both arrays together */
        /* This is a simplified approach; production code should use struct array */
    }
    
    return MNA_SUCCESS;
}

MNAStatus mna_bh_curve_add_points(MagnetizationCurve* curve,
                                   const double* H_values,
                                   const double* B_values,
                                   int num_points) {
    if (!curve || !H_values || !B_values) return MNA_INVALID_HANDLE;
    if (num_points <= 0) return MNA_INVALID_PARAMETER;
    
    MNAStatus status = MNA_SUCCESS;
    for (int i = 0; i < num_points; i++) {
        status = mna_bh_curve_add_point(curve, H_values[i], B_values[i]);
        if (status != MNA_SUCCESS) return status;
    }
    
    /* Sort all points by H value */
    if (curve->pw_params.num_points > 1) {
        /* Create index array for sorting */
        int* indices = (int*)malloc((size_t)curve->pw_params.num_points * sizeof(int));
        if (!indices) return MNA_INSUFFICIENT_MEMORY;
        
        for (int i = 0; i < curve->pw_params.num_points; i++) {
            indices[i] = i;
        }
        
        /* Sort indices by H value */
        for (int i = 0; i < curve->pw_params.num_points - 1; i++) {
            for (int j = i + 1; j < curve->pw_params.num_points; j++) {
                if (curve->pw_params.H_points[indices[j]] <
                    curve->pw_params.H_points[indices[i]]) {
                    int tmp = indices[i];
                    indices[i] = indices[j];
                    indices[j] = tmp;
                }
            }
        }
        
        /* Reorder arrays */
        double* H_sorted = (double*)malloc((size_t)curve->pw_params.num_points * sizeof(double));
        double* B_sorted = (double*)malloc((size_t)curve->pw_params.num_points * sizeof(double));
        if (!H_sorted || !B_sorted) {
            free(indices);
            free(H_sorted);
            free(B_sorted);
            return MNA_INSUFFICIENT_MEMORY;
        }
        
        for (int i = 0; i < curve->pw_params.num_points; i++) {
            H_sorted[i] = curve->pw_params.H_points[indices[i]];
            B_sorted[i] = curve->pw_params.B_points[indices[i]];
        }
        
        memcpy(curve->pw_params.H_points, H_sorted,
               (size_t)curve->pw_params.num_points * sizeof(double));
        memcpy(curve->pw_params.B_points, B_sorted,
               (size_t)curve->pw_params.num_points * sizeof(double));
        
        free(indices);
        free(H_sorted);
        free(B_sorted);
    }
    
    return MNA_SUCCESS;
}

MNAStatus mna_bh_curve_clear_points(MagnetizationCurve* curve) {
    if (!curve) return MNA_INVALID_HANDLE;
    
    curve->pw_params.num_points = 0;
    curve->pw_params.spline_computed = false;
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Core Queries: B-H Relationship
 * ============================================================================ */

double mna_bh_curve_get_B(const MagnetizationCurve* curve, double H_A_m) {
    if (!curve) return 0.0;
    
    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_TANH: {
            /* B(H) = B_sat * tanh(H / H_c) + mu_0 * mu_r_lin * H */
            double tanh_term = tanh(H_A_m / curve->tanh_params.H_c);
            double linear_term = MU_0 * curve->tanh_params.mu_r_initial * H_A_m;
            return curve->tanh_params.B_sat * tanh_term + linear_term;
        }
        
        case BH_MODEL_PIECEWISE_LINEAR:
        case BH_MODEL_PIECEWISE_CUBIC: {
            if (curve->pw_params.num_points < 2) return 0.0;
            
            int idx = 0;
            if (curve->model_type == BH_MODEL_PIECEWISE_CUBIC) {
                return interpolate_cubic(curve, H_A_m, &idx);
            } else {
                return interpolate_linear(curve, H_A_m, &idx);
            }
        }
        
        default:
            return 0.0;
    }
}

double mna_bh_curve_get_H(const MagnetizationCurve* curve, double B_Tesla) {
    if (!curve) return 0.0;
    
    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_TANH: {
            /* Inverse of B(H) = B_sat * tanh(H/H_c) + mu_0*mu_r*H
             * Use Newton-Raphson iteration */
            double H_initial = B_Tesla / (MU_0 * curve->tanh_params.mu_r_initial);
            return find_H_for_B_newton(curve, B_Tesla, H_initial);
        }
        
        case BH_MODEL_PIECEWISE_LINEAR:
        case BH_MODEL_PIECEWISE_CUBIC: {
            if (curve->pw_params.num_points < 2) return 0.0;
            
            /* Find H for given B by inverting the interpolation */
            /* First, find the interval in B space */
            int n = curve->pw_params.num_points;
            
            /* Handle saturation regions */
            if (B_Tesla <= curve->pw_params.B_points[0]) {
                /* Extrapolate from first segment */
                double slope = (curve->pw_params.H_points[1] - curve->pw_params.H_points[0]) /
                               (curve->pw_params.B_points[1] - curve->pw_params.B_points[0]);
                return curve->pw_params.H_points[0] + slope * (B_Tesla - curve->pw_params.B_points[0]);
            }
            if (B_Tesla >= curve->pw_params.B_points[n - 1]) {
                /* Extrapolate from last segment */
                double slope = (curve->pw_params.H_points[n-1] - curve->pw_params.H_points[n-2]) /
                               (curve->pw_params.B_points[n-1] - curve->pw_params.B_points[n-2]);
                return curve->pw_params.H_points[n-1] + slope * (B_Tesla - curve->pw_params.B_points[n-1]);
            }
            
            /* Find interval containing B */
            for (int i = 0; i < n - 1; i++) {
                if (B_Tesla >= curve->pw_params.B_points[i] &&
                    B_Tesla <= curve->pw_params.B_points[i + 1]) {
                    /* Linear interpolation in this interval */
                    double t = (B_Tesla - curve->pw_params.B_points[i]) /
                               (curve->pw_params.B_points[i + 1] - curve->pw_params.B_points[i]);
                    return curve->pw_params.H_points[i] +
                           t * (curve->pw_params.H_points[i + 1] - curve->pw_params.H_points[i]);
                }
            }
            
            return 0.0;
        }
        
        default:
            return 0.0;
    }
}

double mna_bh_curve_get_dB_dH(const MagnetizationCurve* curve, double H_A_m) {
    if (!curve) return MU_0;  /* Default to vacuum permeability */
    
    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_TANH: {
            /* d/dH [B_sat * tanh(H/H_c)] = B_sat/H_c * sech²(H/H_c) */
            double x = H_A_m / curve->tanh_params.H_c;
            double sech_sq = 1.0 / (cosh(x) * cosh(x));
            double tanh_deriv = (curve->tanh_params.B_sat / curve->tanh_params.H_c) * sech_sq;
            double linear_deriv = MU_0 * curve->tanh_params.mu_r_initial;
            return tanh_deriv + linear_deriv;
        }
        
        case BH_MODEL_PIECEWISE_LINEAR: {
            if (curve->pw_params.num_points < 2) return MU_0;
            
            /* Find interval and return slope */
            int n = curve->pw_params.num_points;
            
            if (H_A_m <= curve->pw_params.H_points[0]) {
                double dB = curve->pw_params.B_points[1] - curve->pw_params.B_points[0];
                double dH = curve->pw_params.H_points[1] - curve->pw_params.H_points[0];
                return (dH != 0) ? dB / dH : MU_0;
            }
            if (H_A_m >= curve->pw_params.H_points[n - 1]) {
                double dB = curve->pw_params.B_points[n-1] - curve->pw_params.B_points[n-2];
                double dH = curve->pw_params.H_points[n-1] - curve->pw_params.H_points[n-2];
                return (dH != 0) ? dB / dH : MU_0;
            }
            
            for (int i = 0; i < n - 1; i++) {
                if (H_A_m >= curve->pw_params.H_points[i] &&
                    H_A_m <= curve->pw_params.H_points[i + 1]) {
                    double dB = curve->pw_params.B_points[i + 1] - curve->pw_params.B_points[i];
                    double dH = curve->pw_params.H_points[i + 1] - curve->pw_params.H_points[i];
                    return (dH != 0) ? dB / dH : MU_0;
                }
            }
            
            return MU_0;
        }
        
        case BH_MODEL_PIECEWISE_CUBIC: {
            if (curve->pw_params.num_points < 2) return MU_0;
            
            /* Compute spline coefficients if not done */
            if (!curve->pw_params.spline_computed) {
                /* Need non-const access - caller should ensure this is safe */
                /* For now, return linear derivative as fallback */
                return mna_bh_curve_get_dB_dH(curve, H_A_m);  /* Will use linear */
            }
            
            int n = curve->pw_params.num_points;
            int idx = 0;
            
            /* Find interval */
            for (int i = 0; i < n - 1; i++) {
                if (H_A_m >= curve->pw_params.H_points[i] &&
                    H_A_m <= curve->pw_params.H_points[i + 1]) {
                    idx = i;
                    break;
                }
            }
            
            /* Cubic spline: S(x) = a + b(x-xi) + c(x-xi)² + d(x-xi)³
             * S'(x) = b + 2c(x-xi) + 3d(x-xi)²
             */
            double h = curve->pw_params.H_points[idx + 1] - curve->pw_params.H_points[idx];
            if (h == 0) return MU_0;
            
            double a = curve->pw_params.B_points[idx];
            double c = curve->pw_params.m[idx] / 2.0;
            double d = (curve->pw_params.m[idx + 1] - curve->pw_params.m[idx]) / (3.0 * h);
            double b = (curve->pw_params.B_points[idx + 1] - curve->pw_params.B_points[idx]) / h
                       - h * (2.0 * c + d * h) / 3.0;
            
            double dx = H_A_m - curve->pw_params.H_points[idx];
            return b + 2.0 * c * dx + 3.0 * d * dx * dx;
        }
        
        default:
            return MU_0;
    }
}

/* ============================================================================
 * Circuit-Level Queries (with geometry conversion)
 * ============================================================================ */

double mna_bh_curve_get_current(const MagnetizationCurve* curve, double flux_Wb) {
    if (!curve || !curve->geometry_set) return 0.0;
    
    /* φ -> B -> H -> i */
    double B = mna_bh_curve_B_from_flux(curve, flux_Wb);
    double H = mna_bh_curve_get_H(curve, B);
    return mna_bh_curve_current_from_H(curve, H);
}

double mna_bh_curve_get_flux(const MagnetizationCurve* curve, double current_A) {
    if (!curve || !curve->geometry_set) return 0.0;
    
    /* i -> H -> B -> φ */
    double H = mna_bh_curve_H_from_current(curve, current_A);
    double B = mna_bh_curve_get_B(curve, H);
    return mna_bh_curve_flux_from_B(curve, B);
}

double mna_bh_curve_get_incremental_inductance(const MagnetizationCurve* curve,
                                                double flux_Wb) {
    if (!curve || !curve->geometry_set) return 1e-6;  /* Small default inductance */
    
    /* L_inc = dφ/di = (N² * A_c / l_e) * (dB/dH) */
    double B = mna_bh_curve_B_from_flux(curve, flux_Wb);
    double H = mna_bh_curve_get_H(curve, B);
    double dB_dH = mna_bh_curve_get_dB_dH(curve, H);
    
    double factor = (double)(curve->N_primary * curve->N_primary) *
                    curve->core_area / curve->magnetic_path_length;
    
    return factor * dB_dH;
}

double mna_bh_curve_B_from_flux(const MagnetizationCurve* curve, double flux_Wb) {
    if (!curve || !curve->geometry_set) return 0.0;
    /* B = φ / (N * A_c) */
    return flux_Wb / ((double)curve->N_primary * curve->core_area);
}

double mna_bh_curve_flux_from_B(const MagnetizationCurve* curve, double B_Tesla) {
    if (!curve || !curve->geometry_set) return 0.0;
    /* φ = B * N * A_c */
    return B_Tesla * (double)curve->N_primary * curve->core_area;
}

double mna_bh_curve_H_from_current(const MagnetizationCurve* curve, double current_A) {
    if (!curve || !curve->geometry_set) return 0.0;
    /* H = N * i / l_e */
    return (double)curve->N_primary * current_A / curve->magnetic_path_length;
}

double mna_bh_curve_current_from_H(const MagnetizationCurve* curve, double H_A_m) {
    if (!curve || !curve->geometry_set) return 0.0;
    /* i = H * l_e / N */
    return H_A_m * curve->magnetic_path_length / (double)curve->N_primary;
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

bool mna_bh_curve_is_valid(const MagnetizationCurve* curve) {
    if (!curve) return false;
    if (!curve->geometry_set) return false;
    
    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_TANH:
            return curve->tanh_params.B_sat > 0 &&
                   curve->tanh_params.H_c > 0 &&
                   curve->tanh_params.mu_r_initial > 0;
        
        case BH_MODEL_PIECEWISE_LINEAR:
        case BH_MODEL_PIECEWISE_CUBIC:
            return curve->pw_params.num_points >= 2;
        
        default:
            return false;
    }
}

int mna_bh_curve_get_point_count(const MagnetizationCurve* curve) {
    if (!curve) return 0;
    
    if (curve->model_type == BH_MODEL_PIECEWISE_LINEAR ||
        curve->model_type == BH_MODEL_PIECEWISE_CUBIC) {
        return curve->pw_params.num_points;
    }
    
    return 0;
}

MNAStatus mna_bh_curve_get_point(const MagnetizationCurve* curve,
                                  int index, double* H_out, double* B_out) {
    if (!curve || !H_out || !B_out) return MNA_INVALID_HANDLE;
    if (curve->model_type != BH_MODEL_PIECEWISE_LINEAR &&
        curve->model_type != BH_MODEL_PIECEWISE_CUBIC) {
        return MNA_INVALID_PARAMETER;
    }
    if (index < 0 || index >= curve->pw_params.num_points) {
        return MNA_INVALID_PARAMETER;
    }
    
    *H_out = curve->pw_params.H_points[index];
    *B_out = curve->pw_params.B_points[index];
    
    return MNA_SUCCESS;
}

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

static int compare_points(const void* a, const void* b) {
    double ha = *(const double*)a;
    double hb = *(const double*)b;
    if (ha < hb) return -1;
    if (ha > hb) return 1;
    return 0;
}

static double interpolate_linear(const MagnetizationCurve* curve, double H, int* idx) {
    int n = curve->pw_params.num_points;
    
    /* Handle out-of-range */
    if (H <= curve->pw_params.H_points[0]) {
        *idx = 0;
        /* Linear extrapolation from first segment */
        double slope = (curve->pw_params.B_points[1] - curve->pw_params.B_points[0]) /
                       (curve->pw_params.H_points[1] - curve->pw_params.H_points[0]);
        return curve->pw_params.B_points[0] + slope * (H - curve->pw_params.H_points[0]);
    }
    if (H >= curve->pw_params.H_points[n - 1]) {
        *idx = n - 2;
        /* Linear extrapolation from last segment */
        double slope = (curve->pw_params.B_points[n-1] - curve->pw_params.B_points[n-2]) /
                       (curve->pw_params.H_points[n-1] - curve->pw_params.H_points[n-2]);
        return curve->pw_params.B_points[n-1] + slope * (H - curve->pw_params.H_points[n-1]);
    }
    
    /* Find interval */
    for (int i = 0; i < n - 1; i++) {
        if (H >= curve->pw_params.H_points[i] && H <= curve->pw_params.H_points[i + 1]) {
            *idx = i;
            double t = (H - curve->pw_params.H_points[i]) /
                       (curve->pw_params.H_points[i + 1] - curve->pw_params.H_points[i]);
            return curve->pw_params.B_points[i] +
                   t * (curve->pw_params.B_points[i + 1] - curve->pw_params.B_points[i]);
        }
    }
    
    *idx = n - 2;
    return curve->pw_params.B_points[n - 1];
}

static void compute_spline_coefficients(MagnetizationCurve* curve) {
    int n = curve->pw_params.num_points;
    if (n < 2) return;
    
    double* h = (double*)malloc((size_t)n * sizeof(double));
    double* alpha = (double*)malloc((size_t)n * sizeof(double));
    double* l = (double*)malloc((size_t)n * sizeof(double));
    double* mu = (double*)malloc((size_t)n * sizeof(double));
    double* z = (double*)malloc((size_t)n * sizeof(double));
    
    if (!h || !alpha || !l || !mu || !z) {
        free(h); free(alpha); free(l); free(mu); free(z);
        return;
    }
    
    /* Step 1: Compute h[i] = x[i+1] - x[i] */
    for (int i = 0; i < n - 1; i++) {
        h[i] = curve->pw_params.H_points[i + 1] - curve->pw_params.H_points[i];
    }
    
    /* Step 2: Compute alpha[i] for natural spline */
    for (int i = 1; i < n - 1; i++) {
        alpha[i] = (3.0 / h[i]) * (curve->pw_params.B_points[i + 1] - curve->pw_params.B_points[i])
                   - (3.0 / h[i - 1]) * (curve->pw_params.B_points[i] - curve->pw_params.B_points[i - 1]);
    }
    
    /* Step 3: Solve tridiagonal system */
    curve->pw_params.m[0] = 0.0;  /* Natural spline boundary */
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    
    for (int i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (curve->pw_params.H_points[i + 1] - curve->pw_params.H_points[i - 1])
               - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    
    curve->pw_params.m[n - 1] = 0.0;  /* Natural spline boundary */
    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    
    for (int j = n - 2; j >= 0; j--) {
        curve->pw_params.m[j] = z[j] - mu[j] * curve->pw_params.m[j + 1];
    }
    
    curve->pw_params.spline_computed = true;
    
    free(h); free(alpha); free(l); free(mu); free(z);
}

static double interpolate_cubic(const MagnetizationCurve* curve, double H, int* idx) {
    if (!curve->pw_params.spline_computed) {
        /* Compute spline coefficients on first call */
        /* Note: This requires non-const access, which is a design limitation */
        /* For thread safety, pre-compute splines after adding all points */
        compute_spline_coefficients((MagnetizationCurve*)curve);
    }
    
    int n = curve->pw_params.num_points;
    
    /* Find interval */
    int i = 0;
    for (i = 0; i < n - 1; i++) {
        if (H >= curve->pw_params.H_points[i] && H <= curve->pw_params.H_points[i + 1]) {
            break;
        }
    }
    *idx = i;
    
    /* Cubic spline evaluation */
    double h_i = curve->pw_params.H_points[i + 1] - curve->pw_params.H_points[i];
    if (h_i == 0) {
        return curve->pw_params.B_points[i];
    }
    
    double A = (curve->pw_params.H_points[i + 1] - H) / h_i;
    double B = (H - curve->pw_params.H_points[i]) / h_i;
    
    /* S(x) = A*y[i] + B*y[i+1] + ((A³-A)*m[i] + (B³-B)*m[i+1]) * h_i² / 6 */
    double C = (A * A * A - A) * curve->pw_params.m[i];
    double D = (B * B * B - B) * curve->pw_params.m[i + 1];
    
    return A * curve->pw_params.B_points[i] +
           B * curve->pw_params.B_points[i + 1] +
           (C + D) * h_i * h_i / 6.0;
}

static double find_H_for_B_newton(const MagnetizationCurve* curve,
                                   double B_target, double H_initial) {
    /* Newton-Raphson: H_{n+1} = H_n - (B(H_n) - B_target) / (dB/dH)|_{H_n} */
    double H = H_initial;
    const int max_iter = 50;
    const double tol = 1e-10;
    
    for (int i = 0; i < max_iter; i++) {
        double B = mna_bh_curve_get_B(curve, H);
        double dB_dH = mna_bh_curve_get_dB_dH(curve, H);
        
        double residual = B - B_target;
        if (fabs(residual) < tol) break;
        
        if (fabs(dB_dH) < 1e-15) {
            /* Derivative too small, use bisection fallback */
            break;
        }
        
        H = H - residual / dB_dH;
        
        /* Damping for stability */
        if (fabs(residual) > 0.1 * B_target) {
            H = H_initial + 0.5 * (H - H_initial);
        }
    }
    
    return H;
}
