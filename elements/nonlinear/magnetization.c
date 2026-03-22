#include "magnetization.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* Physical constants */
#define MU_0 (4.0 * M_PI * 1e-7)  /* Vacuum permeability [H/m] */

/* Internal helper functions */
static int compare_points(const void* a, const void* b);
static double interpolate_linear(const MagnetizationCurve* curve, double H, int* idx);
static double find_H_for_B_newton(const MagnetizationCurve* curve,
                                   double B_target, double H_initial);

/* Langevin function: L(x) = coth(x) - 1/x */
static double langevin(double x) {
    if (fabs(x) < 1e-6) {
        /* Taylor series for small x: L(x) ≈ x/3 - x³/45 */
        return x / 3.0 - (x * x * x) / 45.0;
    }
    if (x > 20.0) {
        /* For large positive x: coth(x) → 1, so L(x) → 1 - 1/x */
        return 1.0 - 1.0 / x;
    }
    if (x < -20.0) {
        /* For large negative x: coth(x) → -1, so L(x) → -1 - 1/x */
        return -1.0 - 1.0 / x;
    }
    return 1.0 / tanh(x) - 1.0 / x;
}

/* Derivative of Langevin function: L'(x) = 1/x² - 1/sinh²(x) */
static double langevin_deriv(double x) {
    if (fabs(x) < 1e-6) {
        /* Taylor series for small x: L'(x) ≈ 1/3 - x²/15 */
        return 1.0 / 3.0 - (x * x) / 15.0;
    }
    if (fabs(x) > 20.0) {
        /* For large |x|: L'(x) → 1/x² (always positive) */
        return 1.0 / (x * x);
    }
    double sinh_x = sinh(x);
    return 1.0 / (x * x) - 1.0 / (sinh_x * sinh_x);
}

/* 
 * Complete B-H model: B(H) = B_sat * L(H/H_c) + μ₀ * H
 * 
 * The μ₀*H term represents the vacuum/air path in parallel with the
 * ferromagnetic material. This is physically correct and prevents
 * numerical issues:
 * - B can exceed B_sat (following the air path)
 * - The inverse H(B) is always well-defined
 * - No artificial clamping needed
 */
static double B_from_H_complete(const MagnetizationCurve* curve, double H) {
    double B_sat = curve->langevin_params.B_sat;
    double H_c = curve->langevin_params.H_c;
    double x = H / H_c;
    
    /* Ferromagnetic contribution (Langevin) */
    double B_ferro = B_sat * langevin(x);
    
    /* Vacuum/air contribution */
    double B_vac = MU_0 * H;
    
    return B_ferro + B_vac;
}

/* Derivative: dB/dH = (B_sat/H_c) * L'(H/H_c) + μ₀ */
static double dB_dH_complete(const MagnetizationCurve* curve, double H) {
    double B_sat = curve->langevin_params.B_sat;
    double H_c = curve->langevin_params.H_c;
    double x = H / H_c;
    
    return (B_sat / H_c) * langevin_deriv(x) + MU_0;
}

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

    /* Initialize Langevin params */
    curve->langevin_params.B_sat = 0.0;
    curve->langevin_params.H_c = 0.0;
    curve->langevin_params.B_rem = 0.0;

    /* Initialize piecewise params */
    curve->pw_params.H_points = NULL;
    curve->pw_params.B_points = NULL;
    curve->pw_params.num_points = 0;
    curve->pw_params.allocated_size = 0;

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
 * Analytic Langevin Model Configuration
 * ============================================================================ */

MNAStatus mna_bh_curve_set_langevin_params(MagnetizationCurve* curve,
                                            double B_sat_T,
                                            double mu_r_initial,
                                            double H_c_A_m) {
    if (!curve) return MNA_INVALID_HANDLE;
    if (curve->model_type != BH_MODEL_ANALYTIC_LANGEVIN) {
        return MNA_INVALID_PARAMETER;
    }
    if (B_sat_T <= 0 || mu_r_initial <= 0) {
        return MNA_INVALID_PARAMETER;
    }

    curve->langevin_params.B_sat = B_sat_T;

    /* Auto-compute H_c from initial permeability
     * Initial slope: dB/dH|₀ = B_sat / (3*H_c) = μ₀ * μᵣ
     * So: H_c = B_sat / (3 * μ₀ * μᵣ)
     */
    if (H_c_A_m <= 0) {
        curve->langevin_params.H_c = B_sat_T / (3.0 * MU_0 * mu_r_initial);
    } else {
        curve->langevin_params.H_c = H_c_A_m;
    }

    return MNA_SUCCESS;
}

/* ============================================================================
 * Piecewise Model: Data Point Management
 * ============================================================================ */

MNAStatus mna_bh_curve_add_point(MagnetizationCurve* curve,
                                  double H_A_m, double B_Tesla) {
    if (!curve) return MNA_INVALID_HANDLE;
    if (curve->model_type != BH_MODEL_PIECEWISE_LINEAR) {
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

        if (!new_H || !new_B) {
            return MNA_INSUFFICIENT_MEMORY;
        }

        curve->pw_params.H_points = new_H;
        curve->pw_params.B_points = new_B;
        curve->pw_params.allocated_size = new_size;
    }

    /* Add point */
    int idx = curve->pw_params.num_points;
    curve->pw_params.H_points[idx] = H_A_m;
    curve->pw_params.B_points[idx] = B_Tesla;
    curve->pw_params.num_points++;

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

    return MNA_SUCCESS;
}

/* ============================================================================
 * Core Queries: B-H Relationship
 * ============================================================================ */

double mna_bh_curve_get_B(const MagnetizationCurve* curve, double H_A_m) {
    if (!curve) return 0.0;

    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_LANGEVIN: {
            /* Complete model: B(H) = B_sat * L(H/H_c) + μ₀ * H */
            return B_from_H_complete(curve, H_A_m);
        }

        case BH_MODEL_PIECEWISE_LINEAR: {
            if (curve->pw_params.num_points < 2) return 0.0;

            int idx = 0;
            return interpolate_linear(curve, H_A_m, &idx);
        }

        default:
            return 0.0;
    }
}

double mna_bh_curve_get_H(const MagnetizationCurve* curve, double B_Tesla) {
    if (!curve) return 0.0;

    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_LANGEVIN: {
            /* Inverse of B(H) = B_sat * L(H/H_c)
             * Use Newton-Raphson iteration
             */
            double H_initial = B_Tesla * 3.0 * curve->langevin_params.H_c / 
                               curve->langevin_params.B_sat;
            return find_H_for_B_newton(curve, B_Tesla, H_initial);
        }

        case BH_MODEL_PIECEWISE_LINEAR: {
            if (curve->pw_params.num_points < 2) return 0.0;

            /* Find H for given B by inverting the interpolation */
            int n = curve->pw_params.num_points;

            /* Handle saturation regions */
            if (B_Tesla <= curve->pw_params.B_points[0]) {
                double dB = curve->pw_params.B_points[1] - curve->pw_params.B_points[0];
                double dH = curve->pw_params.H_points[1] - curve->pw_params.H_points[0];
                if (dH == 0) return curve->pw_params.H_points[0];
                double slope = dH / dB;
                return curve->pw_params.H_points[0] + slope * (B_Tesla - curve->pw_params.B_points[0]);
            }
            if (B_Tesla >= curve->pw_params.B_points[n - 1]) {
                double dB = curve->pw_params.B_points[n-1] - curve->pw_params.B_points[n-2];
                double dH = curve->pw_params.H_points[n-1] - curve->pw_params.H_points[n-2];
                if (dH == 0) return curve->pw_params.H_points[n-1];
                double slope = dH / dB;
                return curve->pw_params.H_points[n-1] + slope * (B_Tesla - curve->pw_params.B_points[n-1]);
            }

            /* Find interval containing B */
            for (int i = 0; i < n - 1; i++) {
                if (B_Tesla >= curve->pw_params.B_points[i] &&
                    B_Tesla <= curve->pw_params.B_points[i + 1]) {
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
    if (!curve) return MU_0;

    switch (curve->model_type) {
        case BH_MODEL_ANALYTIC_LANGEVIN: {
            /* Complete model: dB/dH = (B_sat/H_c) * L'(H/H_c) + μ₀ */
            return dB_dH_complete(curve, H_A_m);
        }

        case BH_MODEL_PIECEWISE_LINEAR: {
            if (curve->pw_params.num_points < 2) return MU_0;

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

        default:
            return MU_0;
    }
}

/* ============================================================================
 * Circuit-Level Queries (with geometry conversion)
 * ============================================================================ */

double mna_bh_curve_get_current(const MagnetizationCurve* curve, double flux_Wb) {
    if (!curve || !curve->geometry_set) return 0.0;

    /* φ → B → H → i */
    double B = mna_bh_curve_B_from_flux(curve, flux_Wb);
    double H = mna_bh_curve_get_H(curve, B);
    return mna_bh_curve_current_from_H(curve, H);
}

double mna_bh_curve_get_flux(const MagnetizationCurve* curve, double current_A) {
    if (!curve || !curve->geometry_set) return 0.0;

    /* i → H → B → φ */
    double H = mna_bh_curve_H_from_current(curve, current_A);
    double B = mna_bh_curve_get_B(curve, H);
    return mna_bh_curve_flux_from_B(curve, B);
}

double mna_bh_curve_get_incremental_inductance(const MagnetizationCurve* curve,
                                                double flux_Wb) {
    if (!curve || !curve->geometry_set) return 1e-6;

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
        case BH_MODEL_ANALYTIC_LANGEVIN:
            return curve->langevin_params.B_sat > 0 &&
                   curve->langevin_params.H_c > 0;

        case BH_MODEL_PIECEWISE_LINEAR:
            return curve->pw_params.num_points >= 2;

        default:
            return false;
    }
}

int mna_bh_curve_get_point_count(const MagnetizationCurve* curve) {
    if (!curve) return 0;

    if (curve->model_type == BH_MODEL_PIECEWISE_LINEAR) {
        return curve->pw_params.num_points;
    }

    return 0;
}

MNAStatus mna_bh_curve_get_point(const MagnetizationCurve* curve,
                                  int index, double* H_out, double* B_out) {
    if (!curve || !H_out || !B_out) return MNA_INVALID_HANDLE;
    if (curve->model_type != BH_MODEL_PIECEWISE_LINEAR) {
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
        double slope = (curve->pw_params.B_points[1] - curve->pw_params.B_points[0]) /
                       (curve->pw_params.H_points[1] - curve->pw_params.H_points[0]);
        return curve->pw_params.B_points[0] + slope * (H - curve->pw_params.H_points[0]);
    }
    if (H >= curve->pw_params.H_points[n - 1]) {
        *idx = n - 2;
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

static double find_H_for_B_newton(const MagnetizationCurve* curve,
                                   double B_target, double H_initial) {
    /*
     * Newton-Raphson for complete model: B(H) = B_sat * L(H/H_c) + μ₀ * H
     *
     * Key insight: The μ₀*H term makes dB/dH >= μ₀ always, so:
     * - The function is strictly monotonic
     * - The inverse is always well-defined
     * - No numerical instability even for B >> B_sat
     */

    double B_sat = curve->langevin_params.B_sat;
    double H_c = curve->langevin_params.H_c;

    /* Special case: B = 0 => H = 0 */
    if (fabs(B_target) < 1e-12 * B_sat) {
        return 0.0;
    }

    /* Handle sign */
    int sign = (B_target >= 0) ? 1 : -1;
    double B_abs = fabs(B_target);

    /*
     * Initial guess strategy:
     * - For small B: use linear approximation (initial permeability region)
     * - For large B: use asymptotic formula with μ₀ term
     */
    double H;

    /*
     * Initial slope: dB/dH|₀ = B_sat/(3*H_c) + μ₀
     * For typical values (B_sat=1.6T, H_c=425 A/m, μᵣ=3000):
     *   B_sat/(3*H_c) ≈ 1.26 T·m/A = μ₀ * μᵣ
     *   μ₀ ≈ 1.26e-6 T·m/A (negligible in linear region)
     */
    double initial_slope = B_sat / (3.0 * H_c) + MU_0;

    if (B_abs < B_sat * 0.1) {
        /* Linear region: H ≈ B / initial_slope */
        H = B_abs / initial_slope;
    } else if (B_abs < B_sat * 0.9) {
        /* Transition region: use Langevin inverse approximation */
        double ratio = B_abs / B_sat;
        /* Approximate: L(x) ≈ ratio => x ≈ 3*ratio for small x */
        H = H_c * 3.0 * ratio;
    } else {
        /* Saturation region: B ≈ B_sat + μ₀*H */
        /* H ≈ (B - B_sat) / μ₀ */
        H = (B_abs - B_sat) / MU_0;
        if (H < H_c * 10) H = H_c * 10;  /* Minimum H in saturation */
    }

    /* Ensure H is positive (we're solving for |B| -> |H|) */
    H = fabs(H);
    
    /* Avoid H being exactly zero */
    if (H < 1e-10 * H_c) H = 1e-10 * H_c;

    /* Newton-Raphson iteration */
    const int max_iter = 50;
    const double tol = 1e-10;

    for (int i = 0; i < max_iter; i++) {
        double B = B_from_H_complete(curve, H);  /* Use positive H */
        double dB_dH = dB_dH_complete(curve, H);

        double residual = B - B_abs;  /* Both positive */

        /* Check convergence */
        if (fabs(residual) < tol * B_sat) break;
        if (fabs(residual) < tol * B_abs) break;  /* Relative tolerance */

        /* Newton step with damping */
        double delta = residual / dB_dH;

        /*
         * Damping: limit step to 50% of current H
         * This prevents overshoot in steep regions
         */
        double max_step = fabs(H) * 0.5 + H_c * 0.1;
        if (fabs(delta) > max_step) {
            delta = (delta > 0) ? max_step : -max_step;
        }

        H = H - delta;

        /* Keep H positive */
        if (H < 0) H = fabs(H) * 0.5;
        if (H < 1e-10 * H_c) H = 1e-10 * H_c;
    }

    return sign * H;
}
