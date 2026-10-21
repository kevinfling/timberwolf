// timberwolf.h
// Tail-Integrated Multi-band EMD with Robust Spectral Health Index (ASHI)
// Pure C11 — requires FFTW3 and CHEAPDAN in include/link path
// Compile: -std=c11 -O3 -march=native -ffast-math -I<includes> -lfftw3

#ifndef TIMBERWOLF_H
#define TIMBERWOLF_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// =============================================================================
// CONFIGURATION & CONSTANTS
// =============================================================================

#ifndef TW_MAX_IMFS
#define TW_MAX_IMFS 12          // CHEAPDAN hard-caps 1D CEEMDAN at 12 IMFs
#endif

// Default ensemble size for CEEMDAN. Reduce for large signals or tight latency budgets.
#ifndef TW_CEEMDAN_ENSEMBLE
#define TW_CEEMDAN_ENSEMBLE 20
#endif

// EWMA smoothing factor for inter-window leptokurtosis tracking in tw_ctx_t.
// 0.8 retains 80% of the previous window's estimate and blends in 20% of the new
// batch measurement — equivalent to a half-life of ~3 windows.
#ifndef TW_LEPTO_EWMA_ALPHA
#define TW_LEPTO_EWMA_ALPHA 0.8
#endif

#ifndef TW_ALIGN_BYTES
#define TW_ALIGN_BYTES 64       // Cache line alignment for SIMD
#endif

#define TW_EPSILON 1e-12
#define TW_PI 3.14159265358979323846

// Tukey bisquare constant (typically 4.685 for 95% efficiency at Gaussian)
#define TW_TUKEY_C 4.685

// =============================================================================
// PLATFORM ABSTRACTION
// =============================================================================

#if defined(__GNUC__) || defined(__clang__)
    #define TW_INLINE static inline __attribute__((always_inline))
    #define TW_ALIGN(x) __attribute__((aligned(x)))
#elif defined(_MSC_VER)
    #define TW_INLINE static inline __forceinline
    #define TW_ALIGN(x) __declspec(align(x))
#else
    #define TW_INLINE static inline
    #define TW_ALIGN(x)
#endif

#define TW_RESTRICT restrict

// =============================================================================
// 1. ROBUST STATISTICS CORE (Type 7 Quantiles + Tukey Bisquare)
// =============================================================================

/**
 * @brief Type 7 Quantile (Linear interpolation)
 *
 * Mathematical definition:
 *   Q(p) = (1-γ)·x_{(j)} + γ·x_{(j+1)}
 *   where j = floor((N-1)·p), γ = (N-1)·p - j
 *
 * Input:
 *   @data   : Sorted array of length n (ascending)
 *   @n      : Length of array (n >= 1)
 *   @p      : Quantile probability [0, 1]
 *
 * Output:
 *   Returns the Type 7 quantile value
 */
TW_INLINE double tw_type7_quantile(const double* TW_RESTRICT data, size_t n, double p) {
    if (n == 0) return 0.0;
    if (n == 1) return data[0];
    if (p <= 0.0) return data[0];
    if (p >= 1.0) return data[n-1];

    double h = (n - 1) * p;
    size_t j = (size_t)floor(h);
    double gamma = h - (double)j;

    if (j >= n - 1) return data[n-1];

    return (1.0 - gamma) * data[j] + gamma * data[j + 1];
}

/**
 * @brief QuickSelect (Hoare's algorithm) - O(n) median/quantile
 *
 * Partitions array such that k-th smallest is at index k.
 * Input array is modified (reordered) but not fully sorted.
 *
 * Input:
 *   @arr : Array of doubles (modified in place)
 *   @n   : Length
 *   @k   : Index to select (0-based)
 *
 * Output:
 *   Returns value at index k if array were sorted
 */
TW_INLINE double tw_quickselect(double* TW_RESTRICT arr, size_t n, size_t k) {
    size_t left = 0, right = n - 1;

    while (left < right) {
        /* Median-of-three pivot to avoid O(n²) on sorted input */
        size_t mid = left + (right - left) / 2;
        if (arr[mid] < arr[left])  { double t = arr[left];  arr[left]  = arr[mid];   arr[mid]   = t; }
        if (arr[right] < arr[left]){ double t = arr[left];  arr[left]  = arr[right]; arr[right] = t; }
        if (arr[right] < arr[mid]) { double t = arr[mid];   arr[mid]   = arr[right]; arr[right] = t; }

        double pivot = arr[mid];
        /* Dutch-National-Flag three-way partition (signed indices to avoid size_t underflow) */
        ptrdiff_t lo = (ptrdiff_t)left, eq = (ptrdiff_t)left, hi = (ptrdiff_t)right;
        while (eq <= hi) {
            if      (arr[eq] < pivot) { double t = arr[lo]; arr[lo++] = arr[eq]; arr[eq++] = t; }
            else if (arr[eq] > pivot) { double t = arr[eq]; arr[eq]   = arr[hi]; arr[hi--] = t; }
            else                      { ++eq; }
        }
        /* arr[left..lo-1] < pivot, arr[lo..hi] == pivot, arr[hi+1..right] > pivot */
        if      ((ptrdiff_t)k < lo) right = (size_t)(lo - 1);
        else if ((ptrdiff_t)k > hi) left  = (size_t)(hi + 1);
        else                        break; /* k in the equal-pivot range */
    }
    return arr[k];
}

/**
 * @brief Tukey Bisquare (Biweight) Weight Function
 *
 * w(u) = [1 - (u/c)²]²  for |u| ≤ c
 * w(u) = 0              for |u| > c
 *
 * Input:
 *   @u : Standardized residual (u = x / scale)
 *   @c : Tuning constant (default TW_TUKEY_C ≈ 4.685 for 95% Gaussian efficiency)
 *
 * Output:
 *   Returns weight in [0, 1]
 */
TW_INLINE double tw_tukey_weight(double u, double c) {
    double abs_u = fabs(u);
    if (abs_u >= c) return 0.0;
    double tmp = 1.0 - (abs_u / c) * (abs_u / c);
    return tmp * tmp;
}

/**
 * @brief Robust MAD (Median Absolute Deviation) with Type 7 quantile
 *
 * MAD = median(|x_i - median(x)|)
 *
 * Input:
 *   @x      : Data array (length n)
 *   @n      : Length
 *   @scratch: Pre-allocated scratch buffer (size >= n) for temp storage
 *
 * Output:
 *   Returns MAD, or 0 if n == 0
 */
TW_INLINE double tw_robust_mad(const double* TW_RESTRICT x, size_t n,
                               double* TW_RESTRICT scratch) {
    if (n == 0) return 0.0;

    memcpy(scratch, x, n * sizeof(double));

    double med = tw_quickselect(scratch, n, n / 2);

    for (size_t i = 0; i < n; ++i) scratch[i] = fabs(x[i] - med);

    return tw_quickselect(scratch, n, n / 2);
}

// Comparator for qsort — used by tw_moors_kurtosis and tw_bowley_skew
static inline int tw_compare_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

/**
 * @brief Moors Kurtosis (Octile-based, robust to tails)
 *
 * K_m = [(Q(7/8) - Q(5/8)) + (Q(3/8) - Q(1/8))] / (Q(3/4) - Q(1/4))
 *
 * Input:
 *   @x      : Data array
 *   @n      : Length (n >= 8 recommended)
 *   @scratch: Scratch buffer (size >= n)
 *
 * Output:
 *   Returns Moors kurtosis (Gaussian ≈ 1.233)
 */
TW_INLINE double tw_moors_kurtosis(const double* TW_RESTRICT x, size_t n,
                                   double* TW_RESTRICT scratch) {
    if (n < 8) return 1.233;

    memcpy(scratch, x, n * sizeof(double));
    qsort(scratch, n, sizeof(double), tw_compare_double);

    double q125 = tw_type7_quantile(scratch, n, 0.125);
    double q375 = tw_type7_quantile(scratch, n, 0.375);
    double q625 = tw_type7_quantile(scratch, n, 0.625);
    double q875 = tw_type7_quantile(scratch, n, 0.875);
    double q25  = tw_type7_quantile(scratch, n, 0.25);
    double q75  = tw_type7_quantile(scratch, n, 0.75);

    double denom = q75 - q25;
    return ((q875 - q625) + (q375 - q125)) / (denom + TW_EPSILON);
}

/**
 * @brief Bowley Skewness (Quartile-based, bounded [-1, 1])
 *
 * S_b = (Q(0.25) + Q(0.75) - 2·Q(0.5)) / (Q(0.75) - Q(0.25))
 */
TW_INLINE double tw_bowley_skew(const double* TW_RESTRICT x, size_t n,
                                double* TW_RESTRICT scratch) {
    if (n < 4) return 0.0;

    memcpy(scratch, x, n * sizeof(double));
    qsort(scratch, n, sizeof(double), tw_compare_double);

    double q1 = tw_type7_quantile(scratch, n, 0.25);
    double q2 = tw_type7_quantile(scratch, n, 0.50);
    double q3 = tw_type7_quantile(scratch, n, 0.75);

    return (q3 + q1 - 2.0 * q2) / (q3 - q1 + TW_EPSILON);
}

// =============================================================================
// 2. THREAD-SAFE RNG (xorshift128+)
// =============================================================================

typedef struct {
    uint64_t s[2];
} tw_rng_t;

TW_INLINE void tw_rng_seed(tw_rng_t* rng, uint64_t seed) {
    rng->s[0] = seed;
    rng->s[1] = seed + 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 10; ++i) {
        uint64_t x = rng->s[0];
        uint64_t y = rng->s[1];
        rng->s[0] = y;
        x ^= x << 23;
        rng->s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
    }
}

TW_INLINE uint64_t tw_rng_u64(tw_rng_t* rng) {
    uint64_t x = rng->s[0];
    uint64_t y = rng->s[1];
    rng->s[0] = y;
    x ^= x << 23;
    rng->s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
    return rng->s[1] + y;
}

TW_INLINE double tw_rng_uniform(tw_rng_t* rng, double a, double b) {
    uint64_t r = tw_rng_u64(rng);
    double u = (r >> 11) * (1.0 / (1ULL << 53));
    return a + u * (b - a);
}

// =============================================================================
// 3. SPECTRAL ANALYSIS (DCT-II via FFTW3 C API)
// =============================================================================

#include <fftw3.h>

/**
 * @brief Compute unitary DCT-II via FFTW3 FFTW_REDFT10 plan (O(N log N))
 *
 * FFTW_REDFT10 computes:
 *   Y_k = 2 * sum_{n=0}^{N-1} x_n * cos[pi/N * (n + 0.5) * k]
 *
 * Unitary normalization applied here:
 *   Y_0    *= 1/sqrt(4N)
 *   Y_{k>0} *= 1/sqrt(2N)
 *
 * Input:
 *   @in  : Real input (length n) — contents overwritten
 *   @out : DCT-II coefficients (length n); may alias @in
 *   @n   : Transform length
 */
TW_INLINE void tw_dct_ii(double* TW_RESTRICT in, double* TW_RESTRICT out, size_t n) {
    if (n == 0) return;
    fftw_plan p = fftw_plan_r2r_1d((int)n, in, out, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    double scale_dc = 1.0 / sqrt(4.0 * (double)n);
    double scale_ac = 1.0 / sqrt(2.0 * (double)n);
    out[0] *= scale_dc;
    for (size_t k = 1; k < n; ++k) out[k] *= scale_ac;
}

/**
 * @brief Adaptive Spectral Health Index (ASHI)
 *
 * Computes the ratio of high-frequency to low-frequency energy in DCT domain,
 * with adaptive binning (cumulative energy percentiles) and Tukey robustification.
 *
 * Mathematical definition:
 *   ρ = log(1 + E_high / (E_low + ε))
 *   where E_high = energy above 80th percentile frequency
 *         E_low  = energy below 20th percentile frequency
 *         ε      = η * σ²_noise (regularization)
 *
 * Input:
 *   @signal     : Input signal (length n) — copied internally, not modified
 *   @n          : Length (power of 2 recommended)
 *   @dct_buf    : Pre-allocated buffer (size >= n) used for DCT computation
 *   @scratch    : General scratch buffer (size >= n)
 *   @noise_floor: Estimated noise variance (from MAD²)
 *
 * Output:
 *   Returns ASHI value (0 = healthy/low-freq dominant, >0 = anomalous/high-freq)
 */
TW_INLINE double tw_compute_ashi(const double* TW_RESTRICT signal, size_t n,
                                 double* TW_RESTRICT dct_buf,
                                 double* TW_RESTRICT scratch,
                                 double noise_floor) {
    // Copy signal so tw_dct_ii can operate in-place without modifying caller's buffer
    memcpy(dct_buf, signal, n * sizeof(double));
    tw_dct_ii(dct_buf, dct_buf, n);

    // Compute energy spectrum (squared magnitudes)
    double* energy = scratch;
    double robust_energy_sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        energy[i] = dct_buf[i] * dct_buf[i];
    }

    // Robustify energy with Tukey bisquare (down-weight outlier frequency spikes)
    double med_energy = tw_quickselect(energy, n, n / 2);
    // Recompute deviations (quickselect reorders energy in-place; recalc from dct_buf)
    for (size_t i = 0; i < n; ++i) {
        energy[i] = dct_buf[i] * dct_buf[i];
        dct_buf[i] = fabs(energy[i] - med_energy);
    }
    double mad_energy = tw_quickselect(dct_buf, n, n / 2);

    double c = TW_TUKEY_C * mad_energy;
    for (size_t i = 0; i < n; ++i) {
        double w = tw_tukey_weight(energy[i] - med_energy, c);
        energy[i] *= w;
        robust_energy_sum += energy[i];
    }

    // Adaptive binning via cumulative energy
    // Reuse dct_buf for cumulative sum
    dct_buf[0] = energy[0];
    for (size_t i = 1; i < n; ++i) dct_buf[i] = dct_buf[i-1] + energy[i];
    for (size_t i = 0; i < n; ++i) dct_buf[i] /= (robust_energy_sum + TW_EPSILON);

    size_t idx_20 = 0, idx_80 = n - 1;
    for (size_t i = 0; i < n; ++i) {
        if (dct_buf[i] <= 0.20) idx_20 = i;
        if (dct_buf[i] <= 0.80) idx_80 = i;
    }

    double E_low = 0.0, E_high = 0.0;
    for (size_t i = 0; i <= idx_20; ++i) E_low  += energy[i];
    for (size_t i = idx_80; i < n;  ++i) E_high += energy[i];

    double eps = 1e-6 * robust_energy_sum + 0.01 * noise_floor * n;
    return log1p(E_high / (E_low + eps));
}

// =============================================================================
// 4. WOLF EWMA (Robust Kalman-style Filter)
// =============================================================================

typedef struct {
    double mu;          // State estimate (mean)
    double P;           // Innovation variance estimate
    double Q;           // Process noise
    double R;           // Measurement noise
    double c;           // Tukey cutoff threshold
    int use_hard_limit; // 0 = soft (biweight), 1 = hard (indicator)
} tw_wolf_t;

/**
 * @brief Initialize Wolf EWMA with adaptive threshold calibration
 *
 * Input:
 *   @wf       : Filter state (output)
 *   @mu0      : Initial state estimate
 *   @q        : Process noise variance
 *   @r        : Measurement noise variance
 *   @data     : Calibration data (for threshold estimation)
 *   @n        : Calibration data length
 *   @scratch  : Scratch buffer (size >= n)
 */
TW_INLINE void tw_wolf_init(tw_wolf_t* TW_RESTRICT wf, double mu0, double q, double r,
                            const double* TW_RESTRICT data, size_t n,
                            double* TW_RESTRICT scratch) {
    wf->mu = mu0;
    wf->P = 1.0;
    wf->Q = q;
    wf->R = r;
    wf->use_hard_limit = 0;

    double mad = tw_robust_mad(data, n, scratch);
    double ku = tw_moors_kurtosis(data, n, scratch);
    double sk = fabs(tw_bowley_skew(data, n, scratch));

    wf->c = TW_TUKEY_C * mad * (1.0 + 0.5 * sk) * fmax(1.0, ku / 1.233);
    if (wf->c < 1e-6) wf->c = 1.0;
}

/**
 * @brief Update Wolf EWMA with robust outlier rejection
 *
 * Input:
 *   @wf : Filter state (input/output)
 *   @y  : New measurement
 *
 * Output:
 *   Updates wf->mu (new estimate) and wf->P (new uncertainty)
 */
TW_INLINE void tw_wolf_update(tw_wolf_t* TW_RESTRICT wf, double y) {
    double p_pred = wf->P + wf->Q;
    double e = y - wf->mu;
    double abs_e = fabs(e);

    double w;
    if (wf->use_hard_limit) {
        w = (abs_e <= wf->c) ? 1.0 : 0.0;
    } else {
        w = (abs_e >= wf->c) ? 0.0 : (1.0 - (abs_e / wf->c) * (abs_e / wf->c));
        w = w * w;
    }

    double r_eff = wf->R / (w * w + TW_EPSILON);
    double K = p_pred / (p_pred + r_eff);

    wf->mu += K * e;
    wf->P = (1.0 - K) * p_pred;
}

// =============================================================================
// 5. CEEMDAN via CHEAPDAN
// =============================================================================

#include <cheapdan.h>   // cheap_ctx, cheap_init/destroy, cheap_ceemdan_frac, CHEAP_OK

/**
 * @brief Persistent processing context — caches FFTW plans and heap buffers
 *        across tw_process_ctx calls so plan creation is paid only once.
 *
 * Lifecycle:
 *   tw_ctx_t ctx;
 *   tw_ctx_init(&ctx, n, H_initial);   // allocates + builds FFTW plans
 *   tw_process_ctx(&ctx, raw, ...);    // fast — no plan creation
 *   tw_process_ctx(&ctx, raw2, ...);   // reuses plans if n and H unchanged
 *   tw_ctx_destroy(&ctx);              // frees everything
 *
 * tw_ctx_reinit is called automatically by tw_process_ctx when n or H changes.
 */
typedef struct {
    cheap_ctx  cheap;          // CHEAPDAN context (owns FFTW plans + eigenvalues)
    double*    imfs_flat;      // TW_MAX_IMFS * n doubles, IMF i at imfs_flat + i*n
    double*    energy_scratch; // n doubles, used by tw_compute_ashi
    size_t     n;              // signal length this context was built for
    double     H;              // Hurst exponent this context was built for
    int        initialized;    // 1 after successful tw_ctx_init
    // Inter-window leptokurtosis EWMA — smooths batch Moors/Bowley estimates
    // across successive tw_process_ctx calls rather than recomputing per sample.
    double     lepto_ewma;     // running smoothed leptokurtosis
    double     lepto_alpha;    // smoothing factor (TW_LEPTO_EWMA_ALPHA by default)
    int        n_calls;        // windows processed; first call bootstraps without decay
} tw_ctx_t;

TW_INLINE int tw_ctx_init(tw_ctx_t* ctx, size_t n, double H) {
    if (!ctx || n < 8 || H <= 0.0 || H >= 1.0) return -1;
    ctx->initialized = 0;
    ctx->n = n;
    ctx->H = H;
    ctx->imfs_flat      = (double*)malloc((size_t)TW_MAX_IMFS * n * sizeof(double));
    ctx->energy_scratch = (double*)malloc(n * sizeof(double));
    if (!ctx->imfs_flat || !ctx->energy_scratch) {
        free(ctx->imfs_flat);
        free(ctx->energy_scratch);
        ctx->imfs_flat = ctx->energy_scratch = NULL;
        return -1;
    }
    if (cheap_init(&ctx->cheap, (int)n, H) != CHEAP_OK) {
        free(ctx->imfs_flat);
        free(ctx->energy_scratch);
        ctx->imfs_flat = ctx->energy_scratch = NULL;
        return -1;
    }
    ctx->initialized = 1;
    ctx->lepto_ewma  = 1.233; // Gaussian baseline (Moors kurtosis of normal distribution)
    ctx->lepto_alpha = TW_LEPTO_EWMA_ALPHA;
    ctx->n_calls     = 0;
    return 0;
}

TW_INLINE void tw_ctx_destroy(tw_ctx_t* ctx) {
    if (!ctx) return;
    if (ctx->initialized) cheap_destroy(&ctx->cheap);
    free(ctx->imfs_flat);
    free(ctx->energy_scratch);
    ctx->imfs_flat = ctx->energy_scratch = NULL;
    ctx->initialized = 0;
}

/* Reinitialize when n or H changes — destroys old plans before building new ones. */
TW_INLINE int tw_ctx_reinit(tw_ctx_t* ctx, size_t n, double H) {
    tw_ctx_destroy(ctx);
    return tw_ctx_init(ctx, n, H);
}

/**
 * @brief CEEMDAN decomposition via CHEAPDAN library (context-based).
 *
 * Uses the FFTW plans already built in ctx->cheap; no plan creation overhead.
 * imfs_out must point to a buffer of at least TW_MAX_IMFS * n doubles.
 */
TW_INLINE int tw_ceemdan_ctx(tw_ctx_t* TW_RESTRICT ctx,
                              const double* TW_RESTRICT signal, size_t n,
                              double noise_amp, double noise_ar,
                              double* TW_RESTRICT imfs_out, int* n_imfs_out) {
    if (!ctx || !ctx->initialized || !signal || !imfs_out || !n_imfs_out || n < 8)
        return -1;
    if (noise_amp <= 0.0) noise_amp = 0.01;

    double H = 0.5 + noise_ar * 0.45;
    if (H <= 0.0) H = 0.01;
    if (H >= 1.0) H = 0.99;

    /* Rebuild plans if H drifted enough to matter (>1% relative change). */
    if (fabs(H - ctx->H) > 0.01 * ctx->H) {
        if (cheap_init(&ctx->cheap, (int)n, H) != CHEAP_OK) return -1;
        ctx->H = H;
    }

    int ret = cheap_ceemdan_frac(&ctx->cheap, signal, TW_CEEMDAN_ENSEMBLE,
                                  noise_amp, H, imfs_out, n_imfs_out);
    return (ret == CHEAP_OK) ? 0 : -1;
}

/**
 * @brief CEEMDAN decomposition (stateless convenience — creates/destroys context).
 *
 * Pays FFTW plan creation cost on every call. Prefer tw_ceemdan_ctx for repeated
 * calls on the same signal length.
 */
TW_INLINE int tw_ceemdan(const double* TW_RESTRICT signal, size_t n,
                          double noise_amp, double noise_ar,
                          double* TW_RESTRICT imfs_out, int* n_imfs_out) {
    if (!signal || !imfs_out || !n_imfs_out || n < 8) return -1;
    double H = 0.5 + noise_ar * 0.45;
    if (H <= 0.0) H = 0.01;
    if (H >= 1.0) H = 0.99;
    if (noise_amp <= 0.0) noise_amp = 0.01;
    cheap_ctx cheap;
    if (cheap_init(&cheap, (int)n, H) != CHEAP_OK) return -1;
    int ret = cheap_ceemdan_frac(&cheap, signal, TW_CEEMDAN_ENSEMBLE,
                                  noise_amp, H, imfs_out, n_imfs_out);
    cheap_destroy(&cheap);
    return (ret == CHEAP_OK) ? 0 : -1;
}

// =============================================================================
// 6. TOP-LEVEL API (Timberwolf Process)
// =============================================================================

/* Core pipeline — uses a pre-initialized tw_ctx_t for zero plan overhead. */
TW_INLINE int tw_process_ctx(tw_ctx_t* TW_RESTRICT ctx,
                              const double* TW_RESTRICT raw, size_t n, double temp,
                              double* TW_RESTRICT health_index,
                              double* TW_RESTRICT work_buffers,
                              int use_spectral) {
    if (!ctx || n < 8 || !raw || !health_index || !work_buffers) return -1;

    /* Reinitialize context if signal length changed. H will be set after
     * noise estimation below; use 0.5 (white noise) as a safe initial guess. */
    if (!ctx->initialized || ctx->n != n) {
        if (tw_ctx_reinit(ctx, n, 0.5) != 0) return -1;
    }

    // Buffer layout: log_signal (n) | scratch1 (n) | scratch2 (n)
    double* log_signal = work_buffers;
    double* scratch1   = work_buffers + n;
    double* scratch2   = work_buffers + 2 * n;

    // 1. Log-stabilization
    for (size_t i = 0; i < n; ++i)
        log_signal[i] = log1p(fabs(raw[i]));

    // 2. Adaptive noise parameters
    double mad = tw_robust_mad(log_signal, n, scratch1);
    double ku  = tw_moors_kurtosis(log_signal, n, scratch1);
    double sk  = fabs(tw_bowley_skew(log_signal, n, scratch1));
    double eps = 0.2 * (mad / 0.6745) * (1.0 + fmax(0.0, ku - 1.233)) * (1.0 + 0.5 * sk);
    double ar  = fmin(0.95, 0.7 + 0.25 * sk);

    // 3. CEEMDAN — uses ctx->imfs_flat (no malloc) and cached FFTW plans
    int n_imfs = 0;
    if (tw_ceemdan_ctx(ctx, log_signal, n, eps, ar, ctx->imfs_flat, &n_imfs) != 0)
        return -1;

    // 4. Softmax reconstruction weights (leptokurtosis-based)
    double weights[TW_MAX_IMFS];
    double max_w = -INFINITY;
    for (int i = 0; i < n_imfs; ++i) {
        const double* imf = ctx->imfs_flat + (size_t)i * n;
        double kurt = tw_moors_kurtosis(imf, n, scratch1);
        double skew = tw_bowley_skew(imf, n, scratch1);
        weights[i] = kurt * (1.0 + fabs(skew));
        if (weights[i] > max_w) max_w = weights[i];
    }

    double sum_w = 0.0;
    for (int i = 0; i < n_imfs; ++i) {
        weights[i] = exp((weights[i] - max_w) / temp);
        sum_w += weights[i];
    }

    for (size_t j = 0; j < n; ++j) {
        scratch1[j] = 0.0;
        for (int i = 0; i < n_imfs; ++i)
            scratch1[j] += (weights[i] / sum_w) * ctx->imfs_flat[(size_t)i * n + j];
    }

    // 5. Robust tracking + spectral fusion
    tw_wolf_t wf;
    tw_wolf_init(&wf, scratch1[0], 1e-5, mad * 0.1, scratch1, n, scratch2);

    double lepto_batch = tw_moors_kurtosis(scratch1, n, scratch2) *
                         (1.0 + fabs(tw_bowley_skew(scratch1, n, scratch2)));

    // Blend batch estimate into the inter-window EWMA.
    // First call: seed directly so there's no decayed-from-Gaussian lag.
    if (ctx->n_calls == 0)
        ctx->lepto_ewma = lepto_batch;
    else
        ctx->lepto_ewma = ctx->lepto_alpha * ctx->lepto_ewma +
                          (1.0 - ctx->lepto_alpha) * lepto_batch;
    ctx->n_calls++;

    double ashi_value = 0.0;
    if (use_spectral)
        ashi_value = tw_compute_ashi(scratch1, n, scratch2, ctx->energy_scratch, mad * mad);

    double time_component = tanh(fmax(0.0, ctx->lepto_ewma - 1.233) * 0.5);

    for (size_t i = 0; i < n; ++i) {
        tw_wolf_update(&wf, scratch1[i]);
        double fc = 1.0 - wf.P;
        health_index[i] = time_component * fc;
        if (use_spectral)
            health_index[i] += 0.3 * ashi_value * fc;
    }

    return 0;
}

/**
 * @brief Timberwolf pipeline — stateless convenience wrapper.
 *
 * Creates and destroys a tw_ctx_t on every call (pays FFTW plan cost each time).
 * Prefer tw_process_ctx when processing multiple windows of the same length.
 */
TW_INLINE int tw_process(const double* TW_RESTRICT raw, size_t n, double temp,
                         double* TW_RESTRICT health_index,
                         double* TW_RESTRICT work_buffers,
                         int use_spectral) {
    tw_ctx_t ctx;
    /* H=0.5 is overwritten after noise estimation inside tw_process_ctx */
    if (tw_ctx_init(&ctx, n, 0.5) != 0) return -1;
    int rc = tw_process_ctx(&ctx, raw, n, temp, health_index, work_buffers, use_spectral);
    tw_ctx_destroy(&ctx);
    return rc;
}

// =============================================================================
// 7. UTILITY FUNCTIONS
// =============================================================================

TW_INLINE size_t tw_scratch_size(size_t n) {
    return 3 * n * sizeof(double);
}

#endif // TIMBERWOLF_H
