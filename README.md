# Timberwolf

**Tail-Integrated Multi-band EMD with Robust Spectral Health Index**

A single-header C11 library for detecting signal degradation via high-frequency energy migration. Applies ensemble empirical mode decomposition, leptokurtosis-weighted reconstruction, and a robust Kalman tracker to produce a scalar health index per sample.

## Quick Start

```c
#include "timberwolf.h"
#include <stdlib.h>

int main(void) {
    const size_t n = 1024;
    double* raw  = malloc(n * sizeof(double));
    double* hi   = malloc(n * sizeof(double));
    double* work = malloc(tw_scratch_size(n));

    /* fill raw[] with your signal */

    tw_process(raw, n, /*temp=*/1.0, hi, work, /*use_spectral=*/1);
    /* hi[i] ≈ 0: healthy; hi[i] >> 0: degradation detected */

    free(raw); free(hi); free(work);
}
```

Compile with:
```
gcc -std=c11 -O3 -march=native -ffast-math -I./include -lfftw3 -o my_app my_app.c
```

## Build

Requires CMake ≥ 3.18 and a C11 compiler. Dependencies (FFTW3, CHEAPDAN, CHEAP) are fetched automatically.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/examples/basic_health_index | head -5
./build/benchmarks/bench_process
```

## API

```c
// Main pipeline
int tw_process(const double* raw, size_t n, double temp,
               double* health_index, double* work_buffers, int use_spectral);

// Returns required byte count for work_buffers
size_t tw_scratch_size(size_t n);  // = 3 * n * sizeof(double)
```

`temp` controls softmax sharpness: 0.1 = winner-takes-all IMF selection, 1.0 = uniform blend.
`use_spectral = 1` adds DCT-II spectral analysis (requires FFTW3 at link time).
`use_spectral = 0` uses only time-domain statistics (no FFTW3 needed at runtime).

Lower-level functions (`tw_wolf_init/update`, `tw_compute_ashi`, `tw_ceemdan`, robust statistics) are all exposed and documented inline.

---

## Mathematical Treatment

### 1. Signal Model and Decomposition Rationale

Observed signal:

```
y(t) = x(t) + η(t)
```

where `x(t)` is the true process and `η(t)` is colored noise with AR(1) coefficient `φ ∈ [0, 0.95]`. For impulsive signals (vibration, acoustic emission), the variance can be unstable across time. A log-stabilizing transform is applied first:

```
z(t) = log(1 + |y(t)|)
```

This compresses dynamic range and gives a signal whose variance is approximately stationary, making all subsequent statistics more reliable.

Classical spectral methods (FFT, wavelets) assume stationarity and fixed bases. Real machinery signals are non-stationary and non-linear. EMD decomposes adaptively into **Intrinsic Mode Functions (IMFs)** — narrowband components ordered from highest to lowest frequency — without assuming any basis. This is the motivation for the CEEMDAN stage.

### 2. CEEMDAN

**EMD sifting loop** for a single IMF:

1. Find all local maxima and minima of `z(t)`.
2. Fit cubic spline envelopes `e_upper(t)`, `e_lower(t)` through the extrema.
3. Compute mean envelope: `m(t) = (e_upper(t) + e_lower(t)) / 2`.
4. Candidate IMF: `h(t) = z(t) - m(t)`.
5. Repeat steps 1–4 on `h(t)` until stopping criterion: `||h_j - h_{j-1}||² / ||h_{j-1}||² < δ`.
6. Accept `c₁(t) = h(t)`. Residual: `r₁(t) = z(t) - c₁(t)`. Repeat on residual for subsequent IMFs.

**Mode mixing problem:** pure EMD is sensitive to noise; a single IMF can contain components from widely different frequencies. CEEMDAN (Complete Ensemble EMD with Adaptive Noise) resolves this:

```
c_k(t) = E[ EMD_k( r_{k-1}(t) + ε_{k-1} · E_{k-1}(w^(i))(t) ) ]
```

where `w^(i)` are independent white noise realizations, `E_{k-1}(·)` extracts the `(k-1)`-th IMF of its argument (scaling noise to the IMF energy level), and the expectation is over the ensemble. This ensures the noise cancels in the ensemble mean while preventing mode mixing in each realization.

**CHEAPDAN extension:** uses fractional Gaussian noise (Hurst exponent `H`) instead of white noise. Pink noise (`H → 1`) better models the self-similar noise present in rotating machinery and financial time series. The AR(1) coefficient `φ` maps to `H = 0.5 + 0.45φ`.

### 3. Softmax Reconstruction with Leptokurtosis Weighting

Not all IMFs carry equal diagnostic information. IMFs with leptokurtic (heavy-tailed) distributions and strong asymmetry tend to encode fault-related transients.

**IMF quality score** (asymmetric leptokurtosis):

```
wᵢ = Kₘ(IMFᵢ) · (1 + |Sᵦ(IMFᵢ)|)
```

where `Kₘ` is the **Moors kurtosis** (octile-based, robust to tail contamination):

```
Kₘ = [ (Q(7/8) - Q(5/8)) + (Q(3/8) - Q(1/8)) ] / (Q(3/4) - Q(1/4))
```

and `Sᵦ` is the **Bowley skewness** (quartile-based, bounded in [-1, 1]):

```
Sᵦ = (Q(3/4) + Q(1/4) - 2·Q(1/2)) / (Q(3/4) - Q(1/4))
```

Both use Type 7 linear-interpolation quantiles, making them O(N log N) but robust against outliers — unlike moment-based kurtosis and skewness which are O(N) but heavily distorted by even a single extreme observation.

**Stable softmax** with temperature `T`:

```
ŵᵢ = exp((wᵢ - max_j wⱼ) / T) / Σⱼ exp((wⱼ - max_j wⱼ) / T)
```

Subtracting the maximum prevents overflow while preserving the relative ordering. Temperature `T → 0` selects the highest-quality IMF exclusively; `T → ∞` weights all IMFs equally.

**Reconstruction:**

```
x̂(t) = Σᵢ ŵᵢ · IMFᵢ(t)
```

### 4. ASHI: Adaptive Spectral Health Index

The DCT-II (via FFTW3 `FFTW_REDFT10`) transforms the reconstructed signal to the frequency domain. The unitary form:

```
X_k = c_k · sqrt(2/N) · Σ_{n=0}^{N-1} x_n · cos[π/N · (n + 1/2) · k]

where c_0 = 1/sqrt(2),  c_{k>0} = 1
```

The DCT-II is preferred over the DFT for real signals: it is purely real-valued, has a superior energy compaction property (most energy in the first few coefficients for smooth signals), and avoids the artificial discontinuity at the signal boundaries that causes DFT spectral leakage.

**Energy spectrum:** `E_k = X_k²`

**Tukey robustification** of the energy spectrum down-weights interference spikes before binning:

```
E_k^robust = E_k · w((E_k - median(E)) / MAD(E))

where w(u) = [1 - (u/c)²]² · 𝟙(|u| ≤ c),  c = 4.685
```

**Adaptive binning:** the cumulative energy CDF is computed over frequencies. Frequency indices `f₂₀` and `f₈₀` are found where the CDF crosses 20% and 80% of total robust energy:

```
E_low  = Σ_{k ≤ f₂₀} E_k^robust
E_high = Σ_{k ≥ f₈₀} E_k^robust
```

This is *adaptive*: the band boundaries shift with the signal's spectral content rather than being fixed at e.g. 1 kHz and 10 kHz.

**Log-ratio health index:**

```
ρ = log(1 + E_high / (E_low + ε))

where ε = 10⁻⁶ · ΣE_k^robust + 0.01 · σ²_noise · N
```

`ρ ≈ 0`: energy concentrated at low frequencies (healthy). `ρ >> 0`: energy has migrated to high frequencies (fault onset, bearing defect, cavitation).

### 5. Wolf EWMA: Robust Kalman Derivation

Standard scalar Kalman filter:

```
P⁻ = P + Q                          (prediction)
K  = P⁻ / (P⁻ + R)                  (Kalman gain)
μ⁺ = μ + K·(y - μ)                  (update)
P⁺ = (1 - K)·P⁻                     (covariance update)
```

In the presence of outliers, a single large innovation `e = y - μ` can corrupt the state estimate `μ` permanently. The Wolf EWMA applies Tukey biweight down-weighting to the innovation before the update:

```
w(e) = [1 - (e/c)²]² · 𝟙(|e| ≤ c)
R_eff = R / (w(e)² + ε)
K     = P⁻ / (P⁻ + R_eff)
```

When an outlier is detected (`|e| > c`), `w → 0`, so `R_eff → ∞`, which drives `K → 0` and the state does not update. This is mathematically equivalent to saying: "an outlier provides no information about the true state."

**Adaptive threshold calibration:**

```
c = 4.685 · MAD(data) · (1 + 0.5|Sᵦ|) · max(1, Kₘ / 1.233)
```

The constant 4.685 achieves 95% efficiency relative to the Gaussian MLE at normal distributions. The skewness and kurtosis factors relax the threshold for heavy-tailed or asymmetric data, avoiding over-rejection of legitimate observations.

**Final health index:**

```
lepto         = Kₘ(x̂) · (1 + |Sᵦ(x̂)|)
time_component = tanh(max(0, lepto - 1.233) · 0.5)
confidence     = 1 - P⁺
health_index   = time_component · confidence + 0.3 · ρ · confidence
```

The `tanh` saturates the time-domain contribution; `1 - P` is a proxy for filter confidence (low uncertainty = high weight); the `0.3` spectral fusion weight was chosen empirically.

---

## Pareto Frontier: Honest Assessment

### CEEMDAN vs. Wolf EWMA alone

| Dimension | CEEMDAN (full pipeline) | Wolf EWMA only |
|---|---|---|
| Computational cost | O(N · E · K · log N); E ≈ 100 ensemble members, K ≈ 10 IMFs | O(N), streaming |
| Latency | Batch — requires full signal before producing output | Sample-by-sample, causal |
| SNR detection floor | ~3 dB — ensemble averaging suppresses noise efficiently | ~10 dB — statistical moments degrade with noise |
| Memory | O(N · K) for IMF storage (heap, via CHEAPDAN) | O(1) — two filter state variables |
| Non-stationarity handling | Handles frequency modulation, amplitude modulation, chirps | Assumes slowly-varying mean; struggles with rapid regime changes |
| Implementation complexity | High — depends on CHEAPDAN + FFTW3 | Low — pure C11, no external dependencies |

**When Wolf EWMA alone is sufficient:** impulsive fault detection (bearing defects, gear tooth fracture) where faults appear as amplitude outliers. SNR > 10 dB. Real-time systems with hard latency budgets. Embedded targets without dynamic memory allocation.

**When CEEMDAN adds value:** slowly developing faults that manifest as subtle spectral shifts rather than amplitude spikes. Non-stationary signals (variable speed drives, run-up/coast-down). SNR 3–10 dB where ensemble averaging is necessary to resolve the fault signature.

### ASHI (spectral) vs. time-domain only

| Dimension | ASHI (spectral path) | Time-domain only |
|---|---|---|
| Fault type detected | High-frequency energy migration (e.g., bearing surface fatigue, cavitation) | Amplitude/variance anomalies (impulsive faults, sudden events) |
| SNR where useful | Above ~6 dB (Tukey robustification helps, but DCT is noise-sensitive) | Always — MAD and kurtosis are robust at low SNR |
| Computational cost | O(N log N) per batch + FFTW3 plan creation overhead | Absorbed in Wolf EWMA O(N) pass |
| External dependency | FFTW3 required at link time | libc only |
| Interpretability | Physical interpretation as spectral migration ratio | Direct time-domain statistics |

**FFTW3 plan overhead note:** `tw_dct_ii` creates and destroys an `fftw_plan` on every call. For signals shorter than ~512 samples, plan creation may dominate execution time. For production use, cache the plan keyed on `n` (a `tw_spectral_ctx_t` context type is the natural extension).

### Header-only vs. compiled

The current design embeds everything in a single header with `TW_INLINE`. This means:

- **No separate compilation step** — useful for embedded and single-translation-unit builds.
- **Potential code bloat** — each translation unit that includes the header gets its own copy of all functions. For large projects, a `.c` / `.h` split is preferable.
- **Link-time optimizer friendly** — the compiler sees the full function body and can inline, vectorize, and constant-fold aggressively.

For embedded targets without FFTW3 or dynamic allocation, compile with `use_spectral = 0`. The entire CHEAPDAN + FFTW3 dependency is never referenced at runtime in that path (though it must still be present at link time in the current design — a future compile-time flag `TW_NO_SPECTRAL` could guard the includes).
